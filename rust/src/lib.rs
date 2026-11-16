use ndarray::parallel::prelude::*;
use ndarray::{ArrayView1, ArrayViewMut2, Axis, Zip};
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadwriteArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

const MIN_PARALLEL: usize = 1 << 16;

fn sweep_rows(mut block: ArrayViewMut2<f64>, a: ArrayView1<f64>, b: ArrayView1<f64>) {
    let mut rows = block.rows_mut().into_iter();
    let mut prev = match rows.next() {
        Some(row) => row,
        None => return,
    };
    for (k, mut cur) in rows.enumerate() {
        let (ak, bk) = (a[k], b[k]);
        Zip::from(&mut cur).and(&prev).for_each(|c, &p| *c = ak * p + bk * *c);
        prev = cur;
    }
}

fn filter(mut data: ArrayViewMut2<f64>, a: ArrayView1<f64>, b: ArrayView1<f64>, axis: usize) {
    let (n_rows, n_cols) = data.dim();
    let parallel = n_rows * n_cols >= MIN_PARALLEL;
    if axis == 0 {
        let threads = rayon::current_num_threads().max(1);
        let chunk = n_cols.div_ceil(threads).max(1);
        if parallel && n_cols > 1 {
            data.axis_chunks_iter_mut(Axis(1), chunk)
                .into_par_iter()
                .for_each(|block| sweep_rows(block, a, b));
        } else {
            sweep_rows(data, a, b);
        }
    } else {
        let along = |mut row: ndarray::ArrayViewMut1<f64>| {
            for j in 1..row.len() {
                row[j] = a[j - 1] * row[j - 1] + b[j - 1] * row[j];
            }
        };
        if parallel && n_rows > 1 {
            data.axis_iter_mut(Axis(0)).into_par_iter().for_each(along);
        } else {
            data.axis_iter_mut(Axis(0)).for_each(along);
        }
    }
}

#[pyfunction]
fn ou_filter_inplace(
    py: Python<'_>,
    mut data: PyReadwriteArray2<f64>,
    a: PyReadonlyArray1<f64>,
    b: PyReadonlyArray1<f64>,
    axis: usize,
) -> PyResult<()> {
    if axis > 1 {
        return Err(PyValueError::new_err("axis must be 0 or 1"));
    }
    let view = data.as_array_mut();
    let (a, b) = (a.as_array(), b.as_array());
    let length = view.len_of(Axis(axis));
    if length > 0 && (a.len() != length - 1 || b.len() != length - 1) {
        return Err(PyValueError::new_err("a and b must have one element fewer than the filtered axis"));
    }
    py.detach(|| filter(view, a, b, axis));
    Ok(())
}

fn bin_index(edges: &[f64], value: f64) -> usize {
    let right = edges.partition_point(|&e| e <= value);
    right.saturating_sub(1).min(edges.len() - 2)
}

type Bins = (Vec<f64>, Vec<f64>, Vec<f64>);

fn accumulate(t: &[f64], r: &[f64], lag: usize, edges: &[f64]) -> Bins {
    let n_bins = edges.len() - 1;
    let mut out = (vec![0.0; n_bins], vec![0.0; n_bins], vec![0.0; n_bins]);
    for i in 0..t.len().saturating_sub(lag) {
        let dt = t[i + lag] - t[i];
        if dt > 0.0 {
            let d = r[i + lag] - r[i];
            let k = bin_index(edges, dt);
            out.0[k] += dt;
            out.1[k] += 0.5 * d * d;
            out.2[k] += 1.0;
        }
    }
    out
}

type BinArrays<'py> = (Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>);

#[pyfunction]
fn variogram_bins<'py>(
    py: Python<'py>,
    t: PyReadonlyArray1<f64>,
    r: PyReadonlyArray1<f64>,
    lags: Vec<usize>,
    edges: PyReadonlyArray1<f64>,
) -> PyResult<BinArrays<'py>> {
    let (t, r, edges) = (t.as_slice()?, r.as_slice()?, edges.as_slice()?);
    if t.len() != r.len() {
        return Err(PyValueError::new_err("t and r must have the same length"));
    }
    if edges.len() < 2 {
        return Err(PyValueError::new_err("edges needs at least two values"));
    }
    let n_bins = edges.len() - 1;
    let parts: Vec<Bins> = py.detach(|| lags.par_iter().map(|&lag| accumulate(t, r, lag, edges)).collect());
    let mut total = (vec![0.0; n_bins], vec![0.0; n_bins], vec![0.0; n_bins]);
    for part in parts {
        for k in 0..n_bins {
            total.0[k] += part.0[k];
            total.1[k] += part.1[k];
            total.2[k] += part.2[k];
        }
    }
    Ok((total.0.into_pyarray(py), total.1.into_pyarray(py), total.2.into_pyarray(py)))
}

fn reflect(i: isize, n: isize) -> usize {
    if n == 1 {
        return 0;
    }
    let period = 2 * (n - 1);
    let m = i.rem_euclid(period);
    (if m < n { m } else { period - m }) as usize
}

fn insert_sorted(window: &mut Vec<f64>, value: f64) {
    let at = window.partition_point(|&w| w < value);
    window.insert(at, value);
}

fn remove_sorted(window: &mut Vec<f64>, value: f64) {
    let at = window.partition_point(|&w| w < value);
    window.remove(at);
}

fn median_span(y: &[f64], k: usize, start: usize, stop: usize, out: &mut [f64]) {
    let half = (k / 2) as isize;
    let n = y.len() as isize;
    let mut window: Vec<f64> = Vec::with_capacity(k);
    for j in -half..=half {
        insert_sorted(&mut window, y[reflect(start as isize + j, n)]);
    }
    out[0] = window[k / 2];
    for i in (start + 1)..stop {
        let i = i as isize;
        remove_sorted(&mut window, y[reflect(i - 1 - half, n)]);
        insert_sorted(&mut window, y[reflect(i + half, n)]);
        out[i as usize - start] = window[k / 2];
    }
}

#[pyfunction]
fn running_median<'py>(py: Python<'py>, y: PyReadonlyArray1<f64>, k: usize) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let y = y.as_slice()?;
    let n = y.len();
    if k % 2 == 0 || k < 3 || k >= n {
        return Err(PyValueError::new_err("k must be odd, at least 3 and smaller than the series"));
    }
    if y.iter().any(|v| !v.is_finite()) {
        return Err(PyValueError::new_err("y must be finite"));
    }
    let mut out = vec![0.0; n];
    py.detach(|| {
        let span = (n / rayon::current_num_threads().max(1)).max(4 * k).max(1);
        out.par_chunks_mut(span)
            .enumerate()
            .for_each(|(c, piece)| median_span(y, k, c * span, c * span + piece.len(), piece));
    });
    Ok(out.into_pyarray(py))
}

#[pyfunction]
fn num_threads() -> usize {
    rayon::current_num_threads()
}

#[pymodule(gil_used = false)]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ou_filter_inplace, m)?)?;
    m.add_function(wrap_pyfunction!(variogram_bins, m)?)?;
    m.add_function(wrap_pyfunction!(running_median, m)?)?;
    m.add_function(wrap_pyfunction!(num_threads, m)?)?;
    Ok(())
}
