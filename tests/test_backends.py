import numpy as np
import pytest
import tmpest as tp
from tmpest import _kernels

needs_rust = pytest.mark.skipif(not tp.rust_available(), reason="compiled extension not built")

@pytest.fixture
def both():
    start = tp.backend()

    def run(fn):
        out = {}
        for name in ("numpy", "rust"):
            tp.set_backend(name)
            out[name] = fn()
        return out["numpy"], out["rust"]

    yield run
    tp.set_backend(start)

def test_backend_is_reported_and_switchable():
    start = tp.backend()
    assert start in ("rust", "numpy")
    tp.set_backend("numpy")
    assert tp.backend() == "numpy" and tp.num_threads() == 1
    with pytest.raises(ValueError):
        tp.set_backend("fortran")
    if not tp.rust_available():
        with pytest.raises(ImportError):
            tp.set_backend("rust")
    tp.set_backend(start)

@needs_rust
@pytest.mark.parametrize("shape", [(1, 1), (2, 3), (500, 7), (7, 500), (300, 300)])
@pytest.mark.parametrize("scale", [0.5, 30.0, np.inf])
def test_ou_filter_matches_numpy(both, shape, scale):
    rng = np.random.default_rng(0)
    white = rng.standard_normal(shape)
    for axis in (0, 1):
        x = np.cumsum(rng.uniform(0.0, 2.0, shape[axis]))
        ref, fast = both(lambda: tp.ou_filter(white, x, scale, axis=axis))
        assert np.array_equal(ref, fast)
    x = np.cumsum(rng.uniform(0.0, 2.0, shape[0]))
    ref, fast = both(lambda: tp.ou_filter(white[:, 0], x, scale))
    assert ref.shape == (shape[0],) and np.array_equal(ref, fast)

@needs_rust
def test_ou_filter_leaves_its_input_alone(both):
    white = np.random.default_rng(1).standard_normal((50, 4))
    keep = white.copy()
    both(lambda: tp.ou_filter(white, np.arange(50.0), 3.0))
    assert np.array_equal(white, keep)

@needs_rust
@pytest.mark.parametrize("n", [2, 11, 400, 5000])
def test_variogram_bins_match_numpy(both, n):
    rng = np.random.default_rng(2)
    t = np.sort(np.round(rng.uniform(0.0, 100.0, n), 2))
    r = rng.standard_normal(n)
    lags = np.unique(np.round(np.geomspace(1, max(n // 2, 1), 40)).astype(int))
    edges = np.geomspace(0.01, 100.0, 25)
    ref, fast = both(lambda: _kernels.variogram_bins(t, r, lags, edges))
    for a, b in zip(ref, fast):
        assert np.allclose(a, b, rtol=1e-12, atol=0.0)
    assert np.array_equal(ref[2], fast[2])

@needs_rust
@pytest.mark.parametrize("n,k", [(4, 3), (50, 3), (50, 49), (2000, 101), (30000, 301)])
def test_running_median_matches_numpy(both, n, k):
    y = np.random.default_rng(3).standard_normal(n)
    y[: n // 3] = np.round(y[: n // 3], 1)
    ref, fast = both(lambda: _kernels.running_median(y, k))
    assert np.array_equal(ref, fast)

@needs_rust
def test_whole_pipeline_gives_the_same_numbers_on_both_backends(both):
    def run():
        rng = np.random.default_rng(4)
        t = np.sort(rng.uniform(0.0, 500.0, 800))
        names = ["gaia_G", "gaia_BP", "gaia_RP"]
        red = 4e-3 * tp.ou_filter(rng.standard_normal(800), t, 6.0)
        phot = tp.Photometry({n: t for n in names},
                             {n: 1.0 + red + 1e-3 * rng.standard_normal(800) for n in names})
        prof = tp.catalog(phot, window=80.0)
        spec = tp.Spectra(np.linspace(400.0, 900.0, 120), np.ones((60, 120)), np.sort(rng.uniform(0.0, 500.0, 60)))
        noisy = tp.apply(prof, spec, rng)
        shuffled = tp.Photometry({"tess": t[::-1]}, {"tess": np.ones(800)})
        return prof, noisy.flux, tp.apply(prof, shuffled, rng).flux["tess"]

    (p_ref, s_ref, f_ref), (p_fast, s_fast, f_fast) = both(run)
    assert p_ref.tau == pytest.approx(p_fast.tau, rel=1e-9)
    assert np.allclose(p_ref.red_var, p_fast.red_var, rtol=1e-9)
    assert np.allclose(p_ref.white_var, p_fast.white_var, rtol=1e-9)
    assert np.allclose(s_ref, s_fast, rtol=1e-9) and np.allclose(f_ref, f_fast, rtol=1e-9)