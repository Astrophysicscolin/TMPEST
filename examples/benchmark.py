import sys
import time
import numpy as np
import tmpest as tp
from tmpest import _kernels

SEED = 0
SCALE = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0

def best(fn, repeats=3):
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    return min(times)

def jobs():
    rng = np.random.default_rng(SEED)
    n_long = int(1_000_000 * SCALE)
    n_t, n_lam = int(2000 * SCALE ** 0.5), int(8000 * SCALE ** 0.5)
    t = np.sort(rng.uniform(0.0, 4.0 * n_long / 1000.0, n_long))
    y = 1.0 + 1e-3 * rng.standard_normal(n_long)
    r = y - 1.0
    white = rng.standard_normal((n_t, n_lam))
    x_t = np.sort(rng.uniform(0.0, 900.0, n_t))
    x_lam = np.log(np.linspace(380.0, 1000.0, n_lam))
    lags = np.unique(np.round(np.geomspace(1, n_long // 2, 40)).astype(int))
    edges = np.geomspace(1e-3, t[-1] - t[0], 25)
    long_series = tp.Photometry({"kepler": t}, {"kepler": y})
    names = ["gaia_G", "gaia_BP", "gaia_RP"]
    tt = np.sort(rng.uniform(0.0, 900.0, 3000))
    red = 2e-3 * tp.ou_filter(rng.standard_normal(3000), tt, 5.0)
    profile = tp.catalog(tp.Photometry({n: tt for n in names},
                                       {n: 1.0 + red + 1e-3 * rng.standard_normal(3000) for n in names}))
    spec = tp.Spectra(np.exp(x_lam), np.ones((n_t, n_lam)), x_t)
    return [
        (f"kernel: OU filter along time, {n_t} x {n_lam}", lambda: _kernels.ou_filter_inplace(white.copy(), x_t, 5.0, 0)),
        (f"kernel: OU filter along wavelength, {n_t} x {n_lam}", lambda: _kernels.ou_filter_inplace(white.copy(), x_lam, 0.5, 1)),
        (f"kernel: OU filter, {n_long} x 9", lambda: _kernels.ou_filter_inplace(np.ones((n_long, 9)), t, 5.0, 0)),
        (f"kernel: variogram, {n_long} points", lambda: _kernels.variogram_bins(t, r, lags, edges)),
        (f"kernel: running median, {n_long} points, window 501", lambda: _kernels.running_median(y, 501)),
        (f"catalog: 1 band, {n_long} points", lambda: tp.catalog(long_series)),
        (f"catalog: 1 band, {n_long} points, window detrend", lambda: tp.catalog(long_series, window=2.0)),
        (f"apply: spectra {n_t} x {n_lam}", lambda: tp.apply(profile, spec, np.random.default_rng(1))),
        (f"apply: photometry 1 band, {n_long} points", lambda: tp.apply(profile, long_series, np.random.default_rng(1))),
    ]

def main():
    start = tp.backend()
    names = ["numpy", "rust"] if tp.rust_available() else ["numpy"]
    print(f"tmpest {tp.__version__}, numpy {np.__version__}, rust extension: {tp.rust_available()}")
    table = {}
    for name in names:
        tp.set_backend(name)
        if name == "rust":
            print(f"rust threads: {tp.num_threads()}")
        for label, fn in jobs():
            table.setdefault(label, {})[name] = best(fn)
    tp.set_backend(start)
    print(f"{'job':<58} {'numpy [s]':>10} {'rust [s]':>10} {'speed-up':>9}")
    for label, row in table.items():
        fast = row.get("rust")
        ratio = f"{row['numpy'] / fast:>8.1f}x" if fast else f"{'-':>9}"
        print(f"{label:<58} {row['numpy']:>10.3f} {fast if fast is not None else float('nan'):>10.3f} {ratio}")

if __name__ == "__main__":
    main()