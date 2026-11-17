import csv
import os
import sys
import time

import numpy as np
from datasets import read_hoffman

import tmpest as tp
from tmpest.catalog import _nearest

SEED = 0
STARS = {"LX": ("LX Leo", 0.235247), "MU": ("MU Leo", 0.388442), "V0354": ("V354 UMa", 0.293825)}  # [d]
N_FOURIER = 12
N_PIXELS = 600
LAMBDA = (380.0, 730.0)  # [nm]
REPEATS = 50
N_DRAWS = 200
MATCH_TOL = 0.5
METHODS = ["TMPEST", "white Gaussian", "block bootstrap", "Timmer-Koenig", "dense GP"]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
MARKERS = ["o", "s", "^"]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")


# INPUT
def fourier(t, period, harmonics):
    phase = 2.0 * np.pi * t / period
    cols = [np.ones_like(t)]
    for k in range(1, harmonics + 1):
        cols += [np.cos(k * phase), np.sin(k * phase)]
    return np.stack(cols, axis=1)


def fourier_signal(t, flux, period):
    design = fourier(t, period, N_FOURIER)
    return design @ np.linalg.lstsq(design, flux, rcond=None)[0]


def nearest_band(res, lam):
    return np.argmin(np.abs(np.log(lam)[:, None] - np.log(res.wavelength)[None, :]), axis=1)


# WHITE GAUSSIAN
def white_measure(res):
    return np.array([np.std(r) for r in res.resid])


def white_apply(sigma, res, target, rng):
    order = np.argsort(res.wavelength)
    s = np.exp(np.interp(np.log(target.wavelength), np.log(res.wavelength[order]), np.log(sigma[order])))
    return rng.standard_normal(target.flux.shape) * s[None, :]


# BLOCK BOOTSTRAP
def block_measure(res):
    return [int(np.ceil(len(r) ** (1.0 / 3.0))) for r in res.resid]


def block_apply(blocks, res, target, rng):
    n_t = len(target.time)
    series = np.empty((n_t, len(res.resid)))
    for b, r in enumerate(res.resid):
        size = blocks[b]
        starts = rng.integers(0, len(r) - size + 1, int(np.ceil(n_t / size)))
        series[:, b] = r[(starts[:, None] + np.arange(size)[None, :]).ravel()[:n_t]]
    return series[:, nearest_band(res, target.wavelength)]


# TIMMER KOENIG
def tk_measure(res):
    fits = []
    for t, r in zip(res.time, res.resid):
        step = float(np.median(np.diff(t)))
        grid = np.arange(t[0], t[-1], step)
        power = np.abs(np.fft.rfft(np.interp(grid, t, r))[1:]) ** 2
        freq = np.fft.rfftfreq(len(grid), step)[1:]
        ok = power > 0
        beta = float(np.clip(-np.polyfit(np.log(freq[ok]), np.log(power[ok]), 1)[0], 0.0, 3.0))
        fits.append((beta, step, float(np.std(r))))
    return fits


def tk_apply(fits, res, target, rng):
    t = target.time
    series = np.empty((len(t), len(fits)))
    for b, (beta, step, sigma) in enumerate(fits):
        n = max(int(np.ceil((t[-1] - t[0]) / step)) + 1, 4)
        freq = np.fft.rfftfreq(n, step)[1:]
        amp = np.sqrt(0.5 * freq ** -beta)
        spec = np.concatenate([[0.0], amp * (rng.standard_normal(len(freq)) + 1j * rng.standard_normal(len(freq)))])
        x = np.fft.irfft(spec, n)
        x *= sigma / np.std(x)
        series[:, b] = x[np.clip(np.rint((t - t[0]) / step).astype(int), 0, n - 1)]
    return series[:, nearest_band(res, target.wavelength)]


# DENSE GP
def gp_measure(res):
    fits = []
    for t, r in zip(res.time, res.resid):
        var = float(np.var(r))
        lag = np.abs(t[:, None] - t[None, :])
        eye = np.eye(len(t))
        best = (-np.inf, 0.0, 1.0)
        for tau in np.geomspace(float(np.median(np.diff(t))), float(t[-1] - t[0]), 20):
            kernel = np.exp(-lag / tau)
            for frac in np.linspace(0.0, 0.95, 20):
                chol = np.linalg.cholesky(var * (frac * kernel + (1.0 - frac) * eye))
                alpha = np.linalg.solve(chol, r)
                loglike = -0.5 * float(alpha @ alpha) - float(np.sum(np.log(np.diag(chol))))
                if loglike > best[0]:
                    best = (loglike, float(frac), float(tau))
        fits.append((var, best[1], best[2]))
    return fits


def gp_apply(fits, res, target, rng):
    t = target.time
    lag = np.abs(t[:, None] - t[None, :])
    series = np.empty((len(t), len(fits)))
    for b, (var, frac, tau) in enumerate(fits):
        chol = np.linalg.cholesky(var * (frac * np.exp(-lag / tau) + (1.0 - frac) * np.eye(len(t))))
        series[:, b] = chol @ rng.standard_normal(len(t))
    return series[:, nearest_band(res, target.wavelength)]


# TMPEST
def tmpest_measure(res):
    return tp.catalog(res, match_tol=MATCH_TOL * float(np.median(np.diff(res.time[0]))))


def tmpest_apply(profile, res, target, rng):
    return tp.draw(profile, target, rng)


PAIRS = {"TMPEST": (tmpest_measure, tmpest_apply), "white Gaussian": (white_measure, white_apply),
         "block bootstrap": (block_measure, block_apply), "Timmer-Koenig": (tk_measure, tk_apply),
         "dense GP": (gp_measure, gp_apply)}


# TIMING AND FIDELITY
def clock(fn):
    times = []
    for _ in range(REPEATS):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    return float(np.median(times))


def lag_one(x):
    return float(np.mean(x[1:] * x[:-1]) / np.mean(x * x))


def matched_corr(ta, ra, tb, rb):
    j = _nearest(ta, tb)
    return float(np.mean(ra * rb[j]) / np.sqrt(np.mean(ra * ra) * np.mean(rb[j] * rb[j])))


def band_means(res, target, noise):
    half = 0.5 * res.width
    cols = [noise[:, np.abs(target.wavelength - c) <= h].mean(axis=1) for c, h in zip(res.wavelength, half)]
    return np.stack(cols, axis=1)


def fidelity(apply, model, res, target, rng):
    ratio, auto, cross = [], [], []
    for _ in range(N_DRAWS):
        band = band_means(res, target, apply(model, res, target, rng))
        ratio.append(np.std(band, axis=0))
        auto.append(np.mean([lag_one(band[:, b] - band[:, b].mean()) for b in range(band.shape[1])]))
        cross.append(np.corrcoef(band[:, 0], band[:, -1])[0, 1])
    measured = np.array([np.std(r) for r in res.resid])
    return float(np.mean(np.mean(ratio, axis=0) / measured)), float(np.mean(auto)), float(np.mean(cross))


def run():
    rows = []
    for key, (label, period) in STARS.items():
        phot = read_hoffman(key)
        signal = {b: fourier_signal(phot.time[b], phot.flux[b], period) for b in phot.names}
        res = tp.scrape(phot, signal=signal)
        lam = np.linspace(LAMBDA[0], LAMBDA[1], N_PIXELS)
        epochs = res.time[1]
        target = tp.Spectra(lam, np.ones((len(epochs), N_PIXELS)), epochs)
        real_auto = float(np.mean([lag_one(r) for r in res.resid]))
        real_cross = matched_corr(res.time[0], res.resid[0], res.time[-1], res.resid[-1])
        rows.append({"star": label, "method": "real data", "backend": "", "epochs": len(epochs),
                     "pixels": N_PIXELS, "measure_ms": np.nan, "apply_ms": np.nan, "scatter_ratio": 1.0,
                     "lag_one": real_auto, "cross_band": real_cross})
        for name in METHODS:
            measure, apply = PAIRS[name]
            backends = ["rust", "numpy"] if name == "TMPEST" and tp.rust_available() else [tp.backend()]
            for backend in backends:
                if name == "TMPEST":
                    tp.set_backend(backend)
                rng = np.random.default_rng(SEED)
                model = measure(res)
                t_measure = clock(lambda: measure(res))
                t_apply = clock(lambda: apply(model, res, target, rng))
                ratio, auto, cross = fidelity(apply, model, res, target, rng)
                rows.append({"star": label, "method": name, "backend": backend if name == "TMPEST" else "numpy",
                             "epochs": len(epochs), "pixels": N_PIXELS, "measure_ms": 1e3 * t_measure,
                             "apply_ms": 1e3 * t_apply, "scatter_ratio": ratio, "lag_one": auto,
                             "cross_band": cross})
        if tp.rust_available():
            tp.set_backend("rust")
    return rows


# OUTPUT
def report(rows):
    print(f"{'star':<10}{'method':<17}{'backend':<8}{'measure [ms]':>13}{'apply [ms]':>12}"
          f"{'scatter':>9}{'lag-1':>8}{'B-R corr':>10}")
    for r in rows:
        print(f"{r['star']:<10}{r['method']:<17}{r['backend']:<8}{r['measure_ms']:>13.3f}{r['apply_ms']:>12.3f}"
              f"{r['scatter_ratio']:>9.2f}{r['lag_one']:>8.2f}{r['cross_band']:>10.2f}")
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "latency_hoffman.csv")
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[SAVE] {path}")


def pick(rows, star, method):
    return next(r for r in rows if r["star"] == star and r["method"] == method and r["backend"] in ("rust", "numpy", ""))


def figures(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stars = [v[0] for v in STARS.values()]
    height = 0.8 / len(METHODS)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0), sharey=True)
    for ax, column, title in zip(axes, ("measure_ms", "apply_ms"), ("Measure the noise", "Apply it to synthetic spectra")):
        for m, (name, color) in enumerate(zip(METHODS, COLORS)):
            y = np.arange(len(stars)) + (m - 0.5 * (len(METHODS) - 1)) * height
            x = [pick(rows, s, name)[column] for s in stars]
            ax.barh(y, x, height=0.86 * height, color=color, label=name)
            for yi, xi in zip(y, x):
                ax.text(xi * 1.12, yi, f"{xi:.0f}" if xi >= 100 else f"{xi:.2g}", va="center", ha="left", fontsize=8)
        ax.set_xscale("log")
        ax.set_xlim(right=ax.get_xlim()[1] * 3.0)
        ax.set_xlabel("median time per call [ms]")
        ax.set_title(title, loc="left")
        ax.set_yticks(np.arange(len(stars)))
        ax.set_yticklabels([f"{s}\n{pick(rows, s, 'real data')['epochs']} epochs" for s in stars])
        ax.invert_yaxis()
        ax.grid(axis="x", color="0.9", linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="upper left", bbox_to_anchor=(0.005, 0.95), ncol=len(METHODS), fontsize=9)
    fig.suptitle(f"Photometry to spectra, Hoffman et al. (2022) light curves, {N_PIXELS} pixels", x=0.01, y=0.99, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    path = os.path.join(OUT, "latency_hoffman.png")
    fig.savefig(path, dpi=160)
    print(f"[SAVE] {path}")

    panels = (("scatter_ratio", "Band scatter, synthetic / measured", 1.0),
              ("lag_one", "Lag-one autocorrelation in time", None),
              ("cross_band", "Correlation between B and R", None))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
    x = np.arange(len(METHODS))
    for ax, (column, title, reference) in zip(axes, panels):
        for s, (star, marker) in enumerate(zip(stars, MARKERS)):
            offset = (s - 1) * 0.2
            for m, (name, color) in enumerate(zip(METHODS, COLORS)):
                ax.plot(m + offset, pick(rows, star, name)[column], marker=marker, color=color, markersize=8,
                        markeredgecolor="white", linestyle="none")
            real = reference if reference is not None else pick(rows, star, "real data")[column]
            ax.plot([-0.5, len(METHODS) - 0.5], [real, real], color="0.35", linewidth=1.0,
                    linestyle=("-", "--", ":")[s])
        ax.set_xticks(x)
        ax.set_xticklabels([n.replace(" ", "\n").replace("-", "-\n") for n in METHODS], fontsize=8)
        ax.set_title(title, loc="left", fontsize=10)
        ax.grid(axis="y", color="0.9", linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    handles = [plt.Line2D([], [], marker=mk, color="0.35", linestyle=ls, markersize=7, label=f"{st} (line: real data)")
               for st, mk, ls in zip(stars, MARKERS, ("-", "--", ":"))]
    fig.legend(handles=handles, frameon=False, fontsize=9, loc="upper right", ncol=3)
    fig.suptitle(f"What each method reproduces (mean of {N_DRAWS} draws)", x=0.01, y=0.99, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = os.path.join(OUT, "fidelity_hoffman.png")
    fig.savefig(path, dpi=160)
    print(f"[SAVE] {path}")


if __name__ == "__main__":
    results = run()
    report(results)
    if "--no-figures" not in sys.argv:
        figures(results)
