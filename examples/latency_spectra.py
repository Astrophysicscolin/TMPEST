import csv
import os
import sys
import warnings

import numpy as np
from datasets import DATA, read_lamost_mrs
from latency import COLORS, MARKERS, METHODS, N_DRAWS, OUT, PAIRS, clock, fourier, lag_one

import tmpest as tp
from tmpest.catalog import _chunk_index

SEED = 0
N_CHUNKS = 16
PERIOD = 1.5195284  # [d]
MAX_HARMONICS = 5
BANDS = {"b1": (502.0, 10.0), "b2": (517.0, 20.0), "b3": (530.0, 6.0),
         "r1": (638.0, 12.0), "r2": (656.0, 24.0), "r3": (673.0, 10.0)}  # [nm]
PAIRS_OF_BANDS = (("b1", "b2"), ("b2", "r2"))


def comparison(spectra, folder):
    blue = spectra.wavelength < 580.0
    minutes = np.rint(spectra.time * 1440.0)
    stack = []
    for fiber in sorted({n[6:9] for n in os.listdir(folder) if n.endswith(".fits.gz")}):
        try:
            comp = read_lamost_mrs(folder, fiber=fiber, normalise=False)
        except ValueError:
            continue
        ratio = comp.flux / np.median(comp.flux, axis=0)
        out = np.full(spectra.flux.shape, np.nan)
        for i, j in zip(*np.nonzero(minutes[:, None] == np.rint(comp.time * 1440.0)[None, :])):
            for sel, csel in ((blue, comp.wavelength < 580.0), (~blue, comp.wavelength >= 580.0)):
                lam = spectra.wavelength[sel]
                row = np.interp(lam, comp.wavelength[csel], ratio[j, csel])
                row[(lam < comp.wavelength[csel][0]) | (lam > comp.wavelength[csel][-1])] = np.nan
                out[i, sel] = row
        stack.append(out)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        inst = np.nanmedian(np.array(stack), axis=0)
    return np.where(np.isfinite(inst), inst, 1.0), len(stack)


def chunk_means(spectra, x):
    idx = _chunk_index(spectra.wavelength, N_CHUNKS)
    return np.stack([x[:, idx == c].mean(axis=1) for c in np.unique(idx)], axis=1)


def choose_harmonics(spectra, flux):
    night = np.floor(spectra.time)
    scores = []
    for h in range(MAX_HARMONICS + 1):
        err = []
        for n in np.unique(night):
            fit = night != n
            coef = np.linalg.lstsq(fourier(spectra.time[fit], PERIOD, h), flux[fit], rcond=None)[0]
            pred = fourier(spectra.time[~fit], PERIOD, h) @ coef
            err.append(np.mean((chunk_means(spectra, flux[~fit]) - chunk_means(spectra, pred)) ** 2))
        scores.append(float(np.sqrt(np.mean(err))))
    return int(np.argmin(scores)), scores


def baseline(spectra, folder):
    inst, n_comp = comparison(spectra, folder)
    flux = spectra.flux / inst
    flux = flux / np.median(flux, axis=0)
    harmonics, scores = choose_harmonics(spectra, flux)
    design = fourier(spectra.time, PERIOD, harmonics)
    orbit = design @ np.linalg.lstsq(design, flux, rcond=None)[0]
    print(f"[OK] {n_comp} comparison stars; orbital harmonics {harmonics} (held-out rms by harmonics: "
          + ", ".join(f"{s:.4f}" for s in scores) + ")")
    return np.median(spectra.flux, axis=0)[None, :] * inst * orbit


def real_bands(spectra, signal):
    ratio = spectra.flux / signal
    cols = []
    for centre, width in BANDS.values():
        sel = np.abs(spectra.wavelength - centre) <= 0.5 * width
        cols.append(ratio[:, sel].mean(axis=1) - 1.0)
    band = np.stack(cols, axis=1)
    return band - band.mean(axis=0)


def targets(spectra):
    names = list(BANDS)
    phot = tp.Photometry({n: spectra.time for n in names}, {n: np.ones(len(spectra.time)) for n in names},
                         {n: BANDS[n][0] for n in names}, {n: BANDS[n][1] for n in names})
    grid = tp.Spectra([BANDS[n][0] for n in names], np.ones((len(spectra.time), len(names))), spectra.time)
    return phot, grid


def draw(name, model, res, phot, grid, rng):
    if name == "TMPEST":
        noise = tp.draw(model, phot, rng)
        return np.stack([noise[n] for n in phot.names], axis=1)
    return PAIRS[name][1](model, res, grid, rng)


def summary(band):
    names = list(BANDS)
    scatter = np.std(band, axis=0)
    auto = float(np.mean([lag_one(band[:, b] - band[:, b].mean()) for b in range(band.shape[1])]))
    cross = [float(np.corrcoef(band[:, names.index(a)], band[:, names.index(b)])[0, 1]) for a, b in PAIRS_OF_BANDS]
    return scatter, auto, cross


def run():
    spectra = read_lamost_mrs(os.path.join(DATA, "lamost"), normalise=False)
    signal = baseline(spectra, os.path.join(DATA, "lamost", "comparison"))
    res = tp.scrape(spectra, signal=signal, n_chunks=N_CHUNKS)
    phot, grid = targets(spectra)
    real_scatter, real_auto, real_cross = summary(real_bands(spectra, signal))
    rows = [{"method": "real data", "backend": "", "epochs": len(spectra.time), "pixels": len(spectra.wavelength),
             "chunks": len(res.resid), "measure_ms": np.nan, "apply_ms": np.nan, "scatter_ratio": 1.0,
             "scatter_ratio_min": 1.0, "scatter_ratio_max": 1.0, "lag_one": real_auto,
             "cross_within_arm": real_cross[0], "cross_between_arms": real_cross[1]}]
    for name in METHODS:
        measure = PAIRS[name][0]
        backends = ["rust", "numpy"] if name == "TMPEST" and tp.rust_available() else ["numpy"]
        for backend in backends:
            if name == "TMPEST":
                tp.set_backend(backend)
            rng = np.random.default_rng(SEED)
            model = measure(res)
            t_measure = clock(lambda: measure(res))
            t_apply = clock(lambda: draw(name, model, res, phot, grid, rng))
            ratio, auto, cross = [], [], []
            for _ in range(N_DRAWS):
                s, a, c = summary(draw(name, model, res, phot, grid, rng))
                ratio.append(s)
                auto.append(a)
                cross.append(c)
            ratio = np.mean(ratio, axis=0) / real_scatter
            cross = np.mean(cross, axis=0)
            rows.append({"method": name, "backend": backend, "epochs": len(spectra.time),
                         "pixels": len(spectra.wavelength), "chunks": len(res.resid), "measure_ms": 1e3 * t_measure,
                         "apply_ms": 1e3 * t_apply, "scatter_ratio": float(np.mean(ratio)),
                         "scatter_ratio_min": float(ratio.min()), "scatter_ratio_max": float(ratio.max()),
                         "lag_one": float(np.mean(auto)), "cross_within_arm": float(cross[0]),
                         "cross_between_arms": float(cross[1])})
    if tp.rust_available():
        tp.set_backend("rust")
    return rows


def report(rows):
    print(f"{'method':<17}{'backend':<8}{'measure [ms]':>13}{'apply [ms]':>12}{'scatter':>9}{'min':>7}{'max':>7}"
          f"{'lag-1':>8}{'in arm':>8}{'B-R':>7}")
    for r in rows:
        print(f"{r['method']:<17}{r['backend']:<8}{r['measure_ms']:>13.3f}{r['apply_ms']:>12.3f}"
              f"{r['scatter_ratio']:>9.2f}{r['scatter_ratio_min']:>7.2f}{r['scatter_ratio_max']:>7.2f}"
              f"{r['lag_one']:>8.2f}{r['cross_within_arm']:>8.2f}{r['cross_between_arms']:>7.2f}")
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "latency_lamost.csv")
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[SAVE] {path}")


def pick(rows, method):
    return next(r for r in rows if r["method"] == method and r["backend"] in ("rust", "numpy", ""))


def figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    real = pick(rows, "real data")
    fig, axes = plt.subplots(1, 5, figsize=(16, 4.4))
    y = np.arange(len(METHODS))
    for ax, column, title in zip(axes[:2], ("measure_ms", "apply_ms"), ("Measure the noise [ms]", "Apply to photometry [ms]")):
        x = [pick(rows, m)[column] for m in METHODS]
        ax.barh(y, x, height=0.7, color=COLORS)
        for yi, xi in zip(y, x):
            ax.text(xi * 1.15, yi, f"{xi:.0f}" if xi >= 100 else f"{xi:.2g}", va="center", ha="left", fontsize=8)
        ax.set_xscale("log")
        ax.set_xlim(right=ax.get_xlim()[1] * 4.0)
        ax.set_yticks(y)
        ax.set_yticklabels(METHODS if ax is axes[0] else [])
        ax.invert_yaxis()
        ax.set_title(title, loc="left", fontsize=10)
        ax.grid(axis="x", color="0.9", linewidth=0.8)
    panels = (("scatter_ratio", "Band scatter, synthetic / real", 1.0),
              ("lag_one", "Lag-one autocorrelation", real["lag_one"]),
              ("cross_between_arms", "Correlation, blue and red arm", real["cross_between_arms"]))
    for ax, (column, title, reference) in zip(axes[2:], panels):
        for m, (name, color) in enumerate(zip(METHODS, COLORS)):
            r = pick(rows, name)
            if column == "scatter_ratio":
                ax.plot([m, m], [r["scatter_ratio_min"], r["scatter_ratio_max"]], color=color, linewidth=2)
            ax.plot(m, r[column], marker=MARKERS[0], color=color, markersize=8, markeredgecolor="white",
                    linestyle="none")
        ax.plot([-0.5, len(METHODS) - 0.5], [reference, reference], color="0.35", linewidth=1.0)
        ax.text(len(METHODS) - 0.5, reference, " real", va="center", ha="left", fontsize=8, color="0.35")
        ax.set_xticks(np.arange(len(METHODS)))
        ax.set_xticklabels(METHODS, fontsize=8, rotation=35, ha="right")
        ax.set_title(title, loc="left", fontsize=10)
        ax.grid(axis="y", color="0.9", linewidth=0.8)
    for ax in axes:
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle(f"Spectra to photometry, LAMOST spectra of KIC 5444392: {real['epochs']} epochs, {real['pixels']} pixels, "
                 f"{len(BANDS)} synthetic bands", x=0.01, y=0.99, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = os.path.join(OUT, "latency_lamost.png")
    fig.savefig(path, dpi=160)
    print(f"[SAVE] {path}")


if __name__ == "__main__":
    results = run()
    report(results)
    if "--no-figures" not in sys.argv:
        figure(results)
