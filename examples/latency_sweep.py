import csv
import os
import sys
import time

import numpy as np
from latency import COLORS, METHODS, N_PIXELS, LAMBDA, OUT, PAIRS

import tmpest as tp

SEED = 0
SIZES = [100, 300, 1000, 3000, 10000, 30000]
CADENCE = 259.0 / 86400.0  # [d]
OFFSET = 86.0 / 86400.0  # [d]
SIGMA = {"johnson_B": 0.0083, "johnson_V": 0.0074, "cousins_R": 0.0057}
SHARED = 0.35
BUDGET = 20.0  # [s]
GP_MEASURE_MAX = 1000
GP_APPLY_MAX = 3000
LABEL_GAP = 0.045
RED = "--red" in sys.argv
RED_FRACTION = 0.5
TAU = 20.0 * CADENCE  # [d]
PIXELS = "--pixels" in sys.argv
PIXEL_SIZES = [100, 300, 1000, 3000, 10000]
PIXEL_EPOCHS = 1000
CASES = [(PIXEL_EPOCHS, p) for p in PIXEL_SIZES] if PIXELS else [(n, N_PIXELS) for n in SIZES]
AXIS = "pixels" if PIXELS else "epochs"
SUFFIX = ("_red" if RED else "") + ("_pixels" if PIXELS else "")


def residuals(n, rng):
    names = list(SIGMA)
    common = rng.standard_normal(n)
    times = [CADENCE * np.arange(n) + OFFSET * b for b in range(len(names))]
    resids = [SIGMA[name] * (np.sqrt(SHARED) * common + np.sqrt(1.0 - SHARED) * rng.standard_normal(n))
              for name in names]
    if RED:
        red = tp.ou_filter(rng.standard_normal(n), times[0], TAU)
        resids = [np.sqrt(1.0 - RED_FRACTION) * r + np.sqrt(RED_FRACTION) * SIGMA[name] * red
                  for r, name in zip(resids, names)]
    lam = np.array([tp.BANDS[name][0] for name in names])
    width = np.array([tp.BANDS[name][1] for name in names])
    return tp.Residuals(times, [r - r.mean() for r in resids], lam, width, names, "photometry")


def clock(fn):
    start = time.perf_counter()
    fn()
    first = time.perf_counter() - start
    repeats = int(np.clip(BUDGET / max(first, 1e-9), 0, 20))
    times = [first]
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    return float(np.median(times))


def run():
    rows = []
    for n, pixels in CASES:
        rng = np.random.default_rng(SEED)
        res = residuals(n, rng)
        target = tp.Spectra(np.linspace(LAMBDA[0], LAMBDA[1], pixels), np.ones((n, pixels)), res.time[1])
        for name in METHODS:
            measure, apply = PAIRS[name]
            backends = ["rust", "numpy"] if name == "TMPEST" and tp.rust_available() else ["numpy"]
            for backend in backends:
                if name == "TMPEST":
                    tp.set_backend(backend)
                skip_measure = name == "dense GP" and n > GP_MEASURE_MAX
                skip_apply = name == "dense GP" and n > GP_APPLY_MAX
                model = PAIRS[name][0](res) if not skip_measure else [(SIGMA[b] ** 2, 0.5, 10.0 * CADENCE) for b in res.names]
                t_measure = np.nan if skip_measure else clock(lambda: measure(res))
                t_apply = np.nan if skip_apply else clock(lambda: apply(model, res, target, rng))
                rows.append({"epochs": n, "pixels": pixels, "method": name, "backend": backend,
                             "measure_ms": 1e3 * t_measure, "apply_ms": 1e3 * t_apply})
                print(f"{n:>7d}{pixels:>7d}  {name:<16}{backend:<7}{1e3 * t_measure:>12.3f}{1e3 * t_apply:>12.3f}")
        if tp.rust_available():
            tp.set_backend("rust")
    return rows


def save(rows):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"latency_sweep{SUFFIX}.csv")
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[SAVE] {path}")


def figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))
    for ax, column, title in zip(axes, ("measure_ms", "apply_ms"), ("Measure the noise", "Apply it to synthetic spectra")):
        ends = []
        for name, color in zip(METHODS, COLORS):
            for backend, style in (("rust", "-"), ("numpy", "--" if name == "TMPEST" else "-")):
                sel = [r for r in rows if r["method"] == name and r["backend"] == backend and np.isfinite(r[column])]
                if not sel:
                    continue
                x = [r[AXIS] for r in sel]
                y = [r[column] for r in sel]
                ax.plot(x, y, color=color, linestyle=style, linewidth=2, marker="o", markersize=5)
                ends.append([y[-1], x[-1], f"{name} ({backend})" if name == "TMPEST" else name])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ends.sort()
        low, high = ax.get_ylim()
        gap = (high / low) ** LABEL_GAP
        y_text = 0.0
        for y_end, x_end, label in ends:
            y_text = max(y_end, y_text * gap)
            ax.text(x_end * 1.12, y_text, label, va="center", ha="left", fontsize=8)
        ax.set_xlim(right=ax.get_xlim()[1] * 9.0)
        ax.set_xlabel("pixels" if PIXELS else "epochs per band")
        ax.set_ylabel("median time per call [ms]")
        ax.set_title(title, loc="left")
        ax.grid(color="0.9", linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    scope = f"pixels, three bands, {PIXEL_EPOCHS} epochs" if PIXELS else f"epochs, three bands, {N_PIXELS} pixels"
    fig.suptitle(f"Scaling with the number of {scope} (simulated {'red' if RED else 'white'} residuals)",
                 x=0.01, y=0.99, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = os.path.join(OUT, f"latency_sweep{SUFFIX}.png")
    fig.savefig(path, dpi=160)
    print(f"[SAVE] {path}")


if __name__ == "__main__":
    results = run()
    save(results)
    if "--no-figures" not in sys.argv:
        figure(results)
