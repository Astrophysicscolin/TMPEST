import contextlib
import io
import sys

import numpy as np

from datasets import read_gaia_epoch_csv
from demo import N_EPOCHS, synthetic_spectra

import tmpest as tp

N_SEEDS = 300


def run(path, law, seed):
    rng = np.random.default_rng(seed)
    phot = read_gaia_epoch_csv(path)
    profile = tp.catalog(phot, law=law)
    total = np.sqrt(profile.red_var + profile.white_var)
    t_all = phot.time["gaia_G"]
    epochs = np.sort(rng.choice(t_all, size=min(N_EPOCHS, len(t_all)), replace=False))
    clean = synthetic_spectra(epochs)
    noisy = tp.apply(profile, clean, rng)
    back = tp.catalog(noisy, signal=clean.flux, n_chunks=8, min_points=8, law=law)
    flat = tp.Photometry({b: phot.time[b] for b in phot.names},
                         {b: np.ones(len(phot.time[b])) for b in phot.names})
    synth = tp.apply(back, flat, rng)
    direct = tp.apply(profile, flat, np.random.default_rng(seed + 100000))
    full = np.array([np.std(synth.flux[b]) for b in phot.names]) / total
    one = np.array([np.std(direct.flux[b]) for b in phot.names]) / total
    return full, one, phot.names, back.tau, profile.tau


def q(x):
    lo, med, hi = np.nanpercentile(x, [10, 50, 90])
    return f"{med:.2f} ({lo:.2f}-{hi:.2f})"


def main(paths):
    for path in paths:
        print(path)
        for law in ("powerlaw", "interp"):
            full, one, taus = [], [], []
            with contextlib.redirect_stdout(io.StringIO()):
                for seed in range(N_SEEDS):
                    f, o, names, tau_back, tau = run(path, law, seed)
                    full.append(f)
                    one.append(o)
                    taus.append(tau_back)
            full, one, taus = np.array(full), np.array(one), np.array(taus, dtype=float)
            print(f"  law={law}  profile tau={tau:.2f}  back tau median={np.nanmedian(taus):.2f}  back tau nan frac={np.mean(~np.isfinite(taus)):.2f}")
            for i, b in enumerate(names):
                print(f"    {b:8s} direct {q(one[:, i])}   via spectra {q(full[:, i])}   frac within 30%: direct {np.mean(np.abs(one[:, i]-1)<0.3):.2f} via {np.mean(np.abs(full[:, i]-1)<0.3):.2f}")


if __name__ == "__main__":
    main(sys.argv[1:])
