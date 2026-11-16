import os
import sys

import matplotlib
import numpy as np
from datasets import load

import tmpest as tp
import tmpest.plot as tpp

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 42
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")


def save(fig, name):
    fig.tight_layout()
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"[SAVE] {path}")


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    # GAIA THREE BANDS TO SYNTHETIC SPECTRA
    phot = load("gaia_dr3_4368996662288843264.csv")
    profile = tp.catalog(phot, law="interp")
    print(profile.summary())
    lam = np.linspace(380.0, 1000.0, 300)  # [nm]
    epochs = np.sort(rng.choice(phot.time["gaia_G"], size=30, replace=False))
    clean = tp.Spectra(lam, np.tile((lam / 600.0) ** -1.5, (30, 1)), epochs)
    noisy = tp.apply(profile, clean, rng)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    tpp.plot_photometry(phot, ax=axes[0, 0])
    tpp.plot_residuals(phot, ax=axes[0, 1])
    tpp.plot_variogram(phot, "gaia_G", ax=axes[0, 2])
    tpp.plot_profile(profile, ax=axes[1, 0])
    tpp.plot_spectra(noisy, ax=axes[1, 1], max_lines=5)
    tpp.plot_noise(clean, noisy, ax=axes[1, 2])
    save(fig, "gaia_to_spectra.png")

    # DENSE SINGLE BAND
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for row, name in enumerate(("kepler_kic10002567.csv", "tess_141826334.csv")):
        data = load(name)
        tpp.plot_photometry(data, ax=axes[row, 0])
        tpp.plot_variogram(data, ax=axes[row, 1])
    save(fig, "kepler_tess_time_structure.png")

    # SUPERNOVA 2011FE PHOTOMETRY AND SPECTRA
    sn_phot = load("sn2011fe_photometry.csv")
    sn_spec = load("sn2011fe_spectra.csv")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    tpp.plot_photometry(sn_phot, ax=axes[0], normalize=False)
    axes[0].set_yscale("log")
    tpp.plot_spectra(sn_spec, ax=axes[1], max_lines=6)
    tpp.plot_spectra(sn_spec, ax=axes[2], kind="image", normalize=True)
    save(fig, "sn2011fe_data.png")


if __name__ == "__main__":
    sys.exit(main())
