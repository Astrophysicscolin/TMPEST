import sys

import numpy as np
from datasets import read_gaia_epoch_csv

import tmpest as tp

SEED = 42
N_EPOCHS = 40
N_PIXELS = 600
LINES = [(486.1, 0.6), (656.3, 0.8), (854.2, 0.5)]  # [nm]

def simulated_photometry(rng):
    t = np.sort(rng.uniform(0.0, 900.0, 400))  # [d]
    activity = 4e-3 * tp.ou_filter(rng.standard_normal(400), t, 12.0)
    time, flux = {}, {}
    for band in ("gaia_G", "gaia_BP", "gaia_RP"):
        chroma = (tp.BANDS[band][0] / 600.0) ** -1.0
        time[band] = t
        flux[band] = 1e4 * (1.0 + chroma * activity + 1.5e-3 * rng.standard_normal(400))
    return tp.Photometry(time, flux)

def synthetic_spectra(time):
    lam = np.linspace(380.0, 1000.0, N_PIXELS)  # [nm]
    continuum = (lam / 600.0) ** -1.5
    for centre, depth in LINES:
        continuum = continuum * (1.0 - depth * np.exp(-0.5 * ((lam - centre) / 1.5) ** 2))
    return tp.Spectra(lam, np.tile(continuum, (len(time), 1)), time)

def main():
    rng = np.random.default_rng(SEED)
    phot = read_gaia_epoch_csv(sys.argv[1]) if len(sys.argv) > 1 else simulated_photometry(rng)

    # PHOTOMETRY TO SPECTROSCOPY
    profile = tp.catalog(phot)
    print(profile.summary())
    t_all = phot.time["gaia_G"]
    epochs = np.sort(rng.choice(t_all, size=min(N_EPOCHS, len(t_all)), replace=False))
    clean = synthetic_spectra(epochs)
    noisy = tp.apply(profile, clean, rng)
    scatter = np.std(noisy.flux / clean.flux - 1.0, axis=0)
    print(f"[OK] noisy spectra {noisy.flux.shape}: per-pixel scatter "
          f"{scatter[0]:.2e} at {clean.wavelength[0]:.0f} nm, {scatter[-1]:.2e} at {clean.wavelength[-1]:.0f} nm")

    # SPECTROSCOPY TO PHOTOMETRY
    back = tp.catalog(noisy, signal=clean.flux, n_chunks=8, min_points=8)
    flat = tp.Photometry({b: phot.time[b] for b in phot.names},
                         {b: np.ones(len(phot.time[b])) for b in phot.names})
    synth = tp.apply(back, flat, rng)
    for i, band in enumerate(phot.names):
        print(f"[OK] {band}: synthetic photometry scatter {np.std(synth.flux[band]):.2e}, "
              f"original total {np.sqrt(profile.red_var[i] + profile.white_var[i]):.2e}")

    if len(sys.argv) > 2:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(8, 7))
        for band in phot.names:
            axes[0].plot(phot.time[band], phot.flux[band] / np.median(phot.flux[band]), ".", label=band)
        axes[0].set_xlabel("time")
        axes[0].set_ylabel("relative flux")
        axes[0].legend()
        for k in range(min(3, len(epochs))):
            axes[1].plot(noisy.wavelength, noisy.flux[k], lw=0.6)
        axes[1].plot(clean.wavelength, clean.flux[0], "k", lw=0.8)
        axes[1].set_xlabel("wavelength [nm]")
        axes[1].set_ylabel("flux")
        fig.tight_layout()
        fig.savefig(sys.argv[2], dpi=150)
        print(f"[SAVE] {sys.argv[2]}")

if __name__ == "__main__":
    main()