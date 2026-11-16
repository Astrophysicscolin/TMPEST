# TMPEST

TMPEST: Transfer of Multimodal Perturbations upon Emulated Spectra and Time series.

Measure the chromatic jitter (wavelength-dependent noise) in one kind of data and put it on
synthetic data of the other kind, for the same source:

- real **photometry** (a few broad bands, many epochs) -> synthetic **spectroscopy**
- real **spectroscopy** (many wavelengths, some epochs) -> synthetic **photometry**

The only Python dependency is numpy; matplotlib is optional, for
`tmpest.plot`.

## Install

TMPEST is a Python library with a small compiled core written in Rust (the slow inner loops
only). Installing from source needs Python 3.10+, and the Rust compiler from
https://rustup.rs. Step-by-step guides for each operating system are in `docs/testing/`.
It needs numpy 1.21.6 or newer. The tests pass on Python 3.10 to 3.14 with numpy 1.21.6 to
2.5.5 (16 combinations, Apple Silicon Mac, 2026-11-15); other platforms are listed in
`docs/testing/README.md`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,examples]"     # compiles the Rust core, 1-3 minutes the first time
python examples/system_report.py               # should say: active backend: rust
python -m pytest -q                            # 43 tests, a few seconds
python examples/demo.py                        # simulated source, both directions
python examples/gallery.py                     # figures/ from the files in data/
python examples/benchmark.py                   # numpy vs. Rust timings on this machine
```

`-e` is an editable install: Python imports the code straight from `src/`, so edits to `.py`
files take effect at once. Edits to `rust/src/lib.rs` need the install line run again.
`[dev]` adds pytest; `[plot]` or `[examples]` adds matplotlib.

## Backends

Every compiled routine has a plain-numpy twin that gives bit-for-bit the same numbers, and the
tests check that. The Rust one is used when it is built.

```python
import tmpest as tp
tp.backend()            # "rust" or "numpy"
tp.set_backend("numpy") # switch at run time; "rust" raises ImportError if it is not built
tp.num_threads()        # worker threads the Rust core uses
```

Setting the environment variable `TMPEST_BACKEND=numpy` before starting Python does the
same. `RAYON_NUM_THREADS=4` limits the Rust core to 4 threads (default: all cores). Without
the compiled core (for example running straight from `src/` with `PYTHONPATH`), everything
still works on the numpy backend, more slowly on long series.

What the Rust core covers, measured on a 2-core Linux machine (`examples/benchmark.py`):

| Job | numpy | Rust | Speed-up |
|---|---|---|---|
| Correlated-noise filter, 1,000,000 epochs x 9 | 1.74 s | 0.03 s | 60x |
| Running-median detrend, 1,000,000 points, window 501 | 6.37 s | 0.06 s | 104x |
| Variogram, 1,000,000 points | 0.95 s | 0.13 s | 7x |
| `catalog`, 1,000,000 points with `window=` | 9.2 s | 0.14 s | 66x |
| `apply` to photometry, 1,000,000 points | 2.3 s | 0.35 s | 7x |
| `apply` to spectra, 2,000 x 8,000 | 1.2 s | 0.9 s | 1.3x |

The same jobs on a 14-core Apple Silicon Mac (macOS 15, Python 3.12, 14 threads, one run on
2026-11-15):

| Job | numpy | Rust | Speed-up |
|---|---|---|---|
| Correlated-noise filter, 1,000,000 epochs x 9 | 1.23 s | 0.02 s | 66x |
| Running-median detrend, 1,000,000 points, window 501 | 4.58 s | 0.02 s | 305x |
| Variogram, 1,000,000 points | 0.71 s | 0.02 s | 42x |
| `catalog`, 1,000,000 points with `window=` | 6.6 s | 0.06 s | 113x |
| `apply` to photometry, 1,000,000 points | 1.4 s | 0.15 s | 10x |
| `apply` to spectra, 2,000 x 8,000 | 0.49 s | 0.34 s | 1.4x |

In both tables the last row barely moves because drawing the random numbers (numpy, one thread) is now the
slow part. Ten million points catalogued and applied in under a minute with 1.4 GB of memory.
Everything is still held in memory at once: a spectral cube needs about 3 x 8 bytes per
element, so that, not speed, is the ceiling for very large jobs.

## Quick start

```python
import numpy as np
import tmpest as tp

rng = np.random.default_rng(42)          # every random draw comes from this generator

# real photometry: one time array and one flux array per band (dicts keyed by band name)
phot = tp.Photometry(time={"gaia_G": t_g, "gaia_BP": t_bp, "gaia_RP": t_rp},
                     flux={"gaia_G": f_g, "gaia_BP": f_bp, "gaia_RP": f_rp})

# synthetic, noiseless spectra of the same source: flux has shape (n_epochs, n_wavelengths)
spec = tp.Spectra(wavelength=lam_nm, flux=model_flux, time=t_spec)

noisy_spec, profile = tp.phot_to_spec(phot, spec, rng)     # photometry -> spectroscopy
noisy_phot, profile = tp.spec_to_phot(real_spec, synthetic_phot, rng)   # and the reverse
print(profile.summary())
```

The one-call functions are the three steps below chained together. Use the steps when you want
to look at, save or reuse what was measured.

| Step | Call | Gives you |
|---|---|---|
| Scrape | `res = tp.scrape(data)` | `Residuals`: the noise itself, per channel, as relative flux with the baseline divided out |
| Catalog | `profile = tp.catalog(data)` | `ChromaticProfile`: the measured numbers (calls `scrape` for you) |
| Apply | `noisy = tp.apply(profile, target, rng)` | a new `Photometry` or `Spectra` with noise drawn from the profile |

`tp.draw(profile, target, rng)` returns the noise alone (relative units) without touching the
flux. `profile.save("p.json")` and `tp.ChromaticProfile.load("p.json")` store a profile, so one
catalog pass can feed many synthetic datasets.

## What is measured

Each band (or each wavelength chunk of a spectrum) is a channel. Its residual is modelled as
three parts:

| Part | In time | Across wavelength | Stored as |
|---|---|---|---|
| red jitter | correlated, exp(-dt / `tau`) | coherent over `ell` in ln(wavelength) | `sigma_red(lambda)` |
| coherent white | uncorrelated | coherent over `ell` | `sigma_coherent_white(lambda)` |
| incoherent white | uncorrelated | uncorrelated | `sigma_white(lambda, width)` |

- **Time structure** is fitted per channel from the variogram (half the mean squared
  difference of pairs of points, against their time separation). It works on irregular
  sampling and is in real time units, so a profile from 30-minute data can be applied at any
  other cadence. A channel is only called red if a von Neumann test rejects white noise
  (`significance=2.5`); otherwise all its variance is white.
- **Wavelength structure** comes from the covariance of channels observed at matching times.
  It sets `ell` and the share of the white noise that is common to all channels
  (`coherent_fraction`). This is what decides whether noise survives being averaged over a
  passband, so it matters most in the spectroscopy -> photometry direction.
- **Chromatic law.** Each amplitude is a power law in wavelength,
  `sigma = amp * (lambda / lambda_ref)^alpha`, fitted across channels. `law="interp"` uses
  log-log interpolation between channels instead (flat outside them), which suits a spectrum
  with many chunks.
- **Band integration.** A passband is treated as a top-hat of the given width. Measured band
  variances are converted to per-wavelength values, and synthetic bands are built by averaging
  `n_sub` wavelengths across the band, so the two directions are consistent.
- **Incoherent white noise scales as 1/sqrt(channel width)**, the photon-noise scaling. This
  assumes the target instrument collects light like the source instrument. It usually does
  not, so set `white=` in `apply` (a multiplier; `0` turns white noise off) or add your own
  instrument noise afterwards. `red=` is the same knob for the red part.

## Plotting

`tmpest.plot` is a set of plain matplotlib helpers. It is not imported by `import tmpest`,
so matplotlib stays optional.

```python
import matplotlib.pyplot as plt
import tmpest.plot as tpp

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
tpp.plot_photometry(phot, ax=axes[0])
tpp.plot_profile(profile, ax=axes[1])
axes[1].set_title("my title")        # the return value is an ordinary matplotlib Axes
fig.savefig("check.png")
```

Every function takes `ax=` (an Axes to draw on; it makes a new figure if you leave it out),
returns that Axes, and never touches global style (`rcParams`), so your own style sheets,
subplots and later calls work as usual. Extra keyword arguments go to the underlying
matplotlib call.

| Function | Shows |
|---|---|
| `plot_photometry(phot)` | Each band against time (`normalize=False` for raw flux). |
| `plot_spectra(spec)` | A few spectra as lines, or all of them as an image with `kind="image"`. |
| `plot_noise(clean, noisy)` | What `apply` added: noisy / clean, for either data type. |
| `plot_residuals(data)` | What `scrape` extracted, one offset row per channel. |
| `plot_variogram(data, channel)` | The time-structure fit for one channel: measured variogram, fitted curve, white level. |
| `plot_profile(profile)` | Measured red and white noise per channel, and the fitted laws against wavelength. |
| `plot_correlation(profile)` | Channel-to-channel correlation at matched times. |

## Test data

`data/` holds seven small real CSV files (Gaia, Kepler, TESS, ZTF photometry; SN 2011fe
photometry and spectra); `data/README.md` gives shapes, provenance and what each is good for.
`examples/datasets.py` reads them: `load("kepler_kic10002567.csv")` returns a `Photometry`.

## Options

| Argument | Where | Meaning |
|---|---|---|
| `signal` | `scrape`, `catalog` | The noiseless model of the data (same shapes as the flux). Divided out before measuring. Use it when the data contain a transit or other real signal. |
| `window` | `scrape`, `catalog` | Running-median detrend window, in time units. Default `None` divides by the median only, so slow variability counts as red jitter. |
| `clip` | `scrape`, `catalog` | Outlier rejection in robust sigmas (default 5). |
| `n_chunks` | spectra only | Number of log-spaced wavelength chunks (default 16). |
| `alpha` | `catalog` | Force the red-jitter exponent. Needed with a single band, where it cannot be measured (default then is 0, achromatic). |
| `match_tol` | `catalog` | Largest time offset for two channels to count as simultaneous. Default: a quarter of the shorter median cadence. |
| `min_points` | `catalog` | Channels with fewer usable points are skipped with a `[WARN]` (default 10). |
| `red`, `white` | `apply`, `draw` | Amplitude multipliers. |

## Conventions

- Flux is linear and positive. Convert magnitudes with `tp.mag_to_flux`. Zero-centred (ppm)
  light curves need their mean added back first; a non-positive baseline raises `ValueError`.
- Wavelengths and widths are in nm. Time is in any unit, but the same unit and zero point in
  both datasets, because `tau` is carried in that unit.
- Built-in passbands are in `tp.BANDS` (`gaia_G`, `gaia_BP`, `gaia_RP`, `kepler`, `tess`,
  `ztf_g/r/i`, `sdss_u/g/r/i/z`, `johnson_U/B/V`, `cousins_R/I`, `2mass_J/H/Ks`) as
  (centre, width). They are approximate. For anything else, or for exact values, pass
  `wavelength={"name": nm}` and `width={"name": nm}` to `Photometry`.
- `rng` is always an explicit `np.random.Generator`. Same seed, same noise.
- `tp.scale_depth(flux, wavelength, alpha, lambda_ref)` scales a dip (not the whole flux) by
  a power law in wavelength, for making a transit or other dip chromatic.

## Outside astronomy

Nothing in the maths is astronomical. A "wavelength" is any positive coordinate that orders
the channels (frequency, energy, sensor index + 1), and "bands vs. spectra" is "a few broad
channels vs. many narrow ones". Pass `wavelength=` explicitly and leave `width=None`: the
profile then carries white noise per channel with no width scaling.

## Limits

- The noise is Gaussian and stationary. Flares, gaps in quality and heavy tails are clipped
  out, not reproduced.
- One `tau` and one `ell` per profile. A source with several timescales gets a blend.
- Sparse data cannot establish red noise: with about 60 points and equal red and white
  variance, the red part is detected roughly 4 times in 10 (simulated). Check
  `profile.summary()` before trusting a transfer.
- Noise is assumed to vary smoothly with wavelength. It does not when the channels come from
  different instruments. On real Gaia photometry, where G is several times more precise than
  BP and RP, the default power law makes G about 20% too noisy and BP about 25% too quiet
  (median of 120 sources). `law="interp"` brings BP and RP to within about 5%; G stays about
  10% high.
- Each band's total noise is not guaranteed to survive a transfer when the bands disagree on
  how much of it is red. Red, coherent white and incoherent white each get their own law
  across wavelength, fitted only to the channels where that component was detected, so a
  band can be given a component it never showed and lose one it did. On two real Gaia
  sources with about 40 epochs (`data/`, 2000 seeds each), applying a catalogued profile
  straight back to the same bands gave an rms scatter between 0.77 and 1.14 of the measured
  total under the default law and between 0.85 and 1.27 under `law="interp"`, with the sign
  depending on the source. Neither law is uniformly better.
- With few epochs a single transfer is a noisy draw. For those sources one realisation
  landed within 30% of the measured scatter in 53 to 96% of 300 seeds going directly, and in
  43 to 85% going through 40 synthetic spectral epochs and back. On the second route the
  per-chunk red-noise test passes in only some chunks, and on one of the two sources the
  result came out 15 to 40% too noisy on average. Run several seeds and compare the band scatter with
  `profile.summary()`; `examples/roundtrip_seeds.py` does this for a Gaia file.
- The total a profile reports for a red channel, `sqrt(red_var + white_var)`, comes from
  the variogram fit and is the long-baseline value. It is not the sample scatter of the
  data: on the same two sources the real G-band scatter was 1.18 and 0.89 times it. A
  synthetic series no longer than a few `tau` also shows less scatter than that total.
- Wavelength coherence is only measured from epochs that count as simultaneous. The default
  `match_tol` is a quarter of the shorter cadence, which is too strict for a camera that
  cycles filters: the bands are then offset by a third of the cadence, almost no epochs
  pair up, and the noise is treated as unshared. On BVR light curves taken that way
  (Hoffman et al. 2022) the default found no shared noise where the bands correlate at 0.2
  to 0.5; `match_tol` of half the cadence recovered it. Check `profile.corr` and set
  `match_tol` for interleaved filters.
- `catalog` has a fixed cost of about 12 ms however small the input, almost all of it the
  grid search for the wavelength coherence.
- Three overlapping broad bands constrain `alpha` and `ell` weakly. Spectra need at least
  `min_points` epochs to be catalogued; a single spectrum can be a target but not a source.
- Noise is multiplicative on the flux, so it does not add photon noise in the cores of deep
  absorption lines.
- A fresh realisation is drawn each time. Replaying the actual noise of the real epochs onto
  simultaneous synthetic data is not implemented.

## Citing

The author is Colin Harrison (Florida Institute of Technology). `CITATION.cff` holds the
citation for the software; there is no paper or DOI yet. If you use the files in `data/`, or
your own data from these surveys, cite the surveys, not this repository:

| Data | Cite |
|---|---|
| Gaia DR3 epoch photometry | Gaia Collaboration, Prusti et al. 2016, A&A 595, A1, doi:10.1051/0004-6361/201629272; Gaia Collaboration, Vallenari et al. 2023, A&A 674, A1, doi:10.1051/0004-6361/202243940; Eyer et al. 2023, A&A 674, A13, doi:10.1051/0004-6361/202244242 |
| Kepler | Borucki et al. 2010, Science 327, 977, doi:10.1126/science.1185402 |
| TESS | Ricker et al. 2015, JATIS 1, 014003, doi:10.1117/1.JATIS.1.1.014003 |
| ZTF | Bellm et al. 2019, PASP 131, 018002, doi:10.1088/1538-3873/aaecbe; Masci et al. 2019, PASP 131, 018003, doi:10.1088/1538-3873/aae8ac |
| SN 2011fe (via the Open Supernova Catalog) | Guillochon et al. 2017, ApJ 835, 64, doi:10.3847/1538-4357/835/1/64; spectra: Pereira et al. 2013, A&A 554, A27, doi:10.1051/0004-6361/201221008; photometry: Munari et al. 2013, New Astronomy 20, 30, doi:10.1016/j.newast.2012.09.003, and Richmond & Smith 2012, JAVSO 40, 872 |
| BVR photometry of LX Leo, V354 UMa and MU Leo (`data/hoffman/`; the base dataset of `examples/latency.py`, paper text in `HOFFMAN.md`) | Hoffman, Mazzola & Gokhale 2022, JAAVSO 50, 3; the data are in the AAVSO International Database |

The archives (ESA Gaia Archive, MAST, IRSA) each ask for their own acknowledgement text in
papers; take the current wording from the archive. Provenance of each file is in
`data/README.md`.

### Software and comparison methods

TMPEST is built on NumPy (Harris et al. 2020, Nature 585, 357, doi:10.1038/s41586-020-2649-2),
and its figures on Matplotlib (Hunter 2007, Computing in Science & Engineering 9, 90,
doi:10.1109/MCSE.2007.55).

`examples/latency.py` times TMPEST against four other ways of putting measured noise on
synthetic spectra. No outside package is run: each is a short NumPy version of the method,
written in that script, so its timings are not those of any published code. Cite the method:

| Baseline in `examples/latency.py` | Method | Cite |
|---|---|---|
| white Gaussian | Independent Gaussian noise per pixel, scaled to the band scatter | none; the standard simple case |
| block bootstrap | Resampling the real residuals in contiguous blocks | Efron 1979, Annals of Statistics 7, 1, doi:10.1214/aos/1176344552; Kunsch 1989, Annals of Statistics 17, 1217, doi:10.1214/aos/1176347265 |
| Timmer-Koenig | Simulating a series from a fitted power-law spectrum | Timmer & Koenig 1995, A&A 300, 707 |
| dense GP | Gaussian process with an exponential kernel, by Cholesky factorisation | Rasmussen & Williams 2006, Gaussian Processes for Machine Learning, MIT Press, doi:10.7551/mitpress/3206.001.0001 |

The usual fast Gaussian-process package in astronomy is celerite (Foreman-Mackey et al. 2017,
AJ 154, 220, doi:10.3847/1538-3881/aa9332). It is not benchmarked here and would be much
faster than the dense version above.

## Licence

BSD 3-Clause; see `LICENSE` (copyright Colin Harrison and Florida Institute of Technology).
The compiled core is built from the Rust crates PyO3, rust-numpy, ndarray, rayon and their
dependencies; their licence texts are in `THIRD-PARTY-LICENSES.txt`, which ships in every
wheel and sdist.
