# Test data

Small real datasets for trying TMPEST by hand. No test in `tests/` depends on them. Readers
are in `examples/datasets.py` (`load("<file name>")` returns a `Photometry` or `Spectra`).

| File | What | Shape | Origin |
|---|---|---|---|
| `sn2011fe_photometry.csv` | SN 2011fe, B V Rc Ic magnitudes | 630 rows, long format | Open Supernova Catalog |
| `sn2011fe_spectra.csv` | SN 2011fe, flux-calibrated spectra | 32 epochs x 630 wavelengths (335-964 nm, 1 nm bins) | Open Supernova Catalog |
| `gaia_dr3_1342277958330527360.csv` | Gaia DR3 epoch photometry, G BP RP flux | 49 transits | Gaia archive |
| `gaia_dr3_4368996662288843264.csv` | Gaia DR3 epoch photometry, G BP RP flux | 42 transits | Gaia archive |
| `kepler_kic10002567.csv` | Kepler long-cadence light curve, normalised flux | 1,626 rows | MAST |
| `tess_141826334.csv` | TESS light curve, normalised flux | 1,062 rows | MAST |
| `ztf_1721114100000000.csv` | ZTF g-band magnitudes | 19 rows | IRSA |

## Provenance

The same facts, per file and machine-readable, with a size and SHA-256 for each file, are in
`PROVENANCE.yaml`.

- **SN 2011fe** was downloaded on 2026-11-15 by `examples/fetch_data.py` from the Open
  Supernova Catalog on GitHub (`astrocatalogs/sne-2010-2014`, Guillochon et al. 2017,
  2017ApJ...835...64G). Photometry: Munari et al. 2013 (2013NewA...20...30M, telescope
  "ANS 61" only) and Richmond & Smith 2012 (2012JAVSO..40..872R); the `reference` column says
  which. Spectra: the 32 SNIFS spectra of Pereira et al. 2013 (2013A&A...554A..27P, Nearby
  Supernova Factory), rebinned here from about 0.24 nm to 1 nm by averaging. Re-run the script
  to rebuild both files.
- **Gaia, Kepler, TESS, ZTF** files were not downloaded for this repo. They are copies, with
  unused columns dropped, of files the author pulled from the public archives for earlier work (Gaia DR3 epoch photometry; Kepler and TESS light curves
  through lightkurve; ZTF through IRSA). The number in each file name is the source id used
  there (Gaia DR3 `source_id`, Kepler `kepid`, the TESS and ZTF object ids).
- Cite the original papers and archives, not this folder, if any of this is used in a result.
  The references, with DOIs, are in the Citing section of the top-level `README.md` and in
  `CITATION.cff`.

## Hoffman et al. 2022 light curves

`data/hoffman/` holds nine files: B, V and R differential photometry of LX Leo, MU Leo and
V354 UMa, one night each (81 to 188 points per band), from the project page of Hoffman,
Mazzola & Gokhale 2022 (JAAVSO 50, 3). `examples/fetch_hoffman.py` downloads them again. They
are in the AAVSO extended format and unchanged. The column headed `MAG` appears to hold
relative flux, not magnitudes (its values run from 0.02 to 2.0 and their ratios match the
eclipse depths in the paper), and `examples/datasets.py` reads it that way; this has not been
confirmed with the authors. The Julian dates fall in March 2019. These are eclipsing
binaries, so the light curve is signal: `examples/latency.py` divides out a Fourier fit
before measuring noise.

## LAMOST spectra (not in the repository)

`examples/fetch_lamost.py` downloads ten observation files (about 4 MB) into `data/lamost/`,
which is git-ignored: LAMOST DR9 v2.0 medium-resolution spectra (R about 7500) of the
*Kepler* eclipsing binary KIC 5444392 (period 1.5 d), 41 single 20-minute exposures per arm
on ten nights between May 2019 and June 2020. The blue arm covers 490 to 540 nm and the red
arm 628 to 685 nm, with nothing in between. `examples/datasets.py` reads the single exposures
with its own small FITS reader, drops those with a signal-to-noise ratio under 10 and those
with more than half their pixels flagged, and puts them on a common 0.01 nm grid. By default
it divides each arm of each exposure by its median; `normalise=False` keeps the raw flux. The
observation ids are listed in the fetch script; they came from a 5 arcsecond position search
of the medium-resolution catalogue.

The same script downloads comparison stars (about 67 MB, into `data/lamost/comparison/`):
the stars in fibres 92 to 112 of the same spectrograph, which held the same star on all ten
nights (an observation id is the night's prefix, two digits of spectrograph and three of
fibre; fibres 102, 103, 105, 109 and 112 are absent or the target). Ten of them pass the
signal-to-noise cut in both arms at least once. The spectra are not flux-calibrated, so
`examples/latency_spectra.py` works like differential photometry: it divides the target's
raw flux by the median over comparison stars of each one's flux relative to its own median
spectrum, then by a Fourier fit in orbital phase at the period of Kirk et al. 2016
(1.5195284 d), with the number of harmonics chosen by leave-one-night-out prediction error
(one). What remains is dominated by the throughput of the target's own fibre, which the
comparison stars do not share: 23 per cent scatter, nearly the same at every wavelength.

## What each is good for

- **Gaia**: three bands at matched times, about 40 points each. Sparse; shows the
  different-instrument problem described under Limits in the top-level `README.md`.
- **Kepler, TESS**: one dense band each. Good for the time-structure fit (`plot_variogram`).
- **ZTF**: magnitudes and very few points. Exercises `mag_to_flux` and the sparse limits.
- **SN 2011fe**: the only source here with both photometry and a spectral time series. It is
  a transient, not a steady source: the supernova's rise and decline is signal, not jitter.
  Cataloguing it as-is reports the light curve itself as "red noise" (sigma near 1). Use it to
  test readers, shapes and plots, or with a `signal=` model. It is not a test of noise
  transfer.
