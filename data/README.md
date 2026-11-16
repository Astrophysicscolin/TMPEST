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
