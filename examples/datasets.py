import csv
import gzip
import os

import numpy as np

import tmpest as tp

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
AAVSO_BANDS = {"B": "johnson_B", "V": "johnson_V", "R": "cousins_R"}
FITS_TYPES = {"E": ">f4", "D": ">f8", "I": ">i2", "J": ">i4", "K": ">i8", "B": "u1"}


def _rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def read_gaia_epoch_csv(path):
    cols = {"gaia_G": ("g_transit_time", "g_transit_flux", "variability_flag_g_reject"),
            "gaia_BP": ("bp_obs_time", "bp_flux", "variability_flag_bp_reject"),
            "gaia_RP": ("rp_obs_time", "rp_flux", "variability_flag_rp_reject")}
    rows = _rows(path)
    time, flux = {}, {}
    for band, (tc, fc, rc) in cols.items():
        keep = [r for r in rows if r[tc] and r[fc] and r.get(rc, "false").lower() != "true"]
        time[band] = np.array([float(r[tc]) for r in keep])
        flux[band] = np.array([float(r[fc]) for r in keep])
    return tp.Photometry(time, flux)


def read_single_band_csv(path, band, time_col="time", flux_col="flux"):
    rows = [r for r in _rows(path) if r[time_col] and r[flux_col]]
    return tp.Photometry({band: np.array([float(r[time_col]) for r in rows])},
                         {band: np.array([float(r[flux_col]) for r in rows])})


def read_ztf_csv(path):
    rows = [r for r in _rows(path) if r["catflags"] == "0"]
    time, flux = {}, {}
    for code, band in (("zg", "ztf_g"), ("zr", "ztf_r"), ("zi", "ztf_i")):
        sel = [r for r in rows if r["filtercode"] == code]
        if sel:
            time[band] = np.array([float(r["mjd"]) for r in sel])
            flux[band] = tp.mag_to_flux([float(r["mag"]) for r in sel])
    return tp.Photometry(time, flux)


def read_long_photometry_csv(path, reference=None):
    rows = [r for r in _rows(path) if reference is None or r["reference"] == reference]
    time, flux = {}, {}
    for band in sorted({r["band"] for r in rows}, key=lambda b: tp.BANDS[b][0]):
        sel = [r for r in rows if r["band"] == band]
        time[band] = np.array([float(r["time_mjd"]) for r in sel])
        flux[band] = tp.mag_to_flux([float(r["mag"]) for r in sel])
    return tp.Photometry(time, flux)


def read_wide_spectra_csv(path):
    with open(path) as fh:
        header = fh.readline().strip().split(",")
    table = np.genfromtxt(path, delimiter=",", skip_header=1)
    time = np.array([float(h.split("_")[1]) for h in header[1:]])
    return tp.Spectra(table[:, 0], table[:, 1:].T, time)


def read_aavso_flux_csv(paths):
    time, flux = {}, {}
    for path in paths:
        with open(path) as fh:
            rows = [line.split(",") for line in fh.read().splitlines() if line and not line.startswith("#")]
        band = AAVSO_BANDS[rows[0][4]]
        t = np.array([float(r[1]) for r in rows])
        f = np.array([float(r[2]) for r in rows])
        order = np.argsort(t, kind="stable")
        time[band] = t[order]
        flux[band] = f[order] / f.max()
    names = sorted(time, key=lambda b: tp.BANDS[b][0])
    return tp.Photometry({b: time[b] for b in names}, {b: flux[b] for b in names})


def read_hoffman(star):
    folder = os.path.join(DATA, "hoffman")
    paths = sorted(os.path.join(folder, f) for f in os.listdir(folder) if f.startswith(star) and f.endswith(".csv"))
    if not paths:
        raise FileNotFoundError(f"no files for {star} in {folder}; run examples/fetch_hoffman.py first")
    return read_aavso_flux_csv(paths)


def read_fits_tables(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as fh:
        raw = fh.read()
    tables, pos = [], 0
    while pos < len(raw):
        header = {}
        done = False
        while not done:
            block = raw[pos:pos + 2880].decode("ascii", "replace")
            pos += 2880
            for i in range(0, 2880, 80):
                card = block[i:i + 80]
                if card.startswith("END"):
                    done = True
                    break
                if card[8:10] == "= ":
                    header[card[:8].strip()] = card[10:].split("/")[0].strip().strip("'").strip()
        size = 0
        if int(header.get("NAXIS", 0)) > 0:
            size = abs(int(header["BITPIX"])) // 8
            for k in range(1, int(header["NAXIS"]) + 1):
                size *= int(header[f"NAXIS{k}"])
        columns = {}
        if header.get("XTENSION") == "BINTABLE":
            fields = []
            for k in range(1, int(header["TFIELDS"]) + 1):
                form = header[f"TFORM{k}"]
                count = int(form[:-1]) if form[:-1] else 1
                fields.append((header[f"TTYPE{k}"], FITS_TYPES[form[-1]], (count,)))
            rows = np.frombuffer(raw, dtype=np.dtype(fields), count=int(header["NAXIS2"]), offset=pos)
            columns = {name: rows[name] for name, _, _ in fields}
        tables.append((header, columns))
        pos += (size + int(header.get("PCOUNT", 0)) + 2879) // 2880 * 2880
    return tables


def read_lamost_mrs(folder, min_snr=10.0, step=0.01, fiber=None, normalise=True):
    arms = {"B": {}, "R": {}}
    for name in sorted(os.listdir(folder)):
        if not name.endswith((".fits", ".fits.gz")) or (fiber is not None and name[6:9] != fiber):
            continue
        for header, columns in read_fits_tables(os.path.join(folder, name))[1:]:
            arm = header.get("EXTNAME", "")[:2]
            if arm not in ("B-", "R-") or float(header["SNR"]) < min_snr:
                continue
            lam = 0.1 * columns["WAVELENGTH"][0].astype(float)
            flux = columns["FLUX"][0].astype(float)
            good = (columns["IVAR"][0] > 0) & (columns["PIXMASK"][0] == 0) & np.isfinite(flux)
            if good.sum() < 0.5 * len(flux):
                continue
            arms[arm[0]][int(header["LMJM"])] = (lam[good], flux[good] / np.median(flux[good]) if normalise else flux[good])
    epochs = sorted(set(arms["B"]) & set(arms["R"]))
    grids, cubes = [], []
    for arm in ("B", "R"):
        lo = max(arms[arm][e][0][0] for e in epochs)
        hi = min(arms[arm][e][0][-1] for e in epochs)
        grid = np.arange(lo, hi, step)
        grids.append(grid)
        cubes.append(np.stack([np.interp(grid, *arms[arm][e]) for e in epochs]))
    return tp.Spectra(np.concatenate(grids), np.concatenate(cubes, axis=1), np.array(epochs) / 1440.0)


def load(name):
    path = os.path.join(DATA, name)
    if name.startswith("gaia"):
        return read_gaia_epoch_csv(path)
    if name.startswith("kepler"):
        return read_single_band_csv(path, "kepler")
    if name.startswith("tess"):
        return read_single_band_csv(path, "tess")
    if name.startswith("ztf"):
        return read_ztf_csv(path)
    if name.endswith("photometry.csv"):
        return read_long_photometry_csv(path, reference="Munari2013")
    if name.endswith("spectra.csv"):
        return read_wide_spectra_csv(path)
    raise ValueError(f"no reader for {name}")
