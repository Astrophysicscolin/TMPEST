import csv
import os

import numpy as np

import tmpest as tp

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


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
