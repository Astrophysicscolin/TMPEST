import csv
import json
import os
import urllib.request

import numpy as np

URL = "https://raw.githubusercontent.com/astrocatalogs/sne-2010-2014/master/SN2011fe.json"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
PHOT_SOURCES = {"2013NewA...20...30M": "Munari2013", "2012JAVSO..40..872R": "Richmond2012"}
PHOT_TELESCOPE = {"Munari2013": "ANS 61"}
BANDS = {"B": "johnson_B", "V": "johnson_V", "Rc": "cousins_R", "R": "cousins_R",
         "Ic": "cousins_I", "I": "cousins_I"}
SPEC_INSTRUMENT = "Snifs"
GRID = np.arange(335.0, 965.0, 1.0)  # [nm]

def main():
    print(f"[OK] downloading {URL}")
    with urllib.request.urlopen(URL, timeout=300) as fh:
        event = json.load(fh)["SN2011fe"]
    alias = {s["alias"]: s.get("bibcode") for s in event["sources"]}

    # PHOTOMETRY
    rows = []
    for p in event["photometry"]:
        if "magnitude" not in p or p.get("upperlimit") or p.get("band") not in BANDS:
            continue
        refs = [PHOT_SOURCES[alias[a]] for a in p["source"].split(",") if alias.get(a) in PHOT_SOURCES]
        if not refs:
            continue
        if refs[0] in PHOT_TELESCOPE and p.get("telescope") != PHOT_TELESCOPE[refs[0]]:
            continue
        rows.append((float(p["time"]), BANDS[p["band"]], float(p["magnitude"]),
                     p.get("e_magnitude", ""), refs[0]))
    rows.sort()
    path = os.path.join(OUT, "sn2011fe_photometry.csv")
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["time_mjd", "band", "mag", "mag_err", "reference"])
        writer.writerows(rows)
    print(f"[SAVE] {path}: {len(rows)} points")

    # SPECTRA
    epochs, fluxes = [], []
    for s in event["spectra"]:
        if s.get("instrument") != SPEC_INSTRUMENT:
            continue
        data = np.array([[float(v) for v in row[:2]] for row in s["data"]])
        lam = data[:, 0] / 10.0  # [nm]
        edges = np.concatenate([GRID - 0.5, [GRID[-1] + 0.5]])
        idx = np.digitize(lam, edges) - 1
        ok = (idx >= 0) & (idx < len(GRID))
        total = np.bincount(idx[ok], data[ok, 1], minlength=len(GRID))
        count = np.bincount(idx[ok], minlength=len(GRID))
        epochs.append(float(s["time"]))
        fluxes.append(np.where(count > 0, total / np.maximum(count, 1), np.nan))
    order = np.argsort(epochs)
    path = os.path.join(OUT, "sn2011fe_spectra.csv")
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["wavelength_nm"] + [f"mjd_{epochs[i]:.3f}" for i in order])
        for k, lam in enumerate(GRID):
            writer.writerow([f"{lam:.1f}"] + [f"{fluxes[i][k]:.6e}" for i in order])
    print(f"[SAVE] {path}: {len(order)} epochs x {len(GRID)} wavelengths")

if __name__ == "__main__":
    main()