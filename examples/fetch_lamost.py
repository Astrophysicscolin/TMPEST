import os
import urllib.error
import urllib.request

BASE = "https://www.lamost.org/dr9/v2.0/medspectrum/fits/"
OBSIDS = [744801102, 747001102, 747401102, 824301102, 824501102, 824801102, 825601102, 826101102, 826401102,
          826701102]
COMPARISON_FIBERS = [92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 104, 106, 107, 108, 110, 111]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "lamost")

def fetch(obsid, folder):
    path = os.path.join(folder, f"{obsid}.fits.gz")
    if os.path.exists(path):
        return
    try:
        with urllib.request.urlopen(f"{BASE}{obsid}", timeout=300) as fh:
            raw = fh.read()
    except urllib.error.HTTPError as err:
        print(f"[WARN] {obsid}: HTTP {err.code}, skipped")
        return
    with open(path, "wb") as fh:
        fh.write(raw)
    print(f"[SAVE] {path} ({len(raw)} bytes)")

def main():
    comparison = os.path.join(OUT, "comparison")
    os.makedirs(comparison, exist_ok=True)
    for obsid in OBSIDS:
        fetch(obsid, OUT)
        for fiber in COMPARISON_FIBERS:
            fetch(obsid // 1000 * 1000 + fiber, comparison)

if __name__ == "__main__":
    main()