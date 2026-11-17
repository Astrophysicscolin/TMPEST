import os
import urllib.request

BASE = "https://gokhale.sites.truman.edu/files/2022/05/"
FILES = ["LX-Leo-B-filter-Measurement.csv", "LX-Leo-V-filter-Measurement.csv", "LX-Leo-R-filter-Measurements.csv",
         "MU-Leo-B-filter-Measurements.csv", "MU-Leo-V-filter-measurements.csv", "MU-Leo-R-filter-Measurements.csv",
         "V0354-UMa-B-filter-Measurements.csv", "V0354-UMa-V-filter-Measurements.csv",
         "V0354-UMa-R-filter-Measurements.csv"]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "hoffman")

def main():
    os.makedirs(OUT, exist_ok=True)
    for name in FILES:
        with urllib.request.urlopen(BASE + name, timeout=120) as fh:
            raw = fh.read()
        path = os.path.join(OUT, name)
        with open(path, "wb") as fh:
            fh.write(raw)
        print(f"[SAVE] {path} ({len(raw)} bytes)")

if __name__ == "__main__":
    main()