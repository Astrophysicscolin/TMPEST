import contextlib
import io
import sys

import numpy as np

from roundtrip_via import N_EPOCHS, N_CHUNKS, back_catalog, tp, expected_var, fmt, read_gaia_epoch_csv, synthetic_spectra


def main(paths):
    for path in paths:
        phot = read_gaia_epoch_csv(path)
        t_all = phot.time["gaia_G"]
        print(path)
        for law in ("powerlaw", "interp"):
            with contextlib.redirect_stdout(io.StringIO()):
                p1 = tp.catalog(phot, law=law)
            total = p1.red_var + p1.white_var
            rows = []
            for seed in range(300):
                rng = np.random.default_rng(seed)
                epochs = np.sort(rng.choice(t_all, size=min(N_EPOCHS, len(t_all)), replace=False))
                clean = synthetic_spectra(epochs)
                back = back_catalog(tp.apply(p1, clean, rng), clean, law)
                n_red = int(np.sum(back.red_var > 0))
                rows.append((n_red, back.tau, *expected_var(back, phot)))
            rows = np.array(rows, dtype=float)
            print(f"  LAW {law}")
            for lo, hi in ((0, 0), (1, 2), (3, 5), (6, 8)):
                sel = (rows[:, 0] >= lo) & (rows[:, 0] <= hi)
                if sel.any():
                    print(f"    chunks called red {lo}-{hi}: {sel.mean():.2f} of seeds, rms ratio {fmt(np.sqrt(rows[sel, 2:].mean(axis=0) / total))}")
            fin = np.isfinite(rows[:, 1])
            for lo, hi in ((0, 0.5), (0.5, 2), (2, 5), (5, 1e9)):
                sel = fin & (rows[:, 1] / p1.tau >= lo) & (rows[:, 1] / p1.tau < hi)
                if sel.any():
                    print(f"    back tau / true in [{lo}, {hi}): {sel.mean():.2f} of seeds, rms ratio {fmt(np.sqrt(rows[sel, 2:].mean(axis=0) / total))}")


if __name__ == "__main__":
    main(sys.argv[1:])
