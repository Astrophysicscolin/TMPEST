import os

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

try:
    from . import _core
except ImportError:
    _core = None

BACKENDS = ("rust", "numpy")
STATE = {"rust": _core is not None}
MEDIAN_BLOCK = 65536


def backend():
    return "rust" if STATE["rust"] else "numpy"


def rust_available():
    return _core is not None


def set_backend(name):
    if name not in BACKENDS:
        raise ValueError(f"backend must be one of {BACKENDS}")
    if name == "rust" and _core is None:
        raise ImportError("the compiled tmpest._core extension is not installed; build it with "
                          "'pip install -e .' (needs a Rust toolchain) or use the numpy backend")
    STATE["rust"] = name == "rust"


def num_threads():
    return _core.num_threads() if STATE["rust"] else 1


if os.environ.get("TMPEST_BACKEND"):
    set_backend(os.environ["TMPEST_BACKEND"].lower())


# OU FILTER
def ou_coefficients(x, scale):
    a = np.exp(-np.diff(x) / scale) if np.isfinite(scale) else np.ones(len(x) - 1)
    return a, np.sqrt(np.clip(1.0 - a ** 2, 0.0, None))


def ou_filter_inplace(data, x, scale, axis=0):
    a, b = ou_coefficients(np.asarray(x, dtype=float), scale)
    if STATE["rust"] and data.ndim == 2 and data.flags.c_contiguous and data.dtype == np.float64:
        _core.ou_filter_inplace(data, a, b, axis)
        return data
    view = np.moveaxis(data, axis, 0)
    for i in range(1, view.shape[0]):
        view[i] = a[i - 1] * view[i - 1] + b[i - 1] * view[i]
    return data


# VARIOGRAM
def variogram_bins(t, r, lags, edges):
    if STATE["rust"]:
        return _core.variogram_bins(t, r, [int(k) for k in lags], edges)
    n_bins = len(edges) - 1
    sum_dt = np.zeros(n_bins)
    sum_v = np.zeros(n_bins)
    count = np.zeros(n_bins)
    for k in lags:
        dt = t[k:] - t[:-k]
        keep = dt > 0
        dt = dt[keep]
        v = 0.5 * (r[k:] - r[:-k])[keep] ** 2
        idx = np.clip(np.searchsorted(edges, dt, side="right") - 1, 0, n_bins - 1)
        sum_dt += np.bincount(idx, dt, minlength=n_bins)
        sum_v += np.bincount(idx, v, minlength=n_bins)
        count += np.bincount(idx, minlength=n_bins)
    return sum_dt, sum_v, count


# RUNNING MEDIAN
def running_median(y, k):
    if STATE["rust"]:
        return _core.running_median(y, k)
    windows = sliding_window_view(np.pad(y, k // 2, mode="reflect"), k)
    out = np.empty(len(windows))
    for start in range(0, len(windows), MEDIAN_BLOCK):
        out[start:start + MEDIAN_BLOCK] = np.median(windows[start:start + MEDIAN_BLOCK], axis=-1)
    return out
