import warnings

import numpy as np

from . import _kernels
from .data import Photometry, Residuals, Spectra
from .profile import ChromaticProfile

ELL_GRID = np.concatenate([[np.inf], np.geomspace(100.0, 0.01, 80)])
PHI_GRID = np.linspace(0.0, 1.0, 41)
SHORT_LAG_QUANTILE = 0.06


# SCRAPE
def _running_median(t, y, window):
    cadence = _cadence(t)
    if cadence == 0:
        return np.full(len(y), np.median(y))
    k = max(int(round(window / cadence)) | 1, 3)
    if k >= len(y):
        return np.full(len(y), np.median(y))
    return _kernels.running_median(np.ascontiguousarray(y, dtype=float), k)


def _chunk_index(lam, n):
    edges = np.geomspace(lam[0], lam[-1], n + 1)
    return np.clip(np.searchsorted(edges, lam, side="right") - 1, 0, n - 1)


def _channel_residual(t, y, window, clip, label):
    good = np.isfinite(t) & np.isfinite(y)
    t, y = t[good], y[good]
    order = np.argsort(t, kind="stable")
    t, y = t[order], y[order]
    if len(y) == 0:
        return t, y
    base = np.full(len(y), np.median(y)) if window is None else _running_median(t, y, window)
    if np.any(base <= 0):
        raise ValueError(f"channel {label}: flux baseline is not positive. TMPEST works on linear, "
                         "positive flux; convert magnitudes with mag_to_flux and add the mean back "
                         "to zero-centred (ppm) light curves first")
    r = y / base - 1.0
    if clip is not None:
        dev = np.abs(r - np.median(r))
        mad = 1.4826 * np.median(dev)
        if mad > 0:
            keep = dev <= clip * mad
            t, r = t[keep], r[keep]
    return t, r - np.mean(r)


def scrape(data, signal=None, window=None, clip=5.0, n_chunks=16):
    if isinstance(data, Residuals):
        return data
    if isinstance(data, Photometry):
        times, resids = [], []
        for name in data.names:
            y = data.flux[name]
            if signal is not None:
                y = y / np.asarray(signal[name], dtype=float)
            t, r = _channel_residual(data.time[name], y, window, clip, name)
            times.append(t)
            resids.append(r)
        return Residuals(times, resids, data.wavelength_array(), data.width_array(),
                         data.names, "photometry")
    if isinstance(data, Spectra):
        base = np.nanmedian(data.flux, axis=0) if signal is None else np.asarray(signal, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = data.flux / base
        ratio[~np.isfinite(ratio)] = np.nan
        lam = data.wavelength
        n = int(min(n_chunks, len(lam)))
        idx = _chunk_index(lam, n)
        pix = data.pixel_width()
        times, resids, wavelength, width, names = [], [], [], [], []
        for c in range(n):
            sel = idx == c
            if not sel.any():
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                y = np.nanmean(ratio[:, sel], axis=1)
            lam_c = float(np.exp(np.mean(np.log(lam[sel]))))
            t, r = _channel_residual(data.time, y, window, clip, f"{lam_c:.1f} nm")
            times.append(t)
            resids.append(r)
            wavelength.append(lam_c)
            width.append(float(pix[sel].sum()))
            names.append(f"{lam_c:.1f}nm")
        return Residuals(times, resids, np.array(wavelength), np.array(width), names, "spectra")
    raise TypeError("scrape expects a Photometry, Spectra or Residuals object")


# TIME STRUCTURE
def _variogram(t, r, n_lags=40, n_bins=24):
    n = len(r)
    empty = (np.array([]), np.array([]), np.array([]))
    if n < 2:
        return empty
    t = np.ascontiguousarray(t, dtype=float)
    r = np.ascontiguousarray(r, dtype=float)
    lags = np.unique(np.round(np.geomspace(1, max(n // 2, 1), n_lags)).astype(int))
    steps = np.diff(t)
    steps = steps[steps > 0]
    if len(steps) == 0:
        return empty
    lo = float(np.quantile(steps, SHORT_LAG_QUANTILE))
    hi = float(np.max(t[lags[-1]:] - t[:-lags[-1]]))
    if hi <= lo:
        edges = np.array([0.5 * lo, 2.0 * lo])
    else:
        edges = np.geomspace(lo, hi, n_bins + 1)
    sum_dt, sum_v, count = _kernels.variogram_bins(t, r, lags, edges)
    ok = count >= 5
    return sum_dt[ok] / count[ok], sum_v[ok] / count[ok], np.minimum(count[ok], n)


def _fit_noise(t, r, significance=2.5):
    var = float(np.var(r))
    if var == 0:
        return 0.0, 0.0, np.nan
    von_neumann = float(np.mean(np.diff(r) ** 2) / (2.0 * var))
    if von_neumann > 1.0 - significance / np.sqrt(len(r)):
        return var, 0.0, np.nan
    lag, gamma, count = _variogram(t, r)
    if len(lag) == 0 or np.any(gamma <= 0):
        return var, 0.0, np.nan
    weight = count / gamma ** 2
    cadence = _cadence(t)
    baseline = float(t[-1] - t[0])
    best = None
    if len(lag) >= 4 and baseline > cadence:
        root = np.sqrt(weight)
        for tau in np.geomspace(cadence, baseline, 80):
            g = 1.0 - np.exp(-lag / tau)
            design = np.stack([np.ones_like(g), g], axis=1) * root[:, None]
            w, s = np.linalg.lstsq(design, gamma * root, rcond=None)[0]
            if s <= 0:
                continue
            if w < 0:
                w = 0.0
                s = float(np.sum(weight * g * gamma) / np.sum(weight * g * g))
            sse = float(np.sum(weight * (gamma - w - s * g) ** 2))
            if best is None or sse < best[0]:
                best = (sse, float(w), float(s), float(tau))
    if best is None:
        return var, 0.0, np.nan
    return best[1], min(best[2], 3.0 * var), best[3]


# WAVELENGTH STRUCTURE
def sub_points(wavelength, width, n_sub):
    wavelength = np.asarray(wavelength, dtype=float)
    if width is None:
        return np.log(wavelength)[:, None]
    half = np.minimum(0.5 * np.asarray(width, dtype=float), 0.95 * wavelength)
    frac = (np.arange(n_sub) + 0.5) / n_sub - 0.5
    return np.log(wavelength[:, None] + 2.0 * half[:, None] * frac[None, :])


def _kernel(ua, ub, ell):
    return float(np.mean(np.exp(-np.abs(ua[:, None] - ub[None, :]) / ell)))


def _nearest(ta, tb):
    j = np.clip(np.searchsorted(tb, ta), 1, len(tb) - 1)
    return np.where(np.abs(tb[j - 1] - ta) <= np.abs(tb[j] - ta), j - 1, j)


def _matched_cov(ta, ra, tb, rb, tol):
    j = _nearest(ta, tb)
    ok = np.abs(tb[j] - ta) <= tol
    n = int(ok.sum())
    if n < 8:
        return np.nan, 0
    return float(np.mean(ra[ok] * rb[j[ok]])), n


def _cadence(t):
    steps = np.diff(t)
    steps = steps[steps > 0]
    return float(np.median(steps)) if len(steps) else 0.0


def _fit_coherence(u, cov, weight, red_var, white_var, sample_var):
    total = red_var + white_var
    ia, ib = np.triu_indices(len(u), 1)
    use = (weight[ia, ib] > 0) & (total[ia] > 0) & (total[ib] > 0) & (sample_var[ia] > 0) & (sample_var[ib] > 0)
    if not use.any():
        return np.inf, 0.0, False
    ia, ib = ia[use], ib[use]
    self_dist = [np.abs(u[c][:, None] - u[c][None, :]) for c in range(len(u))]
    pair_dist = [np.abs(u[a][:, None] - u[b][None, :]) for a, b in zip(ia, ib)]
    norm = np.sqrt(total[ia] * total[ib])
    target = cov[ia, ib] / np.sqrt(sample_var[ia] * sample_var[ib])
    wt = weight[ia, ib]
    best = (np.inf, np.inf, 0.0)
    for ell in ELL_GRID:
        diag = np.array([float(np.mean(np.exp(-d / ell))) for d in self_dist])
        kappa = np.array([float(np.mean(np.exp(-d / ell))) for d in pair_dist]) / np.sqrt(diag[ia] * diag[ib])
        for phi in PHI_GRID:
            h = red_var + phi * white_var
            sse = float(np.sum(wt * (target - np.sqrt(h[ia] * h[ib]) / norm * kappa) ** 2))
            if sse < best[0] - 1e-12:
                best = (sse, float(ell), float(phi))
    return best[1], best[2], True


def _power_fit(lam, var, lam_ref, alpha=None):
    ok = var > 0
    if not ok.any():
        return 0.0, 0.0
    x = np.log(lam[ok] / lam_ref)
    y = 0.5 * np.log(var[ok])
    if alpha is None:
        alpha = float(np.polyfit(x, y, 1)[0]) if ok.sum() >= 2 and np.ptp(x) > 0 else 0.0
    return float(np.exp(np.mean(y - alpha * x))), float(alpha)


# CATALOG
def catalog(data, alpha=None, law="powerlaw", min_points=10, match_tol=None, significance=2.5,
            n_sub=9, **scrape_kw):
    if law not in ("powerlaw", "interp"):
        raise ValueError("law must be 'powerlaw' or 'interp'")
    res = scrape(data, **scrape_kw)
    keep = []
    for i in range(len(res.resid)):
        if len(res.resid[i]) >= min_points:
            keep.append(i)
        else:
            print(f"[WARN] tmpest: channel {res.names[i]} has {len(res.resid[i])} usable points "
                  f"(< {min_points}); skipped")
    if not keep:
        raise ValueError(f"no channel has at least {min_points} usable points; a noise profile "
                         "needs a time series (many epochs) in at least one channel")
    times = [res.time[i] for i in keep]
    resids = [res.resid[i] for i in keep]
    names = [res.names[i] for i in keep]
    lam = res.wavelength[keep]
    width = None if res.width is None else res.width[keep]
    n_ch = len(keep)

    fits = np.array([_fit_noise(t, r, significance) for t, r in zip(times, resids)])
    white_var, red_var, tau_channel = fits[:, 0], fits[:, 1], fits[:, 2]
    red = red_var > 0
    if red.any():
        weight = red_var[red] * np.array([len(r) for r in resids])[red]
        tau = float(np.exp(np.sum(weight * np.log(tau_channel[red])) / np.sum(weight)))
    else:
        tau = np.nan

    sample_var = np.array([np.var(r) for r in resids])
    cov = np.full((n_ch, n_ch), np.nan)
    weight = np.zeros((n_ch, n_ch))
    cadence = [_cadence(t) for t in times]
    for a in range(n_ch):
        for b in range(a + 1, n_ch):
            tol = match_tol if match_tol is not None else 0.25 * min(cadence[a], cadence[b])
            c, n = _matched_cov(times[a], resids[a], times[b], resids[b], tol)
            if n:
                cov[a, b] = cov[b, a] = c
                weight[a, b] = n
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.clip(cov / np.sqrt(np.outer(sample_var, sample_var)), -1.0, 1.0)
    u = sub_points(lam, width, n_sub)
    ell, phi, measured = _fit_coherence(u, cov, weight, red_var, white_var, sample_var)

    dilution = np.array([_kernel(u[c], u[c], ell) for c in range(n_ch)])
    red_point_var = red_var / dilution
    cw_point_var = phi * white_var / dilution
    white_density = (1.0 - phi) * white_var * (1.0 if width is None else width)
    lambda_ref = float(np.exp(np.mean(np.log(lam))))
    red_amp, red_alpha = _power_fit(lam, red_point_var, lambda_ref, alpha)
    cw_amp, cw_alpha = _power_fit(lam, cw_point_var, lambda_ref)
    white_amp, white_alpha = _power_fit(lam, white_density, lambda_ref)
    return ChromaticProfile(kind=res.kind, names=names, wavelength=lam, width=width,
                            n_points=np.array([len(r) for r in resids], dtype=float),
                            red_var=red_var, white_var=white_var, tau_channel=tau_channel,
                            corr=corr, tau=tau, ell=ell, coherent_fraction=phi,
                            coherence_measured=measured, red_point_var=red_point_var,
                            cw_point_var=cw_point_var, white_density=white_density,
                            lambda_ref=lambda_ref, red_amp=red_amp, red_alpha=red_alpha,
                            cw_amp=cw_amp, cw_alpha=cw_alpha,
                            white_amp=white_amp, white_alpha=white_alpha, law=law)
