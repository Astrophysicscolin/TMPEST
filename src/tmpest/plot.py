import numpy as np

from .catalog import _fit_noise, _variogram, scrape
from .data import Photometry, Spectra


def _axes(ax):
    if ax is not None:
        return ax
    try:
        import matplotlib.pyplot as plt
    except ImportError as err:
        raise ImportError("tmpest.plot needs matplotlib: pip install 'tmpest[plot]'") from err
    return plt.subplots()[1]


# DATA
def plot_photometry(photometry, ax=None, normalize=True, **kw):
    ax = _axes(ax)
    kw.setdefault("marker", ".")
    kw.setdefault("linestyle", "none")
    for name in photometry.names:
        f = photometry.flux[name]
        if normalize:
            f = f / np.nanmedian(f)
        ax.plot(photometry.time[name], f, label=f"{name} ({photometry.wavelength[name]:.0f} nm)", **kw)
    ax.set_xlabel("time")
    ax.set_ylabel("flux / median" if normalize else "flux")
    ax.legend()
    return ax


def plot_spectra(spectra, ax=None, kind="lines", max_lines=8, normalize=False, **kw):
    ax = _axes(ax)
    flux = spectra.flux
    if normalize:
        flux = flux / np.nanmedian(flux, axis=0)
    if kind == "image":
        order = np.argsort(spectra.time, kind="stable")
        kw.setdefault("shading", "auto")
        mesh = ax.pcolormesh(spectra.wavelength, np.arange(len(order)), flux[order], **kw)
        ax.figure.colorbar(mesh, ax=ax, label="flux / median" if normalize else "flux")
        ax.set_ylabel("epoch index (time order)")
    elif kind == "lines":
        rows = np.unique(np.linspace(0, len(spectra.time) - 1, min(max_lines, len(spectra.time))).astype(int))
        kw.setdefault("linewidth", 0.8)
        for i in rows:
            ax.plot(spectra.wavelength, flux[i], label=f"t = {spectra.time[i]:.4g}", **kw)
        ax.set_ylabel("flux / median" if normalize else "flux")
        ax.legend()
    else:
        raise ValueError("kind must be 'lines' or 'image'")
    ax.set_xlabel("wavelength [nm]")
    return ax


def plot_noise(clean, noisy, ax=None, **kw):
    if type(clean) is not type(noisy):
        raise TypeError("clean and noisy must both be Photometry or both be Spectra")
    if isinstance(clean, Spectra):
        ratio = Spectra(clean.wavelength, noisy.flux / clean.flux - 1.0, clean.time)
        kw.setdefault("cmap", "RdBu_r")
        ax = plot_spectra(ratio, ax=ax, kind="image", **kw)
        ax.collections[-1].colorbar.set_label("noisy / clean - 1")
        return ax
    ratio = Photometry(clean.time, {n: noisy.flux[n] / clean.flux[n] for n in clean.names},
                       clean.wavelength, clean.width)
    ax = plot_photometry(ratio, ax=ax, normalize=False, **kw)
    ax.set_ylabel("noisy / clean")
    return ax


# DIAGNOSTICS
def plot_residuals(data, ax=None, offset=None, **scrape_kw):
    ax = _axes(ax)
    res = scrape(data, **scrape_kw)
    if offset is None:
        spread = [np.std(r) for r in res.resid if len(r)]
        offset = 6.0 * max(spread) if spread else 1.0
    for i, name in enumerate(res.names):
        ax.plot(res.time[i], res.resid[i] + i * offset, ".", label=str(name))
        ax.axhline(i * offset, color="0.7", linewidth=0.5)
    ax.set_xlabel("time")
    ax.set_ylabel("relative residual + offset")
    ax.legend()
    return ax


def plot_variogram(data, channel=0, ax=None, significance=2.5, **scrape_kw):
    ax = _axes(ax)
    res = scrape(data, **scrape_kw)
    i = res.names.index(channel) if channel in res.names else int(channel)
    t, r = res.time[i], res.resid[i]
    lag, gamma, _ = _variogram(t, r)
    white, red, tau = _fit_noise(t, r, significance)
    ax.loglog(lag, gamma, "o", label=f"{res.names[i]} measured")
    grid = np.geomspace(lag.min(), lag.max(), 200)
    if red > 0:
        ax.loglog(grid, white + red * (1.0 - np.exp(-grid / tau)), label=f"fit: tau = {tau:.3g}")
    ax.axhline(white, color="0.5", linestyle=":", label="white level")
    ax.set_xlabel("time separation")
    ax.set_ylabel("variogram (relative flux squared)")
    ax.legend()
    return ax


def plot_profile(profile, ax=None, n=200):
    ax = _axes(ax)
    lam = profile.wavelength
    span = np.geomspace(0.8 * lam.min(), 1.25 * lam.max(), n)
    red = np.sqrt(profile.red_var)
    white = np.sqrt(profile.white_var)
    ax.plot(lam[red > 0], red[red > 0], "o", color="C3", label="red jitter, measured per channel")
    ax.plot(lam[white > 0], white[white > 0], "s", color="C0", label="time-white, measured per channel")
    if profile.has_red:
        ax.plot(span, profile.sigma_red(span), color="C3", label=f"red law (per wavelength), tau = {profile.tau:.3g}")
    if profile.coherent_fraction > 0:
        ax.plot(span, profile.sigma_coherent_white(span), color="C0", linestyle="--",
                label=f"coherent white law, share = {profile.coherent_fraction:.2f}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("wavelength [nm]")
    ax.set_ylabel("relative noise (sigma)")
    ax.legend(fontsize="small")
    return ax


def plot_correlation(profile, ax=None, **kw):
    ax = _axes(ax)
    corr = np.array(profile.corr, dtype=float)
    np.fill_diagonal(corr, 1.0)
    kw.setdefault("vmin", -1.0)
    kw.setdefault("vmax", 1.0)
    kw.setdefault("cmap", "RdBu_r")
    image = ax.imshow(corr, **kw)
    ax.figure.colorbar(image, ax=ax, label="correlation at matched times")
    ticks = np.arange(len(profile.names))
    ax.set_xticks(ticks, [str(n) for n in profile.names], rotation=90)
    ax.set_yticks(ticks, [str(n) for n in profile.names])
    return ax
