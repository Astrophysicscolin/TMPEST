import numpy as np

from . import _kernels
from .catalog import catalog, sub_points
from .data import Photometry, Spectra
from .profile import ChromaticProfile


# GENERATION
def _check_rng(rng):
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy Generator, e.g. np.random.default_rng(42)")


def ou_filter(white, x, scale, axis=0):
    data = np.array(white, dtype=float, order="C", copy=True)
    if data.ndim == 1:
        return _kernels.ou_filter_inplace(data[:, None], x, scale, axis=0)[:, 0]
    return _kernels.ou_filter_inplace(data, x, scale, axis=axis)


def _is_sorted(x):
    return bool(np.all(x[1:] >= x[:-1]))


def _coherent_field(profile, time, log_lam, rng, tau):
    z = rng.standard_normal((len(time), len(log_lam)))
    ordered = _is_sorted(time) and _is_sorted(log_lam)
    it = slice(None) if ordered else np.argsort(time, kind="stable")
    il = slice(None) if ordered else np.argsort(log_lam, kind="stable")
    if tau is not None:
        _kernels.ou_filter_inplace(z, time[it], tau, axis=0)
    _kernels.ou_filter_inplace(z, log_lam[il], profile.ell, axis=1)
    if ordered:
        return z
    out = np.empty_like(z)
    out[np.ix_(it, il)] = z
    return out


def _coherent_noise(profile, time, lam, rng, red, white):
    noise = None
    log_lam = np.log(lam)
    if red > 0 and profile.has_red:
        noise = _coherent_field(profile, time, log_lam, rng, profile.tau)
        noise *= red
        noise *= profile.sigma_red(lam)[None, :]
    if white > 0 and profile.coherent_fraction > 0:
        field = _coherent_field(profile, time, log_lam, rng, None)
        field *= white
        field *= profile.sigma_coherent_white(lam)[None, :]
        noise = field if noise is None else np.add(noise, field, out=noise)
    return np.zeros((len(time), len(lam))) if noise is None else noise


def _draw_spectra(profile, target, rng, red, white):
    lam = target.wavelength
    noise = _coherent_noise(profile, target.time, lam, rng, red, white)
    if white > 0:
        z = rng.standard_normal(noise.shape)
        z *= white
        z *= profile.sigma_white(lam, target.pixel_width())[None, :]
        noise += z
    return noise


def _draw_photometry(profile, target, rng, red, white, n_sub):
    names = target.names
    lam = target.wavelength_array()
    width = target.width_array()
    sizes = [len(target.time[n]) for n in names]
    union, inverse = np.unique(np.concatenate([target.time[n] for n in names]), return_inverse=True)
    rows = np.split(inverse, np.cumsum(sizes)[:-1])
    u = sub_points(lam, width, n_sub)
    field = _coherent_noise(profile, union, np.exp(u.ravel()), rng, red, white)
    band = field.reshape(len(union), *u.shape).mean(axis=2)
    noise = {n: band[rows[b], b].copy() for b, n in enumerate(names)}
    if white > 0:
        sigma = profile.sigma_white(lam, width)
        for b, n in enumerate(names):
            noise[n] += white * sigma[b] * rng.standard_normal(sizes[b])
    return noise


def draw(profile, target, rng, red=1.0, white=1.0, n_sub=9):
    _check_rng(rng)
    if isinstance(target, Spectra):
        return _draw_spectra(profile, target, rng, red, white)
    if isinstance(target, Photometry):
        return _draw_photometry(profile, target, rng, red, white, n_sub)
    raise TypeError("target must be a Photometry or Spectra object")


# APPLY
def apply(profile, target, rng, red=1.0, white=1.0, n_sub=9):
    noise = draw(profile, target, rng, red=red, white=white, n_sub=n_sub)
    if isinstance(target, Spectra):
        noise += 1.0
        noise *= target.flux
        return Spectra(target.wavelength.copy(), noise, target.time.copy())
    for n in target.names:
        noise[n] += 1.0
        noise[n] *= target.flux[n]
    flux = noise
    time = {n: target.time[n].copy() for n in target.names}
    width = None if target.width is None else dict(target.width)
    return Photometry(time, flux, dict(target.wavelength), width)


def transfer(source, target, rng, red=1.0, white=1.0, n_sub=9, **catalog_kw):
    profile = source if isinstance(source, ChromaticProfile) else catalog(source, n_sub=n_sub, **catalog_kw)
    return apply(profile, target, rng, red=red, white=white, n_sub=n_sub), profile


def phot_to_spec(photometry, spectra, rng, **kw):
    if not isinstance(photometry, Photometry) or not isinstance(spectra, Spectra):
        raise TypeError("phot_to_spec(photometry, spectra, rng): real Photometry first, synthetic Spectra second")
    return transfer(photometry, spectra, rng, **kw)


def spec_to_phot(spectra, photometry, rng, **kw):
    if not isinstance(spectra, Spectra) or not isinstance(photometry, Photometry):
        raise TypeError("spec_to_phot(spectra, photometry, rng): real Spectra first, synthetic Photometry second")
    return transfer(spectra, photometry, rng, **kw)


# CHROMATIC SIGNAL
def scale_depth(flux, wavelength, alpha, lambda_ref, baseline=1.0):
    factor = (np.asarray(wavelength, dtype=float) / lambda_ref) ** alpha
    return baseline - (baseline - np.asarray(flux, dtype=float)) * factor
