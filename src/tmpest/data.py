from dataclasses import dataclass

import numpy as np

from .bands import BANDS


def mag_to_flux(mag, zero_point=0.0):
    return 10.0 ** (-0.4 * (np.asarray(mag, dtype=float) - zero_point))


@dataclass
class Photometry:
    time: dict
    flux: dict
    wavelength: dict = None  # [nm]
    width: dict = None  # [nm]

    def __post_init__(self):
        names = list(self.flux)
        if not names:
            raise ValueError("Photometry needs at least one band")
        if self.wavelength is None:
            unknown = [n for n in names if n not in BANDS]
            if unknown:
                raise KeyError(f"no built-in passband for {unknown}; pass wavelength= (and width=) "
                               f"or name bands from {sorted(BANDS)}")
            self.wavelength = {n: BANDS[n][0] for n in names}
            if self.width is None:
                self.width = {n: BANDS[n][1] for n in names}
        self.time = {n: np.asarray(self.time[n], dtype=float) for n in names}
        self.flux = {n: np.asarray(self.flux[n], dtype=float) for n in names}
        for n in names:
            if self.time[n].shape != self.flux[n].shape or self.time[n].ndim != 1:
                raise ValueError(f"band {n}: time and flux must be 1-D arrays of equal length")

    @property
    def names(self):
        return list(self.flux)

    def wavelength_array(self):
        return np.array([self.wavelength[n] for n in self.names], dtype=float)

    def width_array(self):
        if self.width is None:
            return None
        return np.array([self.width[n] for n in self.names], dtype=float)


@dataclass
class Spectra:
    wavelength: np.ndarray  # [nm]
    flux: np.ndarray
    time: np.ndarray = None

    def __post_init__(self):
        self.wavelength = np.asarray(self.wavelength, dtype=float)
        self.flux = np.atleast_2d(np.asarray(self.flux, dtype=float))
        if self.time is None:
            self.time = np.arange(self.flux.shape[0], dtype=float)
        self.time = np.atleast_1d(np.asarray(self.time, dtype=float))
        if self.flux.shape != (len(self.time), len(self.wavelength)):
            raise ValueError(f"flux must have shape (n_epochs, n_wavelengths) = "
                             f"({len(self.time)}, {len(self.wavelength)}), got {self.flux.shape}")
        if len(self.wavelength) < 2:
            raise ValueError("Spectra needs at least 2 wavelength samples")
        order = np.argsort(self.wavelength, kind="stable")
        self.wavelength = self.wavelength[order]
        self.flux = self.flux[:, order]

    def pixel_width(self):
        return np.abs(np.gradient(self.wavelength))


@dataclass
class Residuals:
    time: list
    resid: list
    wavelength: np.ndarray  # [nm]
    width: np.ndarray  # [nm]
    names: list
    kind: str
