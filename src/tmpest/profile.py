import json
from dataclasses import asdict, dataclass

import numpy as np

ARRAY_FIELDS = ("wavelength", "width", "red_var", "white_var", "tau_channel", "n_points",
                "red_point_var", "cw_point_var", "white_density", "corr")


@dataclass
class ChromaticProfile:
    kind: str
    names: list
    wavelength: np.ndarray  # [nm]
    width: np.ndarray  # [nm]
    n_points: np.ndarray
    red_var: np.ndarray
    white_var: np.ndarray
    tau_channel: np.ndarray
    corr: np.ndarray
    tau: float
    ell: float
    coherent_fraction: float
    coherence_measured: bool
    red_point_var: np.ndarray
    cw_point_var: np.ndarray
    white_density: np.ndarray
    lambda_ref: float  # [nm]
    red_amp: float
    red_alpha: float
    cw_amp: float
    cw_alpha: float
    white_amp: float
    white_alpha: float
    law: str = "powerlaw"

    @property
    def has_red(self):
        return bool(np.isfinite(self.tau) and np.any(self.red_point_var > 0))

    def _interp(self, lam, var):
        ok = var > 0
        if not ok.any():
            return np.zeros_like(lam)
        order = np.argsort(self.wavelength[ok])
        x = np.log(self.wavelength[ok][order])
        y = 0.5 * np.log(var[ok][order])
        return np.exp(np.interp(np.log(lam), x, y))

    def sigma_red(self, wavelength):
        lam = np.asarray(wavelength, dtype=float)
        if self.law == "interp":
            return self._interp(lam, self.red_point_var)
        return self.red_amp * (lam / self.lambda_ref) ** self.red_alpha

    def sigma_coherent_white(self, wavelength):
        lam = np.asarray(wavelength, dtype=float)
        if self.law == "interp":
            return self._interp(lam, self.cw_point_var)
        return self.cw_amp * (lam / self.lambda_ref) ** self.cw_alpha

    def sigma_white(self, wavelength, width=None):
        lam = np.asarray(wavelength, dtype=float)
        if self.law == "interp":
            density = self._interp(lam, self.white_density)
        else:
            density = self.white_amp * (lam / self.lambda_ref) ** self.white_alpha
        if self.width is None:
            return density
        if width is None:
            raise ValueError("this profile stores white noise per unit wavelength, so the target "
                             "needs channel widths (Photometry width= or a Spectra wavelength grid)")
        return density / np.sqrt(np.asarray(width, dtype=float))

    def summary(self):
        lines = [f"ChromaticProfile from {self.kind}: {len(self.names)} channels",
                 f"  red timescale tau = {self.tau:.4g} (time units of the input)",
                 f"  wavelength coherence ell = {self.ell:.4g} in ln(wavelength)"
                 f" ({'measured' if self.coherence_measured else 'assumed, not measured'})",
                 f"  red law: sigma = {self.red_amp:.4g} * (lambda / {self.lambda_ref:.1f} nm)^{self.red_alpha:.3f}",
                 f"  wavelength-coherent share of the time-white noise = {self.coherent_fraction:.2f}",
                 f"  white law exponent = {self.white_alpha:.3f}",
                 "  channel      lambda[nm]  n_points  sigma_red   sigma_white  tau"]
        for i, name in enumerate(self.names):
            lines.append(f"  {str(name):<12} {self.wavelength[i]:>9.1f}  {int(self.n_points[i]):>8d}  "
                         f"{np.sqrt(self.red_var[i]):.3e}   {np.sqrt(self.white_var[i]):.3e}    "
                         f"{self.tau_channel[i]:.4g}")
        return "\n".join(lines)

    def save(self, path):
        out = asdict(self)
        for key in ARRAY_FIELDS:
            out[key] = None if out[key] is None else np.asarray(out[key]).tolist()
        out["names"] = [str(n) for n in self.names]
        with open(path, "w") as fh:
            json.dump(out, fh, indent=2)

    @classmethod
    def load(cls, path):
        with open(path) as fh:
            raw = json.load(fh)
        for key in ARRAY_FIELDS:
            raw[key] = None if raw[key] is None else np.asarray(raw[key], dtype=float)
        return cls(**raw)
