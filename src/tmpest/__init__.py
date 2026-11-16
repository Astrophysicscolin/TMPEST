from ._kernels import backend, num_threads, rust_available, set_backend
from .apply import apply, draw, ou_filter, phot_to_spec, scale_depth, spec_to_phot, transfer
from .bands import BANDS
from .catalog import catalog, scrape
from .data import Photometry, Residuals, Spectra, mag_to_flux
from .profile import ChromaticProfile

__version__ = "1.4.0"
__all__ = ["apply", "draw", "ou_filter", "phot_to_spec", "scale_depth", "spec_to_phot", "transfer",
           "BANDS", "catalog", "scrape", "Photometry", "Residuals", "Spectra", "mag_to_flux",
           "ChromaticProfile", "backend", "set_backend", "rust_available", "num_threads"]
