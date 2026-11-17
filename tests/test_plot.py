import numpy as np
import pytest
import tmpest as tp

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import tmpest.plot as tpp

@pytest.fixture
def source():
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0.0, 300.0, 300))
    red = 5e-3 * tp.ou_filter(rng.standard_normal(300), t, 8.0)
    names = ["gaia_G", "gaia_BP", "gaia_RP"]
    phot = tp.Photometry({n: t for n in names},
                         {n: 100.0 * (1.0 + red + 1e-3 * rng.standard_normal(300)) for n in names})
    spec = tp.Spectra(np.linspace(400.0, 900.0, 50), np.ones((12, 50)), np.linspace(0.0, 300.0, 12))
    yield phot, spec, tp.catalog(phot), rng
    plt.close("all")

def test_every_plot_draws_on_the_given_axes_and_returns_it(source):
    phot, spec, prof, rng = source
    noisy = tp.apply(prof, spec, rng)
    calls = [lambda ax: tpp.plot_photometry(phot, ax=ax),
             lambda ax: tpp.plot_spectra(noisy, ax=ax),
             lambda ax: tpp.plot_spectra(noisy, ax=ax, kind="image", normalize=True),
             lambda ax: tpp.plot_noise(spec, noisy, ax=ax),
             lambda ax: tpp.plot_noise(phot, tp.apply(prof, phot, rng), ax=ax),
             lambda ax: tpp.plot_residuals(phot, ax=ax),
             lambda ax: tpp.plot_variogram(phot, "gaia_BP", ax=ax),
             lambda ax: tpp.plot_profile(prof, ax=ax),
             lambda ax: tpp.plot_correlation(prof, ax=ax)]
    for call in calls:
        fig, ax = plt.subplots()
        assert call(ax) is ax
        assert ax.lines or ax.collections or ax.images
        plt.close(fig)

def test_plots_make_their_own_axes_and_leave_style_alone(source):
    phot, _, prof, _ = source
    before = dict(matplotlib.rcParams)
    ax = tpp.plot_profile(prof)
    assert ax.figure is not None and len(ax.get_legend().get_texts()) >= 2
    assert tpp.plot_photometry(phot, color="k").lines[0].get_color() == "k"
    assert dict(matplotlib.rcParams) == before


def test_plot_errors_are_clear(source):
    phot, spec, _, _ = source
    with pytest.raises(ValueError):
        tpp.plot_spectra(spec, kind="surface")
    with pytest.raises(TypeError):
        tpp.plot_noise(phot, spec)