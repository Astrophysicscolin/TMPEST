import numpy as np
import pytest
import tmpest as tp

def truth_profile(ell=np.inf, red_amp=5e-3, cw_amp=0.0):
    lam = np.array([500.0, 600.0, 800.0])
    return tp.ChromaticProfile(kind="truth", names=["a", "b", "c"], wavelength=lam,
                               width=np.full(3, 100.0), n_points=np.zeros(3), red_var=np.zeros(3),
                               white_var=np.zeros(3), tau_channel=np.full(3, 5.0),
                               corr=np.full((3, 3), np.nan), tau=5.0, ell=ell,
                               coherent_fraction=float(cw_amp > 0), coherence_measured=True,
                               red_point_var=np.full(3, float(red_amp > 0)), cw_point_var=np.ones(3),
                               white_density=np.ones(3), lambda_ref=600.0, red_amp=red_amp,
                               red_alpha=-1.0, cw_amp=cw_amp, cw_alpha=-1.0, white_amp=2e-2,
                               white_alpha=0.0)

def flat_photometry(rng, n=1500, baseline=2000.0):
    t = np.sort(rng.uniform(0.0, baseline, n))
    names = ["gaia_G", "gaia_BP", "gaia_RP"]
    return tp.Photometry({k: t + i * 1e-4 for i, k in enumerate(names)},
                         {k: np.full(n, 1000.0) for k in names})

def flat_spectra(rng, n_t=1500, n_lam=200, baseline=2000.0):
    t = np.sort(rng.uniform(0.0, baseline, n_t))
    lam = np.linspace(350.0, 1050.0, n_lam)
    return tp.Spectra(lam, np.full((n_t, n_lam), 50.0), t)

def band_average(spec, name):
    centre, width = tp.BANDS[name]
    sel = np.abs(spec.wavelength - centre) <= 0.5 * width
    return spec.flux[:, sel].mean(axis=1)

def test_ou_filter_has_unit_variance_and_exponential_correlation():
    rng = np.random.default_rng(0)
    x = np.cumsum(rng.uniform(0.5, 1.5, 4000))
    z = tp.ou_filter(rng.standard_normal((4000, 50)), x, 10.0, axis=0)
    assert abs(z.var() - 1.0) < 0.1
    lag1 = np.mean(z[1:] * z[:-1])
    assert abs(lag1 - np.mean(np.exp(-np.diff(x) / 10.0))) < 0.05

def test_catalog_recovers_known_photometric_noise():
    rng = np.random.default_rng(1)
    truth = truth_profile()
    noisy = tp.apply(truth, flat_photometry(rng), rng)
    prof = tp.catalog(noisy)
    assert prof.kind == "photometry"
    assert 2.5 < prof.tau < 10.0
    assert abs(prof.red_alpha - (-1.0)) < 0.35
    assert abs(prof.red_amp / truth.sigma_red(prof.lambda_ref) - 1.0) < 0.25
    expect_white = truth.sigma_white(prof.wavelength, prof.width) ** 2
    assert np.all(np.abs(prof.white_var / expect_white - 1.0) < 0.35)
    assert prof.coherence_measured and prof.ell > 5.0

def test_white_series_is_not_called_red():
    rng = np.random.default_rng(2)
    t = np.sort(rng.uniform(0.0, 100.0, 800))
    phot = tp.Photometry({"x": t}, {"x": 1.0 + 1e-3 * rng.standard_normal(800)},
                         wavelength={"x": 600.0})
    prof = tp.catalog(phot)
    assert prof.red_var[0] == 0.0 and not prof.has_red
    assert abs(np.sqrt(prof.white_var[0]) / 1e-3 - 1.0) < 0.1
    assert prof.width is None

def test_phot_to_spec_preserves_band_integrated_red_noise():
    rng = np.random.default_rng(3)
    noisy_phot = tp.apply(truth_profile(), flat_photometry(rng), rng)
    spec, prof = tp.phot_to_spec(noisy_phot, flat_spectra(rng), rng, white=0.0)
    assert spec.flux.shape == (1500, 200)
    for i, name in enumerate(prof.names):
        got = np.std(band_average(spec, name) / 50.0 - 1.0)
        assert abs(got / np.sqrt(prof.red_var[i]) - 1.0) < 0.3

def test_spec_to_phot_matches_direct_synthetic_photometry():
    rng = np.random.default_rng(4)
    noisy_spec = tp.apply(truth_profile(ell=1.0), flat_spectra(rng), rng)
    phot, prof = tp.spec_to_phot(noisy_spec, flat_photometry(rng), rng, n_chunks=8)
    assert prof.kind == "spectra" and len(prof.names) == 8
    assert abs(prof.red_alpha - (-1.0)) < 0.35
    for name in phot.names:
        direct = np.std(band_average(noisy_spec, name) / 50.0 - 1.0)
        got = np.std(phot.flux[name] / 1000.0 - 1.0)
        assert abs(got / direct - 1.0) < 0.3

def test_time_white_but_wavelength_coherent_noise_survives_band_integration():
    rng = np.random.default_rng(10)
    truth = truth_profile(red_amp=0.0, cw_amp=5e-3)
    noisy_spec = tp.apply(truth, flat_spectra(rng, n_t=60), rng)
    phot, prof = tp.spec_to_phot(noisy_spec, flat_photometry(rng), rng, n_chunks=8)
    assert not prof.has_red and prof.coherent_fraction > 0.5
    for name in phot.names:
        direct = np.std(band_average(noisy_spec, name) / 50.0 - 1.0)
        got = np.std(phot.flux[name] / 1000.0 - 1.0)
        assert abs(got / direct - 1.0) < 0.3

def test_shared_red_noise_stays_coherent_when_tau_reaches_the_baseline():
    names = ["gaia_G", "gaia_BP", "gaia_RP"]
    t = np.linspace(0.0, 40.0, 40)
    clean = tp.Photometry({n: t.copy() for n in names}, {n: np.full(40, 1000.0) for n in names})
    truth = truth_profile(red_amp=5e-3)
    truth.tau = 40.0
    ratios = []
    for seed in range(10):
        rng = np.random.default_rng(seed)
        noisy = tp.apply(truth, clean, rng)
        prof = tp.catalog(noisy)
        draws = [tp.draw(prof, clean, rng) for _ in range(20)]
        real = np.corrcoef(noisy.flux["gaia_BP"], noisy.flux["gaia_RP"])[0, 1]
        ratios.append(np.mean([np.corrcoef(d["gaia_BP"], d["gaia_RP"])[0, 1] for d in draws]) / real)
    assert np.median(ratios) > 0.7

def test_independent_bands_are_not_called_coherent():
    rng = np.random.default_rng(11)
    truth = truth_profile(red_amp=0.0)
    prof = tp.catalog(tp.apply(truth, flat_photometry(rng), rng))
    assert prof.coherent_fraction < 0.1

def test_zero_noise_is_identity_and_seed_is_reproducible():
    rng = np.random.default_rng(5)
    spec = flat_spectra(rng, n_t=20, n_lam=30)
    prof = truth_profile()
    same = tp.apply(prof, spec, rng, red=0.0, white=0.0)
    assert np.array_equal(same.flux, spec.flux)
    a = tp.apply(prof, spec, np.random.default_rng(9))
    b = tp.apply(prof, spec, np.random.default_rng(9))
    assert np.array_equal(a.flux, b.flux) and not np.array_equal(a.flux, spec.flux)

def test_single_spectrum_target_works():
    prof = truth_profile()
    spec = tp.Spectra(np.linspace(400.0, 900.0, 64), np.ones(64))
    out = tp.apply(prof, spec, np.random.default_rng(6))
    assert out.flux.shape == (1, 64) and np.all(np.isfinite(out.flux))

def test_profile_save_load_round_trip(tmp_path):
    rng = np.random.default_rng(7)
    prof = tp.catalog(tp.apply(truth_profile(), flat_photometry(rng, n=300), rng))
    path = tmp_path / "profile.json"
    prof.save(path)
    back = tp.ChromaticProfile.load(path)
    assert back.names == prof.names and back.tau == prof.tau and back.ell == prof.ell
    assert np.allclose(back.sigma_red([450.0, 900.0]), prof.sigma_red([450.0, 900.0]))
    assert "tau" in back.summary()

def test_interp_law_passes_through_measured_channels():
    rng = np.random.default_rng(8)
    prof = tp.catalog(tp.apply(truth_profile(), flat_photometry(rng, n=600), rng), law="interp")
    assert np.allclose(prof.sigma_red(prof.wavelength) ** 2, prof.red_point_var)

def test_scale_depth_scales_the_drop_only():
    flux = np.array([[1.0, 1.0], [0.9, 0.9]])
    out = tp.scale_depth(flux, [500.0, 1000.0], alpha=-1.0, lambda_ref=500.0)
    assert np.allclose(out, [[1.0, 1.0], [0.9, 0.95]])

def test_clear_errors():
    t = np.arange(20.0)
    with pytest.raises(KeyError):
        tp.Photometry({"mystery": t}, {"mystery": t + 1.0})
    with pytest.raises(ValueError):
        tp.catalog(tp.Photometry({"tess": t}, {"tess": np.sin(t) - 2.0}))
    with pytest.raises(TypeError):
        tp.apply(truth_profile(), tp.Spectra([400.0, 500.0], [1.0, 1.0]), None)
    with pytest.raises(ValueError):
        tp.catalog(tp.Spectra(np.linspace(400.0, 500.0, 10), np.ones((3, 10))))
    with pytest.raises(TypeError):
        tp.phot_to_spec(tp.Spectra([400.0, 500.0], [1.0, 1.0]), None, np.random.default_rng(0))