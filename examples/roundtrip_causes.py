import contextlib
import io
import sys
import numpy as np
from datasets import read_gaia_epoch_csv
import tmpest as tp

cat = sys.modules["tmpest.catalog"]
ORIG_FIT_NOISE = cat._fit_noise
ORIG_FIT_COH = cat._fit_coherence
N_MC = 2000
N_SUB = 9

def time_factor(t, tau):
    t = t[np.isfinite(t)]
    if tau is None:
        return 1.0 - 1.0 / len(t)
    return 1.0 - float(np.mean(np.exp(-np.abs(t[:, None] - t[None, :]) / tau)))

def predict(profile, phot, taus=None):
    lam = phot.wavelength_array()
    width = phot.width_array()
    u = cat.sub_points(lam, width, N_SUB)
    out = []
    for b, name in enumerate(phot.names):
        kw = np.exp(-np.abs(u[b][:, None] - u[b][None, :]) / profile.ell)
        l = np.exp(u[b])
        sr = profile.sigma_red(l) if profile.has_red else np.zeros_like(l)
        sc = profile.sigma_coherent_white(l) if profile.coherent_fraction > 0 else np.zeros_like(l)
        red = float(np.mean(np.outer(sr, sr) * kw))
        cw = float(np.mean(np.outer(sc, sc) * kw))
        inc = float(profile.sigma_white(lam, width)[b] ** 2)
        tau = profile.tau if taus is None else taus[b]
        t = phot.time[name]
        f_red = time_factor(t, tau) if np.isfinite(tau) else 0.0
        out.append((red, cw, inc, f_red, time_factor(t, None)))
    return np.array(out)

def monte_carlo(profile, phot):
    flat = tp.Photometry({b: phot.time[b] for b in phot.names},
                         {b: np.ones(len(phot.time[b])) for b in phot.names})
    var = np.zeros(len(phot.names))
    for seed in range(N_MC):
        synth = tp.apply(profile, flat, np.random.default_rng(seed))
        var += [np.var(synth.flux[b]) for b in phot.names]
    return var / N_MC

def quiet_catalog(phot, law):
    with contextlib.redirect_stdout(io.StringIO()):
        return tp.catalog(phot, law=law)

def ratio(profile, phot, taus=None):
    p = predict(profile, phot, taus)
    total = profile.red_var + profile.white_var
    return np.sqrt((p[:, 0] * p[:, 3] + (p[:, 1] + p[:, 2]) * p[:, 4]) / total)

def fmt(x):
    return "  ".join(f"{v:5.2f}" for v in x)

def main(paths):
    for path in paths:
        phot = read_gaia_epoch_csv(path)
        res = tp.scrape(phot)
        sample_var = np.array([np.var(r) for r in res.resid])
        print(f"\n{path}\n  bands: {phot.names}")
        base = quiet_catalog(phot, "powerlaw")
        total = base.red_var + base.white_var
        print(f"  ell={base.ell:.3g} phi={base.coherent_fraction:.2f} tau={base.tau:.4g} tau_channel={np.round(base.tau_channel, 2)}")
        print(f"  red fraction of total per band:            {fmt(base.red_var / total)}")
        print(f"  sqrt(sample var of real resid / total):    {fmt(np.sqrt(sample_var / total))}")
        for law in ("powerlaw", "interp"):
            prof = quiet_catalog(phot, law)
            p = predict(prof, phot)
            print(f"  LAW {law}")
            print(f"    population red  generated / measured (var): {fmt(p[:, 0] / np.where(prof.red_var > 0, prof.red_var, np.nan))}   generated red var where measured 0: {['%.2e' % v for v in p[prof.red_var == 0, 0]]}")
            print(f"    population white generated / measured (var): {fmt((p[:, 1] + p[:, 2]) / prof.white_var)}")
            print(f"    population total sigma ratio (law only):     {fmt(np.sqrt((p[:, 0] + p[:, 1] + p[:, 2]) / total))}")
            print(f"    red finite-sample var factor, shared tau:    {fmt(p[:, 3])}")
            print(f"    red finite-sample var factor, channel tau:   {fmt(predict(prof, phot, prof.tau_channel)[:, 3])}")
            print(f"    predicted rms ratio (as catalogued):         {fmt(ratio(prof, phot))}")
            print(f"    Monte Carlo rms ratio ({N_MC} seeds):          {fmt(np.sqrt(monte_carlo(prof, phot) / total))}")

            fits = [ORIG_FIT_NOISE(t, r) for t, r in zip(res.time, res.resid)]
            cat._fit_coherence = lambda *a: (prof.ell, 0.0, True)
            print(f"    counterfactual phi=0 (ell kept):             {fmt(ratio(quiet_catalog(phot, law), phot))}")
            cat._fit_coherence = lambda *a: (np.inf, prof.coherent_fraction, True)
            print(f"    counterfactual ell=inf (phi kept):           {fmt(ratio(quiet_catalog(phot, law), phot))}")
            cat._fit_coherence = ORIG_FIT_COH
            tot = np.array([f[0] + f[1] for f in fits])
            frac = sum(f[1] for f in fits) / tot.sum()
            queue = iter([((1 - frac) * v, frac * v, prof.tau) for v in tot])
            cat._fit_noise = lambda *a: next(queue)
            same = quiet_catalog(phot, law)
            cat._fit_noise = ORIG_FIT_NOISE
            print(f"    counterfactual same red fraction in all bands (refit ell={same.ell:.3g} phi={same.coherent_fraction:.2f}): {fmt(ratio(same, phot))}")
            print(f"      same, population only:                     {fmt(np.sqrt(predict(same, phot)[:, :3].sum(axis=1) / tot))}")

if __name__ == "__main__":
    main(sys.argv[1:])