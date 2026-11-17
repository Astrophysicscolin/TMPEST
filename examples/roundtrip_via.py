import contextlib
import dataclasses
import io
import sys

import numpy as np

from datasets import read_gaia_epoch_csv
from demo import N_EPOCHS, synthetic_spectra
from roundtrip_causes import fmt, time_factor

import tmpest as tp

cat = sys.modules["tmpest.catalog"]
ORIG = (cat._fit_noise, cat._matched_cov, cat._fit_coherence)
N_SEEDS = 300
N_SUB = 9
N_CHUNKS = 8


def restore():
    cat._fit_noise, cat._matched_cov, cat._fit_coherence = ORIG


def expected_var(profile, phot):
    lam = phot.wavelength_array()
    width = phot.width_array()
    u = cat.sub_points(lam, width, N_SUB)
    out = np.zeros(len(phot.names))
    for b, name in enumerate(phot.names):
        kw = np.exp(-np.abs(u[b][:, None] - u[b][None, :]) / profile.ell)
        l = np.exp(u[b])
        sr = profile.sigma_red(l) if profile.has_red else np.zeros_like(l)
        sc = profile.sigma_coherent_white(l) if profile.coherent_fraction > 0 else np.zeros_like(l)
        t = phot.time[name]
        red = float(np.mean(np.outer(sr, sr) * kw)) * (time_factor(t, profile.tau) if profile.has_red else 0.0)
        white = float(np.mean(np.outer(sc, sc) * kw)) + float(profile.sigma_white(lam, width)[b] ** 2)
        out[b] = red + white * time_factor(t, None)
    return out


def chunk_population(profile, clean):
    lam = clean.wavelength
    edges = np.geomspace(lam[0], lam[-1], N_CHUNKS + 1)
    idx = np.clip(np.searchsorted(edges, lam, side="right") - 1, 0, N_CHUNKS - 1)
    logl = np.log(lam)
    kw = np.exp(-np.abs(logl[:, None] - logl[None, :]) / profile.ell)
    sr = profile.sigma_red(lam) if profile.has_red else np.zeros_like(lam)
    sc = profile.sigma_coherent_white(lam) if profile.coherent_fraction > 0 else np.zeros_like(lam)
    sw = profile.sigma_white(lam, clean.pixel_width())
    red = np.zeros((N_CHUNKS, N_CHUNKS))
    cw = np.zeros((N_CHUNKS, N_CHUNKS))
    inc = np.zeros(N_CHUNKS)
    for a in range(N_CHUNKS):
        ia = idx == a
        inc[a] = np.sum(sw[ia] ** 2) / ia.sum() ** 2
        for b in range(N_CHUNKS):
            ib = idx == b
            k = kw[np.ix_(ia, ib)]
            red[a, b] = np.mean(np.outer(sr[ia], sr[ib]) * k)
            cw[a, b] = np.mean(np.outer(sc[ia], sc[ib]) * k)
    return red, cw, inc


def back_catalog(noisy, clean, law, fits=None, covs=None, coherence=None):
    if fits is not None:
        q = iter(fits)
        cat._fit_noise = lambda *a: next(q)
    if covs is not None:
        c = iter(covs)
        cat._matched_cov = lambda *a: next(c)
    if coherence is not None:
        cat._fit_coherence = lambda *a: coherence
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            return tp.catalog(noisy, signal=clean.flux, n_chunks=N_CHUNKS, min_points=8, law=law)
    finally:
        restore()


def main(paths):
    for path in paths:
        phot = read_gaia_epoch_csv(path)
        flat = tp.Photometry({b: phot.time[b] for b in phot.names},
                             {b: np.ones(len(phot.time[b])) for b in phot.names})
        t_all = phot.time["gaia_G"]
        print(f"\n{path}")
        print(f"  band centre, width [nm]: {[(round(float(a)), round(float(w))) for a, w in zip(phot.wavelength_array(), phot.width_array())]}; spectra cover 380-1000 nm")
        for law in ("powerlaw", "interp"):
            with contextlib.redirect_stdout(io.StringIO()):
                p1 = tp.catalog(phot, law=law)
            total = p1.red_var + p1.white_var
            ref = synthetic_spectra(np.sort(t_all)[:N_EPOCHS])
            red, cw, inc = chunk_population(p1, ref)
            pop_fits = [(cw[c, c] + inc[c], red[c, c], p1.tau if p1.has_red else np.nan) for c in range(N_CHUNKS)]
            pop_covs = [(red[a, b] + cw[a, b], N_EPOCHS) for a in range(N_CHUNKS) for b in range(a + 1, N_CHUNKS)]
            ideal = back_catalog(ref, ref, law, pop_fits, pop_covs)
            pop_tot = np.array([f[0] + f[1] for f in pop_fits])
            pop_red = np.array([f[1] for f in pop_fits])
            acc = {k: np.zeros(3) for k in ("all", "time_est", "coh_est", "all_tau_true", "all_coh_ideal_tau_true", "mc")}
            est_tot, est_red, gate, taus, ells, phis, samp = [], [], [], [], [], [], []
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(seed)
                epochs = np.sort(rng.choice(t_all, size=min(N_EPOCHS, len(t_all)), replace=False))
                clean = synthetic_spectra(epochs)
                noisy = tp.apply(p1, clean, rng)
                back = back_catalog(noisy, clean, law)
                acc["all"] += expected_var(back, phot)
                synth = tp.apply(back, flat, rng)
                acc["mc"] += [np.var(synth.flux[b]) for b in phot.names]
                acc["time_est"] += expected_var(back_catalog(noisy, clean, law, coherence=(ideal.ell, ideal.coherent_fraction, True)), phot)
                acc["coh_est"] += expected_var(back_catalog(noisy, clean, law, fits=pop_fits), phot)
                if p1.has_red and back.has_red:
                    acc["all_tau_true"] += expected_var(dataclasses.replace(back, tau=p1.tau), phot)
                else:
                    acc["all_tau_true"] += expected_var(back, phot)
                est_tot.append(back.red_var + back.white_var)
                est_red.append(back.red_var)
                gate.append(back.red_var > 0)
                taus.append(back.tau)
                ells.append(back.ell)
                phis.append(back.coherent_fraction)
                res = tp.scrape(noisy, signal=clean.flux, n_chunks=N_CHUNKS)
                samp.append([np.var(r) for r in res.resid])
            est_tot, est_red, gate, samp = map(np.array, (est_tot, est_red, gate, samp))
            taus = np.array(taus, dtype=float)
            r = lambda v: fmt(np.sqrt(v / total))
            print(f"  LAW {law}   (rms synthetic / measured total; G BP RP)")
            print(f"    direct, as catalogued:                         {r(expected_var(p1, phot))}")
            print(f"    via, re-catalogued from population values:     {r(expected_var(ideal, phot))}   [ideal ell={ideal.ell:.3g} phi={ideal.coherent_fraction:.2f}; p1 ell={p1.ell:.3g} phi={p1.coherent_fraction:.2f}]")
            print(f"    via, coherence estimated, time fits exact:     {r(acc['coh_est'] / N_SEEDS)}")
            print(f"    via, time fits estimated, coherence exact:     {r(acc['time_est'] / N_SEEDS)}")
            print(f"    via, everything estimated:                     {r(acc['all'] / N_SEEDS)}")
            print(f"    via, everything estimated, tau reset to truth: {r(acc['all_tau_true'] / N_SEEDS)}")
            print(f"    via, full Monte Carlo with second draw:        {r(acc['mc'] / N_SEEDS)}")
            print(f"    chunks: red share of population variance:      {fmt(pop_red / pop_tot)}")
            print(f"    chunks: mean sample var / population:          {fmt(samp.mean(axis=0) / pop_tot)}")
            print(f"    chunks: mean catalogued total / population:    {fmt(est_tot.mean(axis=0) / pop_tot)}")
            print(f"    chunks: mean catalogued red / population red:  {fmt(est_red.mean(axis=0) / np.where(pop_red > 0, pop_red, np.nan))}")
            print(f"    chunks: red-noise gate pass rate:              {fmt(gate.mean(axis=0))}")
            print(f"    expected sample-var factor of red at 40 epochs: {time_factor(np.sort(t_all)[:N_EPOCHS], p1.tau):.2f}")
            print(f"    back tau: true {p1.tau:.1f}, median {np.nanmedian(taus):.1f}, 10-90% {np.nanpercentile(taus, 10):.1f}-{np.nanpercentile(taus, 90):.1f}, no red in {np.mean(~np.isfinite(taus)):.2f}")
            print(f"    back ell median {np.median(ells):.3g}, phi median {np.median(phis):.2f} (10-90% {np.percentile(phis, 10):.2f}-{np.percentile(phis, 90):.2f})")


if __name__ == "__main__":
    main(sys.argv[1:])
