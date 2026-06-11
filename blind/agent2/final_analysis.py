"""Final characterization at the critical point estimates.

For each (construction, p) group in phase4.npz:
  - disorder-averaged rho(t) with per-sample spread
  - alpha_eff(t) in sliding log-windows (+/- 0.4 decades), with jackknife error
    bars over disorder samples
  - global power-law fit  ln rho = -alpha ln t + c   over t in [1e3, T]
  - log-decay fit         ln rho = -x ln(ln t) + c   over the same range
  - comparison of fit residuals (RMS) for the two forms
"""

import sys
from collections import defaultdict

import numpy as np

WIN = 0.4  # half-window in ln t ... NOTE: in decades use 0.4*ln(10)


def load(path):
    d = np.load(path)
    rec = d["rec"]
    runs = defaultdict(list)
    for k in d.files:
        if k == "rec":
            continue
        constr, p, offset, rep = k.split("_")
        runs[(constr, float(p))].append((k, d[k]))
    return rec, runs


def alpha_at(rec, rho, tcen, half=0.4):
    """Slope of ln rho vs ln t in window log10(t) in [log10(tcen)-half, +half]."""
    m = (rho > 0) & (rec >= tcen / 10**half) & (rec <= tcen * 10**half)
    if m.sum() < 5:
        return np.nan
    A = np.polyfit(np.log(rec[m].astype(float)), np.log(rho[m]), 1)
    return -A[0]


def main(path):
    rec, runs = load(path)
    tmax = float(rec[-1])
    probe_ts = [1e3, 1e4, 1e5, 3e5, tmax]
    for (constr, p) in sorted(runs):
        keys = [k for k, _ in runs[(constr, p)]]
        rhos = np.array([r for _, r in runs[(constr, p)]])
        n = rhos.shape[0]
        rho = rhos.mean(axis=0)
        print(f"\n=== {constr}  p={p:.4f}  (n={n} disorder samples, T={int(tmax)}) ===")
        print("samples:", ", ".join(keys))
        # rho at probe times
        for tv in [1e2, 1e3, 1e4, 1e5, 3e5, tmax]:
            i = min(np.searchsorted(rec, tv), rec.size - 1)
            vals = rhos[:, i]
            print(f"  rho(t={int(rec[i])}) = {rho[i]:.4e}   per-sample: "
                  + " ".join(f"{v:.3e}" for v in vals))
        # alpha_eff at probe times: mean curve + jackknife over samples
        print("  alpha_eff (from averaged rho; +/- = jackknife over samples):")
        for tv in probe_ts:
            tc = min(tv, tmax / 10**0.4)  # keep window inside data
            a_full = alpha_at(rec, rho, tc)
            if n > 1:
                jk = []
                for i in range(n):
                    sub = np.delete(rhos, i, axis=0).mean(axis=0)
                    jk.append(alpha_at(rec, sub, tc))
                jk = np.array(jk)
                err = np.sqrt((n - 1) * np.mean((jk - jk.mean()) ** 2))
            else:
                err = np.nan
            print(f"    t~{tv:.0e} (window center {tc:.2e}): "
                  f"alpha_eff = {a_full:.4f} +/- {err:.4f}")
        # global fits over [1e3, tmax]
        m = (rho > 0) & (rec >= 1e3)
        lt = np.log(rec[m].astype(float))
        lr = np.log(rho[m])
        llt = np.log(lt)
        A_pl = np.polyfit(lt, lr, 1)
        res_pl = lr - np.polyval(A_pl, lt)
        A_lg = np.polyfit(llt, lr, 1)
        res_lg = lr - np.polyval(A_lg, llt)
        print(f"  power-law fit t in [1e3,{int(tmax)}]: alpha = {-A_pl[0]:.4f}, "
              f"RMS resid = {res_pl.std():.4f}")
        print(f"  (ln t)^-x  fit t in [1e3,{int(tmax)}]: x     = {-A_lg[0]:.4f}, "
              f"RMS resid = {res_lg.std():.4f}")
        # curvature check: alpha in first vs second half of the fit range
        mid = np.exp(0.5 * (lt[0] + lt[-1]))
        a_lo = alpha_at(rec, rho, np.exp(0.5 * (lt[0] + np.log(mid))), half=0.5)
        a_hi = alpha_at(rec, rho, np.exp(0.5 * (np.log(mid) + lt[-1])), half=0.5)
        print(f"  curvature: alpha(first half) = {a_lo:.4f}, alpha(second half) = {a_hi:.4f}")


if __name__ == "__main__":
    main(sys.argv[1])
