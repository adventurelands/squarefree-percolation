"""Analysis helpers: print rho(t) tables and effective exponents from scan .npz files.

alpha_eff(t) = -d ln rho / d ln t, estimated by least-squares slope of ln rho vs ln t
in a sliding window of width `win` decades.
"""

import sys
from collections import defaultdict

import numpy as np


def load(path):
    d = np.load(path)
    rec = d["rec"]
    runs = defaultdict(list)  # (constr, p) -> list of rho arrays
    for k in d.files:
        if k == "rec":
            continue
        constr, p, offset, rep = k.split("_")
        runs[(constr, float(p))].append(d[k])
    return rec, runs


def alpha_eff_curve(rec, rho, win=0.5, npts=40):
    """Return (t_mid, alpha_eff) over sliding log windows. Ignores rho<=0 points."""
    m = rho > 0
    lt, lr = np.log(rec[m].astype(float)), np.log(rho[m])
    tmids, alphas = [], []
    lo, hi = lt[0], lt[-1]
    for cen in np.linspace(lo + win, hi - 0.01, npts):
        sel = (lt >= cen - win) & (lt <= cen + win)
        if sel.sum() < 5:
            continue
        A = np.polyfit(lt[sel], lr[sel], 1)
        tmids.append(np.exp(cen))
        alphas.append(-A[0])
    return np.array(tmids), np.array(alphas)


def summarize(path, win=0.5):
    rec, runs = load(path)
    pick = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000]
    for (constr, p) in sorted(runs):
        rhos = np.array(runs[(constr, p)])
        rho = rhos.mean(axis=0)
        n = rhos.shape[0]
        line = [f"{constr} p={p:.4f} (n={n}):"]
        for tv in pick:
            i = np.searchsorted(rec, tv)
            if i < rec.size:
                line.append(f"t={rec[i]}:{rho[i]:.3e}")
        print("  ".join(line))
        # alpha_eff at a few times
        tm, al = alpha_eff_curve(rec, rho, win=win)
        if tm.size:
            sel_pts = [1e2, 1e3, 1e4, 3e4, 1e5, 3e5, 1e6]
            parts = []
            for tv in sel_pts:
                if tm[0] <= tv <= tm[-1]:
                    j = np.argmin(np.abs(np.log(tm) - np.log(tv)))
                    parts.append(f"a({tv:.0e})={al[j]:.3f}")
            print("    alpha_eff: " + "  ".join(parts))


if __name__ == "__main__":
    win = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    summarize(sys.argv[1], win)
