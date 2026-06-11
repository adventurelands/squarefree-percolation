"""Final analysis: disorder-averaged rho(t), alpha_eff(t) drift, power-law tests."""
import json
import sys
import numpy as np
from analyze import load, alpha_eff


def collect(files, cons, p):
    runs = []
    for f in files:
        for r in load(f):
            if r["construction"] == cons and abs(r["p"] - p) < 1e-9:
                runs.append(r)
    t = np.array(runs[0]["t"], dtype=float)
    rhos = np.array([r["rho"] for r in runs])
    return t, rhos


def report(files, cons, p):
    t, rhos = collect(files, cons, p)
    n = len(rhos)
    rho = rhos.mean(axis=0)
    print(f"--- {cons} p={p} ({n} disorder samples)")
    # alpha_eff at log-spaced anchor times, with sample-to-sample spread
    anchors = [1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]
    anchors = [a for a in anchors if a <= t[-1] * 1.01]
    ta_m, aa_m = alpha_eff(t, rho)
    per_sample = []
    for r in rhos:
        ta, aa = alpha_eff(t, r)
        per_sample.append((ta, aa))
    print(f"{'t':>9} {'rho(t)':>10} {'a_eff(avg)':>10} {'a_eff mean±sem over samples':>28}")
    rows = []
    for a in anchors:
        i = np.argmin(np.abs(np.log(t) - np.log(a)))
        j = np.argmin(np.abs(np.log(ta_m) - np.log(a)))
        vals = []
        for ta, aa in per_sample:
            if len(ta) and ta[-1] >= a * 0.9:
                k = np.argmin(np.abs(np.log(ta) - np.log(a)))
                vals.append(aa[k])
        vals = np.array(vals)
        sem = vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else np.nan
        print(f"{t[i]:9.0f} {rho[i]:10.3e} {aa_m[j]:10.3f} "
              f"{vals.mean():14.3f} +/- {sem:.3f}")
        rows.append((t[i], rho[i], aa_m[j], vals.mean(), sem))
    # global power-law fit over t in [1e3, T]
    good = (t >= 1e3) & (rho > 0)
    A = np.polyfit(np.log(t[good]), np.log(rho[good]), 1)
    resid = np.log(rho[good]) - np.polyval(A, np.log(t[good]))
    print(f"global fit t>=1e3: alpha={-A[0]:.4f}, rms log-resid={resid.std():.4f}")
    # drift test: fit alpha_eff vs ln t linearly over ta>=1e3
    sel = ta_m >= 1e3
    if sel.sum() > 5:
        D = np.polyfit(np.log(ta_m[sel]), aa_m[sel], 1)
        print(f"alpha_eff drift slope d(a_eff)/d(ln t) over t>=1e3: {D[0]:+.4f}")
    # ln-power test: fit ln rho vs ln ln t over t>=1e3
    sel2 = good & (t >= 1e3)
    B = np.polyfit(np.log(np.log(t[sel2])), np.log(rho[sel2]), 1)
    resid2 = np.log(rho[sel2]) - np.polyval(B, np.log(np.log(t[sel2])))
    print(f"log-power fit rho~(ln t)^-x, t>=1e3: x={-B[0]:.3f}, rms log-resid={resid2.std():.4f}")
    return rows


if __name__ == "__main__":
    files = sys.argv[3:]
    report(files, sys.argv[1], float(sys.argv[2]))
