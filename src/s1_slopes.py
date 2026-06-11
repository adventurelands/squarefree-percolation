#!/usr/bin/env python3.11
"""Local-slope diagnostics for deep runs: alpha_eff(t) over sliding log windows,
averaged across disorder samples, for each model/p. Discriminates:
  stable power law: alpha_eff(t) -> const > 0
  activated (IRFP): alpha_eff(t) -> 0, roughly linear in 1/ln(t)
"""
import sys, glob, re
import numpy as np
from collections import defaultdict

outdir = sys.argv[1]
groups = defaultdict(list)
for f in glob.glob(outdir + "/*.dat"):
    m = re.match(r".*/(sf|sh|clean)_?([0-9.]*)_?(\d*)\.dat", f)
    tag = m.group(1); p = m.group(2) or "ref"
    d = np.loadtxt(f)
    if d.ndim == 1: d = d.reshape(1, -1)
    groups[(tag, p)].append(d)

for (tag, p), runs in sorted(groups.items()):
    # align on common times
    tsets = [set(r[:, 0].astype(int)) for r in runs]
    common = sorted(set.intersection(*tsets))
    t = np.array(common, float)
    rho = np.mean([[dict(zip(r[:, 0].astype(int), r[:, 1]))[c] for c in common] for r in runs], axis=0)
    ok = rho > 0
    t, rho = t[ok], np.array(rho)[ok]
    lt, lr = np.log(t), np.log(rho)
    print(f"\n== {tag} p={p}  ({len(runs)} samples)  rho(T)={rho[-1]:.4e}")
    print(f"{'t_mid':>10} {'alpha_eff':>10} {'1/ln(t)':>9}")
    W = 0.6  # half-window in ln t
    for center in np.arange(lt[0] + 1.5, lt[-1] - 0.3, 0.55):
        m = (lt > center - W) & (lt < center + W)
        if m.sum() < 4: continue
        a = -np.polyfit(lt[m], lr[m], 1)[0]
        tm = np.exp(center)
        print(f"{tm:>10.0f} {a:>10.4f} {1/np.log(tm):>9.4f}")
