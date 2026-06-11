#!/usr/bin/env python3.11
"""Summarize S1 scan: for each (model,p) average rho over samples, print late-time
effective exponent alpha_eff over last decade and curvature sign."""
import json, sys
import numpy as np
from collections import defaultdict

res = json.load(open(sys.argv[1]))
groups = defaultdict(list)
for r in res:
    groups[(r["model"], r["p"])].append(r)

print(f"{'model':<11}{'p':<7}{'rho(T)':<10}{'a_eff[1e2-1e3]':<16}{'a_eff[2e3-2e4]':<16}{'curv'}")
for (model, p), rs in sorted(groups.items()):
    # common time grid
    t0 = [pt[0] for pt in rs[0]["rho"]]
    if any([pt[0] for pt in r["rho"]] != t0 for r in rs):
        tmin = min(len(r["rho"]) for r in rs)
        t0 = [pt[0] for pt in rs[0]["rho"][:tmin]]
        rho = np.mean([[pt[1] for pt in r["rho"][:tmin]] for r in rs], axis=0)
    else:
        rho = np.mean([[pt[1] for pt in r["rho"]] for r in rs], axis=0)
    t = np.array(t0, float)
    ok = rho > 0
    t, rho = t[ok], rho[ok]
    def slope(t1, t2):
        m = (t >= t1) & (t <= t2)
        if m.sum() < 3: return np.nan
        return -np.polyfit(np.log(t[m]), np.log(rho[m]), 1)[0]
    a1, a2 = slope(100, 1000), slope(2000, 20000)
    curv = "decel" if a2 < a1 else "accel"
    print(f"{model:<11}{p:<7}{rho[-1]:<10.5f}{a1:<16.4f}{a2:<16.4f}{curv}")
