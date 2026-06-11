#!/usr/bin/env python3.11
"""Figures for the arXiv paper, from results/deep1/*.dat."""
import glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8,
                     "lines.linewidth": 1.4, "figure.dpi": 200})
DD = "../data/deep1/"

def avg(pattern):
    runs = []
    for f in sorted(glob.glob(DD + pattern)):
        d = np.loadtxt(f)
        if d.ndim == 1: d = d.reshape(1, -1)
        runs.append(dict(zip(d[:, 0].astype(int), d[:, 1])))
    common = sorted(set.intersection(*[set(r) for r in runs]))
    t = np.array(common, float)
    rho = np.array([[r[c] for c in common] for r in runs])
    return t, rho.mean(0), rho.std(0) / np.sqrt(len(runs)), len(runs)

def alpha_eff(t, rho, W=0.6, step=0.45):
    lt, lr = np.log(t), np.log(rho)
    out = []
    for c in np.arange(lt[0] + 1.2, lt[-1] - 0.2, step):
        m = (lt > c - W) & (lt < c + W)
        if m.sum() >= 4:
            out.append((np.exp(c), -np.polyfit(lt[m], lr[m], 1)[0]))
    return np.array(out).T

curves = {
    "clean ($p_c=0.70549$)": ("clean_*.dat", "0.4", "-"),
    "squarefree ($p=0.8207$)": ("sf_0.8207_*.dat", "tab:blue", "-"),
    "shuffled ($p=0.8465$)": ("sh_0.8465_*.dat", "tab:red", "-"),
    "shuffled ($p=0.846$)": ("sh_0.846_*.dat", "tab:red", ":"),
    "shuffled ($p=0.848$)": ("sh_0.848_*.dat", "tab:red", "--"),
}

# ---- Fig 1: rho(t)
fig, ax = plt.subplots(figsize=(3.5, 2.9))
for label, (pat, color, ls) in curves.items():
    if label.endswith("0.846$)") or "0.848" in label:
        continue
    t, r, e, n = avg(pat)
    ax.loglog(t, r, color=color, ls=ls, label=label + f"  [n={n}]")
tg = np.logspace(1.0, 6, 50)
ax.loglog(tg, 0.62 * tg ** -0.1595, "k--", lw=0.9, label=r"$t^{-0.1595}$ (DP)")
ax.set_xlabel(r"$t$"); ax.set_ylabel(r"$\rho(t)$")
ax.legend(frameon=False, loc="lower left")
fig.tight_layout(); fig.savefig("fig1_rho.pdf")

# ---- Fig 2: alpha_eff(t)
fig, ax = plt.subplots(figsize=(3.5, 2.9))
for label, (pat, color, ls) in curves.items():
    t, r, e, n = avg(pat)
    ta, aa = alpha_eff(t, r)
    ax.semilogx(ta, aa, color=color, ls=ls, label=label)
ax.axhline(0.1595, color="k", ls="--", lw=0.9)
ax.text(1.5e0, 0.168, "DP: 0.1595", fontsize=8)
ax.set_xlabel(r"$t$"); ax.set_ylabel(r"$\alpha_{\rm eff}(t)$")
ax.set_ylim(-0.02, 0.30)
ax.legend(frameon=False, loc="lower center", ncol=2, fontsize=6.5, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(); fig.savefig("fig2_alpha.pdf")

# ---- Fig 3: disorder statistics (run lengths + wandering)
import importlib.util
spec = importlib.util.spec_from_file_location("s1", "../src/s1_dp.py")
s1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s1)
L = 2**20
sf = s1.squarefree_mask(1000003, L); w = ~sf
rng = np.random.default_rng(7)
ws = rng.permutation(w)
def runlens(v):
    d = np.diff(np.concatenate([[0], v.view(np.int8), [0]]))
    return np.nonzero(d == -1)[0] - np.nonzero(d == 1)[0]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(5.2, 2.5))
for v, lab, c in [(w, "squarefree positions", "tab:blue"), (ws, "shuffled", "tab:red")]:
    r = runlens(v)
    ks = np.arange(1, r.max() + 1)
    pk = np.array([(r >= k).mean() for k in ks])
    a1.semilogy(ks, pk, "o-", ms=3, color=c, label=lab)
a1.set_xlabel(r"weak-run length $k$"); a1.set_ylabel(r"$P(\mathrm{run}\geq k)$")
a1.legend(frameon=False, fontsize=7)
x = sf - sf.mean()
ls_, st_ = [], []
for ell in [2**k for k in range(6, 19, 2)]:
    n = len(x) // ell
    ls_.append(ell); st_.append(x[:n*ell].reshape(n, ell).sum(1).std())
xs = sf.astype(float).copy(); rng.shuffle(xs); xs -= xs.mean()
ls2, st2 = [], []
for ell in [2**k for k in range(6, 19, 2)]:
    n = len(xs) // ell
    ls2.append(ell); st2.append(xs[:n*ell].reshape(n, ell).sum(1).std())
a2.loglog(ls_, st_, "o-", ms=3, color="tab:blue", label=r"squarefree, $\omega=0.25$")
a2.loglog(ls2, st2, "o-", ms=3, color="tab:red", label=r"shuffled, $\omega=0.5$")
a2.set_xlabel(r"window $\ell$"); a2.set_ylabel(r"std of weak-site count")
a2.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig("fig3_disorder.pdf")
print("figures written")
