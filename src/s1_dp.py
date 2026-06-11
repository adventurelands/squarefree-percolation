#!/usr/bin/env python3.11
"""S1: 1+1D directed site percolation with arithmetic (squarefree) columnar disorder.

Model: diagonal square lattice. Site i active at t+1 iff (site i or i+1 active at t)
AND rand < p_i. Quenched column probabilities:
  clean      : p_i = p
  squarefree : p_i = p if (i + offset) is squarefree else c*p
  shuffled   : same multiset of good/bad columns as squarefree window, randomly permuted
Start fully active, measure density rho(t) at log-spaced times.
"""
import sys, json
import numpy as np

def squarefree_mask(offset, L):
    n0 = offset
    mask = np.ones(L, dtype=bool)
    top = n0 + L
    d = 2
    while d * d < top:
        q = d * d
        start = ((n0 + q - 1) // q) * q - n0
        mask[start::q] = False
        d += 1
    return mask

def run(model, p, L, T, offset, seed, c):
    rng = np.random.default_rng(seed)
    if model == "clean":
        parr = np.full(L, p)
    else:
        sf = squarefree_mask(offset, L)
        if model == "shuffled":
            rng2 = np.random.default_rng(seed + 777)
            sf = rng2.permutation(sf)
        parr = np.where(sf, p, c * p)
    state = np.ones(L, dtype=bool)
    tsamp = np.unique(np.round(np.logspace(0, np.log10(T), 60)).astype(int))
    out = []
    u = np.empty(L)
    j = 0
    for t in range(1, T + 1):
        rng.random(out=u)
        state = (state | np.roll(state, -1)) & (u < parr)
        if j < len(tsamp) and t == tsamp[j]:
            out.append((t, state.mean()))
            j += 1
            if out[-1][1] == 0.0:
                break
    return out

if __name__ == "__main__":
    model, p, L, T, offset, seed, c = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), float(sys.argv[7])
    res = run(model, p, L, T, offset, seed, c)
    print(json.dumps({"model": model, "p": p, "L": L, "T": T, "offset": offset,
                      "seed": seed, "c": c, "rho": res}))
