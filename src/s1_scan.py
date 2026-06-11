#!/usr/bin/env python3.11
"""Coarse p-scan driver for S1, multiprocessing across 8 cores."""
import json, sys, itertools
from multiprocessing import Pool
import importlib.util
spec = importlib.util.spec_from_file_location("s1", "s1_dp.py")
s1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s1)

def job(a):
    model, p, L, T, off, seed, c = a
    return {"model": model, "p": p, "offset": off, "seed": seed, "c": c,
            "rho": s1.run(model, p, L, T, off, seed, c)}

if __name__ == "__main__":
    out = sys.argv[1]
    L, T, c = 131072, 20000, 0.6
    jobs = []
    for p in [0.705489]:
        jobs.append(("clean", p, L, T, 1000003, 11, c))
    for model in ["squarefree", "shuffled"]:
        for p in [0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 0.90]:
            for off, seed in [(1000003, 11), (50000017, 23)]:
                jobs.append((model, p, L, T, off, seed, c))
    with Pool(7) as pool:
        res = pool.map(job, jobs)
    with open(out, "w") as f:
        json.dump(res, f)
    print("done", len(res))
