#!/usr/bin/env python3.11
"""Refined S1: bisect p_c per model at T=1e5, L=2^18, 4 disorder samples."""
import json, sys
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
    L, T, c = 2**18, 100000, 0.6
    samples = [(1000003, 11), (50000017, 23), (360000001, 37), (900000007, 53)]
    jobs = [("clean", 0.705489, L, T, 1000003, 11, c),
            ("clean", 0.705489, L, T, 1000003, 12, c)]
    for p in [0.822, 0.828, 0.834]:
        for off, seed in samples:
            jobs.append(("squarefree", p, L, T, off, seed, c))
    for p in [0.845, 0.851, 0.857]:
        for off, seed in samples:
            jobs.append(("shuffled", p, L, T, off, seed, c))
    with Pool(7) as pool:
        res = pool.map(job, jobs)
    json.dump(res, open(out, "w"))
    print("done", len(res))
