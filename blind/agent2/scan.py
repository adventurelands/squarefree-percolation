"""Driver: run a set of (construction, p, sample) simulations with a 2-process pool.

Usage:
  python3.11 scan.py --out results.npz --T 30000 \
      --jobs "A:0.80,0.81:200000077:0" "B:0.80:200000077:1"

Each --jobs item is constr:p1,p2,...:offset:tag  (tag = rep & shuffle seed).
Results saved to an .npz with one entry per run, key = "constr_p_offset_tag".
"""

import argparse
import sys
import time
from multiprocessing import Pool

import numpy as np

from sim import run_sample, rec_times

L = 131072


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--T", type=int, required=True)
    ap.add_argument("--jobs", nargs="+", required=True)
    args = ap.parse_args()

    tasks = []
    for spec in args.jobs:
        constr, ps, offset, tag = spec.split(":")
        offset, tag = int(offset), int(tag)
        for p in ps.split(","):
            p = float(p)
            tasks.append((constr, offset, L, p, args.T, tag, tag))

    print(f"{len(tasks)} tasks, T={args.T}", flush=True)
    rec = rec_times(args.T)
    out = {"rec": rec}
    t0 = time.time()
    with Pool(2) as pool:
        for k, res in enumerate(pool.imap_unordered(run_sample, tasks)):
            constr, offset, p, rep, sseed, rec_r, rho = res
            key = f"{constr}_{p:.4f}_{offset}_{rep}"
            out[key] = rho
            print(f"[{k+1}/{len(tasks)}] {key} done, elapsed {time.time()-t0:.0f}s, "
                  f"final rho={rho[-1]:.3e}", flush=True)
            np.savez(args.out, **out)  # checkpoint as we go
    print("all done", flush=True)


if __name__ == "__main__":
    main()
