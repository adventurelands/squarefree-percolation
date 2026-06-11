"""Stochastic process simulation per SPEC.md (blind measurement task).

Process: L binary cells, periodic boundary.
  s_i(t+1) = Bernoulli(q_i) if (s_i(t) or s_{i+1}(t)) else 0, synchronous.
Construction A: q_i = p if (OFFSET+i) squarefree else c*p, c=0.6.
Construction B: same multiset of q values, uniformly shuffled.
Initial condition: all ones. Observable rho(t) at log-spaced times.
"""
import json
import math
import os
import sys
import time
import numpy as np
from multiprocessing import Pool

C = 0.6
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")


def squarefree_mask(offset, L):
    """Boolean mask: True where OFFSET+i is squarefree, i=0..L-1."""
    hi = offset + L
    lim = math.isqrt(hi - 1)
    sieve = np.ones(lim + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, math.isqrt(lim) + 1):
        if sieve[p]:
            sieve[p * p:: p] = False
    primes = np.nonzero(sieve)[0]
    mask = np.ones(L, dtype=bool)
    for p in primes:
        p2 = int(p) * int(p)
        start = (-offset) % p2
        if start < L:
            mask[start::p2] = False
    return mask


def make_q(construction, p, offset, L, shuffle_seed):
    sf = squarefree_mask(offset, L)
    q = np.where(sf, p, C * p).astype(np.float32)
    if construction == "B":
        rng = np.random.default_rng(shuffle_seed)
        rng.shuffle(q)
    return q


def record_times(T, pts_per_decade=120):
    n = int(pts_per_decade * math.log10(T)) + 1
    ts = np.unique(np.round(np.logspace(0, math.log10(T), n)).astype(np.int64))
    return ts[ts >= 1]


def run_one(task):
    """task: dict with construction, p, offset, L, T, sim_seed, shuffle_seed"""
    cons = task["construction"]
    p = task["p"]
    offset = task["offset"]
    L = task["L"]
    T = task["T"]
    q = make_q(cons, np.float32(p), offset, L, task["shuffle_seed"])
    rng = np.random.default_rng(task["sim_seed"])
    s = np.ones(L, dtype=bool)
    ts = record_times(T)
    rhos = np.zeros(len(ts), dtype=np.float64)
    k = 0
    t0 = time.time()
    dead = False
    for t in range(1, T + 1):
        if not dead:
            active = s | np.roll(s, -1)
            s = active & (rng.random(L, dtype=np.float32) < q)
        if k < len(ts) and t == ts[k]:
            if not dead:
                m = int(s.sum())
                rhos[k] = m / L
                if m == 0:
                    dead = True
            k += 1
        if dead and k < len(ts):
            # fast-forward: remaining records are zero
            rhos[k:] = 0.0
            k = len(ts)
            break
    out = dict(task)
    out["t"] = ts.tolist()
    out["rho"] = rhos.tolist()
    out["walltime"] = time.time() - t0
    return out


def main():
    jobfile = sys.argv[1]
    outfile = sys.argv[2]
    with open(jobfile) as f:
        tasks = json.load(f)
    os.makedirs(RESULTS, exist_ok=True)
    if outfile.endswith(".jsonl"):
        # incremental: skip tasks already done, append as each finishes
        done = set()
        if os.path.exists(outfile):
            with open(outfile) as f:
                for line in f:
                    r = json.loads(line)
                    done.add((r["construction"], round(r["p"], 6),
                              r["offset"], r["sim_seed"]))
        todo = [t for t in tasks
                if (t["construction"], round(t["p"], 6), t["offset"],
                    t["sim_seed"]) not in done]
        print(f"{len(done)} done, {len(todo)} to run", flush=True)
        with Pool(2) as pool:
            for r in pool.imap_unordered(run_one, todo):
                with open(outfile, "a") as f:
                    f.write(json.dumps(r) + "\n")
                print(f"  {r['construction']} p={r['p']:.4f} "
                      f"off={r['offset']} seed={r['sim_seed']} "
                      f"wall={r['walltime']:.1f}s "
                      f"rho_final={r['rho'][-1]:.3e}", flush=True)
        return
    with Pool(2) as pool:
        results = pool.map(run_one, tasks)
    with open(outfile, "w") as f:
        json.dump(results, f)
    print(f"wrote {outfile} ({len(results)} runs)")
    for r in results:
        print(f"  {r['construction']} p={r['p']:.4f} off={r['offset']} "
              f"seed={r['sim_seed']} T={r['T']} wall={r['walltime']:.1f}s "
              f"rho_final={r['rho'][-1]:.3e}")


if __name__ == "__main__":
    main()
