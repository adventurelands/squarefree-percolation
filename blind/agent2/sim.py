"""Core simulation for the blind stochastic-process measurement task.

Process: L binary cells, periodic boundary.
  s_i(t+1) = 1 with prob q_i  if (s_i(t)=1 or s_{i+1}(t)=1), else 0.
Initial condition: all ones. Observable: rho(t) at log-spaced times.

Construction A: q_i = p if (OFFSET+i) squarefree else c*p, c=0.6.
Construction B: same multiset of q values, uniformly shuffled.
"""

import numpy as np

C = 0.6
SEED_BASE = 0x5EEDFAB1E  # base entropy for this agent's seeding scheme


def squarefree_mask(offset: int, L: int) -> np.ndarray:
    """Boolean mask m[i] = True iff (offset+i) is squarefree, for i in 0..L-1."""
    hi = offset + L  # exclusive upper bound on n
    mask = np.ones(L, dtype=bool)
    limit = int(np.sqrt(hi)) + 2
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(np.sqrt(limit)) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    primes = np.nonzero(sieve)[0]
    for p in primes:
        p2 = int(p) * int(p)
        if p2 >= hi:
            break
        start = ((offset + p2 - 1) // p2) * p2  # first multiple of p2 >= offset
        if start < hi:
            mask[start - offset : L : p2] = False
    return mask


def build_q(constr: str, offset: int, L: int, p: float, shuffle_seed: int = 0) -> np.ndarray:
    mask = squarefree_mask(offset, L)
    q = np.where(mask, p, C * p)
    if constr == "B":
        ss = np.random.SeedSequence([SEED_BASE, 0xB, offset & 0x3FFFFFFF, shuffle_seed])
        rng = np.random.default_rng(ss)
        rng.shuffle(q)
    return q


def rec_times(T: int, n: int = 260) -> np.ndarray:
    ts = np.unique(np.round(np.logspace(0, np.log10(T), n)).astype(np.int64))
    return ts[ts >= 1]


def dyn_rng(constr: str, offset: int, p: float, rep: int) -> np.random.Generator:
    ss = np.random.SeedSequence(
        [SEED_BASE, 0xD, ord(constr[0]), offset & 0x3FFFFFFF, int(round(p * 1_000_000)), rep]
    )
    return np.random.default_rng(ss)


def evolve(q: np.ndarray, T: int, rec: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Run the process for T steps; return rho at the times listed in rec (sorted)."""
    L = q.size
    q = q.astype(np.float32)          # float32 draws: ~2x faster, precision ample
    out = np.zeros(rec.size)
    s = np.ones(L, dtype=bool)        # dense representation
    idx = None                        # sparse representation (sorted active indices)
    sparse = False
    thresh = L // 8                   # switch to index-based updates below this count
    j = 0
    for t in range(1, int(rec[-1]) + 1):
        if not sparse:
            active = s.copy()
            active[:-1] |= s[1:]
            active[-1] |= s[0]
            s = active & (rng.random(L, dtype=np.float32) < q)
            cnt = int(np.count_nonzero(s))
            if cnt == 0:
                break
            if cnt < thresh:
                idx = np.nonzero(s)[0]
                sparse = True
        else:
            left = idx - 1
            if left[0] < 0:
                left[0] = L - 1
            cand = np.unique(np.concatenate((idx, left)))
            idx = cand[rng.random(cand.size, dtype=np.float32) < q[cand]]
            cnt = idx.size
            if cnt == 0:
                break
        while j < rec.size and rec[j] == t:
            out[j] = cnt / L
            j += 1
    # times after extinction stay 0 (out already zero-initialised)
    return out


def run_sample(args):
    """Worker entry: args = (constr, offset, L, p, T, rep, shuffle_seed). Returns rho array."""
    constr, offset, L, p, T, rep, shuffle_seed = args
    q = build_q(constr, offset, L, p, shuffle_seed)
    rec = rec_times(T)
    rng = dyn_rng(constr, offset, p, rep)
    rho = evolve(q, T, rec, rng)
    return constr, offset, p, rep, shuffle_seed, rec, rho
