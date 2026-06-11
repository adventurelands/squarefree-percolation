# Blind measurement report — agent2

(Provenance: agent2's tool access was revoked before it could write this file itself; this is its
final report text, preserved verbatim from its completion message by the orchestrating session.
All supporting code and raw data written by agent2 itself are in this directory: sim.py, scan.py,
analyze.py, final_analysis.py, phase1-4.npz, phase*.log.)

**Setup measured/used:** L = 131072, synchronous update s_i(t+1)=1 w.p. q_i if s_i(t) or s_{i+1}(t),
all-ones start, c = 0.6. Squarefree sieve verified against brute force (squarefree fraction 0.6078,
weak fraction 0.3922). OFFSETs: 200000077, 700000559, 481729333, 314159311. B shuffle seeds
11/22/33; all RNG via SeedSequence([0x5EEDFAB1E, ...]). Campaign: coarse scan p=0.80-0.88 step 0.01
(T=3e4) -> refine step 0.002 (T=1e5, n=3) -> fine step 0.0005 (T=3e5, n=3) -> long runs T=1e6
(A: n=4; B: n=3 at each of two p). Max 2 worker processes. alpha_eff from LSQ slope in
±0.4-decade log windows; errors = jackknife over disorder samples.

## Q3a: p* estimates

| construction | subcritical bracket | supercritical bracket | p* |
|---|---|---|---|
| A (squarefree) | 0.8205 (rho(3e5)=0.059 falling; alpha_eff rising 0.16->0.37) | 0.8210 (plateau rho=0.132) | **p*_A = 0.8207 ± 0.0003** |
| B (shuffled) | 0.8450 (rho(3e5)=0.126 falling; alpha_eff upturn to 0.275 at t=1e5) | 0.8460 (plateau rho≈0.19-0.20 at t=1e6) | **p*_B = 0.8455 ± 0.0005** |

Delta p* = 0.0248 ± 0.0006 — unambiguous.

## Q1/Q2: long-time behavior at p*

**Construction A at p=0.8207 (n=4, T=1e6).** rho(t): 0.250 (1e3), 0.170 (1e4), 0.114 (1e5), 0.069 (1e6).
alpha_eff: **0.163 ± 0.009 (t=1e3)**, 0.183 ± 0.013 (1e4), 0.177 ± 0.019 (1e5), 0.202 ± 0.079 (3e5),
**0.214 ± 0.094 (last window, data to t=1e6)**.
Global power-law fit [1e3,1e6]: **alpha = 0.180, RMS resid 0.035**; (ln t)^-x fit: x=1.81, RMS 0.059
(clearly worse).
-> **A is consistent with a single stable power law, alpha = 0.18 ± 0.03**; alpha_eff flat within
errors over 3 decades, slight late UPWARD tendency (consistent with ±0.0003 p-tuning), no drift
toward zero slope.

**Construction B near p* (n=3 each).**
- p=0.8450 (just below p*): alpha_eff = **0.104 (t=1e3)**, 0.114 (1e4), 0.108 (3e4), then upturn
  **0.275 (1e5)** — subcritical peel-off.
- p=0.8460 (just above): alpha_eff = **0.108 ± 0.004 (t=1e3)** -> 0.082 ± 0.018 (1e4) ->
  0.080 ± 0.038 (3e5) -> **0.063 ± 0.029 (largest window, data to 1e6)**; rho plateaus ≈0.20 by 1e6.
  p=0.8470: 0.082 -> 0.043 ± 0.055.
- Fits at 0.8460 over [1e3,1e6]: power law alpha=0.096 (RMS 0.029) vs (ln t)^-x x=0.96 (RMS 0.035)
  — comparable, neither clean.
- Strong sample heterogeneity: rho(1e6) spans 0.157-0.298 at p=0.847; at 0.845 one of three samples
  plateaus while two decay (smeared transition).
-> **B shows no stable power law**: alpha_eff drifts systematically downward toward zero on the
critical/supercritical side (0.108 -> 0.063) and flips to sharp upward drift just below p*;
consistent with a logarithmic-type rho ~ (ln t)^(-x), x ≈ 1 (suggestive, not conclusive at T=1e6 —
the robust result is the drift and absence of any alpha plateau).

## Q3b verdict: A and B behave DIFFERENTLY

1. Different p*: 0.8207(3) vs 0.8455(5) despite identical q multisets.
2. Different critical dynamics: A = stable power law alpha ≈ 0.18; B = drifting alpha_eff with no
   plateau, slower decay (alpha_eff(1e3) ≈ 0.10, about half of A's), and large rare-sample scatter.
3. Measured structural cause (from my own q arrays): same 39.2% weak fraction, but A's weak-site
   runs are capped at length 7 (squarefree gaps short), while B's shuffle has an exponential
   run-length tail (runs up to 13 observed) — B contains rare long weak regions A lacks.
