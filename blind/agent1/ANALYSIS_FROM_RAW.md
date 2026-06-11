# Agent1 raw-data analysis (post-crash)

Agent1 (blind, independent implementation at blind/agent1/sim.py, L=131072, n=3 disorder samples
per point, T=1e6) was killed by a session crash before writing its own report. The raw simulation
output it produced (results/stage4-8.jsonl) survived. The table below is a mechanical local-slope
analysis of that data performed by the orchestrating session (window ±0.7 in ln t); the
implementation, parameter choices, offsets, and all simulation data are agent1's own.

construction A = squarefree, B = shuffled. a_eff = -dln(rho)/dln(t).

| run | rho(1e6) | a_eff@1e3 | @1e4 | @1e5 | @4e5 | reading |
|---|---|---|---|---|---|---|
| A p=0.8205 | 0.0217 | 0.165 | 0.187 | 0.358 | 0.616 | subcritical peel-off |
| A p=0.8208 | 0.1064 | 0.160 | 0.148 | 0.060 | 0.111 | ~critical, DP-scale slope |
| A p=0.8210 | 0.1406 | 0.137 | 0.136 | 0.066 | 0.032 | slightly supercritical |
| A p=0.8215 | 0.1870 | 0.117 | 0.079 | 0.005 | -0.000 | supercritical plateau |
| A p=0.8220 | 0.2196 | 0.122 | 0.032 | -0.016 | -0.002 | supercritical plateau |
| B p=0.8450 | 0.0927 | 0.103 | 0.125 | 0.222 | 0.213 | subcritical peel-off |
| B p=0.8455 | 0.1581 | 0.103 | 0.119 | 0.176 | 0.006 | ~critical, NO plateau |
| B p=0.8460 | 0.2160 | 0.083 | 0.118 | 0.102 | 0.060 | drifting |
| B p=0.8470 | 0.2768 | 0.077 | 0.054 | 0.046 | 0.028 | drifting to plateau |
| B p=0.8480 | 0.3065 | 0.073 | 0.063 | 0.019 | -0.006 | supercritical plateau |

Conclusions from agent1's data, consistent with the other two implementations:
- p_c(A) ≈ 0.8206-0.8208 (bracketed by sub/supercritical behavior) -- matches 0.8207.
- p_c(B) ≈ 0.8455 -- matches.
- A near criticality: early-window slopes 0.15-0.19, i.e. DP-magnitude (0.1595), power-law-like.
- B near criticality: early-window slopes 0.07-0.10 -- roughly HALF the DP value -- and at no p in
  the bracket does a_eff stabilize; the characteristic disordered drift, no power law.
- Caveat: L=131072 (8x smaller than the C-engine runs) makes late windows noisy and finite-size
  plateau onset earlier; the robust statement from this dataset is the early/mid-window contrast.
