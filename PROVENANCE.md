# Provenance and audit trail

This repository was produced on 2026-06-10/11 in one extended session by Anthropic's Claude
models, Claude Fable 5 and Claude Opus 4.8, operating under human direction, on a consumer machine
(8 core M2, 16 GB). The session followed a fixed protocol set by the human director before
any results existed:

1. **Original work only.** Both engines (`src/dp.c`, `src/s1_dp.py`), all analysis, and all
   harnesses were written from scratch in the session. Prior-art searches were performed
   before pursuing the idea; the closest published work found is the substitution-sequence
   aperiodic contact process literature, which is cited in the paper.
2. **Idea funnel.** Several candidate experiments were screened and killed on the same day
   (a prime-gap driven quantum kicked rotor produced a null result; gcd-matrix spectral
   statistics were judged artifact prone; a cavitation calculation was judged textbook).
   The squarefree directed percolation experiment was the survivor.
3. **Blind replication.** The headline measurements were reproduced by two further
   implementations written by context-isolated agents from `blind/SPEC.md` alone: no
   hypothesis, no expected values, independent RNGs, offsets, and analysis. Their code,
   raw data, and reports are included verbatim under `blind/`. Agent 2 could not write its
   report file due to a sandbox restriction at the end of its run; its report text was
   preserved from its final message and is marked as such. Agent 1 was interrupted by a
   session crash after producing its raw data; the mechanical slope analysis of that data
   (`blind/agent1/ANALYSIS_FROM_RAW.md`) was performed afterward and is labeled accordingly.
4. **Cross validation.** The C and NumPy engines agree on clean and disordered models within
   sampling error. The clean model reproduces the literature critical point and decay
   exponent before any disordered result was trusted.

Known weaknesses of this provenance: the blind agents are instances of the same underlying
model family as the orchestrating assistant, so they are independent in context and code but
not in training; and the time horizon (t = 10^6) bounds every claim, as stated in the paper.
We believe the protocol meaningfully reduces self-confirmation risk; we do not claim it
eliminates it. Independent replication by others is the real test, and the engines plus
spec make that a one-afternoon exercise.
