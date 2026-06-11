# Stochastic process measurement task

You are given a stochastic process. Your job is to implement it yourself, from this spec only,
and measure what it does. Do not search the literature. Do not guess what answer is "expected."
Report only what your own simulation shows.

## Process definition

State: a row of L binary cells, s_i(t) ∈ {0,1}, periodic boundary, i = 0..L-1.
Each cell i has a fixed probability q_i assigned before the run (see constructions below).

Update (synchronous, all cells at once):
  s_i(t+1) = 1 with probability q_i if s_i(t) = 1 OR s_{i+1 mod L}(t) = 1, else s_i(t+1) = 0.
  (Independent random draw per cell per step.)

Initial condition: all cells 1.
Observable: density rho(t) = mean of s_i(t), recorded at log-spaced times.

## Constructions of q_i (parameter p, constant c = 0.6)

Let n_i = OFFSET + i for a given integer OFFSET.
- Construction A: q_i = p if n_i is squarefree (not divisible by any perfect square > 1), else q_i = c*p.
- Construction B: take the exact multiset of q values from construction A (same OFFSET, same L) and
  randomly permute their positions (uniform shuffle, your own RNG seed).

## Task (single directive)

For each construction, find the critical value p* (the p separating decay of rho to 0 from
saturation at rho > 0; search within p ∈ [0.80, 0.88]), then at your best estimate of p*
characterize the long-time behavior of rho(t):
1. Is the decay consistent with a single stable power law rho ~ t^(-alpha)? If yes, report alpha
   with an error bar.
2. Or does the local log-log slope alpha_eff(t) = -d ln rho / d ln t drift systematically with time
   (e.g., consistent with rho ~ (ln t)^(-x), slope heading toward zero)? If so, report the drift
   direction and the alpha_eff values at t ~ 10^3 and t ~ your largest time.
3. Report p* for each construction, and whether the two constructions behave the same or differently.

Practical guidance: L = 131072 is adequate; run to at least T = 3*10^4 (longer if feasible);
average 3+ independent disorder samples (different OFFSETs around 10^6..10^9 for A; different
shuffle seeds for B); use multiprocessing if helpful. Write your code and run it on this machine
(python3.11 with numpy is available). Report NUMBERS, not vibes.
