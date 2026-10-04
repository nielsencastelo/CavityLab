# EXP-0003: Global energy balance in static cavities

**Run:** `python experiments/EXP-0003_static_energy_balance/run.py`
**Outputs:** `results/metrics.json`, `figures/exp0003_energy_residuals.png`, `figures/exp0003_naive_residual_scaling.png`

**Hypothesis.** With the discrete-exact ledger, `E_net = E_out + E_loss + ΔW − E_in`
vanishes to round-off in every static scenario. Naive (continuum-formula) ledgers
leave residuals that are pure discretization artifacts.

## Results (N = 400, S = 0.5)

| Scenario | exact max \|E_net\|/scale | naive max \|E_net\|/scale | naive final E_net/scale |
|---|---|---|---|
| A closed, lossless pulse (100 round trips) | 4.2e-15 | 1.5e-4 | −8.1e-5 |
| B free decay, Q = 50 | 6.0e-13 | 1.6e-4 | **+1.6e-4** |
| C driven on resonance, Q = 100 (≈4.8e5 steps) | 1.9e-12 | 2.3e-5 | −2.2e-5 |
| D open cavity: partial mirror + absorber port | 3.5e-15 | 2.7e-4 | −4.0e-6 |

- **Q measured from the decay:** 50.001, against the analytic `ωε/σ = 50`.
- **Naive residual vs N (asymptotic order):** A ≈ 2.1 (time-averaged energy,
  O(Δt²)); C ≈ 0.9 (left-point source and dissipation rules, O(Δt)); D ≈ 1.2 (mixed).
- **Exact residual:** grows only with step count, at the round-off level (~n_steps·ε).

**Result.** **PASS.** All exact ledgers are below tolerance 1e-10, and Q is within 1%.

## Interpretation

- The exact ledger closes in every case. The residual is floating-point round-off,
  which accumulates with the number of steps (4.8e5 steps → 1.9e-12).
- The naive ledger, which uses formulas found in textbooks and typical
  post-processing scripts, leaves 10⁻⁵ to 10⁻⁴ residuals at a standard resolution
  (λ/400 for mode 1). They can be **positive**: in scenario B, a purely passive
  decaying cavity, the naive ledger reports apparent net generation of +1.6e-4 of
  the energy scale.
- These residuals vanish under refinement at the predicted orders. They are
  artifacts, and they are exactly the "candidate anomalies" the falsification
  checklist must catch.

## Limitations

- The absorber is graded electric conductivity, not a PML (some reflection), which
  does not affect the ledger.
- Mirror and absorber parameters are not optimized.

## Evidence for or against H1

There is none for H1. Every apparent gain is a naive-ledger artifact that
disappears under refinement and under the exact ledger.

## Next action

EXP-0004: add a time-varying permittivity and book pump work exactly.
