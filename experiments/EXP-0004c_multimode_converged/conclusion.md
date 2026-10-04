# EXP-0004c: Converged multimode parametric benchmark against an independent solver

**Run:** `python experiments/EXP-0004c_multimode_converged/run.py` (about 10 min CPU, mostly the MOL reference)
**Outputs:** `results/metrics.json`, `figures/exp0004c_multimode_converged.png`

| Field | Value |
|---|---|
| Hypothesis | With a non-equidistant spectrum, multimode parametric pumping has a resolvable continuum limit. The native FDTD converges to it at 2nd order, as confirmed by an independent solver. |
| Cavity | Smooth dielectric step, ε_s from 1 to 3 (tanh, center 0.6L, width 0.05L). Spectrum ω_k/ω₁ = 1, 2.07, 3.01, 4.11, 5.04, … |
| Modulation | Smooth Gaussian profile (center 0.3L, width 0.05L), δ = 0.2, Ω = 2ω₁, phase π/4. Lossless. Initial E = sin(πx/L), H = 0 (multimode). |
| Solver A | Native Yee/leapfrog FDTD, N = 100…1600, S = 0.5, exact ledger |
| Solver B (independent) | Method of lines: 4th-order staggered FD with PEC images + adaptive DOP853 (rtol 1e-11). N = 200/400/800. Verified 4th order (4.00) on an analytic mode. |
| ω₁ | From the lowest non-null eigenvalue of the fine-grid (N = 4000) operator |
| Result W(20 periods)/W(0) | MOL reference **1.248046**. FDTD 1.27394 (N = 100) → 1.24809 (N = 1600) |
| FDTD vs MOL | Relative error 2.1e-2, 5.1e-3, 1.3e-3, 3.2e-4, **8.1e-5**. **Order 1.995.** |
| MOL self-convergence | 3.3e-6 (N = 200) and 1.9e-7 (N = 400), both vs N = 800, which is more than 400× below the finest FDTD error |
| Exact ledger | PASS in every FDTD run (< 1e-10) |
| Result | **PASS**, all four gates |

## Interpretation

- This closes the "independent validation" requirement of the project rules for the
  1D parametric problem. A different spatial discretization and a different time
  integrator agree with the native solver at its formal order.
- Together with EXP-0004 part B it isolates the reason the equidistant case fails
  to converge: the equidistant spectrum, not the solver.
- The growth is modest here (×1.25 in 20 periods) because the spectrum is detuned:
  ω₃ − ω₁ = 2.01ω₁ is close to, but not on, the pump frequency.

## Bug found and fixed during this experiment

The first run used an ω₁ that was about 10⁴ times too small. The null mode of the
rank-deficient operator `G diag(1/ε) Gᵀ` (H = const) has a numerical eigenvalue of
about 1e-8, which passed a `> 1e-9` filter. The fix is to skip index 0 explicitly.
The symptom was a run that never finished (3.5 h of CPU before it was stopped),
not a wrong result, because the simulated time was enormous.

## Limitations

- 1D only.
- The MOL reference has no energy ledger; it is used only as the field and energy
  reference.

## Next action

- Use this as the converged-multimode figure of Preprint 01 (Fig. 5, right panel).
