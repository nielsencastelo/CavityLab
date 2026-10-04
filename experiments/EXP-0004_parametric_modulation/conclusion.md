# EXP-0004: Classical parametric modulation with explicit pump work

**Run:** `python experiments/EXP-0004_parametric_modulation/run.py`
**Outputs:** `results/metrics.json`, `figures/exp0004a_uniform_modulation.png`, `figures/exp0004b_slab_modulation.png`

**Hypothesis.** Field-energy growth under permittivity modulation is fully accounted
for by the discrete pump work `½ DⁿD^{n+1}(1/ε^{n+1} − 1/εⁿ)Δx`, so the exact ledger
gives `E_net = 0`. Uniform modulation reproduces Hill/Floquet theory, including the
threshold `δ_th = 2/Q`.

## Part A: uniform ε(t), Q = 200, Ω = 2ω₁ (benchmark B-004)

| δ (δ/δ_th) | Floquet rate from FDTD monodromy [1/s] | Hill theory [1/s] | exact ledger | naive final E_net/scale |
|---|---|---|---|---|
| 0.005 (0.5) | −2.3547e7 | −2.3546e7 | 2.7e-14 | −2.7e-3 |
| 0.01 (1.0) | −5.7e3 (≈ 0) | −1.2e3 (≈ 0) | 1.6e-14 | −4.8e-3 |
| 0.02 (2.0) | 4.7097e7 | 4.7087e7 | 8.0e-15 | −1.5e-4 |
| 0.04 (4.0) | 1.4127e8 | 1.4127e8 | 6.5e-15 | **+3.7e-4** |

- **Threshold.** The threshold is bracketed exactly at δ = 2/Q. FDTD and theory agree
  to within 0.1% of the loss rate ω/Q (gate `fdtd_floquet_matches_theory`).
- **Convergence.** The FDTD stored energy converges to the Hill ODE at order
  **2.00** (relative error 3.3e-2 at N = 50, down to 1.3e-4 at N = 800). The naive-
  ledger residual decays only at first order (9.8e-4 down to 6.0e-5), and its sign is
  positive at coarse grids.
- **Pump work and loss.** The FDTD values agree with the ODE to within 0.2–1.8% at
  N = 200. That is O(Δx²) error amplified by exponential growth.

**Methodological finding.** A start in one field quadrature with pump phase 0 lies
*exactly* on the decaying Floquet branch (rate −γ − δω/2). The growing branch is then
seeded only by discretization differences, so naive "growth-rate" fits can disagree
with theory by orders of magnitude. CavityLab therefore uses:

- the full monodromy, with both quadratures measured from FDTD, for Floquet
  comparisons;
- a generic pump phase (π/4) for time-series comparisons.

## Part B: localized modulated slab (multimode), lossless

The slab spans [0.05, 0.075] m with δ = 0.3 and Ω = 2ω₁. The field is seeded with
30 random modes (seed 20261004).

- **Ledger.** Stored energy grows ×7.8e3 in 40 periods, and the cumulative pump work
  matches it to **4.4e-15** (exact ledger). Energy cascades to odd modes (1, 3, 5, 7, …).
- **Not a converged benchmark** (documented negative result). W(t) depends on N:
  - At 20 periods, a 3-mode seed gives W/W₀ = 117 (N = 200) and 247 (N = 1600).
  - Cause: with an equidistant spectrum the resonant cascade keeps feeding higher
    modes, so the continuum problem has no resolvable limit. Yee dispersion also
    detunes high modes in a resolution-dependent way.
  - The exact ledger still closes (as it must), while the naive one reports apparent
    gains of **+0.35% (N = 100)** and **+0.13% (N = 200)** of the pump input.

## Interpretation

- No energy appears without pump work. Every apparent gain is a naive-ledger
  artifact.
- **Ledger closure is necessary but not sufficient.** Part B closes to round-off yet
  is unresolved. Convergence and an independent reference (here Hill/Floquet) are
  required in addition.
- Physically, a classical parametric amplifier needs a seed. Growth from "nothing"
  requires the quantum model (EXP-0006).

## Limitations

- The modulation is prescribed. The pump's internal efficiency is not modelled.
- Materials are non-dispersive.
- Part B is not converged (see above).

## Evidence for or against H1

None for H1. Pump work accounts for 100% of field-energy growth to round-off.

## Next action

- EXP-0005: map the (Ω/ω, δ) plane and the naive-ledger artifact.
- EXP-0006: quantum DCE from vacuum.
- Future: a non-equidistant spectrum (dielectric loading) for a *converged*
  multimode benchmark.
