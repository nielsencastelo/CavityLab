# EXP-0006 (part 1): Minimal quantum DCE: Gaussian model, work audit, classical bridge

**Run:** `python experiments/EXP-0006_quantum_dce_gaussian/run.py`
**Outputs:** `results/metrics.json`, `figures/exp0006_quantum_dce.png`

**Hypothesis.**

- A modulated cavity mode creates photons from vacuum at the parametric rate, and
  the pump supplies all of their energy.
- The transfer matrix of the classical FDTD solver predicts the quantum photon number
  from vacuum. This is exact for linear lossless media, because a quadratic
  Hamiltonian maps covariances as V → MVMᵀ.
- Ergotropy never exceeds pump work.

**Model.** `H = p²/(2f) + q²/2` with `f = 1 + δ sin(2τ)`, exact Gaussian dynamics,
Lindblad amplitude damping (κ = 1/Q) into a bath with occupation n_th. Units:
ħ = 1, τ = ω₀t.

## Results (δ = 0.02)

**Lossless, from vacuum.**

- N(τ = 191π ≈ 600) = 100.45, against the RWA prediction `sinh²(δτ/4)` = 100.40 (0.05%).
- Pump work = 100.45 ħω₀ = N·ħω₀, since the zero-point energy is unchanged.
- Ergotropy/N = 0.99999998 (pure state).
- Ledger residual: 2.5e-9 (ODE tolerance).

**Classical FDTD → quantum bridge** (200 pump periods, N_quantum = 133.4).

| N (cells) | 50 | 100 | 200 | 400 |
|---|---|---|---|---|
| \|N_FDTD/N_q − 1\| | 5.3e-4 | 3.8e-4 | 1.1e-4 | 2.9e-5 |

The convergence order is about 1.9, and det M = 1 to 2e-15 (symplectic).

**Lossy work audit at T = 0** (τ = 637π ≈ 2001, a stroboscopic f = 1 instant; δ_th = 2/Q).

| Q (δ/δ_th) | photons | ergotropy | passive energy | pump work | to bath | ergotropy / pump work |
|---|---|---|---|---|---|---|
| 400 (4) | 1.10e6 | 1.10e6 | 469 | 1.47e6 | 3.67e5 | **0.750** |
| 200 (2) | 1.11e4 | 1.10e4 | 60 | 2.22e4 | 1.11e4 | **0.497** |
| 100 (1, threshold) | 4.88 | 3.76 | 1.12 | 52.5 | 47.6 | **0.072** |

**Emergent law.** Above threshold, ergotropy / pump work → **η = 1 − δ_th/δ** (0.75,
0.50 and 0 for δ/δ_th = 4, 2, 1):

- The net energy growth rate is `δ/2 − κ` and the pump power per unit energy is
  `δ/2`, so the stored fraction is `1 − 2κ/δ`. The state is nearly pure (passive
  energy is much smaller than the photon number), so this stored energy is ergotropy.
- At threshold the efficiency vanishes: all pump work goes to the bath.

With n_th = 1, every quantity scales by (2n_th + 1) = 3, and the efficiency ratios are
unchanged (linear Gaussian dynamics). The ledger closes to below 1e-8 relative in all
cases.

**Result.** **PASS**, all gates.

## Interpretation

1. DCE photons are paid for by pump work, quantum by quantum. A T = 0 bath takes no
   energy from the vacuum state: zero-point energy is not a flow.
2. Only part of the pump work becomes extractable work:
   - The fraction is η = 1 − δ_th/δ: 75% at four times the threshold, 7% at threshold.
   - The rest is dissipated or locked in passive (entropic) energy.
   - This is the first quantitative quantum-work-audit result (RQ8, EXP-0013 preview).
3. For linear media a classical FDTD solver plus vacuum-covariance propagation gives
   exact quantum photon numbers. This justifies the classical-solver-first strategy
   for DCE benchmarks. Loss needs the fluctuation–dissipation noise term, which is
   not present in the classical FDTD.

## Limitations

- Single mode, with a non-RWA Hamiltonian but RWA (Lindblad) damping.
- The modulation is prescribed. There is no Kerr or pump depletion, so growth above
  threshold is unbounded.
- No comparison yet with published circuit-DCE parameters.

## Next action

- **Part 2:** reproduce a published DCE model, using the Wilson et al. 2011 SQUID
  parameters (two-mode squeezing) or the PRA 96, 033851 optimal-control case.
- Extend the bridge to multimode FDTD transfer matrices with Wigner-sampled seeding.
