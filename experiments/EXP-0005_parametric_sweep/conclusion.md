# EXP-0005: (Ω/ω, δ) resonance map, naive-ledger artifact map, audited dataset

**Run:** `python experiments/EXP-0005_parametric_sweep/run.py` (about 15 min on a laptop CPU)
**Outputs:** `results/metrics.json`, `results/dataset.csv` (688 audited simulations), `figures/exp0005_resonance_map.png`

| Field | Value |
|---|---|
| Hypothesis | (i) The Floquet growth rate measured from FDTD reproduces Hill theory over the whole (ν, δ) plane. (ii) The exact ledger closes everywhere. (iii) Naive-ledger "gains" are discretization artifacts. |
| Setup | Uniformly modulated lossy cavity, Q = 100, N = 64 (32 cells per mode-1 wavelength), S = 0.5. Grid ν ∈ [0.5, 2.6] (43) × δ ∈ [0, 0.3] (16). Floquet from a two-quadrature monodromy over 40 pump periods. Ledger runs of 30 mode periods with pump phase π/4. |
| Floquet agreement | Stable/unstable classification **100%** (688/688); median \|μ_FDTD − μ_theory\| = 1.4e-9 per ω₁t; max 8.3e-4 (tongue edge, N = 64) |
| Threshold | Tongue tip at δ = 2/Q = 0.02 in both maps |
| Exact ledger | max \|E_net\|/scale = **1.7e-13** over all points |
| Naive ledger | Apparent E_net/scale **> 0 at 92% of points**; median +0.78%; max **+26.5%** at ν = 2.55, δ = 0.3; min −10.7% (left flank of the tongue) |
| Artifact check | At the worst point the naive value halves with every refinement: 0.265 → 0.132 → 0.064 → 0.032 for N = 64 to 512. **Order 1.02.** The exact ledger stays at about 1e-14. |
| Dataset | `dataset.csv` follows the CavityLab schema (ids, config hash, solver, dt, precision, E_in/E_out/E_loss/ΔE, exact and naive E_net, residual, Floquet exponents, audit status) |
| Result | **PASS** (all three gates) |

## Interpretation

1. The native solver with a two-quadrature monodromy is a quantitative Floquet
   analyzer for modulated cavities, even at coarse resolution.
2. The left-point pump-work rule is **biased toward false gain**. It reports net
   energy generation over most of the parameter space, with a first-order bias, at
   a resolution (32 points per wavelength) that many practitioners would accept.
   An automated pipeline or a human reading that ledger would "discover" over-unity
   behavior. EXP-0005b shows that an optimizer does exactly that.
3. The higher-order tongues (ν = 1, 2/3) are absent at Q = 100 and δ ≤ 0.3. Their
   widths scale as δ² and δ³, so their thresholds lie beyond this δ range at this Q.
   Mapping them needs a larger Q or δ (future run).

## Limitations

- Single resolution for the map (N = 64). Floquet errors concentrate at the tongue
  boundaries.
- The dataset is small (688 rows) and single-mode. It is a seed for the surrogate
  and PINN benchmarks (EXP-0007+).

## Evidence for or against H1

None for H1.

## Next action

- EXP-0004c (independent solver, converged multimode case).
- Extend the dataset with geometry, multimode and loss variations for EXP-0007.
