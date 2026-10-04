# EXP-0005: (Ω/ω, δ) resonance maps, naive-ledger artifact map, audited dataset

**Run (GPU):** `.venv-gpu/Scripts/python.exe experiments/EXP-0005_parametric_sweep/run.py`
**Wall time:** 387 s on an RTX A2000 (torch, float64, batched). The previous CPU version was estimated at more than 1 h.
**Outputs:** `results/metrics.json`, `results/dataset.csv` (688 audited simulations), `figures/exp0005_resonance_map.png`

**Hypothesis.**

- (i) The Floquet growth rate measured from FDTD reproduces Hill theory over the
  whole (ν, δ) plane.
- (ii) The exact ledger closes everywhere.
- (iii) Naive-ledger "gains" are discretization artifacts.

## Map A: Q = 100, N = 128, ν ∈ [0.5, 2.6] × δ ∈ [0, 0.3] (688 points)

| Quantity | Value |
|---|---|
| Stable/unstable agreement, FDTD monodromy vs theory | **100%** |
| Median / max \|μ_FDTD − μ_theory\| per ω₁t | 5.7e-10 / 3.8e-4 (tongue edges) |
| Exact ledger, max \|E_net\|/scale | **1.5e-13** |
| Naive ledger | **Apparent gain at 92% of points**; median +0.36%; max +13.7% (ν = 2.5, δ = 0.3); min −5.1% |
| Naive gain under refinement (worst point, N = 64 → 512) | Decays at **order 0.97**; exact stays ~1e-14 |

## Map B: Q = 2000, ν ∈ [0.5, 2.6] (85) × δ ∈ [0, 0.5] (26), 2210 points

- Agreement is **100%**, with a median error of 7e-11.
- The main tongue (ν ≈ 2) has 281 unstable points in both FDTD and theory.

## Finding: the ν = 1 tongue is absent for permittivity modulation (a known exact result, reproduced)

The single-mode equation is a Hill equation with `ω²(t) = ω₀²/(1 + δ sin Ωt)`. It is
not of Mathieu form. Fine scans (lossless, 200 pump periods, GPU) give:

| Region | Hill/Floquet theory | FDTD (N = 128) |
|---|---|---|
| ν ∈ [0.9, 1.2], δ = 0.5 | **No instability** (max μ = −6e-10, the residual damping) | No instability (max μ = 3.7e-7, at discretization level) |
| ν ∈ [0.672, 0.684], δ = 0.3 (3rd tongue) | Unstable, μ = 2.9e-5 at ν = 0.678 | Unstable, μ = 1.9e-5 at ν = 0.678 |

**Literature check: this is known.** Multiplying by `f = 1 + δ sin` turns the
equation into Ince's equation `(1 + a cos 2x) y'' + λ y = 0`, the
"frequency-modulation equation" of Carson and Cambi:

- By Magnus–Winkler coexistence theory, every *even* instability interval
  (Ω ≈ ω₀/k) has exactly zero width for all 0 < δ < 1.
- The *odd* intervals (Ω ≈ 2ω₀/(2k−1)) stay open, with width ∝ δⁿ.

Sources:

- Cambi, Proc. IRE 36, 42 (1948).
- Magnus & Winkler, *Hill's Equation* (1966).
- Figotin, arXiv:2606.02893 (2026), for an LC circuit with modulated capacitance:
  "every even resonance is exactly stable at all modulation amplitudes".

With a sinusoidal *impermittivity* (1/ε) the equation is exactly Mathieu, and all
tongues, including the even ones, open.

**Use in the paper.** This is a validation result, not a new one. The batched FDTD
reproduces an exact coexistence theorem: the ν = 1 tongue is absent and the narrow
ν = 2/3 tongue is present at the right position. A one-line remark may be useful,
because the closure of even momentum gaps for sinusoidal ε(t) is rarely stated in
the photonic-time-crystal literature (none found in this check, which was not
exhaustive).

## Interpretation

1. The batched FDTD with a two-quadrature monodromy is a quantitative Floquet
   analyzer. It reproduces an exact coexistence theorem (closed even tongues, an open
   narrow third tongue) and needs a GPU for dense maps.
2. The left-point pump-work rule is biased toward false gain. It reports net
   generation at 92% of the points at a resolution (256 points per mode-1 wavelength)
   that many practitioners would accept, with a first-order bias. EXP-0005b shows an
   optimizer exploiting exactly this.

## Limitations

- Single mode, uniform modulation.
- The dataset is small (688 rows); it is a seed for EXP-0007 surrogates.

## Evidence for or against H1

None for H1. All 688 ledger runs close the exact ledger (max 1.5e-13). Map B computes
Floquet exponents only.

## Next action

- ✅ Literature check on the missing ν = 1 tongue: known (Ince/Cambi/Magnus–Winkler).
- Optional: compare the measured third-tongue width with Figotin's closed form (§16.6).
- Generate the larger GPU dataset for EXP-0007 (multimode, loss and geometry
  variations).
