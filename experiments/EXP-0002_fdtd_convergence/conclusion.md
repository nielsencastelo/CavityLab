# EXP-0002: FDTD convergence

**Run:** `python experiments/EXP-0002_fdtd_convergence/run.py`
**Outputs:** `results/metrics.json`, `figures/exp0002_convergence.png`

| Field | Value |
|---|---|
| Hypothesis | The solver converges at second order in space and time. |
| Setup | Gaussian pulse (center 0.4L, width 0.05L) in a PEC cavity, evaluated at t = 1.37 round trips. Reference: exact d'Alembert/image solution, including the consistent initial H at t = −Δt/2. |
| Resolutions | N = 50…1600; S ∈ {0.25, 0.5, 0.9, 1.0} |
| Observed order | **1.998 (S = 0.25), 1.999 (S = 0.5), 2.001 (S = 0.9)** |
| Magic step S = 1 | Error at round-off (≤ 3.9e-14). The 1D Yee scheme is exact at S = 1. |
| Energy residual | Discrete W drift at round-off in all runs |
| Result | **PASS** (`gate_second_order`) |

**Interpretation.** The spatial/temporal discretization behaves as theory predicts.
Errors at fixed resolution grow as S decreases, because dispersion vanishes as S → 1
in 1D. Experiments therefore fix S = 0.5, a representative non-magic value, so that
dispersion effects are visible and honestly reported.

**Limitations.** 1D only. The magic step has no analogue in 2D/3D.

**Next action.** EXP-0003 (global energy balance with losses, sources and ports).
