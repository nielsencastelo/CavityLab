# EXP-0001: Analytical modes of a 1D PEC cavity

**Run:** `python experiments/EXP-0001_analytic_modes/run.py`
**Outputs:** `results/metrics.json`, `results/manifest.json`, `figures/exp0001_mode_frequencies.png`

| Field | Value |
|---|---|
| Hypothesis | The native Yee solver reproduces the analytical modes `ω_n = nπc/L`. Its frequency error equals the Yee dispersion relation, and the discrete energy is conserved. |
| Solver | `cavitylab.native_fdtd1d`, float64, S = 0.5, N ∈ {50, 100, 200, 400}, modes 1–10 |
| Geometry / BC | L = 0.1 m, vacuum, PEC walls |
| Initial condition | Exact discrete eigenmode (amplitude 1 V/m) |
| Energy in / out / loss | 0 / 0 / 0 |
| Energy residual | max drift of discrete W: **3.2e-15** (relative) |
| Convergence | Frequency error vs analytic follows `-(kΔx)²(1-S²)/24` (e.g. −1.928e-6 measured vs −1.928e-6 predicted, N = 400, n = 1) |
| Independent validation | Analytic dispersion relation; frequency measured from the time series agrees with Yee dispersion to **2.7e-12** |
| Result | **PASS.** Both gates (`frequency_vs_yee < 1e-10`, `energy_drift < 1e-12`) met. |

**Interpretation.** The solver is exactly the Yee scheme. The only frequency error is
the known numerical dispersion, which is second order in Δx. The discrete stored
energy differs from the continuum value `εA²L/4` by `≈ 2×` the frequency error.
That is expected: it is the discrete norm of a discrete mode.

**Limitations.** Homogeneous cavity only.

**Evidence for H0 / H1.** Not applicable (a passive static cavity).

**Next action.** EXP-0002 (convergence for non-modal initial data).
