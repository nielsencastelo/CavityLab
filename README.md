# CavityLab

**Reproducible, energy-audited simulation of electromagnetic and quantum cavities.**

CavityLab is a scientific Python library developed alongside the **ZPE-AI** research
program. ZPE-AI is a falsifiable computational investigation of energy transfer in
dynamically modulated cavities: dynamical Casimir effect, parametric resonance,
and later toroidal geometries, scientific ML and inverse design.

> **Scientific note.** CavityLab does not assume net energy generation from the
> vacuum. It provides the infrastructure to simulate, reproduce, **audit** and
> **falsify** such hypotheses. No project experiment has produced net energy: all
> field-energy growth observed so far is accounted for by pump work to round-off.
> See [Scientific claims and limitations](docs/scientific_claims_and_limitations.md).

## Why

In time-modulated systems the field energy is not conserved. Telling physical pump
work apart from numerical error is the central problem. CavityLab makes that
distinction a default:

- **Discrete-exact energy ledger.** The native FDTD solver books stored energy, pump
  work, losses, port output and source work consistently with the discretization.
  `E_net = E_out + E_loss + ΔE_stored − E_in` closes to floating-point round-off
  (≈1e-15 to 1e-12). See [theory](docs/theory/discrete_energy_identity.md).
- **Naive ledger alongside it.** Textbook estimators are booked too, to show how
  they create false "gains". EXP-0005 finds apparent net generation at most points
  of a parameter map, and it vanishes at first order under refinement.
- **Benchmarks against theory.** Yee dispersion, d'Alembert pulses, Hill/Floquet
  parametric resonance, and an exact Gaussian quantum model of the DCE.
- **Reproducibility.** Every experiment writes a manifest with the git commit,
  environment and configuration hash, plus metrics and figures.

## Quickstart

```bash
pip install -e ".[dev]"
python -m pytest
python experiments/EXP-0004_parametric_modulation/run.py
```

**GPU (optional).** Parameter sweeps use the batched solver `FDTD1DBatch`, which has
numpy, torch (CUDA) and cupy backends. To use an NVIDIA GPU, create an isolated
environment with `pip install numpy scipy matplotlib pytest` and
`pip install torch --index-url https://download.pytorch.org/whl/cu126`, then
`pip install -e .`. Experiments with `backend="auto"` pick the GPU automatically.
Audits stay in float64 (see EXP-0003b).

```python
import math
from cavitylab import FDTD1D, Modulation, EnergyAudit
from cavitylab.benchmarks.analytic1d import discrete_mode_fields

L = 0.1                                   # 10 cm PEC cavity
w1 = math.pi * 299792458.0 / L            # fundamental mode
sim = FDTD1D(L, n_cells=200, courant=0.5, sigma=0.02,
             modulation=Modulation(depth=0.04, angular_frequency=2 * w1))
E, H, _ = discrete_mode_fields(sim, 1)
sim.set_fields(E, H)
result = sim.run(50_000)

print(EnergyAudit.from_fdtd(result).report())              # exact ledger: PASS, ~1e-14
print(EnergyAudit.from_fdtd(result, naive=True).report())  # naive ledger: spurious residual
```

## Experiments

| ID | Topic | Status | Key result |
|---|---|---|---|
| [EXP-0001](experiments/EXP-0001_analytic_modes/conclusion.md) | Analytic 1D modes | ✅ | ω_FDTD equals Yee dispersion to 3e-12 |
| [EXP-0002](experiments/EXP-0002_fdtd_convergence/conclusion.md) | FDTD convergence | ✅ | Order 2.00; magic step exact |
| [EXP-0003](experiments/EXP-0003_static_energy_balance/conclusion.md) | Static energy balance | ✅ | Exact ≤ 2e-12; naive false gain +1.6e-4 in a passive cavity |
| [EXP-0003b](experiments/EXP-0003b_precision_and_acceleration/conclusion.md) | Ledger precision + GPU throughput | ✅ | float64 floor 1e-13 (CPU = GPU); float32 unfit for audits; GPU ~6× faster for large batches |
| [EXP-0004](experiments/EXP-0004_parametric_modulation/conclusion.md) | Parametric modulation | ✅ | Floquet from FDTD equals Hill theory; δ_th = 2/Q; multimode stress test |
| [EXP-0004c](experiments/EXP-0004c_multimode_converged/conclusion.md) | Converged multimode vs independent solver | ✅ | FDTD → 4th-order MOL at order 1.995 |
| [EXP-0005](experiments/EXP-0005_parametric_sweep/conclusion.md) | (Ω/ω, δ) resonance maps + dataset (GPU) | ✅ | 100% stability agreement; naive "gain" at 92% of points; ν = 1 tongue absent for ε(t) |
| [EXP-0005b](experiments/EXP-0005b_optimizer_exploit/conclusion.md) | Optimizer exploiting a naive ledger | ✅ | Apparent +154% "gain" on coarse grids; 5/5 candidates rejected by the audit |
| [EXP-0006](experiments/EXP-0006_quantum_dce_gaussian/conclusion.md) | Quantum DCE (Gaussian) | ✅ part 1 | Ergotropy / pump work = 1 − δ_th/δ; classical FDTD predicts vacuum photons |

The roadmap (EXP-0007 to EXP-0013: PINNs, neural operators, 2D/3D, toroids, inverse
design, quantum work audit) is in [docs/context/project_context.md](docs/context/project_context.md).
The plan for the first preprint is in
[docs/research-plan/preprint-01-action-plan.md](docs/research-plan/preprint-01-action-plan.md).

## Repository layout

```
src/cavitylab/     library: core, geometry, solvers, energy, benchmarks, quantum, experiments
tests/             pytest: energy invariants, analytical benchmarks, quantum model
experiments/       EXP-XXXX_*/run.py → results/{manifest,metrics}.json, figures/, conclusion.md
docs/              context (incl. Portuguese source documents), theory, literature, research plan
papers/            preprint outlines and (later) LaTeX sources
```

## Citation and license

See [CITATION.cff](CITATION.cff). MIT License.
