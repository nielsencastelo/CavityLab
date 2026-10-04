# Architecture

CavityLab keeps a **physics-neutral core**: ZPE/DCE questions live in experiments and,
later, in an optional `zpe` profile. Heavy backends (Meep, FEniCSx, QuTiP, PyTorch)
enter as optional adapters, so they are never core dependencies.

## Package layout

| Path | Status | Responsibility |
|---|---|---|
| `src/cavitylab/core/` | v0.1 | SI constants, run manifests (git commit, environment, config hash), JSON I/O |
| `src/cavitylab/geometry/` | v0.1 | `Cavity1D` (analytic modes, Q). Rectangular, annular and toroidal are planned |
| `src/cavitylab/solvers/` | v0.1 | `FDTD1D`: native Yee solver with a discrete-exact energy ledger, ε(x,t), σ(x), ports, soft sources. `mol1d`: independent reference (4th-order staggered FD + adaptive DOP853) |
| `src/cavitylab/energy/` | v0.1 | `EnergyAudit`: E_in / E_out / E_loss / ΔE, E_net, PASS/FAIL, reports |
| `src/cavitylab/benchmarks/` | v0.1 | Analytic 1D references (B-001..003), Hill ODE and Floquet theory (B-004), FDTD monodromy (`parametric.py`) |
| `src/cavitylab/experiments/` | v0.1 | `ExperimentContext`: manifest + metrics + figures per run |
| `src/cavitylab/quantum/` | v0.1 (Gaussian) | Exact Gaussian moments of a parametric mode with Lindblad loss: photons, energy ledger, ergotropy, passive energy, entropy; photons from a classical transfer matrix. QuTiP adapter planned for nonlinear models |
| `ml/` | planned v0.4–0.5 | PINN, energy-loss ablation, neural operators, surrogates |
| `optimization/` | planned v0.8 | Auditable objectives `J = E_out - E_in - E_loss - penalties`, revalidation loop |
| `datasets/` | started (EXP-0005 CSV) | Schema-conformant audited datasets |

## Data flow

```
config (dict) ──► solver.run() ──► FDTDResult (time series + exact & naive ledgers)
                                        │
                                        ├─► EnergyAudit ──► PASS/FAIL, E_net, report
                                        ├─► benchmark comparison (analytic / ODE / Floquet)
                                        └─► ExperimentContext ──► manifest.json, metrics.json, figures/
```

## Conventions

- SI units everywhere. In 1D, energies are per unit transverse area [J/m²].
- `float64` by default. Precision is recorded in each manifest.
- `E_net = E_out + E_loss + ΔE_stored - E_in`. See
  [theory/discrete_energy_identity.md](theory/discrete_energy_identity.md).
- Every experiment is a script `experiments/EXP-XXXX_*/run.py`. Its outputs are
  committed: small JSON/CSV files and PNG figures. Raw arrays are git-ignored.
