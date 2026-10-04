# ZPE-AI / CavityLab: consolidated project context (English)

> This is a condensed English version of the internal documents in
> [`source-documents/`](source-documents/) (Portuguese originals, 26–27 Sep 2026):
> - *ZPE-AI: Contexto Consolidado da Discussão*
> - *ZPE-AI: Revisão aprofundada da literatura científica, v2 (with the Stargate/ZPM extension)*
> - *CavityLab: Proposta de Biblioteca Científica*
>
> Where this summary and the originals disagree, the originals are authoritative
> for intent. This file is authoritative for the conventions adopted in code
> (for example the E_net sign; see [theory](../theory/discrete_energy_identity.md)).

## 1. Project in one paragraph

**ZPE-AI** (*Computational Investigation of Zero-Point Energy for Autonomous Clean
Energy Generation*) asks, in a falsifiable way, whether dynamically modulated
electromagnetic and quantum cavities, eventually including toroidal geometries, can
show **positive net energy transfer over a complete cycle once all inputs, losses
and stored-energy changes are booked**. **CavityLab** is the open-source,
physics-neutral Python library built alongside it. Every experiment adds a validated,
reusable component, and the library stays useful even if every positive hypothesis
is refuted.

The project does **not** assume that net energy can be extracted from the vacuum.
It also does not declare that impossible before testing. Established physics is the
baseline:

- Zero-point energy exists in the formalism. Its existence does not imply a
  reservoir that can be cycled.
- The static Casimir effect is established, but Jaffe (2005) shows it does not
  require the absolute vacuum energy.
- The dynamical Casimir effect (DCE, Wilson et al. 2011) creates real photons, and
  the modulation supplies their energy.
- Moddel (2019) found no reliably demonstrated vacuum-energy extraction proposal.

## 2. Central question and hypotheses

Exact energy bookkeeping with the CavityLab sign convention:

```
E_net = E_out + E_loss + ΔE_stored - E_in        (energy that appeared without a booked input)
H0: E_net ≤ 0      H1 (exploratory): E_net > 0 after full audit
```

Any positive result is first a **candidate computational anomaly**. It must survive
the falsification checklist below before any physical interpretation.

**Extension from the Stargate/ZPM discussion.** Fiction supplies questions, never
evidence. The added axis is quantum thermodynamics: how much of the field energy
is *extractable work*? The relevant quantities are ergotropy, passive energy,
entropy and control cost. New research question:

> **RQ8**: In dynamically modulated cavities, what is the maximum extractable work
> once the initial quantum state, modulation/control energy, dissipation, stored-
> energy change, entropy and passivity constraints are all accounted for?

### Research questions (RQ1–RQ8)

1. **RQ1.** Which time-dependent cavity formulations reproduce published DCE results
   quantitatively with an explicit energy balance?
2. **RQ2.** What is the energy error of FDTD/FEM, PINNs and neural operators versus
   resolution, precision, frequency and modulation amplitude?
3. **RQ3.** Does an explicit conservation loss reduce the energy residual without
   degrading the fields?
4. **RQ4.** How do rectangular, annular and toroidal geometries change modal density,
   Q, field concentration, pump power and conversion efficiency?
5. **RQ5.** Can an optimizer find an apparent E_net > 0, and does it survive
   refinement, a second solver, float64 and full accounting?
6. **RQ6.** Do surrogates and neural operators preserve conservation out of
   distribution?
7. **RQ7.** Which negative results set useful bounds on future vacuum-energy
   proposals?
8. **RQ8.** The extractable-work question above.

## 3. Falsification checklist (mandatory for any E_net > 0)

Apply these in order:

1. numerical precision
2. time step
3. mesh
4. number of modes
5. domain size
6. boundary implementation
7. signs and units
8. initial stored energy
9. work done by the modulation and the boundaries
10. material and resistive losses
11. Poynting flux
12. independent solver
13. sensitivity analysis
14. repetition
15. adversarial code review

Only after all of these may an unmodelled physical hypothesis be discussed. Hardware
claims additionally need traceable calibration, simultaneous measurement at every
port, thermal and RF control, an uncertainty budget and blind tests.

## 4. Methodological principles

- Do not presuppose free energy. Do not presuppose its impossibility as an
  experimental outcome.
- Reproduce published results before proposing extensions.
- Close the energy balance of every experiment.
- Separate stored energy from produced energy, and modulation energy from output.
- Study convergence. Use more than one numerical method.
- Version code, parameters and results.
- Preserve negative results. Never adjust results to fit a hypothesis.
- An LLM is never a physical authority. Validate claims against primary literature
  with DOIs.

## 5. Roadmap (gates)

| Phase | Content | Exit gate |
|---|---|---|
| 0 | Governance: repo, units, seeds, environment, literature matrix | Reproducible run on a clean machine |
| 1 | Maxwell 1D static: analytic + FDTD | Convergence and energy residual characterized |
| 2 | Classical boundary/permittivity modulation | Balance closes within a convergent tolerance |
| 3 | Minimal DCE: reproduce a published model | Observables and scaling match the literature |
| 4 | Scientific ML: PINN, hard constraints/operator, surrogate | Benchmark vs classical solver, ≥5 seeds |
| 5 | 2D/3D reference: rectangular, cylindrical, annular | A second solver confirms |
| 6 | Toroid: eigenmodes, Q, fields | Validation against the toroidal literature or experiment |
| 7 | Inverse design with the full energy objective | Top candidates survive revalidation |
| 8 | Hardware | Only hypotheses that survived; metrology with an uncertainty budget |

### Experiment list

| ID | Content |
|---|---|
| EXP-0001 | Analytic 1D modes |
| EXP-0002 | FDTD convergence |
| EXP-0003 | Static energy balance |
| EXP-0004 | Classical parametric modulation with pump work |
| EXP-0005 | Ω/ω × δ sweep |
| EXP-0006 | Reproduction of a published DCE model |
| EXP-0007 | PINN Maxwell |
| EXP-0008 | PINN energy-loss ablation |
| EXP-0009 | Neural operator / hard constraints |
| EXP-0010 | Rectangular vs annular |
| EXP-0011 | Toroidal eigenmodes |
| EXP-0012 | Auditable inverse design |
| EXP-0013 | Quantum work audit (ergotropy, passivity, entropy, control cost) |

## 6. Standard experiment record

Each experiment folder contains `README.md`, `run.py`, `results/manifest.json`,
`results/metrics.json`, `figures/` and `conclusion.md`. The manifest records the git
commit, environment, configuration hash and date. The conclusion records:

- hypothesis
- energy in, out, stored and lost
- residual
- convergence
- independent validation
- result
- interpretation
- limitations
- the next action with the highest information gain
- evidence for and against

## 7. Publication strategy (from the source documents, refined)

- **Methodology / software papers first.** These are publishable regardless of the
  E_net outcome.
- **Anomaly claims only after independent validation.** For example: "Investigation
  of anomalous net energy transfer…".
- See [research-plan](../research-plan/preprint-01-action-plan.md) for the concrete
  plan for the first preprint.
