# Action plan: first experiments to an arXiv preprint

**Status date:** 2026-10-04. **Owner:** N. C. D. Dantas. **Repository:** CavityLab v0.1.0.dev0.

## 1. Decision: what the first preprint is

Of the publication lines in the source documents, the first preprint is the
**methodology paper**. It is publishable whatever the E_net outcome, and it builds
directly on EXP-0001 to EXP-0005, which are already implemented.

> **Working title.** *Discrete-exact energy accounting for FDTD simulations of
> time-modulated cavities: separating pump work from numerical artifacts*
>
> **Target.** arXiv physics.comp-ph, cross-listed to physics.optics. After feedback,
> submit to Journal of Computational Physics, Computer Physics Communications or
> IEEE Trans. Antennas Propag.
>
> **Tone.** Verification methodology. Vacuum-energy claims are motivation and
> context; the paper does not argue for or against them. See the positioning in
> [prior_art_preprint01.md](../literature/prior_art_preprint01.md).

### Why this paper first

- It answers the most basic question the project must solve: can we tell physical
  energy from numerical artifacts? That is the gate for every later step: DCE,
  ML, toroids and inverse design.
- The prior-art search found the Yee energy identity in the literature but **not**:
  - the exact discrete pump-work term for ε(t);
  - the quantification of naive-ledger false gains;
  - an audit protocol for time-modulated FDTD.
- Everything runs on a laptop in minutes. There are no external solvers or GPUs on
  the critical path.

The second paper is already seeded by EXP-0006: *Quantum work audit of the parametric
dynamical Casimir effect*. Its content is the efficiency law η = 1 − δ_th/δ and a
classical FDTD transfer matrix that predicts vacuum photon numbers. See §7.

## 2. Claims and evidence status

| # | Claim | Evidence | Status |
|---|---|---|---|
| C1 | The Yee scheme with ε(x,t), σ, soft sources and ports satisfies an exact discrete energy balance. The pump term is `½DⁿD^{n+1}(1/ε^{n+1}−1/εⁿ)Δx`. | Derivation ([theory](../theory/discrete_energy_identity.md)); EXP-0003/0004/0005 residuals of 1e-15 to 1e-12, growing like n_steps·ε | ✅ done |
| C2 | Naive ledgers produce residuals of O(Δt²) (time-averaged energy) or O(Δt) (left-point work and dissipation), at the 10⁻⁵ to 10⁻² level, **with either sign**. Positive values read as false net gain, even in a passive decaying cavity. | EXP-0003 (orders 2.1 / 0.9 / 1.2; B: +1.6e-4), EXP-0004 (+3.7e-4 at δ = 0.04; +0.35% of input in the multimode case), EXP-0005 map | ✅ done (map pending final run) |
| C3 | The solver and ledger reproduce Hill/Floquet theory for uniform modulation: tongues, threshold δ_th = 2/Q, both Floquet branches. FDTD converges to the ODE at second order. | EXP-0004A (rates within 0.1% of ω/Q; order 2.00), EXP-0005 map | ✅ done |
| C4 | Ledger closure is necessary, not sufficient. In an equidistant-spectrum multimode case the ledger closes while the solution is unresolved. Convergence plus an independent reference are required. | EXP-0004B (W/W₀ diverges with N; centroid of modal energy) | ✅ done (negative result) |
| C5 | A naive ledger inside an optimization loop is exploited: the optimizer finds "over-unity" designs that the exact ledger and refinement reject. | EXP-0005b: naive objective > 0 in 93% of 300 trials; best +154% of the energy scale; optimizer picks the coarsest grids; 5/5 top candidates rejected (first-order decay); exact ≤ 1.5e-14 | ✅ done |
| C6 | A widely used open-source FDTD code's energy diagnostics behave like a naive ledger under ε(t). | Planned: Meep check | ⏳ to do (R4). Include only if verified in source code. |
| C7 | A converged multimode benchmark (non-equidistant spectrum) agrees with an independent coupled-mode solver. | Planned | ⏳ to do (R2) |

## 3. Remaining work before submission

Items are in order of information gain per effort.

| ID | Task | Output | Effort | Gate |
|---|---|---|---|---|
| R1 | Finish EXP-0005 (running) and write its conclusion; make figures publication-grade | Fig. resonance map + naive map | 0.5 d | Stability agreement ≥ 98%, exact ledger PASS |
| R2 | **EXP-0004c, converged multimode benchmark.** Dielectric-loaded cavity (non-equidistant spectrum) with a modulated slab, compared against an independent coupled-mode ODE solver built from the discrete eigenmodes | Second-solver validation (the "independent validation" rule) | 2–3 d | W(t) converges at 2nd order; agrees with the coupled-mode solver to < 1% |
| R3 ✅ | **EXP-0005b, optimizer exploit demo (done 2026-10-04).** Random search or Optuna over (ν, δ, N, S, slab position), maximizing the naive `E_net/E_in`. Re-audit the top-k candidates with the exact ledger, refinement and float64 | The "how artifacts become discoveries" figure; motivates C5 | 1–2 d | Every naive "gain" rejected; report the rejection rate |
| R4 | **Meep check.** Implement uniform ε(t) in Meep; compare `field_energy_in_box` with the exact ledger. Read the Meep source to learn the H time alignment | C6, or drop it | 2 d | Stated only if reproducible; otherwise omitted |
| R5 | Precision study: float32 vs float64 (GPU relevance), and Kahan summation of the ledger | Table: residual floor vs precision and steps | 0.5 d | Documented floor |
| R6 | Prior-art full-text checks: Sarkar (2022), Taravati et al. (2019, 2024), the MOTL Poynting-FDTD paper, Dodonov & Klimov (1996) on the 1D cascade | Updated novelty statement | 1 d | No overlooked prior derivation of the pump term |
| R7 | Reproducibility package: CI (GitHub Actions, pytest + all EXP runs), `v0.1.0` tag, Zenodo DOI, `environment.yml`, `make reproduce` | DOI cited in the paper | 1 d | A clean-machine run reproduces all figures |
| R8 | Writing: draft from [papers/preprint-01/outline.md](../../papers/preprint-01/outline.md); internal adversarial review (checklist §5) | `papers/preprint-01/main.tex` | 4–5 d | All numbers come from committed `metrics.json` |

**Critical path:** R1 → R2 → R6 → R8 → R7 (R3 done). About **2–3 weeks** at a steady part-time
pace. R4 is optional and must not block submission.

## 4. Planned figures

1. **Ledger closure** across static scenarios A–D: exact vs naive residual time
   series. *(EXP-0003; exists.)*
2. **Naive residual vs resolution**, with orders 2 and 1. *(EXP-0003 and EXP-0004A
   convergence; exists.)*
3. **Uniform modulation vs Hill ODE**: energy growth, pump work, loss and residual.
   *(EXP-0004A; exists.)*
4. **Resonance map**: Floquet rate from FDTD monodromy vs theory, plus the naive
   apparent-gain map. *(EXP-0005; running.)*
5. **Multimode**: the equidistant cascade (non-convergence, ledger closes) vs the
   non-equidistant case (converged, matches the coupled-mode solver).
   *(EXP-0004B exists; R2 pending.)*
6. **Optimizer exploit**: naive objective landscape, top candidates, re-audit
   outcome. *(EXP-0005b; exists.)*
7. *(Optional)* Meep comparison. *(R4.)*

## 5. Internal adversarial review checklist (before upload)

- [ ] Every number in the text is traceable to a `metrics.json` and a git commit.
- [ ] Signs: E_net convention stated once and used consistently. The pump work sign
      is checked against the continuum `-½E²∂ε/∂t`.
- [ ] Units: 1D energies are per unit area [J/m²]; normalized units are stated.
- [ ] Every claimed order of accuracy comes with the fitted range.
- [ ] The negative result (EXP-0004B) is reported, not hidden.
- [ ] No claim of novelty for the Yee energy identity itself.
- [ ] The vacuum-energy motivation is phrased neutrally, with Jaffe (2005) and
      Moddel (2019) cited.
- [ ] Code and data DOI (Zenodo) included; license checked (MIT).
- [ ] Every reference DOI is resolved through doi.org. No unverified references.

## 6. Experiment log so far (all on a laptop CPU, float64)

| EXP | Topic | Key result | Gate |
|---|---|---|---|
| 0001 | Analytic modes | ω_FDTD = Yee dispersion to 2.7e-12; energy drift 3e-15 | ✅ |
| 0002 | Convergence | Order 2.00 (S = 0.25/0.5/0.9); magic step exact | ✅ |
| 0003 | Static balance | Exact ≤ 1.9e-12; naive 2e-5 to 3e-4, including a **+1.6e-4 false gain** in a passive cavity; Q = 50.001 | ✅ |
| 0004 | Parametric modulation | Floquet from FDTD = theory (≤ 0.1% of ω/Q); threshold 2/Q; order 2.00 vs ODE; multimode cascade ledger 4e-15 but not converged (negative result) | ✅ |
| 0005 | (ν, δ) sweep + dataset | *Running*: tongue map, naive false-gain map, 688-row audited CSV | ⏳ |
| 0005b | Optimizer exploit | The naive objective finds +154% false gain on the coarsest grids; all top candidates rejected by the audit | ✅ |
| 0006 | Quantum DCE (Gaussian) | N = sinh²(δτ/4); FDTD transfer matrix predicts vacuum photons (order ~1.9); **ergotropy/pump work = 1 − δ_th/δ** | ✅ |

## 7. Pipeline after Preprint 01

| Preprint | Content | Depends on |
|---|---|---|
| 02 | *Quantum work audit of the parametric DCE.* η = 1 − δ_th/δ, a classical-FDTD prediction of vacuum photon creation, a published-model reproduction (Wilson 2011 / PRA 96 033851) and a multimode extension | EXP-0006 part 2, EXP-0013 |
| 03 | *Benchmarking PINNs and neural operators under energy-conservation constraints.* Uses the EXP-0005 audited dataset and the exact ledger as the metric | EXP-0007 to EXP-0009 |
| 04 | *Toroidal resonators under strict energy accounting.* 2D/3D, second solver (Meep/FEM), inverse design with the exact objective | EXP-0010 to EXP-0012 |

## 8. Decisions taken autonomously on 2026-10-04 (for review)

1. **Repository language and layout.** English. `src/` layout following the
   CavityLab proposal; Portuguese originals kept in `docs/context/source-documents/`.
2. **E_net sign convention.** `E_net = E_out + E_loss + ΔE − E_in`. The literal
   formula in the source documents does not vanish under conservation. The exit
   criterion of the context document is unchanged.
3. **D-formulation Yee with semi-implicit σ.** It handles ε(t) exactly and gives a
   closed-form discrete pump term.
4. **Naive ledger alongside the exact one.** It quantifies artifacts and becomes a
   central result of the paper.
5. **Generic pump phase (π/4)** for time-series comparisons, and full monodromy for
   Floquet comparisons, after the degenerate-quadrature finding.
6. **Multimode slab case reclassified** as a documented stress test (not converged),
   instead of being tuned until it "looked converged".
7. **Quantum model: exact Gaussian moments in numpy** instead of QuTiP. It is exact
   for this linear problem, has no new dependency, and enables the classical-
   quantum bridge. QuTiP stays an optional extra for nonlinear models.
8. **License kept as MIT** (already in the repository). The proposal suggested
   evaluating BSD/Apache. MIT is permissive and compatible; revisit if Meep (GPL)
   code is ever vendored, though calling it as an external package is fine.
9. **Nothing pushed to GitHub.** Commits are local only. Publishing is the owner's
   call.
