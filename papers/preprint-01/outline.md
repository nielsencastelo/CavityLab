# Preprint 01: outline and draft abstract

**Working title.** Discrete-exact energy accounting for FDTD simulations of
time-modulated cavities: separating pump work from numerical artifacts

**Authors.** N. C. D. Dantas (and collaborators to be defined)

## Draft abstract (≈200 words)

Time-modulated electromagnetic media and cavities are central to parametric
amplification, temporal metamaterials and analogues of the dynamical Casimir effect.
In these systems the field energy is not conserved: the modulation does work on the
field. Telling that physical pump work apart from discretization error is essential.
It matters most when simulations feed automated optimization or test extraordinary
energy claims.

We show that the Yee FDTD scheme in D-formulation, with time-varying permittivity,
semi-implicit conduction, soft sources and absorbing ports, satisfies an exact
discrete energy balance. The modulation contributes a closed-form pump-work term,
½DⁿD^{n+1}(1/ε^{n+1} − 1/εⁿ)Δx, so the full ledger closes to floating-point round-off.

Textbook energy estimators leave residuals at the 10⁻⁵ to 10⁻² level at practical
resolutions. They converge at first or second order and have either sign, including
apparent net generation in a passive decaying cavity.

The audited solver reproduces Hill–Floquet theory for uniformly modulated lossy
cavities: instability tongues, the threshold δ_th = 2/Q and both Floquet branches,
with second-order convergence. A multimode stress test shows that ledger closure is
necessary but not sufficient. We release CavityLab, an open-source library that
makes the audit a default part of every simulation.

## Structure

1. **Introduction.** Time-varying media, DCE analogues, verification, and the risk
   of artifacts that look like discoveries. Positioning: the Yee energy identity is
   known (Joly; Chen–Li–Liang; Bekmambetova et al.); the pump term, the audit
   protocol and the artifact taxonomy are the contribution.
2. **Discrete energy identity.** Lossless case (known); pump, conduction, source and
   port terms; round-off floor; continuum limit of the pump term.
3. **Naive estimators and their errors.** Taxonomy table and order analysis.
4. **Verification.**
   1. Modes and dispersion (EXP-0001).
   2. Convergence (EXP-0002).
   3. Static ledgers A–D (EXP-0003).
5. **Parametric cavities.**
   1. Hill reduction and Floquet monodromy from FDTD (EXP-0004A).
   2. Resonance map (EXP-0005).
   3. The degenerate quadrature start.
6. **Multimode cases.** Equidistant cascade (stress test) vs a converged
   non-equidistant benchmark with an independent coupled-mode solver (EXP-0004B and
   R2).
7. **Artifacts under optimization.** Naive-objective exploit and its rejection by
   the audit (R3).
8. **Discussion.**
   - Necessary and sufficient checks: the falsification checklist.
   - Implications for automated discovery pipelines and for evaluating energy
     claims.
   - Extension to 2D/3D, PML and dispersive media.
9. **Code and data availability.** CavityLab version, Zenodo DOI, `run.py` per
   figure.

## Figure → source mapping

| Figure | Script | Metrics |
|---|---|---|
| 1 | `experiments/EXP-0003_static_energy_balance/run.py` | `scenarios.*` |
| 2 | EXP-0003 + EXP-0004 | `naive_scaling`, `A.convergence` |
| 3 | `experiments/EXP-0004_parametric_modulation/run.py` | `A.rows` |
| 4 | `experiments/EXP-0005_parametric_sweep/run.py` | all |
| 5 | EXP-0004B + R2 | `B`, `C_naive_under_pumping` |
| 6 | `experiments/EXP-0005b_optimizer_exploit/run.py` | `top_candidates` |
