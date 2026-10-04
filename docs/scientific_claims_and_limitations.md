# Scientific claims and limitations

This page states what CavityLab results do and do not support.

## What CavityLab claims

- **Energy accounting, not energy generation.** The library books every energy
  term (stored, pump/control, source, loss, port) consistently with the
  discretization. A closed ledger shows that the *numerical* bookkeeping is
  consistent.
- **Verified reference solvers.** Native solvers are checked against analytical
  solutions: Yee dispersion, d'Alembert pulses and Hill/Floquet theory. They are
  also checked for convergence order. The status of each check is recorded in the
  experiment `conclusion.md` files.

## What CavityLab does not claim

- **No net energy generation from the vacuum has been observed.** No project
  experiment has produced `E_net > 0` in the exact ledger. All field-energy growth
  seen so far (EXP-0004, EXP-0005) is fully paid for by pump work.
- **Classical simulations say nothing about vacuum fluctuations.** The classical
  parametric runs need a nonzero seed field. Photon creation from vacuum (DCE)
  needs a quantum model (EXP-0006, planned).
- **A closed ledger is necessary but not sufficient.** The exact ledger closes to
  round-off even for unresolved solutions. EXP-0004 part B is an example: the
  equidistant-spectrum cascade reaches grid scale. Accuracy needs convergence
  studies and independent references as well.

## Rules for interpreting a positive residual

1. Report it as a **candidate computational anomaly**, never as an energy source.
2. Apply the falsification checklist in
   [project_context.md §3](context/project_context.md) in order.
3. Check the residual against the expected ledger error:
   - The exact ledger has an expected residual of about `n_steps · 1e-16`.
   - Naive ledgers have O(Δt) to O(Δt²) residuals at the 10⁻⁵ to 10⁻² level.
   - A residual of that size is an artifact by construction.
4. Results that disappear under refinement are artifacts. Keep and publish them as
   negative results.

## Known limitations (v0.1)

- The solver is 1D, with Ez/Hy polarization and PEC walls. Absorbers are graded
  electric conductivity, not a PML, so they reflect somewhat. Reflection does not
  affect the ledger.
- Materials have no dispersion yet. That excludes Lorentz/Drude media and the
  physical cost of modulating a dispersive medium.
- Modulation is imposed: `ε(t)` is prescribed. The pump's own dynamics and
  efficiency (the energy drawn from the wall plug) are not modelled. Pump work is
  the work done *on the field*.
- There is no independent second solver yet. Meep and FEM integration are planned,
  as are the quantum models (QuTiP/Gaussian), PINNs and neural operators.
