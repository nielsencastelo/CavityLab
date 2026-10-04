# Literature review (condensed English version, v2)

This is the English summary of *ZPE-AI: Revisão aprofundada da literatura
científica, v2* (26 Sep 2026); the Portuguese original is in
[`../context/source-documents/`](../context/source-documents/). It is a scoping
review driven by the hypothesis. It is **not** a PRISMA systematic review (see
[search_protocol.md](search_protocol.md)). The machine-readable matrix is
[literature_matrix.csv](literature_matrix.csv).

## Key conclusions

1. **Zero-point energy is part of standard QFT**, but its existence is not a
   demonstrated mechanism for continuous work extraction. Observables come from
   regularized differences between configurations. Jaffe (2005) shows that the
   Casimir force can be computed as a relativistic quantum force between charges and
   currents, without invoking the absolute vacuum energy.
2. **Extraction proposals exist historically, but none is demonstrated.**
   - Forward (1984) proposed a vacuum-fluctuation battery; its cycle needs re-expansion.
   - Cole & Puthoff (1993) gave a thermodynamic analysis.
   - Moddel (2019) found no extraction approach reliably demonstrated, and several
     incompatible with detailed balance.
3. **The static Casimir effect is established and depends strongly on geometry**
   (Bordag et al. 2001; Maclay 2000 for rectangular cavities). More confinement does
   not mean more usable energy, so geometries must be compared normalized: volume,
   Q, pump power, losses and stored energy.
4. **The DCE is the key experimental bridge.**
   - Wilson et al. (2011) observed photons and two-mode squeezing from a SQUID-
     modulated transmission line.
   - Dodonov (2020, 2025) reviews circuit, optomechanical and parametric platforms.
   - The pump supplies the energy, so `E_pump(t)` must be a first-class variable.
5. **Toroids.** Physical toroidal cavities, compactified toroidal topologies and
   toroidal dipole modes are three different things.
   - Physical toroidal cavities: Janaki & Dasgupta (1990) solved the eigenmodes;
     Giraldez (USP) measured a copper toroidal cavity.
   - High-Q microtoroids: Armani 2003 (Q > 10⁸); Spillane 2005 (cavity QED).
   - Toroidal metamaterial modes: Tasolamprou 2016.
   - Casimir energy on T^n: Kirsten & Elizalde 1995.
   - The gap is the *combination* toroid + DCE-type modulation + full energy
     accounting + ML + inverse design. The toroid on its own is not new.
6. **Scientific ML for Maxwell is growing but immature.**
   - Schmeing & Pioch (2026) reviewed 139 PINN studies, mostly 2D and trained on
     physics only, and found a lack of benchmarks and of reproducibility.
   - Hard-constraint and operator approaches are promising, e.g. the Fourier–
     Helmholtz–Maxwell neural operator of Leon & Scheinker 2024.
   - AI must come *after* a classical ground truth.
7. **Inverse design is mature in photonics** (Schul et al. 2026). The ZPE-AI
   differentiator is the objective `J = E_out - E_pump - E_other - E_loss - ΔE_stored
   - λR_num - λR_phys`. Without it, an optimizer exploits omitted pump energy or
   numerical error. Extreme candidates must be re-run at higher resolution, in
   float64, with a second solver and with parameter perturbations.
8. **Quantum thermodynamics (v2 extension).**
   - Ergotropy is the maximum work unitarily extractable from a state.
   - Passive states allow no further unitary extraction.
   - Key references: Francica et al. 2020 (coherence and ergotropy); Šafránek, Rosa,
     Binder 2023 (unknown sources); Hokkyo & Ueda 2025 (ETH no-go bound); Kamin,
     Salimi, Santos 2021 (exergy of passive states).
   - Optimal control can amplify the parametric DCE (PRA 96, 033851, 2017).
   - Quantum reservoir engineering: Tóth et al. 2017.
   - This motivates EXP-0013 (quantum work audit) as a second audit layer.

## Most defensible novelty (to be confirmed by a registered systematic review)

No located work implements the full pipeline: reproducible DCE and dynamic cavities,
an explicit and testable energy-accounting module, FDTD/FEM vs PINN/operator
benchmarks, controlled progression to toroids, and inverse design under a net-
energy objective with adversarial revalidation. Every component has precedents; the
integration is the contribution. Even with `E_net ≤ 0` there are at least three
publishable products:

- an energy-accounting methodology for dynamic cavities
- a scientific-ML benchmark for time-dependent Maxwell
- toroidal inverse design under strict conservation

Preprint 01 is the first of these, narrowed to FDTD verification (see the
[prior-art note](prior_art_preprint01.md)).

## Confidence

- **High** for the Casimir, DCE, toroid and PINN state-of-the-art claims.
- **Moderate** for the claim that the integrated pipeline is absent, because proving
  absence needs an exhaustive registered review.
