# Prior art for Preprint 01: energy accounting in FDTD of time-modulated cavities

The search ran on 2026-10-04 with web search, done by a research agent. Every
identifier below was seen on a publisher, arXiv, PMC, ADS or Springer page. Items
marked *(DOI not seen)* are confirmed citations whose identifier still needs
checking. The table is a starting point and is **not** a systematic review (see
[search_protocol.md](search_protocol.md)).

## Verified references

| # | Reference | Venue | DOI / arXiv | Relevance | Cat. |
|---|---|---|---|---|---|
| 1 | Taflove & Hagness (2005), *Computational Electrodynamics: The FDTD Method*, 3rd ed. | Artech House | ISBN 978-1580538329 | Standard FDTD reference | a |
| 2 | Gedney (2011), *Introduction to the FDTD Method for Electromagnetics* | Morgan & Claypool | 10.1007/978-3-031-01712-4 | Textbook | a |
| 3 | Joly (2003), Variational methods for time-dependent wave propagation problems | Springer LNCSE | 10.1007/978-3-642-55483-4_6 | Energy-method stability; Yee leapfrog energy | a |
| 4 | Chen, Li, Liang (2008), Energy-conserved splitting FDTD methods for Maxwell's equations | Numer. Math. 108, 445 | 10.1007/s00211-007-0123-9 | Discrete energy identities | a |
| 5 | Bekmambetova, Zhang, Triverio (2016), A dissipative systems theory for FDTD | arXiv | 1606.08761 | Storage energy and dissipation inequality for FDTD | a |
| 6 | Xiao et al. (2017), Lattice Maxwell system with discrete space-time symmetry and local energy-momentum conservation | arXiv | 1709.09593 | Exact discrete local conservation | a |
| 7 | Gruninger & Griffith (2026), Composite B-spline current deposition, thin-wire FDTD | arXiv | 2605.21450 | States the leapfrog modified energy explicitly | a |
| 8 | Galiffi et al. (2022), Photonics of time-varying media | Adv. Photonics 4, 014002 | 10.1117/1.AP.4.1.014002 | Main review of time-varying media | b |
| 9 | Galiffi, Solís, Yin, Engheta, Alù (2025), Electrodynamics of photonic temporal interfaces | Light Sci. Appl. | 10.1038/s41377-025-01947-2 | Energy at temporal interfaces (D- vs B-conserving); analytic | b |
| 10 | Caloz & Deck-Léger (2020), Spacetime metamaterials, Part I | IEEE TAP 68, 1569 | 10.1109/TAP.2019.2944225 | Framework | b |
| 11 | Pacheco-Peña & Engheta (2020), Antireflection temporal coatings | Optica 7, 323 | *(DOI not seen)* | Temporal metamaterials | b |
| 12 | Solís, Kastner, Engheta (2021), Time-varying materials in presence of dispersion | Photonics Res. 9 | 10.1364/PRJ.427368 | Dispersive temporal discontinuity | b |
| 13 | Hayran, Khurgin, Monticone (2022), ħω versus ħk: dispersion and energy constraints on time-varying photonic materials | Opt. Mater. Express 12, 3904 | 10.1364/OME.471672 | Pump-energy realism | b |
| 14 | Hayran & Monticone (2023), Using time-varying systems to challenge fundamental limitations in electromagnetics | IEEE AP Mag. 65, 29 | arXiv 2205.07142 | Energy budget of claims to beat bounds | b |
| 15 | Mirmoosa, Ptitcyn, Asadchy, Tretyakov (2019), Time-varying reactive elements for extreme accumulation of EM energy | Phys. Rev. Applied 11, 014024 | 10.1103/PhysRevApplied.11.014024 | Lumped pump-work accounting | b |
| 16 | Deshmukh & Milton (2022), An energy conserving mechanism for temporal metasurfaces | — | arXiv 2205.09030 | Energy-conserving parameter switching | b |
| 17 | Zurita-Sánchez, Halevi, Cervantes-González (2009), Reflection and transmission of a wave incident on a slab with a time-periodic dielectric function | Phys. Rev. A 79, 053821 | *(DOI not seen)* | Mathieu/Hill analysis of an ε(t) slab | b/c |
| 18 | Stewart, Smy, Gupta (2016), FDTD modelling of space-time modulated metasurfaces | arXiv | 1612.02087 | FDTD with ε(x,t), no energy audit | b |
| 19 | Taravati & Kishk (2019), Space-time modulation: principles and applications | arXiv | 1903.01272 | Review | b |
| 20 | Taravati, Kishk, Eleftheriades (2024), FDTD simulation of wave transmission through space-time-varying media | arXiv | 2409.19923 | General ε, μ, σ(x,t) FDTD | b |
| 21 | Sarkar (2022), FDTD analysis of guided EM wave interaction with time-modulated dielectric medium | SpringerBriefs | 10.1007/978-981-19-1630-4 | 1D FDTD with ε(t), parametric amplification, code | b |
| 22 | Attiya & Eldesouki (2024), Numerical analysis for temporal and spectral responses of time varying medium | Sci. Rep. | 10.1038/s41598-024-64874-z | ODE solver vs FDTD | b |
| 23 | Ruser (2006), Numerical approach to the dynamical Casimir effect | J. Phys. A 39, 6711 | 10.1088/0305-4470/39/21/S72 | Numerical DCE, moving walls | c |
| 24 | Crocce, Dalvit, Mazzitelli (2001), Resonant photon creation in a 3D oscillating cavity | Phys. Rev. A 64, 013808 | 10.1103/PhysRevA.64.013808 | Parametric resonance conditions | c |
| 25 | Dodonov (2010), Current status of the dynamical Casimir effect | Phys. Scr. 82, 038105 | 10.1088/0031-8949/82/03/038105 | Review | c |
| 26 | Dodonov (2020), Fifty years of the dynamical Casimir effect | Physics 2, 67 | 10.3390/physics2010007 | Review | c |
| 27 | Dalvit, Maia Neto, Mazzitelli (2011), Fluctuations, dissipation and the dynamical Casimir effect | Lect. Notes Phys. | 10.1007/978-3-642-20288-9_13 | Review including analog setups | c |
| 28 | Faccio & Carusotto (2011), Dynamical Casimir effect in optically modulated cavities | EPL 96, 24006 | 10.1209/0295-5075/96/24006 | Closest classical FDTD analogue of the DCE | c |
| 29 | Antunes (2003), Numerical simulation of vacuum particle production with time-dependent non-homogeneous dielectrics | arXiv | hep-ph/0310131 | Lattice DCE with ε(x,t) | c |
| 30 | Wilson et al. (2011), Observation of the DCE in a superconducting circuit | Nature 479, 376 | 10.1038/nature10561 | Effective-length modulation experiment | c |
| 31 | Oskooi et al. (2010), Meep | Comput. Phys. Commun. 181, 687 | 10.1016/j.cpc.2009.11.008 | Open-source FDTD | e |
| 32 | Lehman et al. (2018), The surprising creativity of digital evolution | arXiv | 1803.03453 | Optimizers exploiting simulator bugs | d |

## Novelty assessment

**Already known.** Do not claim these as new:

- The leapfrog modified energy and its use in stability proofs (refs 3, 4, 5, 7).
- The physics of energy at temporal interfaces and of pump work in time-varying
  media (refs 9, 13, 15, 16).
- FDTD of ε(t) and ε(x,t) media (refs 18 to 22, 28).
- Hill/Mathieu tongues of modulated cavities and slabs (refs 17, 24 to 26).

**Gap found.** No paper located in the search does the following:

- Derives the exact discrete pump-work term `½ DⁿD^{n+1}(1/ε^{n+1} - 1/εⁿ)Δx` for Yee
  with ε(t) and closes the FDTD energy budget to machine precision.
- Quantifies the O(Δt) and O(Δt²) spurious residuals of naive ledgers and frames
  them as false "gain" signals.
- Provides an audit protocol (exact ledger, convergence, an independent reference)
  for time-modulated FDTD.

Absence from a web search is not proof. Before claiming priority, check the full
texts of Sarkar (2022) and Taravati et al. (2019, 2024), and the source code of
Meep's `*_energy_in_box`.

**Positioning.** The paper is a verification-methodology contribution in
physics.comp-ph, cross-listed to physics.optics, with a neutral tone. The energy
identity is cited as known, and the contribution is a complete audit protocol plus
benchmarks and an artifact taxonomy.

## Still to verify

- The authors of the MOTL paper "Poynting's theorem for the FDTD method"
  (10.1002/mop.4650080512) and of the URSI GA 2008 paper "Energy conservation and
  Poynting's theorem in the staggered FDTD grid".
- The volume of Monk & Süli (SIAM J. Numer. Anal., 1994).
- The DOIs of Pacheco-Peña & Engheta (2020) and Zurita-Sánchez et al. (2009).
- Whether Meep's energy-in-box functions time-average H.
- Dodonov & Klimov (1996), on the exponential energy growth and packet formation in
  a resonantly vibrating 1D cavity. This is needed to discuss EXP-0004 part B and
  is not yet verified.
