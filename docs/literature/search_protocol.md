# Systematic review protocol (draft v0.1, to be registered on OSF)

The goal is to turn the scoping review into a PRISMA 2020-compliant systematic map
that confirms or refutes the "integrated pipeline gap" hypothesis and builds
`literature_matrix.csv` up to 100–300 prioritized studies.

## Clusters and base search strings

| Cluster | Search string |
|---|---|
| DCE | ("dynamical Casimir" OR "time-dependent cavity" OR "moving boundary quantum field") AND (energy OR photons OR modulation OR superconducting) |
| Geometry | (Casimir OR "electromagnetic cavity" OR resonator) AND (toroidal OR torus OR annular OR "complex geometry") |
| Scientific ML | ("physics-informed" OR "neural operator" OR surrogate) AND (Maxwell OR electromagnetics) |
| Inverse design | ("inverse design" OR "topology optimization" OR "Bayesian optimization") AND ("electromagnetic cavity" OR resonator OR photonics) |
| Time-varying FDTD *(added for Preprint 01)* | (FDTD OR "finite-difference time-domain" OR Yee) AND ("time-varying" OR "time-modulated" OR "temporal" OR parametric) AND (energy OR conservation OR Poynting) |
| Quantum thermodynamics *(added in v2)* | (ergotropy OR "passive state" OR "quantum battery" OR "work extraction") AND (cavity OR resonator OR parametric OR Casimir) |

## Databases

- Web of Science, Scopus and IEEE Xplore.
- APS (PROLA) and INSPIRE-HEP.
- arXiv for preprints.
- Crossref and Semantic Scholar for metadata and snowballing.
- Google Scholar for backward and forward citation snowballing only.

## Extraction fields per study

The fields match `literature_matrix.csv`:

- id
- authors
- year
- title
- venue
- doi_or_arxiv
- cluster
- model type (classical / quantum / ML)
- geometry
- equations
- boundary conditions
- energy source
- energy accounting (none / partial / explicit / discrete-exact)
- numerical method
- experimental validation
- code/data availability
- limitations
- relevance to RQ1–RQ8
- verified (yes / no)

## Screening

There are two stages: title/abstract, then full text. A second screener (a human or
an independent agent run) checks a 20% sample, and agreement is reported as Cohen's
κ. An LLM may help screening but never decides inclusion alone. Every identifier
must be resolved through doi.org or arXiv before inclusion.

## Status

- **Done.** Scoping review v2 and the prior-art search for Preprint 01
  ([prior_art_preprint01.md](prior_art_preprint01.md)).
- **Pending.** OSF registration, database exports and a PRISMA flow diagram.
