# Instructions for AI agents working on CavityLab / ZPE-AI

Read [docs/context/project_context.md](docs/context/project_context.md) before you
change code or hypotheses. The Portuguese originals are in
`docs/context/source-documents/`.

## Repository

- The repository is in English: code, docs, commit messages and figures.
- Python environment: `C:/Users/niels/anaconda3/python.exe` (3.13) with numpy,
  scipy, matplotlib and pytest. The package is installed editable with
  `pip install -e .`.
- GPU environment: `.venv-gpu/Scripts/python.exe` (torch 2.14 + CUDA 12.6 for the RTX
  A2000 12 GB). Use it for batched sweeps (`FDTD1DBatch(..., xp="torch")` or
  `backend="auto"`). Do not import torch inside Anaconda (OpenMP clash). CuPy is
  blocked by Windows application control.
- Prefer the batched solver and the GPU for any sweep with ≥ 1000 simulations.
  Audit-grade results must be float64.
- Run the tests with `python -m pytest`.
- Run an experiment with `python experiments/EXP-XXXX_*/run.py`. It writes
  `results/manifest.json`, `results/metrics.json` and `figures/`.

## Git

- The owner makes every commit and push. Do not run `git commit` or `git push`.
  Prepare the changes and suggest a commit message.
- Never add "Co-Authored-By" or any AI attribution to commits or PRs.

## Scientific rules (non-negotiable)

1. Never treat zero-point-energy extraction as a demonstrated fact.
2. Keep an explicit energy balance in every experiment. Use `EnergyAudit` and the
   sign convention `E_net = E_out + E_loss + ΔE_stored - E_in`.
3. A positive E_net is a *candidate computational anomaly*. Run the falsification
   checklist (refinement, precision, signs/units, pump work, losses, independent
   solver) before anything else.
4. Never accept a positive result without a convergence study.
5. Do not discard results just because they contradict expectations. Do not change
   results to fit a hypothesis. Preserve negative results.
6. Distinguish explicitly between scientific fact, hypothesis, numerical result and
   speculation.
7. Cite primary literature with a DOI. Never invent references. Mark unverified
   identifiers in `docs/literature/literature_matrix.csv`.
8. Write a `conclusion.md` after each experiment: evidence for and against, and the
   next experiment with the highest information gain.
9. An LLM is never a physical authority.
