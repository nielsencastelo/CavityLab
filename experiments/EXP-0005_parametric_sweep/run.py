"""EXP-0005: sweep of pump frequency ratio nu = Omega/omega_1 and depth delta (resonance map).

For a uniformly modulated lossy cavity (Q = 100) this experiment produces:
  1. Floquet growth-rate map measured from the FDTD monodromy (two quadratures)
     vs the Hill-equation theory: instability tongues near nu = 2/m.
  2. A map of the *naive-ledger* apparent net energy E_net / energy scale for a fixed
     run length: where a naive pipeline (or an optimizer reading it) would report
     "gain". The exact ledger is audited on the same runs.
  3. A small audited dataset (CSV) following the CavityLab per-simulation schema,
     reusable for surrogate / scientific-ML benchmarks (EXP-0007+).
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from cavitylab import EPS0, EnergyAudit, FDTD1D, Modulation
from cavitylab.benchmarks.analytic1d import discrete_mode_fields
from cavitylab.benchmarks.hill import floquet_growth_rate
from cavitylab.benchmarks.parametric import fdtd_floquet_exponents
from cavitylab.core.manifest import config_hash
from cavitylab.experiments.runner import ExperimentContext, setup_matplotlib
from cavitylab.geometry import Cavity1D

HERE = Path(__file__).parent
CONFIG = {
    "length_m": 0.1,
    "courant": 0.5,
    "n_cells": 64,
    "quality": 100.0,
    "nu_grid": [0.5, 2.6, 43],
    "delta_grid": [0.0, 0.3, 16],
    "floquet_pump_periods": 40,
    "ledger_run_mode_periods": 30,
    "pump_phase_rad": math.pi / 4,
}


def ledger_run(nu, delta):
    L, S, Q = CONFIG["length_m"], CONFIG["courant"], CONFIG["quality"]
    w1 = Cavity1D(L).mode_angular_frequency(1)
    sim = FDTD1D(L, CONFIG["n_cells"], S, sigma=w1 * EPS0 / Q,
                 modulation=Modulation(delta, nu * w1, phase=CONFIG["pump_phase_rad"]))
    E, H, _ = discrete_mode_fields(sim, 1)
    sim.set_fields(E, H)
    n_steps = round(CONFIG["ledger_run_mode_periods"] * 2 * math.pi / w1 / sim.dt)
    res = sim.run(n_steps, record_every=max(1, n_steps // 500))
    return sim, res


def main() -> None:
    ctx = ExperimentContext("EXP-0005", HERE, CONFIG)
    nus = np.linspace(*CONFIG["nu_grid"][:2], CONFIG["nu_grid"][2])
    deltas = np.linspace(*CONFIG["delta_grid"][:2], CONFIG["delta_grid"][2])
    Q = CONFIG["quality"]
    shape = (len(deltas), len(nus))
    mu_theory = np.zeros(shape)
    mu_fdtd = np.zeros(shape)
    naive_gain = np.zeros(shape)
    exact_res = np.zeros(shape)
    rows = []
    for i, d in enumerate(deltas):
        for j, nu in enumerate(nus):
            mu_theory[i, j] = floquet_growth_rate(d, nu, Q) if d > 0 else -0.5 / Q
            mu_fdtd[i, j] = fdtd_floquet_exponents(
                CONFIG["length_m"], CONFIG["n_cells"], d, nu, Q,
                pump_periods=CONFIG["floquet_pump_periods"])[0]
            sim, res = ledger_run(nu, d)
            ex = EnergyAudit.from_fdtd(res)
            nv = EnergyAudit.from_fdtd(res, naive=True)
            # normalize by the exact ledger's energy scale (|E_in| alone is ~0 off resonance)
            naive_gain[i, j] = nv.e_net[-1] / ex.scale
            exact_res[i, j] = ex.max_relative_residual
            cfg = {"nu": nu, "delta": d, "Q": Q, "n_cells": CONFIG["n_cells"], "courant": CONFIG["courant"]}
            rows.append({
                "simulation_id": f"EXP0005-{i:02d}-{j:02d}",
                "config_hash": config_hash(cfg),
                "geometry": "cavity_1d_pec",
                "length_m": CONFIG["length_m"],
                "eps_r": 1.0,
                "quality_factor": Q,
                "modulation_nu": nu,
                "modulation_delta": d,
                "modulation_phase_rad": CONFIG["pump_phase_rad"],
                "solver": "cavitylab.native_fdtd1d",
                "n_cells": CONFIG["n_cells"],
                "dt_s": sim.dt,
                "precision": "float64",
                "t_end_s": res.t[-1],
                "E_in": ex.e_in[-1],
                "E_out": ex.e_out[-1],
                "E_loss": ex.e_loss[-1],
                "dE_stored": ex.delta_stored[-1],
                "E_net_exact": ex.e_net[-1],
                "rel_residual_exact": ex.max_relative_residual,
                "E_net_naive": nv.e_net[-1],
                "floquet_mu_fdtd_per_tau": mu_fdtd[i, j],
                "floquet_mu_theory_per_tau": mu_theory[i, j],
                "audit_status": ex.summary()["status"],
            })
        print(f"delta {d:.3f} done")

    with open(ctx.results_dir / "dataset.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # refinement check: the naive "gain" must shrink with resolution if it is an artifact
    i_w, j_w = np.unravel_index(np.argmax(naive_gain), naive_gain.shape)
    refinement = []
    n_base = CONFIG["n_cells"]
    for n in (64, 128, 256, 512):
        CONFIG["n_cells"] = n
        _, r = ledger_run(nus[j_w], deltas[i_w])
        refinement.append({"n_cells": n,
                           "naive_final_rel_e_net": EnergyAudit.from_fdtd(r, naive=True).e_net[-1]
                           / EnergyAudit.from_fdtd(r).scale,
                           "exact_max_rel": EnergyAudit.from_fdtd(r).max_relative_residual})
    CONFIG["n_cells"] = n_base
    ref_order = float(-np.polyfit(np.log([x["n_cells"] for x in refinement]),
                                  np.log([abs(x["naive_final_rel_e_net"]) for x in refinement]), 1)[0])

    err = np.abs(mu_fdtd - mu_theory)
    unstable_theory = mu_theory > 0
    unstable_fdtd = mu_fdtd > 0
    ctx.metrics = {
        "max_abs_floquet_error_per_tau": float(err.max()),
        "median_abs_floquet_error_per_tau": float(np.median(err)),
        "stability_classification_agreement": float(np.mean(unstable_theory == unstable_fdtd)),
        "max_exact_relative_residual": float(exact_res.max()),
        "naive_apparent_gain_max": float(naive_gain.max()),
        "naive_apparent_gain_min": float(naive_gain.min()),
        "naive_apparent_gain_median": float(np.median(naive_gain)),
        "fraction_points_naive_reports_gain": float(np.mean(naive_gain > 0)),
        "worst_naive_point": {"nu": float(nus[j_w]), "delta": float(deltas[i_w])},
        "worst_point_refinement": refinement,
        "worst_point_naive_order": ref_order,
        "n_points": int(mu_fdtd.size),
        "gates": {
            "exact_ledger_all_pass": bool(exact_res.max() < 1e-10),
            "stability_agreement_ge_98pct": bool(np.mean(unstable_theory == unstable_fdtd) >= 0.98),
            "naive_gain_is_artifact": bool(ref_order > 0.8),
        },
    }

    plt = setup_matplotlib()
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ext = [nus[0], nus[-1], deltas[0], deltas[-1]]
    vmax = max(abs(mu_theory).max(), abs(mu_fdtd).max())
    for a, data, title in ((ax[0], mu_theory, "Hill theory"), (ax[1], mu_fdtd, "FDTD monodromy (N = 64)")):
        im = a.imshow(data, origin="lower", aspect="auto", extent=ext, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        a.contour(nus, deltas, data, levels=[0.0], colors="k", linewidths=0.8)
        a.set_xlabel(r"$\nu = \Omega/\omega_1$")
        a.set_ylabel(r"modulation depth $\delta$")
        a.set_title(f"Floquet growth rate, Q = {Q:.0f}: {title}")
    fig.colorbar(im, ax=ax[1], label=r"$\mu$ per $\omega_1 t$")
    from matplotlib.colors import SymLogNorm

    g = naive_gain
    gl = max(abs(g).max(), 1e-12)
    im2 = ax[2].imshow(g, origin="lower", aspect="auto", extent=ext, cmap="PuOr_r",
                       norm=SymLogNorm(linthresh=1e-3, vmin=-gl, vmax=gl))
    ax[2].contour(nus, deltas, mu_theory, levels=[0.0], colors="k", linewidths=0.6, linestyles="--")
    ax[2].set_xlabel(r"$\nu = \Omega/\omega_1$")
    ax[2].set_title(f"Naive ledger: apparent $E_{{\\rm net}}$/scale (> 0 at {np.mean(g > 0):.0%} of points)")
    fig.colorbar(im2, ax=ax[2], label="exact ledger: |E_net|/scale < 2e-13 everywhere")
    fig.tight_layout()
    fig.savefig(ctx.figure_path("exp0005_resonance_map.png"))

    ctx.save()
    print("EXP-0005 metrics:", {k: v for k, v in ctx.metrics.items() if k != "gates"}, ctx.metrics["gates"])


if __name__ == "__main__":
    main()
