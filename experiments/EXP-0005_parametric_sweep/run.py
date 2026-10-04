"""EXP-0005: sweep of pump ratio nu = Omega/omega_1 and depth delta (resonance maps).

Runs on the batched FDTD solver (cavitylab.solvers.fdtd1d_batch): every map point
and both Floquet quadratures advance together (about 15-20x faster than serial runs).

Map A (Q = 100, delta <= 0.3):
  1. Floquet growth rate from the FDTD monodromy (two quadratures) vs Hill theory.
  2. Naive-ledger apparent net energy E_net / energy scale after 30 mode periods,
     with the exact ledger audited on the same runs.
  3. Audited dataset (CSV) following the CavityLab per-simulation schema.
Map B (Q = 2000, delta <= 0.5): higher-order tongues (nu = 1, 2/3) whose widths scale
  as delta^2 and delta^3, so they appear only at high Q.
Refinement check: the worst naive "gain" must shrink with resolution.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from cavitylab import EPS0
from cavitylab.benchmarks.analytic1d import discrete_mode_fields
from cavitylab.benchmarks.hill import floquet_growth_rate
from cavitylab.benchmarks.parametric import batch_floquet_exponents
from cavitylab.core.manifest import config_hash
from cavitylab.experiments.runner import ExperimentContext, setup_matplotlib
from cavitylab.geometry import Cavity1D
from cavitylab.solvers.fdtd1d import FDTD1D
from cavitylab.solvers.fdtd1d_batch import FDTD1DBatch

HERE = Path(__file__).parent
CONFIG = {
    "length_m": 0.1,
    "courant": 0.5,
    "n_cells": 128,
    "map_A": {"quality": 100.0, "nu_grid": [0.5, 2.6, 43], "delta_grid": [0.0, 0.3, 16]},
    "map_B": {"quality": 2000.0, "nu_grid": [0.5, 2.6, 85], "delta_grid": [0.0, 0.5, 26]},
    "floquet_pump_periods": 40,
    "ledger_run_mode_periods": 30,
    "pump_phase_rad": math.pi / 4,
    "refinement_n_cells": [64, 128, 256, 512],
    "backend": "numpy",
}


def grid(spec):
    nus = np.linspace(*spec["nu_grid"][:2], spec["nu_grid"][2])
    deltas = np.linspace(*spec["delta_grid"][:2], spec["delta_grid"][2])
    return nus, deltas


def floquet_maps(spec, n_cells):
    nus, deltas = grid(spec)
    Q = spec["quality"]
    D, NU = np.meshgrid(deltas, nus, indexing="ij")
    mu_fdtd = batch_floquet_exponents(CONFIG["length_m"], n_cells, D.ravel(), NU.ravel(), Q,
                                      pump_periods=CONFIG["floquet_pump_periods"],
                                      courant=CONFIG["courant"], xp=CONFIG["backend"])[:, 0].reshape(D.shape)
    mu_theory = np.array([[floquet_growth_rate(d, nu, Q) if d > 0 else -0.5 / Q for nu in nus] for d in deltas])
    return nus, deltas, mu_fdtd, mu_theory


def ledger_batch(deltas_flat, nus_flat, Q, n_cells):
    L = CONFIG["length_m"]
    w1 = Cavity1D(L).mode_angular_frequency(1)
    sim = FDTD1DBatch(L, n_cells, len(deltas_flat), CONFIG["courant"], sigma=w1 * EPS0 / Q,
                      depth=deltas_flat, angular_frequency=nus_flat * w1, phase=CONFIG["pump_phase_rad"],
                      xp=CONFIG["backend"])
    E, H, _ = discrete_mode_fields(FDTD1D(L, n_cells, CONFIG["courant"]), 1)
    # H^{-1/2} must use the batch dt
    k = math.pi / L
    w_d = 2.0 / sim.dt * math.asin(299792458.0 * sim.dt / sim.dx * math.sin(0.5 * k * sim.dx))
    H = (1.0 / 376.730313668) * np.cos(k * sim.x_h) * math.sin(-0.5 * w_d * sim.dt)
    sim.set_fields(E, H)
    n_steps = round(CONFIG["ledger_run_mode_periods"] * 2 * math.pi / w1 / sim.dt)
    res = sim.run(n_steps, record_every=max(1, n_steps // 300))
    return sim, res


def main() -> None:
    ctx = ExperimentContext("EXP-0005", HERE, CONFIG)
    N = CONFIG["n_cells"]

    # ---- map A: Floquet + ledger + dataset
    A = CONFIG["map_A"]
    nus, deltas, mu_fdtd, mu_theory = floquet_maps(A, N)
    D, NU = np.meshgrid(deltas, nus, indexing="ij")
    sim, res = ledger_batch(D.ravel(), NU.ravel(), A["quality"], N)
    scale = res.scale()
    e_net = res.e_net()
    e_net_naive = res.e_net(naive=True)
    exact_res = (np.abs(e_net).max(0) / scale).reshape(D.shape)
    naive_gain = (e_net_naive[-1] / scale).reshape(D.shape)
    print("map A done")

    rows = []
    for idx, (d, nu) in enumerate(zip(D.ravel(), NU.ravel())):
        i, j = np.unravel_index(idx, D.shape)
        cfg = {"nu": nu, "delta": d, "Q": A["quality"], "n_cells": N, "courant": CONFIG["courant"]}
        rows.append({
            "simulation_id": f"EXP0005-{i:02d}-{j:02d}", "config_hash": config_hash(cfg),
            "geometry": "cavity_1d_pec", "length_m": CONFIG["length_m"], "eps_r": 1.0,
            "quality_factor": A["quality"], "modulation_nu": nu, "modulation_delta": d,
            "modulation_phase_rad": CONFIG["pump_phase_rad"], "solver": res.meta["solver"],
            "n_cells": N, "dt_s": sim.dt, "precision": "float64", "t_end_s": res.t[-1],
            "E_in": res.pump[-1, idx], "E_out": res.port[-1, idx], "E_loss": res.loss[-1, idx],
            "dE_stored": res.W[-1, idx] - res.W[0, idx], "E_net_exact": e_net[-1, idx],
            "rel_residual_exact": exact_res[i, j], "E_net_naive": e_net_naive[-1, idx],
            "floquet_mu_fdtd_per_tau": mu_fdtd[i, j], "floquet_mu_theory_per_tau": mu_theory[i, j],
            "audit_status": "PASS" if exact_res[i, j] < 1e-10 else "FAIL",
        })
    with open(ctx.results_dir / "dataset.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- refinement of the worst naive point
    i_w, j_w = np.unravel_index(np.argmax(naive_gain), naive_gain.shape)
    refinement = []
    for n in CONFIG["refinement_n_cells"]:
        _, r = ledger_batch(np.array([deltas[i_w]]), np.array([nus[j_w]]), A["quality"], n)
        refinement.append({"n_cells": n, "naive_final_rel_e_net": float(r.e_net(True)[-1, 0] / r.scale()[0]),
                           "exact_max_rel": float(np.abs(r.e_net()[:, 0]).max() / r.scale()[0])})
    ref_order = float(-np.polyfit(np.log([x["n_cells"] for x in refinement]),
                                  np.log([abs(x["naive_final_rel_e_net"]) for x in refinement]), 1)[0])

    # ---- map B: high-Q Floquet map (higher-order tongues)
    B = CONFIG["map_B"]
    nus_b, deltas_b, mu_fdtd_b, mu_theory_b = floquet_maps(B, N)
    print("map B done")

    def agreement(a, b):
        return float(np.mean((a > 0) == (b > 0)))

    def tongue_points(mu, nus_, lo, hi):
        sel = (nus_ >= lo) & (nus_ <= hi)
        return int(np.sum(mu[:, sel] > 0))

    ctx.metrics = {
        "map_A": {
            "n_points": int(mu_fdtd.size),
            "stability_classification_agreement": agreement(mu_fdtd, mu_theory),
            "median_abs_floquet_error_per_tau": float(np.median(np.abs(mu_fdtd - mu_theory))),
            "max_abs_floquet_error_per_tau": float(np.abs(mu_fdtd - mu_theory).max()),
            "max_exact_relative_residual": float(exact_res.max()),
            "naive_apparent_gain_max": float(naive_gain.max()),
            "naive_apparent_gain_min": float(naive_gain.min()),
            "naive_apparent_gain_median": float(np.median(naive_gain)),
            "fraction_points_naive_reports_gain": float(np.mean(naive_gain > 0)),
            "worst_naive_point": {"nu": float(nus[j_w]), "delta": float(deltas[i_w])},
            "worst_point_refinement": refinement,
            "worst_point_naive_order": ref_order,
        },
        "map_B": {
            "n_points": int(mu_fdtd_b.size),
            "stability_classification_agreement": agreement(mu_fdtd_b, mu_theory_b),
            "median_abs_floquet_error_per_tau": float(np.median(np.abs(mu_fdtd_b - mu_theory_b))),
            "unstable_points_nu2_tongue_fdtd_theory": [tongue_points(mu_fdtd_b, nus_b, 1.7, 2.4),
                                                       tongue_points(mu_theory_b, nus_b, 1.7, 2.4)],
            "unstable_points_nu1_tongue_fdtd_theory": [tongue_points(mu_fdtd_b, nus_b, 0.85, 1.15),
                                                       tongue_points(mu_theory_b, nus_b, 0.85, 1.15)],
            "unstable_points_nu2over3_tongue_fdtd_theory": [tongue_points(mu_fdtd_b, nus_b, 0.6, 0.72),
                                                            tongue_points(mu_theory_b, nus_b, 0.6, 0.72)],
        },
        "gates": {
            "exact_ledger_all_pass": bool(exact_res.max() < 1e-10),
            "stability_agreement_A_ge_98pct": agreement(mu_fdtd, mu_theory) >= 0.98,
            "stability_agreement_B_ge_98pct": agreement(mu_fdtd_b, mu_theory_b) >= 0.98,
            "naive_gain_is_artifact": ref_order > 0.8,
        },
    }

    plt = setup_matplotlib()
    from matplotlib.colors import SymLogNorm

    fig, ax = plt.subplots(2, 3, figsize=(15, 8.4))
    for row, (nu_, d_, mf, mt, Q) in enumerate(((nus, deltas, mu_fdtd, mu_theory, A["quality"]),
                                                 (nus_b, deltas_b, mu_fdtd_b, mu_theory_b, B["quality"]))):
        ext = [nu_[0], nu_[-1], d_[0], d_[-1]]
        vmax = max(abs(mt).max(), abs(mf).max())
        for col, (data, title) in enumerate(((mt, "Hill theory"), (mf, f"FDTD monodromy (N = {N})"))):
            a = ax[row, col]
            im = a.imshow(data, origin="lower", aspect="auto", extent=ext, cmap="RdBu_r",
                          norm=SymLogNorm(linthresh=1e-4, vmin=-vmax, vmax=vmax))
            a.contour(nu_, d_, data, levels=[0.0], colors="k", linewidths=0.8)
            a.set_xlabel(r"$\nu = \Omega/\omega_1$")
            a.set_ylabel(r"$\delta$")
            a.set_title(f"Floquet rate, Q = {Q:.0f}: {title}")
        fig.colorbar(im, ax=ax[row, 1], label=r"$\mu$ per $\omega_1 t$")
    g = naive_gain
    gl = max(abs(g).max(), 1e-12)
    ext = [nus[0], nus[-1], deltas[0], deltas[-1]]
    im2 = ax[0, 2].imshow(g, origin="lower", aspect="auto", extent=ext, cmap="PuOr_r",
                          norm=SymLogNorm(linthresh=1e-3, vmin=-gl, vmax=gl))
    ax[0, 2].contour(nus, deltas, mu_theory, levels=[0.0], colors="k", linewidths=0.6, linestyles="--")
    ax[0, 2].set_xlabel(r"$\nu = \Omega/\omega_1$")
    ax[0, 2].set_title(f"Naive ledger: apparent $E_{{\\rm net}}$/scale (> 0 at {np.mean(g > 0):.0%})")
    fig.colorbar(im2, ax=ax[0, 2], label=f"exact ledger: max |E_net|/scale = {exact_res.max():.0e}")
    ns = [x["n_cells"] for x in refinement]
    ax[1, 2].loglog(ns, [x["naive_final_rel_e_net"] for x in refinement], "o-",
                    label=f"naive (order {ref_order:.2f})")
    ax[1, 2].loglog(ns, [max(x["exact_max_rel"], 1e-17) for x in refinement], "x:", label="exact")
    ax[1, 2].set_xlabel("cells N")
    ax[1, 2].set_ylabel(r"$E_{\rm net}$ / scale")
    ax[1, 2].set_title(fr"Worst naive point ($\nu$={nus[j_w]:.2f}, $\delta$={deltas[i_w]:.2f}) under refinement")
    ax[1, 2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ctx.figure_path("exp0005_resonance_map.png"))
    ctx.save()
    print({k: v for k, v in ctx.metrics["map_A"].items() if k != "worst_point_refinement"})
    print(ctx.metrics["map_B"])
    print(ctx.metrics["gates"])


if __name__ == "__main__":
    main()
