"""EXP-0001: analytical modes of a 1D PEC cavity vs native FDTD.

Checks (benchmark B-001):
  1. Exact Yee eigenmodes oscillate at the Yee-dispersion frequency to round-off.
  2. Continuum-initialized modes: relative frequency error vs analytical omega_n,
     compared with the leading-order prediction -(k dx)^2 (1 - S^2) / 24.
  3. Discrete stored energy of a mode vs analytical eps A^2 L / 4.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from cavitylab import EPS0, FDTD1D
from cavitylab.benchmarks.analytic1d import discrete_mode_fields, estimate_frequency
from cavitylab.experiments.runner import ExperimentContext, setup_matplotlib
from cavitylab.geometry import Cavity1D

HERE = Path(__file__).parent
CONFIG = {
    "length_m": 0.1,
    "eps_r": 1.0,
    "modes": list(range(1, 11)),
    "n_cells": [50, 100, 200, 400],
    "courant": 0.5,
    "periods_measured": 20,
}


def main() -> None:
    ctx = ExperimentContext("EXP-0001", HERE, CONFIG)
    cav = Cavity1D(CONFIG["length_m"], CONFIG["eps_r"])
    S = CONFIG["courant"]

    rows = []
    for n_cells in CONFIG["n_cells"]:
        for n in CONFIG["modes"]:
            sim = FDTD1D(cav.length, n_cells, S)
            w_exact = cav.mode_angular_frequency(n)
            k = cav.wavenumber(n)

            # (1) exact discrete eigenmode
            E, H, w_yee = discrete_mode_fields(sim, n)
            sim.set_fields(E, H)
            n_steps = round(CONFIG["periods_measured"] * 2 * math.pi / w_yee / sim.dt)
            res = sim.run(n_steps, modes=[n])
            w_meas_discrete = estimate_frequency(res.modes[n], sim.dt)
            w_analytic_energy = EPS0 * cav.eps_r * cav.length / 4.0  # A = 1

            rows.append({
                "n_cells": n_cells,
                "mode": n,
                "points_per_wavelength": 2 * n_cells / n,
                "omega_exact": w_exact,
                "omega_yee": w_yee,
                "omega_measured": w_meas_discrete,
                "measured_vs_yee_rel": w_meas_discrete / w_yee - 1.0,
                "yee_vs_exact_rel": w_yee / w_exact - 1.0,
                "leading_order_prediction_rel": -((k * sim.dx) ** 2) * (1 - S**2) / 24.0,
                "stored_energy_discrete": float(np.mean(res.W)),
                "stored_energy_analytic": w_analytic_energy,
                "stored_energy_rel_err": float(np.mean(res.W)) / w_analytic_energy - 1.0,
                "energy_drift_rel": float(np.ptp(res.W) / res.W[0]),
            })

    worst_meas = max(abs(r["measured_vs_yee_rel"]) for r in rows)
    worst_drift = max(r["energy_drift_rel"] for r in rows)
    ctx.metrics = {
        "rows": rows,
        "max_abs_measured_vs_yee_rel": worst_meas,
        "max_energy_drift_rel": worst_drift,
        "gate_frequency_vs_yee": worst_meas < 1e-10,
        "gate_energy_drift": worst_drift < 1e-12,
    }

    plt = setup_matplotlib()
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for n_cells in CONFIG["n_cells"]:
        sel = [r for r in rows if r["n_cells"] == n_cells]
        ppw = [r["points_per_wavelength"] for r in sel]
        ax[0].loglog(ppw, [abs(r["yee_vs_exact_rel"]) for r in sel], "o", label=f"N={n_cells}")
    ppw_line = np.logspace(1, 3.2, 50)
    pred = (2 * math.pi / ppw_line) ** 2 * (1 - S**2) / 24
    ax[0].loglog(ppw_line, pred, "k--", lw=1, label=r"$(k\Delta x)^2(1-S^2)/24$")
    ax[0].set_xlabel("points per wavelength")
    ax[0].set_ylabel(r"$|\omega_{\rm FDTD}/\omega_n - 1|$")
    ax[0].set_title("Mode frequency error (S = 0.5)")
    ax[0].legend(fontsize=8)

    ax[1].semilogy([r["mode"] for r in rows if r["n_cells"] == 200],
                   [max(abs(r["measured_vs_yee_rel"]), 1e-17) for r in rows if r["n_cells"] == 200],
                   "s-", label="measured vs Yee dispersion")
    ax[1].semilogy([r["mode"] for r in rows if r["n_cells"] == 200],
                   [max(r["energy_drift_rel"], 1e-17) for r in rows if r["n_cells"] == 200],
                   "^-", label="discrete energy drift")
    ax[1].set_xlabel("mode number n (N = 200)")
    ax[1].set_ylabel("relative value")
    ax[1].set_title("Round-off-level agreement")
    ax[1].legend(fontsize=8)
    fig.savefig(ctx.figure_path("exp0001_mode_frequencies.png"))

    ctx.save()
    print(f"EXP-0001 done: max |w_meas/w_yee-1| = {worst_meas:.2e}, max drift = {worst_drift:.2e}")


if __name__ == "__main__":
    main()
