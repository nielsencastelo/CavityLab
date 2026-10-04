"""EXP-0002: spatial/temporal convergence of the native 1D FDTD solver.

Gaussian pulse in a PEC cavity, exact d'Alembert/image solution as reference
(benchmark B-002). Measures the L2 field error at a fixed physical time for a
range of resolutions and Courant numbers and fits the observed order.
S = 1 is the 1D "magic time step" where Yee is dispersionless.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from cavitylab import C0, FDTD1D
from cavitylab.benchmarks.analytic1d import gaussian, pulse_fields
from cavitylab.experiments.runner import ExperimentContext, setup_matplotlib

HERE = Path(__file__).parent
CONFIG = {
    "length_m": 0.1,
    "pulse_center_frac": 0.4,
    "pulse_width_frac": 0.05,
    "t_end_round_trips": 1.37,
    "n_cells": [50, 100, 200, 400, 800, 1600],
    "courants": [0.25, 0.5, 0.9, 1.0],
}


def run_case(n_cells: int, courant: float) -> dict:
    L = CONFIG["length_m"]
    f = gaussian(CONFIG["pulse_center_frac"] * L, CONFIG["pulse_width_frac"] * L)
    sim = FDTD1D(L, n_cells, courant)
    sim.set_fields(pulse_fields(f, L, sim.x_e, 0.0)[0], pulse_fields(f, L, sim.x_h, -0.5 * sim.dt)[1])
    t_target = CONFIG["t_end_round_trips"] * 2 * L / C0
    n_steps = round(t_target / sim.dt)
    res = sim.run(n_steps, record_every=max(1, n_steps // 200))
    exact = pulse_fields(f, L, sim.x_e, sim.time)[0]
    err = math.sqrt(np.mean((sim.E - exact) ** 2)) / np.max(np.abs(exact))
    return {
        "n_cells": n_cells,
        "courant": courant,
        "dx": sim.dx,
        "dt": sim.dt,
        "n_steps": n_steps,
        "l2_rel_error": err,
        "energy_drift_rel": float(np.ptp(res.W) / res.W[0]),
    }


def main() -> None:
    ctx = ExperimentContext("EXP-0002", HERE, CONFIG)
    rows = [run_case(n, s) for s in CONFIG["courants"] for n in CONFIG["n_cells"]]

    orders = {}
    for s in CONFIG["courants"]:
        sel = [r for r in rows if r["courant"] == s]
        if s == 1.0:
            continue
        dx = np.log([r["dx"] for r in sel[2:]])
        er = np.log([r["l2_rel_error"] for r in sel[2:]])
        orders[str(s)] = float(np.polyfit(dx, er, 1)[0])

    ctx.metrics = {
        "rows": rows,
        "observed_order_by_courant": orders,
        "max_error_magic_step": max(r["l2_rel_error"] for r in rows if r["courant"] == 1.0),
        "max_energy_drift_rel": max(r["energy_drift_rel"] for r in rows),
        "gate_second_order": all(1.9 < p < 2.1 for p in orders.values()),
    }

    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(5.2, 4))
    for s in CONFIG["courants"]:
        sel = [r for r in rows if r["courant"] == s]
        lab = f"S = {s}" + (f" (p = {orders[str(s)]:.2f})" if str(s) in orders else " (magic step)")
        ax.loglog([2 * r["n_cells"] for r in sel], [max(r["l2_rel_error"], 1e-17) for r in sel], "o-", label=lab)
    x = np.array([100, 3200])
    ax.loglog(x, 3e1 * (x / 100.0) ** -2, "k--", lw=1, label="slope -2")
    ax.set_xlabel("cells per cavity round trip (2N)")
    ax.set_ylabel("relative L2 error of $E_z$")
    ax.set_title("EXP-0002: FDTD convergence (Gaussian pulse)")
    ax.legend(fontsize=8)
    fig.savefig(ctx.figure_path("exp0002_convergence.png"))

    ctx.save()
    print("EXP-0002 done: orders", orders)


if __name__ == "__main__":
    main()
