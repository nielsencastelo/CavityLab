"""EXP-0003: global energy balance in static (unmodulated) cavities.

Scenarios (benchmark B-003):
  A  closed, lossless, multi-mode pulse, long run          -> stored energy constant
  B  uniform conductivity, free decay                       -> Q matches eps*omega/sigma
  C  lossy cavity driven on resonance by a soft source      -> E_in = E_loss + dW
  D  open cavity: partial mirror + graded absorber (port)   -> E_out booked as output

For every scenario both ledgers are audited:
  exact  discrete-consistent energy and dissipation (cavitylab default)
  naive  textbook estimators (time-averaged H, left-point dissipation/source work)
and the naive residual is measured vs resolution to show that it is a pure
discretization artifact (O(dx^2) at fixed Courant number).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from cavitylab import C0, EPS0, EnergyAudit, FDTD1D, SoftSource
from cavitylab.benchmarks.analytic1d import discrete_mode_fields, gaussian, pulse_fields
from cavitylab.experiments.runner import ExperimentContext, setup_matplotlib
from cavitylab.geometry import Cavity1D

HERE = Path(__file__).parent
CONFIG = {
    "length_m": 0.1,
    "courant": 0.5,
    "n_cells": 400,
    "A_round_trips": 100,
    "B_quality": 50.0,
    "B_mode": 2,
    "B_periods": 60,
    "C_quality": 100.0,
    "C_source_node_frac": 0.23,
    "C_drive_amplitude_A_m2": 1.0,
    "C_periods": 300,
    "D_domain_m": 0.25,
    "D_mirror_eps_r": 60.0,
    "D_mirror_cells": 4,
    "D_absorber_start_m": 0.16,
    "D_absorber_sigma_max": 0.8,
    "D_round_trips": 40,
    "naive_scaling_n_cells": [100, 200, 400, 800, 1600],
    "tolerance_exact": 1e-12,
}


def pulse_init(sim, length, center=0.4, width=1 / 30):
    f = gaussian(center * length, width * length)
    return pulse_fields(f, length, sim.x_e, 0.0)[0], pulse_fields(f, length, sim.x_h, -0.5 * sim.dt)[1]


def scenario_a(n_cells):
    L, S = CONFIG["length_m"], CONFIG["courant"]
    sim = FDTD1D(L, n_cells, S)
    sim.set_fields(*pulse_init(sim, L))
    n_steps = round(CONFIG["A_round_trips"] * 2 * L / C0 / sim.dt)
    return sim, sim.run(n_steps, record_every=max(1, n_steps // 4000))


def scenario_b(n_cells):
    L, S = CONFIG["length_m"], CONFIG["courant"]
    cav = Cavity1D(L)
    n = CONFIG["B_mode"]
    w = cav.mode_angular_frequency(n)
    sigma = w * EPS0 / CONFIG["B_quality"]
    sim = FDTD1D(L, n_cells, S, sigma=sigma)
    E, H, _ = discrete_mode_fields(sim, n)
    sim.set_fields(E, H)
    n_steps = round(CONFIG["B_periods"] * 2 * math.pi / w / sim.dt)
    return sim, sim.run(n_steps, record_every=max(1, n_steps // 4000))


def scenario_c(n_cells):
    L, S = CONFIG["length_m"], CONFIG["courant"]
    w = Cavity1D(L).mode_angular_frequency(1)
    sigma = w * EPS0 / CONFIG["C_quality"]
    node = round(CONFIG["C_source_node_frac"] * n_cells)
    amp = CONFIG["C_drive_amplitude_A_m2"]
    sim = FDTD1D(L, n_cells, S, sigma=sigma, sources=[SoftSource(node, lambda t: amp * math.sin(w * t))])
    n_steps = round(CONFIG["C_periods"] * 2 * math.pi / w / sim.dt)
    return sim, sim.run(n_steps, record_every=max(1, n_steps // 4000))


def scenario_d(n_cells_cavity):
    L, S = CONFIG["length_m"], CONFIG["courant"]
    Ld = CONFIG["D_domain_m"]
    n_cells = round(n_cells_cavity * Ld / L)
    probe = FDTD1D(Ld, n_cells, S)
    x = probe.x_e
    eps_r = np.ones_like(x)
    mirror = (x > L) & (x <= L + CONFIG["D_mirror_cells"] * probe.dx + 1e-12)
    eps_r[mirror] = CONFIG["D_mirror_eps_r"]
    xa = CONFIG["D_absorber_start_m"]
    port = x >= xa
    sigma = np.where(port, CONFIG["D_absorber_sigma_max"] * ((x - xa) / (Ld - xa)) ** 3, 0.0)
    sim = FDTD1D(Ld, n_cells, S, eps_r=eps_r, sigma=sigma, port_mask=port, cavity_length=L)
    E0, Hm = pulse_init(probe, L)
    E0[x > L] = 0.0
    Hm[probe.x_h > L] = 0.0
    sim.set_fields(E0, Hm)
    n_steps = round(CONFIG["D_round_trips"] * 2 * L / C0 / sim.dt)
    return sim, sim.run(n_steps, record_every=max(1, n_steps // 4000))


SCENARIOS = {"A_closed_lossless": scenario_a, "B_free_decay": scenario_b,
             "C_driven_lossy": scenario_c, "D_open_port": scenario_d}


def main() -> None:
    ctx = ExperimentContext("EXP-0003", HERE, CONFIG)
    tol = CONFIG["tolerance_exact"]
    out, results = {}, {}
    for name, fn in SCENARIOS.items():
        sim, res = fn(CONFIG["n_cells"])
        results[name] = res
        exact = EnergyAudit.from_fdtd(res, tolerance=tol)
        naive = EnergyAudit.from_fdtd(res, naive=True, tolerance=tol)
        out[name] = {"exact": exact.summary(), "naive": naive.summary(), "dt": sim.dt, "dx": sim.dx}
        print(name, "\n" + exact.report())

    # Q measured from the free-decay envelope of the exact stored energy
    rb = results["B_free_decay"]
    gamma_fit = -np.polyfit(rb.t, np.log(rb.W), 1)[0]
    w_b = Cavity1D(CONFIG["length_m"]).mode_angular_frequency(CONFIG["B_mode"])
    out["B_free_decay"]["Q_measured"] = w_b / gamma_fit
    out["B_free_decay"]["Q_analytic"] = CONFIG["B_quality"]

    # naive residual vs resolution
    scaling = {}
    for name in ("A_closed_lossless", "C_driven_lossy", "D_open_port"):
        scaling[name] = []
        for n in CONFIG["naive_scaling_n_cells"]:
            _, res = SCENARIOS[name](n)
            scaling[name].append({
                "n_cells": n,
                "exact_max_rel": EnergyAudit.from_fdtd(res).max_relative_residual,
                "naive_max_rel": EnergyAudit.from_fdtd(res, naive=True).max_relative_residual,
                "naive_final_rel_e_net": EnergyAudit.from_fdtd(res, naive=True).final_relative_e_net,
            })
    orders = {}
    for name, rows in scaling.items():
        n = np.log([r["n_cells"] for r in rows])
        e = np.log([r["naive_max_rel"] for r in rows])
        orders[name] = float(-np.polyfit(n, e, 1)[0])

    ctx.metrics = {
        "scenarios": out,
        "naive_scaling": scaling,
        "naive_residual_order_in_n": orders,
        "gate_all_exact_pass": all(v["exact"]["status"] == "PASS" for v in out.values()),
        "gate_Q_within_1pct": abs(out["B_free_decay"]["Q_measured"] / CONFIG["B_quality"] - 1) < 0.01,
    }

    plt = setup_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for ax, (name, res) in zip(axes.flat, results.items()):
        ex = EnergyAudit.from_fdtd(res)
        nv = EnergyAudit.from_fdtd(res, naive=True)
        tn = res.t * 1e9
        ax.plot(tn, nv.e_net / nv.scale, color="tab:red", lw=0.8, label="naive ledger")
        ax.plot(tn, ex.e_net / ex.scale, color="tab:blue", lw=1.2, label="exact discrete ledger")
        ax.set_title(name.replace("_", " "))
        ax.set_xlabel("t [ns]")
        ax.set_ylabel(r"$E_{\rm net}$ / energy scale")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("EXP-0003: energy residual E_net = E_out + E_loss + dW - E_in (N = %d, S = %.2f)"
                 % (CONFIG["n_cells"], CONFIG["courant"]))
    fig.tight_layout()
    fig.savefig(ctx.figure_path("exp0003_energy_residuals.png"))

    fig, ax = plt.subplots(figsize=(5.2, 4))
    for name, rows in scaling.items():
        ns = [r["n_cells"] for r in rows]
        ax.loglog(ns, [r["naive_max_rel"] for r in rows], "o-", label=f"naive, {name} (p={orders[name]:.2f})")
        ax.loglog(ns, [max(r["exact_max_rel"], 1e-17) for r in rows], "x:", color="gray")
    ax.loglog([], [], "x:", color="gray", label="exact (all scenarios)")
    ax.set_xlabel("cells in cavity N")
    ax.set_ylabel(r"max $|E_{\rm net}|$ / energy scale")
    ax.set_title("Naive-ledger residual is a discretization artifact")
    ax.legend(fontsize=7)
    fig.savefig(ctx.figure_path("exp0003_naive_residual_scaling.png"))

    ctx.save()
    print("EXP-0003 done: naive orders", orders, "Q", out["B_free_decay"]["Q_measured"])


if __name__ == "__main__":
    main()
