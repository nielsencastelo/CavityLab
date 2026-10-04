"""EXP-0004: classical parametric modulation of a cavity with explicit pump-work accounting.

Part A  uniform eps(t) = eps0 (1 + delta sin(2 omega_1 t)), lossy (Q = 200).
        Reduces exactly to a damped Hill equation per mode (benchmark B-004):
        FDTD stored energy, pump work and loss are compared with the ODE, and the
        measured growth rate with the Floquet exponent. Threshold delta_th = 2/Q.
Part B  localized modulated slab (multimode parametric coupling), lossless,
        seeded with a deterministic multimode classical field (NOT vacuum
        fluctuations: classical fields need a seed). Pump work must equal the
        growth of stored energy exactly.
Part C  naive-ledger artifact under pumping: apparent E_net vs resolution.

Energy that appears in the field is always booked against pump work. Any
E_net != 0 in the exact ledger would be a bug, not physics.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from cavitylab import EPS0, MU0, EnergyAudit, FDTD1D, Modulation
from cavitylab.benchmarks.analytic1d import discrete_mode_fields
from cavitylab.benchmarks.hill import floquet_growth_rate, integrate_mode
from cavitylab.experiments.runner import ExperimentContext, setup_matplotlib
from cavitylab.geometry import Cavity1D

HERE = Path(__file__).parent
CONFIG = {
    "length_m": 0.1,
    "courant": 0.5,
    "A_quality": 200.0,
    "A_deltas": [0.005, 0.01, 0.02, 0.04],
    "A_nu": 2.0,
    "A_n_cells": 200,
    "A_periods": 150,
    "A_convergence_n_cells": [50, 100, 200, 400, 800],
    "A_convergence_delta": 0.04,
    "A_convergence_periods": 60,
    "B_n_cells": 400,
    "B_slab_m": [0.05, 0.075],
    "B_delta": 0.3,
    "B_periods": 40,
    "B_seed": 20261004,
    "B_seed_modes": 30,
    "B_seed_amplitude_V_m": 1e-3,
    "B_tracked_modes": 40,
    "C_n_cells": [100, 200, 400, 800],
}


def mode_omega(n=1):
    return Cavity1D(CONFIG["length_m"]).mode_angular_frequency(n)


def run_uniform(n_cells, delta, periods, record_target=3000):
    L, S = CONFIG["length_m"], CONFIG["courant"]
    w0 = mode_omega(1)
    sigma = w0 * EPS0 / CONFIG["A_quality"]
    sim = FDTD1D(L, n_cells, S, sigma=sigma, modulation=Modulation(delta, CONFIG["A_nu"] * w0))
    E, H, _ = discrete_mode_fields(sim, 1)
    sim.set_fields(E, H)
    n_steps = round(periods * 2 * math.pi / w0 / sim.dt)
    rec = max(1, n_steps // record_target)
    n_steps -= n_steps % rec
    return sim, sim.run(n_steps, record_every=rec)


def energy_growth_rate(t, W):
    half = len(t) // 2
    return float(np.polyfit(t[half:], np.log(W[half:]), 1)[0])


def part_a(ctx, plt):
    w0 = mode_omega(1)
    Q = CONFIG["A_quality"]
    rows, curves = [], {}
    for delta in CONFIG["A_deltas"]:
        sim, res = run_uniform(CONFIG["A_n_cells"], delta, CONFIG["A_periods"])
        sol = integrate_mode(w0 * res.t[-1], delta, CONFIG["A_nu"], Q, n_eval=len(res.t))
        scale_w = res.W[0] / sol["w"][0]
        exact = EnergyAudit.from_fdtd(res)
        naive = EnergyAudit.from_fdtd(res, naive=True)
        rate_floquet = 2.0 * w0 * floquet_growth_rate(delta, CONFIG["A_nu"], Q)
        rows.append({
            "delta": delta,
            "delta_over_threshold": delta / (2.0 / Q),
            "energy_rate_fdtd_per_s": energy_growth_rate(res.t, res.W),
            "energy_rate_floquet_per_s": rate_floquet,
            "W_final_over_W0_fdtd": res.W[-1] / res.W[0],
            "W_final_over_W0_ode": sol["w"][-1] / sol["w"][0],
            "pump_work_fdtd": res.pump[-1],
            "pump_work_ode": sol["w_pump"][-1] * scale_w,
            "loss_fdtd": res.loss[-1],
            "loss_ode": sol["w_loss"][-1] * scale_w,
            "exact": exact.summary(),
            "naive": naive.summary(),
        })
        curves[delta] = (res.t, res.W / res.W[0], sol["w"] / sol["w"][0], res.pump / res.W[0], res.loss / res.W[0])

    conv = []
    for n in CONFIG["A_convergence_n_cells"]:
        sim, res = run_uniform(n, CONFIG["A_convergence_delta"], CONFIG["A_convergence_periods"])
        sol = integrate_mode(w0 * res.t[-1], CONFIG["A_convergence_delta"], CONFIG["A_nu"], Q, n_eval=len(res.t))
        ratio_f, ratio_o = res.W[-1] / res.W[0], sol["w"][-1] / sol["w"][0]
        conv.append({
            "n_cells": n,
            "rel_err_final_energy_vs_ode": abs(ratio_f / ratio_o - 1.0),
            "exact_max_rel": EnergyAudit.from_fdtd(res).max_relative_residual,
            "naive_max_rel": EnergyAudit.from_fdtd(res, naive=True).max_relative_residual,
            "naive_final_rel_e_net": EnergyAudit.from_fdtd(res, naive=True).final_relative_e_net,
        })
    order = float(-np.polyfit(np.log([c["n_cells"] for c in conv[-3:]]),
                              np.log([c["rel_err_final_energy_vs_ode"] for c in conv[-3:]]), 1)[0])

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for delta, (t, wf, wo, pump, loss) in curves.items():
        line, = ax[0].semilogy(t * 1e9, wf, lw=1, label=fr"$\delta$ = {delta} ({delta / (2 / Q):.1f}$\delta_{{th}}$)")
        ax[0].semilogy(t * 1e9, wo, "--", color=line.get_color(), lw=0.8)
    ax[0].set_xlabel("t [ns]")
    ax[0].set_ylabel("W(t) / W(0)")
    ax[0].set_title("A: uniform modulation, Q = 200 (solid FDTD, dashed Hill ODE)")
    ax[0].legend(fontsize=7)

    t, wf, wo, pump, loss = curves[max(CONFIG["A_deltas"])]
    ax[1].plot(t * 1e9, wf - 1.0, label=r"$\Delta W$")
    ax[1].plot(t * 1e9, pump, label=r"$E_{\rm pump}$ (work by modulation)")
    ax[1].plot(t * 1e9, loss, label=r"$E_{\rm loss}$")
    ax[1].plot(t * 1e9, pump - loss - (wf - 1.0), "k", lw=0.8, label=r"$E_{\rm in}-E_{\rm loss}-\Delta W$")
    ax[1].set_yscale("symlog", linthresh=1.0)
    ax[1].set_xlabel("t [ns]")
    ax[1].set_ylabel("energy / W(0)")
    ax[1].set_title(fr"Ledger at $\delta$ = {max(CONFIG['A_deltas'])}: field energy is pump work")
    ax[1].legend(fontsize=7)

    ns = [c["n_cells"] for c in conv]
    ax[2].loglog(ns, [c["rel_err_final_energy_vs_ode"] for c in conv], "o-", label=f"FDTD vs Hill ODE (p = {order:.2f})")
    ax[2].loglog(ns, [c["naive_max_rel"] for c in conv], "s-", label="naive ledger max |E_net|/scale")
    ax[2].loglog(ns, [max(c["exact_max_rel"], 1e-17) for c in conv], "x:", label="exact ledger max |E_net|/scale")
    ax[2].set_xlabel("cells N")
    ax[2].set_title(fr"Convergence, $\delta$ = {CONFIG['A_convergence_delta']}")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(ctx.figure_path("exp0004a_uniform_modulation.png"))

    ctx.metrics["A"] = {"rows": rows, "convergence": conv, "fdtd_vs_ode_order": order,
                        "delta_threshold_theory": 2.0 / Q}


def seeded_multimode_fields(sim, n_modes, amplitude, seed):
    rng = np.random.default_rng(seed)
    E = np.zeros_like(sim.x_e)
    H = np.zeros_like(sim.x_h)
    for m in range(1, n_modes + 1):
        a, b = amplitude * rng.standard_normal(2)
        e_m, h_m, w_m = discrete_mode_fields(sim, m)
        # quadrature partner: E ~ sin(kx) sin(w t) at t=0 is zero, H ~ cos(kx) cos(w (-dt/2))
        k = m * math.pi / sim.length
        h_q = -(1.0 / math.sqrt(MU0 / EPS0)) * np.cos(k * sim.x_h) * math.cos(-0.5 * w_m * sim.dt)
        E += a * e_m
        H += a * h_m + b * h_q
    return E, H


def run_slab(n_cells, periods, modes=()):
    L, S = CONFIG["length_m"], CONFIG["courant"]
    w1 = mode_omega(1)
    probe = FDTD1D(L, n_cells, S)
    x0, x1 = CONFIG["B_slab_m"]
    profile = ((probe.x_e >= x0) & (probe.x_e <= x1)).astype(float)
    sim = FDTD1D(L, n_cells, S, modulation=Modulation(CONFIG["B_delta"], 2 * w1, profile=profile))
    sim.set_fields(*seeded_multimode_fields(sim, CONFIG["B_seed_modes"], CONFIG["B_seed_amplitude_V_m"], CONFIG["B_seed"]))
    n_steps = round(periods * 2 * math.pi / w1 / sim.dt)
    rec = max(1, n_steps // 3000)
    n_steps -= n_steps % rec
    return sim, sim.run(n_steps, record_every=rec, modes=modes)


def part_b(ctx, plt):
    L = CONFIG["length_m"]
    modes = list(range(1, CONFIG["B_tracked_modes"] + 1))
    sim, res = run_slab(CONFIG["B_n_cells"], CONFIG["B_periods"], modes)
    exact = EnergyAudit.from_fdtd(res)
    naive = EnergyAudit.from_fdtd(res, naive=True)
    Wm = np.array([(L / 4) * (EPS0 * res.modes[m] ** 2 + MU0 * res.modes_h[m] ** 2) for m in modes])
    # period-average to remove the 2-omega breathing of E/H split
    win = max(1, int(round(2 * math.pi / mode_omega(1) / (res.t[1] - res.t[0]))))
    kernel = np.ones(win) / win
    Wm_avg = np.array([np.convolve(w, kernel, mode="same") for w in Wm])
    final_share = Wm_avg[:, -win] / Wm_avg[:, -win].sum()

    conv = []
    for n in CONFIG["C_n_cells"]:
        _, r = run_slab(n, CONFIG["B_periods"])
        ex, nv = EnergyAudit.from_fdtd(r), EnergyAudit.from_fdtd(r, naive=True)
        conv.append({"n_cells": n, "exact_max_rel": ex.max_relative_residual,
                     "naive_max_rel": nv.max_relative_residual,
                     "naive_final_rel_e_net": nv.final_relative_e_net,
                     "naive_apparent_gain_over_input": float(nv.e_net[-1] / nv.e_in[-1]),
                     "W_final_over_W0": float(r.W[-1] / r.W[0])})

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].semilogy(res.t * 1e9, res.W / res.W[0], label="stored energy W")
    ax[0].semilogy(res.t * 1e9, np.maximum(res.pump, 1e-30) / res.W[0], "--", label="pump work (cumulative)")
    ax[0].set_xlabel("t [ns]")
    ax[0].set_ylabel("energy / W(0)")
    ax[0].set_title(r"B: localized slab, $\delta$ = 0.3, $\Omega = 2\omega_1$, lossless")
    ax[0].legend(fontsize=8)
    for m in (1, 2, 3, 5, 7, 9, 11):
        ax[1].semilogy(res.t * 1e9, np.maximum(Wm_avg[m - 1], 1e-40) / res.W[0], lw=1, label=f"mode {m}")
    ax[1].set_xlabel("t [ns]")
    ax[1].set_ylabel("period-averaged modal energy / W(0)")
    ax[1].set_title("Parametric cascade between modes")
    ax[1].legend(fontsize=7, ncol=2)
    ns = [c["n_cells"] for c in conv]
    ax[2].loglog(ns, [c["naive_max_rel"] for c in conv], "s-", label="naive ledger")
    ax[2].loglog(ns, [max(c["exact_max_rel"], 1e-17) for c in conv], "x:", label="exact ledger")
    ax[2].set_xlabel("cells N")
    ax[2].set_ylabel(r"max $|E_{\rm net}|$ / energy scale")
    ax[2].set_title("C: residual under pumping")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ctx.figure_path("exp0004b_slab_modulation.png"))

    ctx.metrics["B"] = {"exact": exact.summary(), "naive": naive.summary(),
                        "W_final_over_W0": float(res.W[-1] / res.W[0]),
                        "final_modal_energy_share_top10": sorted(
                            [(int(m), float(s)) for m, s in zip(modes, final_share)],
                            key=lambda p: -p[1])[:10]}
    ctx.metrics["C_naive_under_pumping"] = conv


def main() -> None:
    ctx = ExperimentContext("EXP-0004", HERE, CONFIG)
    plt = setup_matplotlib()
    part_a(ctx, plt)
    part_b(ctx, plt)
    A = ctx.metrics["A"]
    ctx.metrics["gates"] = {
        "exact_ledger_all_pass": all(r["exact"]["status"] == "PASS" for r in A["rows"])
        and ctx.metrics["B"]["exact"]["status"] == "PASS",
        "threshold_bracketed": A["rows"][0]["energy_rate_fdtd_per_s"] < 0 < A["rows"][2]["energy_rate_fdtd_per_s"],
        "growth_rate_matches_floquet_5pct": all(
            abs(r["energy_rate_fdtd_per_s"] - r["energy_rate_floquet_per_s"])
            <= 0.05 * abs(r["energy_rate_floquet_per_s"]) + 0.02 * mode_omega(1) / CONFIG["A_quality"]
            for r in A["rows"]),
    }
    ctx.save()
    for r in A["rows"]:
        print(f"delta={r['delta']}: rate FDTD {r['energy_rate_fdtd_per_s']:.4e} vs Floquet "
              f"{r['energy_rate_floquet_per_s']:.4e} 1/s; exact {r['exact']['max_relative_residual']:.1e}, "
              f"naive {r['naive']['max_relative_residual']:.1e}")
    print("FDTD vs ODE order:", A["fdtd_vs_ode_order"])
    print("B:", ctx.metrics["B"]["exact"]["status"], ctx.metrics["B"]["W_final_over_W0"])
    print("gates:", ctx.metrics["gates"])


if __name__ == "__main__":
    main()
