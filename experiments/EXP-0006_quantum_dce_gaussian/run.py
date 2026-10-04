"""EXP-0006 (part 1): minimal quantum DCE model with an energy ledger and a work audit.

Model: one cavity mode with modulated permittivity, quantized exactly as a Gaussian
state (hbar = 1, tau = omega_0 t). See cavitylab.quantum.gaussian.

Checks:
  1. Lossless resonant pumping from vacuum: photons N(tau) vs the RWA prediction
     sinh^2(delta tau / 4); the state stays pure, so all excitation is ergotropy.
  2. Quantum-classical bridge: the transfer matrix measured from the *classical*
     FDTD solver predicts the quantum photon number from vacuum (exact for linear
     lossless media) and converges with resolution.
  3. Lossy cavity (Q) at T = 0 and T > 0: exact energy ledger (pump work = dE +
     energy to bath), photons, ergotropy, passive energy, entropy, and the ratio
     ergotropy / pump work. This is the first quantum-work-audit result (EXP-0013
     preview).

This is the minimal reproduction target. Part 2 (to do) reproduces a published
circuit-DCE model (Wilson et al. 2011 parameters, two-mode squeezing).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from cavitylab.benchmarks.parametric import fdtd_monodromy
from cavitylab.experiments.runner import ExperimentContext, setup_matplotlib
from cavitylab.quantum import photons_from_monodromy, simulate_parametric_mode

HERE = Path(__file__).parent
CONFIG = {
    "delta": 0.02,
    "nu": 2.0,
    "tau_end_lossless": 191 * math.pi,  # stroboscopic: f = 1 at tau = k pi/2
    "bridge_length_m": 0.1,
    "bridge_pump_periods": 200,
    "bridge_n_cells": [50, 100, 200, 400],
    "lossy_qualities": [400.0, 200.0, 100.0],
    "lossy_n_th": [0.0, 1.0],
    "tau_end_lossy": 637 * math.pi,
}


def main() -> None:
    ctx = ExperimentContext("EXP-0006", HERE, CONFIG)
    plt = setup_matplotlib()
    d, nu = CONFIG["delta"], CONFIG["nu"]

    # 1. lossless from vacuum
    r0 = simulate_parametric_mode(CONFIG["tau_end_lossless"], d, nu, n_eval=1201)
    rwa = np.sinh(d * r0.tau / 4) ** 2
    lossless = {
        "photons_final": float(r0.photons[-1]),
        "rwa_sinh2_final": float(rwa[-1]),
        "rel_diff_vs_rwa": float(r0.photons[-1] / rwa[-1] - 1),
        "entropy_max": float(r0.entropy.max()),
        "ergotropy_over_photons_final": float(r0.ergotropy[-1] / r0.photons[-1]),
        "ledger_residual_max": float(np.abs(r0.ledger_residual).max()),
        "pump_work_final": float(r0.pump_work[-1]),
    }

    # 2. classical FDTD transfer matrix -> quantum photon number
    bridge = []
    for n in CONFIG["bridge_n_cells"]:
        M, tau = fdtd_monodromy(CONFIG["bridge_length_m"], n, d, nu, math.inf,
                                pump_periods=CONFIG["bridge_pump_periods"])
        M_qp = M[::-1, ::-1]  # FDTD state is (p, q); quantum ordering is (q, p)
        n_fdtd = photons_from_monodromy(M_qp)
        n_q = simulate_parametric_mode(tau, d, nu, n_eval=3).photons[-1]
        bridge.append({"n_cells": n, "tau": tau, "photons_fdtd_map": n_fdtd, "photons_quantum": float(n_q),
                       "rel_err": float(abs(n_fdtd / n_q - 1)), "det_M": float(np.linalg.det(M))})
    order = float(-np.polyfit(np.log([b["n_cells"] for b in bridge[-3:]]),
                              np.log([b["rel_err"] for b in bridge[-3:]]), 1)[0])

    # 3. lossy work audit
    lossy = []
    curves = {}
    for n_th in CONFIG["lossy_n_th"]:
        for Q in CONFIG["lossy_qualities"]:
            r = simulate_parametric_mode(CONFIG["tau_end_lossy"], d, nu, quality=Q, n_th=n_th,
                                         V0=(n_th + 0.5) * np.eye(2), n_eval=2001)
            k = 1.0 / Q
            lossy.append({
                "Q": Q, "n_th": n_th, "delta_over_threshold": d / (2 / Q),
                "photons_final": float(r.photons[-1]),
                "photons_initial": float(r.photons[0]),
                "ergotropy_final": float(r.ergotropy[-1]),
                "passive_energy_final": float(r.passive_energy[-1]),
                "entropy_final": float(r.entropy[-1]),
                "pump_work": float(r.pump_work[-1]),
                "energy_to_bath": float(r.bath_energy_out[-1]),
                "delta_energy": float(r.energy[-1] - r.energy[0]),
                "ergotropy_over_pump_work": float(r.ergotropy[-1] / r.pump_work[-1]),
                "ledger_residual_max_rel": float(np.abs(r.ledger_residual).max() / max(r.pump_work.max(), 1e-30)),
                "kappa_per_tau": k,
            })
            curves[(n_th, Q)] = r

    ctx.metrics = {
        "lossless": lossless,
        "classical_fdtd_bridge": bridge,
        "bridge_convergence_order": order,
        "lossy_work_audit": lossy,
        "gates": {
            "rwa_agreement_1pct": abs(lossless["rel_diff_vs_rwa"]) < 0.01,
            # purity: all excitation is extractable (entropy ~1e-4 here is det(V) integration error)
            "pure_state_lossless": lossless["ergotropy_over_photons_final"] > 0.9999,
            "ledger_closes": lossless["ledger_residual_max"] < 1e-8
            and all(x["ledger_residual_max_rel"] < 1e-8 for x in lossy),
            "bridge_second_order": 1.8 < order < 2.2,
            "ergotropy_never_exceeds_pump_work": all(x["ergotropy_final"] <= x["pump_work"] for x in lossy),
        },
    }

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].semilogy(r0.tau, np.maximum(r0.photons, 1e-12), label="Gaussian (exact)")
    ax[0].semilogy(r0.tau, np.maximum(rwa, 1e-12), "--", label=r"RWA $\sinh^2(\delta\tau/4)$")
    ax[0].semilogy(r0.tau, np.maximum(r0.pump_work, 1e-12), ":", label="pump work [$\\hbar\\omega_0$]")
    ax[0].set_xlabel(r"$\tau = \omega_0 t$")
    ax[0].set_title(fr"Lossless DCE from vacuum, $\delta$ = {d}, $\Omega = 2\omega_0$")
    ax[0].legend(fontsize=8)

    ns = [b["n_cells"] for b in bridge]
    ax[1].loglog(ns, [b["rel_err"] for b in bridge], "o-", label=f"|N_FDTD/N_quantum - 1| (p = {order:.2f})")
    ax[1].set_xlabel("FDTD cells N")
    ax[1].set_title(f"Classical FDTD map predicts photons from vacuum\n(N_quantum = {bridge[-1]['photons_quantum']:.1f})")
    ax[1].legend(fontsize=8)

    for (n_th, Q), r in curves.items():
        if n_th != 0.0:
            continue
        line, = ax[2].plot(r.tau, r.pump_work, lw=1, label=f"pump work, Q={Q:.0f}")
        ax[2].plot(r.tau, r.ergotropy, "--", color=line.get_color(), lw=1)
        ax[2].plot(r.tau, r.bath_energy_out, ":", color=line.get_color(), lw=1)
    ax[2].set_yscale("log")
    ax[2].set_xlabel(r"$\tau$")
    ax[2].set_ylabel(r"energy [$\hbar\omega_0$]")
    ax[2].set_title("T = 0 work audit: solid pump work, dashed ergotropy, dotted to bath")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(ctx.figure_path("exp0006_quantum_dce.png"))

    ctx.save()
    print("lossless:", lossless)
    print("bridge:", [(b["n_cells"], round(b["rel_err"], 8)) for b in bridge], "order", order)
    for x in lossy:
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in x.items()
               if k in ("Q", "n_th", "photons_final", "ergotropy_final", "pump_work", "ergotropy_over_pump_work")})
    print("gates:", ctx.metrics["gates"])


if __name__ == "__main__":
    main()
