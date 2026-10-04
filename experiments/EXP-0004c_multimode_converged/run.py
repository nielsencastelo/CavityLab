"""EXP-0004c: converged multimode parametric benchmark against an independent solver.

EXP-0004 part B showed that an equidistant spectrum (empty 1D cavity) gives a
non-convergent parametric cascade. Here a smooth dielectric loading breaks the
equidistance, and a smooth localized modulation pumps at Omega = 2 omega_1:

    eps_s(x) = 1 + 2 * (1 + tanh((x/L - 0.6)/0.05)) / 2        (1 -> 3)
    eps(x,t) = eps_s(x) * (1 + delta * exp(-((x/L - 0.3)/0.05)^2) * sin(Omega t + phi))

Initial field: E = sin(pi x / L), H = 0 (multimode content). Lossless.

Solvers:
  A. native Yee/leapfrog FDTD (2nd order), exact energy ledger
  B. independent method of lines: 4th-order staggered FD + adaptive DOP853
     (cavitylab.solvers.mol1d)
omega_1 comes from the eigenvalues of a fine-grid (N = 4000) operator.

Gate: FDTD stored energy W(t)/W(0) converges (2nd order) to the MOL reference,
and the MOL reference is self-converged far below the FDTD error.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

from cavitylab import C0, MU0, EnergyAudit, FDTD1D, Modulation
from cavitylab.experiments.runner import ExperimentContext, setup_matplotlib
from cavitylab.solvers.mol1d import run_mol

HERE = Path(__file__).parent
CONFIG = {
    "length_m": 0.1,
    "courant": 0.5,
    "eps_contrast": 2.0,
    "eps_step_center": 0.6,
    "eps_step_width": 0.05,
    "mod_center": 0.3,
    "mod_width": 0.05,
    "delta": 0.2,
    "pump_phase_rad": math.pi / 4,
    "checkpoints_periods": [5, 10, 15, 20],
    "fdtd_n_cells": [100, 200, 400, 800, 1600],
    "mol_n_cells": [200, 400, 800],
    "eig_n_cells": 4000,
    "n_spectrum": 8,
}


def eps_static(xn):
    c, w, a = CONFIG["eps_step_center"], CONFIG["eps_step_width"], CONFIG["eps_contrast"]
    return 1.0 + a * 0.5 * (1.0 + np.tanh((xn - c) / w))


def mod_profile(xn):
    return np.exp(-(((xn - CONFIG["mod_center"]) / CONFIG["mod_width"]) ** 2))


def spectrum_normalized(n_cells, n_modes):
    """Lowest angular frequencies (units c/L) of the static cavity, Yee operator on a fine grid."""
    dx = 1.0 / n_cells
    xi = np.arange(1, n_cells) * dx
    inv_eps = 1.0 / eps_static(xi)
    # omega^2 H = G diag(1/eps) G^T H / dx^2, G: interior nodes -> half nodes (tridiagonal N x N)
    diag = np.zeros(n_cells)
    diag[:-1] += inv_eps
    diag[1:] += inv_eps
    off = -inv_eps
    # index 0 is the static null mode (H = const, rank-deficient G diag G^T); skip it explicitly
    w2 = eigh_tridiagonal(diag / dx**2, off / dx**2, select="i", select_range=(1, n_modes))[0]
    return np.sqrt(w2)


def main() -> None:
    ctx = ExperimentContext("EXP-0004c", HERE, CONFIG)
    L = CONFIG["length_m"]
    spec = spectrum_normalized(CONFIG["eig_n_cells"], CONFIG["n_spectrum"])
    w1n = spec[0]
    big_omega_n = 2.0 * w1n
    period_n = 2 * math.pi / w1n
    ctx.metrics["spectrum_ratios_omega_k_over_omega_1"] = (spec / w1n).tolist()

    def eps_of(xn, tn):
        return eps_static(xn) * (1.0 + CONFIG["delta"] * mod_profile(xn) * math.sin(big_omega_n * tn + CONFIG["pump_phase_rad"]))

    # ---- FDTD runs, segment by segment so checkpoints fall exactly on time steps
    fdtd = []
    for n in CONFIG["fdtd_n_cells"]:
        probe = FDTD1D(L, n, CONFIG["courant"])
        xn = probe.x_e / L
        sim = FDTD1D(L, n, CONFIG["courant"], eps_r=eps_static(xn),
                     modulation=Modulation(CONFIG["delta"], big_omega_n * C0 / L,
                                           phase=CONFIG["pump_phase_rad"], profile=mod_profile(xn)))
        E0 = np.sin(math.pi * xn)
        H_half = -(0.5 * sim.dt / MU0) * (math.pi / L) * np.cos(math.pi * sim.x_h / L)  # H(-dt/2), H(0) = 0
        sim.set_fields(E0, H_half)
        W0 = None
        times, ratios, pump, resid = [], [], 0.0, 0.0
        done = 0
        for k in CONFIG["checkpoints_periods"]:
            target = round(k * period_n * L / C0 / sim.dt)
            res = sim.run(target - done, record_every=max(1, (target - done) // 2000))
            if W0 is None:
                W0 = res.W[0]
            done = target
            times.append(sim.time * C0 / L)
            ratios.append(res.W[-1] / W0)
            pump += res.pump[-1]
            resid = max(resid, EnergyAudit.from_fdtd(res).max_relative_residual)
        fdtd.append({"n_cells": n, "t_norm": times, "W_ratio": ratios, "exact_max_rel": resid})
        print(f"FDTD N={n} done: W/W0 = {ratios}")

    # ---- MOL reference at the exact FDTD checkpoint times
    all_t = sorted({t for f in fdtd for t in f["t_norm"]})
    mol = []
    for n in CONFIG["mol_n_cells"]:
        x = np.arange(n + 1) / n
        xh = (np.arange(n) + 0.5) / n
        r = run_mol(n, eps_of, np.sin(math.pi * x), np.zeros(n), np.array([0.0] + all_t), rtol=1e-11, atol=1e-14)
        mol.append({"n_cells": n, "t": r["t"], "W": r["W"]})
        print(f"MOL N={n} done ({r['nfev']} rhs evaluations)")
    ref = mol[-1]

    def ref_ratio(t):
        i = int(np.argmin(np.abs(ref["t"] - t)))
        return ref["W"][i] / ref["W"][0]

    for f in fdtd:
        f["ref_ratio"] = [ref_ratio(t) for t in f["t_norm"]]
        f["rel_err_vs_mol"] = [abs(a / b - 1) for a, b in zip(f["W_ratio"], f["ref_ratio"])]
    mol_self = [abs((m["W"][-1] / m["W"][0]) / (ref["W"][-1] / ref["W"][0]) - 1) for m in mol[:-1]]
    errs_final = [f["rel_err_vs_mol"][-1] for f in fdtd]
    order = float(-np.polyfit(np.log(CONFIG["fdtd_n_cells"][-3:]), np.log(errs_final[-3:]), 1)[0])

    ctx.metrics.update({
        "omega1_norm": w1n,
        "W_ratio_final_reference": ref["W"][-1] / ref["W"][0],
        "fdtd": [{k: v for k, v in f.items()} for f in fdtd],
        "mol_self_convergence_final": dict(zip([m["n_cells"] for m in mol[:-1]], mol_self)),
        "fdtd_vs_mol_order_final": order,
        "gates": {
            "fdtd_second_order_vs_independent_solver": 1.7 < order < 2.3,
            "reference_self_converged": mol_self[-1] < 0.1 * min(errs_final),
            "fdtd_finest_within_1e-3": errs_final[-1] < 1e-3,
            "exact_ledger_all_pass": all(f["exact_max_rel"] < 1e-10 for f in fdtd),
        },
    })

    plt = setup_matplotlib()
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    xs = np.linspace(0, 1, 400)
    ax[0].plot(xs, eps_static(xs), label=r"$\varepsilon_s(x)$")
    ax[0].plot(xs, 1 + CONFIG["delta"] * mod_profile(xs), "--", label=r"$1+\delta\,m(x)$ (modulation envelope)")
    ax[0].set_xlabel("x / L")
    ax[0].set_title(r"Non-equidistant cavity: $\omega_k/\omega_1$ = " + ", ".join(f"{r:.2f}" for r in (spec / w1n)[:5]))
    ax[0].legend(fontsize=8)
    for f in fdtd:
        ax[1].plot(CONFIG["checkpoints_periods"], f["W_ratio"], "o-", ms=3, lw=0.8, label=f"FDTD N={f['n_cells']}")
    ax[1].plot(CONFIG["checkpoints_periods"], fdtd[-1]["ref_ratio"], "k*", ms=9, label="MOL 4th-order (N=800)")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("mode-1 periods")
    ax[1].set_ylabel("W(t)/W(0)")
    ax[1].set_title(r"Parametric growth, $\Omega = 2\omega_1$, $\delta$ = 0.2")
    ax[1].legend(fontsize=7)
    ax[2].loglog(CONFIG["fdtd_n_cells"], errs_final, "o-", label=f"FDTD vs MOL (p = {order:.2f})")
    ax[2].axhline(mol_self[-1], color="gray", ls=":", label="MOL self-convergence (N=400 vs 800)")
    ax[2].set_xlabel("FDTD cells N")
    ax[2].set_ylabel("relative error in W(20 periods)")
    ax[2].set_title("Independent-solver validation")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ctx.figure_path("exp0004c_multimode_converged.png"))
    ctx.save()
    print("order", order, "errors", errs_final, "mol self", mol_self, "gates", ctx.metrics["gates"])


if __name__ == "__main__":
    main()
