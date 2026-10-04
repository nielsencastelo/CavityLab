"""EXP-0003b: floating-point precision of the energy ledger, and CPU/GPU acceleration.

Questions:
  1. How does the exact-ledger residual grow with step count in float64 and float32?
     (The ledger is exact algebraically, so its residual is pure round-off.)
  2. Is float32 still able to separate physics from artifacts? (Compare its floor
     with the naive-ledger artifact level, about 1e-5 to 1e-2.)
  3. Does a float32 run reproduce the Floquet growth rate?
  4. Throughput of the batched solver: numpy (CPU) vs torch (CUDA GPU), float64 vs float32.

Runs on any backend available; GPU parts are skipped without CUDA.
Recommended interpreter: .venv-gpu/Scripts/python.exe (torch + CUDA, no MKL/OpenMP clash).
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np

from cavitylab import EPS0, ETA0
from cavitylab.benchmarks.hill import floquet_growth_rate
from cavitylab.experiments.runner import ExperimentContext, setup_matplotlib
from cavitylab.solvers.fdtd1d_batch import FDTD1DBatch, gpu_available

HERE = Path(__file__).parent
CONFIG = {
    "length_m": 0.1,
    "courant": 0.5,
    "n_cells": 200,
    "quality": 200.0,
    "deltas": [0.0, 0.005, 0.02, 0.04],
    "nu": 2.0,
    "mode_periods": 400,
    "record_points": 400,
    "throughput_batches": [16, 256, 2048, 8192],
    "throughput_n_cells": 128,
    "throughput_steps": 2000,
}


def mode1_fields(sim):
    k = math.pi / sim.length
    eta = ETA0
    w_d = 2.0 / sim.dt * math.asin(299792458.0 * sim.dt / sim.dx * math.sin(0.5 * k * sim.dx))
    E = np.sin(k * sim.x_e)
    H = (1.0 / eta) * np.cos(k * sim.x_h) * math.sin(-0.5 * w_d * sim.dt)
    return E, H, w_d


def accuracy_run(xp, dtype):
    L, Q = CONFIG["length_m"], CONFIG["quality"]
    w1 = math.pi * 299792458.0 / L
    deltas = np.array(CONFIG["deltas"])
    sim = FDTD1DBatch(L, CONFIG["n_cells"], len(deltas), CONFIG["courant"], sigma=w1 * EPS0 / Q,
                      depth=deltas, angular_frequency=CONFIG["nu"] * w1, phase=math.pi / 4, xp=xp, dtype=dtype)
    E, H, _ = mode1_fields(sim)
    # amplitude 1e3 V/m keeps energies O(1e-7) J/m^2: far from float32 under/overflow
    sim.set_fields(1e3 * E, 1e3 * H)
    n_steps = round(CONFIG["mode_periods"] * 2 * math.pi / w1 / sim.dt)
    rec = max(1, n_steps // CONFIG["record_points"])
    n_steps -= n_steps % rec
    res = sim.run(n_steps, record_every=rec)
    rel = np.abs(res.e_net()) / res.scale()
    rel_naive = np.abs(res.e_net(naive=True)) / res.scale()
    half = len(res.t) // 2
    rates = [float(np.polyfit(res.t[half:], np.log(res.W[half:, i].astype(float)), 1)[0]) for i in range(len(deltas))]
    theory = [2 * w1 * floquet_growth_rate(d, CONFIG["nu"], Q) if d > 0 else -w1 / Q for d in deltas]
    return {
        "steps": res.t / sim.dt,
        "rel_exact": rel,
        "rel_naive": rel_naive,
        "max_rel_exact": rel.max(0).tolist(),
        "max_rel_naive": rel_naive.max(0).tolist(),
        "energy_rate_series": rates,
        "energy_rate_theory_dominant": theory,
        "n_steps": n_steps,
    }


def throughput(xp, dtype, batch):
    L = CONFIG["length_m"]
    w1 = math.pi * 299792458.0 / L
    sim = FDTD1DBatch(L, CONFIG["throughput_n_cells"], batch, CONFIG["courant"], sigma=w1 * EPS0 / 100,
                      depth=np.linspace(0, 0.3, batch), angular_frequency=2 * w1, xp=xp, dtype=dtype)
    E, H, _ = mode1_fields(sim)
    sim.set_fields(E, H)
    sim.run(50, record_every=50)  # warm-up (CUDA context, kernels)
    t0 = time.perf_counter()
    sim.run(CONFIG["throughput_steps"], record_every=CONFIG["throughput_steps"])
    if xp == "torch":
        import torch

        torch.cuda.synchronize()
    dt = time.perf_counter() - t0
    cell_updates = batch * CONFIG["throughput_n_cells"] * CONFIG["throughput_steps"]
    return {"backend": xp, "dtype": dtype, "batch": batch, "seconds": dt, "Mcell_updates_per_s": cell_updates / dt / 1e6}


def main() -> None:
    ctx = ExperimentContext("EXP-0003b", HERE, CONFIG)
    backends = ["numpy"] + (["torch"] if gpu_available() else [])
    acc = {}
    for xp in backends:
        for dtype in ("float64", "float32"):
            acc[f"{xp}-{dtype}"] = accuracy_run(xp, dtype)
            print(xp, dtype, "max exact", ["%.1e" % v for v in acc[f"{xp}-{dtype}"]["max_rel_exact"]])
    thr = [throughput(xp, dt, b) for xp in backends for dt in ("float64", "float32") for b in CONFIG["throughput_batches"]]
    for t in thr:
        print(t)

    ctx.metrics = {
        "backends": backends,
        "accuracy": {k: {kk: vv for kk, vv in v.items() if kk not in ("steps", "rel_exact", "rel_naive")}
                     for k, v in acc.items()},
        "throughput": thr,
        "gates": {
            "float64_floor_below_1e-12": max(acc["numpy-float64"]["max_rel_exact"]) < 1e-12,
        },
        # expected negative result, recorded but not a pass/fail gate: float32 cannot audit
        "findings": {
            "float32_floor_below_naive_artifacts": max(acc["numpy-float32"]["max_rel_exact"])
            < 0.1 * min(v for v in acc["numpy-float64"]["max_rel_naive"] if v > 0),
        },
    }

    plt = setup_matplotlib()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for key, v in acc.items():
        ls = "-" if "float64" in key else "--"
        ax[0].loglog(v["steps"][1:], np.maximum(v["rel_exact"][1:, 3], 1e-18), ls, lw=1, label=f"exact, {key}")
    ax[0].loglog(acc["numpy-float64"]["steps"][1:], acc["numpy-float64"]["rel_naive"][1:, 3], ":", color="k",
                 label="naive ledger (float64)")
    n = acc["numpy-float64"]["steps"][1:]
    ax[0].loglog(n, 1.1e-16 * np.sqrt(n), color="gray", lw=0.6, label=r"$\epsilon_{64}\sqrt{n}$")
    ax[0].loglog(n, 6e-8 * np.sqrt(n), color="gray", lw=0.6, ls="--", label=r"$\epsilon_{32}\sqrt{n}$")
    ax[0].set_xlabel("time steps n")
    ax[0].set_ylabel(r"$|E_{\rm net}|$ / energy scale")
    ax[0].set_title(r"Ledger round-off floor ($\delta$ = 0.04, Q = 200)")
    ax[0].legend(fontsize=7)
    for xp in backends:
        for dtype, mk in (("float64", "o-"), ("float32", "s--")):
            sel = [t for t in thr if t["backend"] == xp and t["dtype"] == dtype]
            ax[1].loglog([t["batch"] for t in sel], [t["Mcell_updates_per_s"] for t in sel], mk, label=f"{xp} {dtype}")
    ax[1].set_xlabel("batch size (simultaneous simulations)")
    ax[1].set_ylabel("Mcell-updates / s")
    ax[1].set_title(f"Batched FDTD throughput (N = {CONFIG['throughput_n_cells']})")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ctx.figure_path("exp0003b_precision_throughput.png"))
    ctx.save()
    print(ctx.metrics["gates"])


if __name__ == "__main__":
    main()
