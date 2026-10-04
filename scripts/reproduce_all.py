"""Re-run the experiments behind Preprint 01 and check every gate.

Usage:
    python scripts/reproduce_all.py            # full set (GPU recommended for EXP-0005)
    python scripts/reproduce_all.py --quick    # fast subset (CI smoke test)
    python scripts/reproduce_all.py --only EXP-0003 EXP-0004c

Each experiment writes results/manifest.json (git commit, environment, config hash),
results/metrics.json and figures/. This script reads the gates from metrics.json and
exits non-zero if any gate fails, so a clean-machine run certifies the paper numbers.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = [
    # (id, folder, quick?)
    ("EXP-0001", "EXP-0001_analytic_modes", True),
    ("EXP-0002", "EXP-0002_fdtd_convergence", True),
    ("EXP-0003", "EXP-0003_static_energy_balance", False),
    ("EXP-0003b", "EXP-0003b_precision_and_acceleration", False),
    ("EXP-0004", "EXP-0004_parametric_modulation", False),
    ("EXP-0004c", "EXP-0004c_multimode_converged", False),
    ("EXP-0005", "EXP-0005_parametric_sweep", False),
    ("EXP-0005b", "EXP-0005b_optimizer_exploit", False),
    ("EXP-0006", "EXP-0006_quantum_dce_gaussian", True),
]


def gates_of(metrics: dict) -> dict:
    """Collect boolean gates wherever an experiment stores them."""
    out = {}
    for key, val in metrics.items():
        if key == "gates" and isinstance(val, dict):
            out.update(val)
        elif key.startswith("gate_") and isinstance(val, bool):
            out[key] = val
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args()

    selected = [e for e in EXPERIMENTS
                if (not args.quick or e[2]) and (not args.only or e[0] in args.only)]
    failures = []
    for exp_id, folder, _ in selected:
        script = ROOT / "experiments" / folder / "run.py"
        t0 = time.perf_counter()
        print(f"=== {exp_id}: {script.relative_to(ROOT)}", flush=True)
        proc = subprocess.run([sys.executable, "-u", str(script)], cwd=ROOT)
        dt = time.perf_counter() - t0
        if proc.returncode != 0:
            failures.append(f"{exp_id}: exit code {proc.returncode}")
            continue
        metrics = json.loads((script.parent / "results" / "metrics.json").read_text(encoding="utf-8"))
        gates = gates_of(metrics)
        bad = [k for k, v in gates.items() if v is not True]
        status = "PASS" if not bad else "FAIL " + ", ".join(bad)
        print(f"--- {exp_id}: {status} ({len(gates)} gates, {dt:.0f} s)", flush=True)
        if bad:
            failures.append(f"{exp_id}: {bad}")

    print("\nSUMMARY:", "all gates passed" if not failures else failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
