"""Minimal experiment context: every run writes manifest.json + metrics.json + figures."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from cavitylab.core.manifest import build_manifest, write_json


class ExperimentContext:
    """Collects metrics for one experiment and writes them with a reproducibility manifest.

    Layout written under ``exp_dir``::

        results/manifest.json   git commit, environment, config hash
        results/metrics.json    numbers quoted in conclusion.md / papers
        figures/*.png
    """

    def __init__(self, experiment_id: str, exp_dir: Path, config: dict[str, Any]) -> None:
        self.experiment_id = experiment_id
        self.exp_dir = Path(exp_dir)
        self.config = config
        self.metrics: dict[str, Any] = {}
        self.results_dir = self.exp_dir / "results"
        self.figures_dir = self.exp_dir / "figures"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        self._t0 = time.perf_counter()

    def figure_path(self, name: str) -> Path:
        return self.figures_dir / name

    def save(self) -> None:
        manifest = build_manifest(self.experiment_id, self.config)
        manifest["wall_time_s"] = round(time.perf_counter() - self._t0, 2)
        write_json(self.results_dir / "manifest.json", manifest)
        write_json(self.results_dir / "metrics.json", self.metrics)


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "axes.grid": True,
        "grid.alpha": 0.3,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    return plt
