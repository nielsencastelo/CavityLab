"""Run manifests: everything needed to reproduce a numerical result.

Every experiment writes a manifest next to its metrics (git commit, environment,
configuration hash, timestamp). Results without a manifest are not citable.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np


def _git(*args: str) -> str | None:
    try:
        out = subprocess.run(
            ["git", *args], capture_output=True, text=True, check=True, timeout=10
        )
        return out.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


def config_hash(config: dict[str, Any]) -> str:
    """Stable SHA-256 of a JSON-serializable configuration."""
    blob = json.dumps(config, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def build_manifest(experiment_id: str, config: dict[str, Any]) -> dict[str, Any]:
    import scipy

    import cavitylab

    dirty = _git("status", "--porcelain")
    return {
        "experiment_id": experiment_id,
        "date_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "git_commit": _git("rev-parse", "HEAD"),
        "git_dirty": bool(dirty) if dirty is not None else None,
        "cavitylab_version": cavitylab.__version__,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "precision": "float64",
        "config_hash": config_hash(config),
        "config": config,
    }


def to_jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, float) and not np.isfinite(obj):
        return str(obj)
    return obj


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(to_jsonable(payload), indent=2) + "\n", encoding="utf-8")
