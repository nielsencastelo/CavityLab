r"""EnergyAudit: first-class energy accounting for every simulation.

Sign convention used throughout CavityLab (project-wide, see docs/theory):

    E_in   = work delivered to the field by pumps/modulation and sources (signed)
    E_out  = energy leaving through ports
    E_loss = energy dissipated in materials
    dE     = W(t) - W(0), change of stored energy

    E_net(t) = E_out + E_loss + dE - E_in

E_net is the energy that *appeared without a booked input*. Exact conservation
means E_net = 0. The ZPE-AI null hypothesis is H0: E_net <= 0 (within
uncertainty); any E_net > 0 is first treated as a *candidate computational
anomaly* and must survive the falsification checklist.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class EnergyAudit:
    t: np.ndarray
    stored: np.ndarray
    e_in: np.ndarray
    e_out: np.ndarray
    e_loss: np.ndarray
    tolerance: float = 1e-10
    label: str = "exact"

    @classmethod
    def from_fdtd(cls, result, naive: bool = False, tolerance: float = 1e-10) -> "EnergyAudit":
        """Build an audit from an :class:`~cavitylab.solvers.fdtd1d.FDTDResult`.

        ``naive=True`` audits the textbook estimators instead of the discrete-exact ones.
        """
        if naive:
            return cls(
                t=result.t,
                stored=result.naive_W,
                e_in=result.naive_pump + result.naive_source,
                e_out=result.naive_port,
                e_loss=result.naive_loss,
                tolerance=tolerance,
                label="naive",
            )
        return cls(
            t=result.t,
            stored=result.W,
            e_in=result.pump + result.source,
            e_out=result.port,
            e_loss=result.loss,
            tolerance=tolerance,
            label="exact",
        )

    # ------------------------------------------------------------ quantities
    @property
    def delta_stored(self) -> np.ndarray:
        return self.stored - self.stored[0]

    @property
    def e_net(self) -> np.ndarray:
        return self.e_out + self.e_loss + self.delta_stored - self.e_in

    @property
    def scale(self) -> float:
        """Energy scale used to normalize residuals: the largest booked quantity."""
        cands = [np.abs(self.stored), np.abs(self.e_in), np.abs(self.e_out), np.abs(self.e_loss)]
        s = max(float(np.max(c)) for c in cands)
        return s if s > 0 else 1.0

    @property
    def max_relative_residual(self) -> float:
        return float(np.max(np.abs(self.e_net))) / self.scale

    @property
    def final_relative_e_net(self) -> float:
        return float(self.e_net[-1]) / self.scale

    def passed(self) -> bool:
        return self.max_relative_residual <= self.tolerance

    # --------------------------------------------------------------- reports
    def summary(self) -> dict:
        return {
            "ledger": self.label,
            "stored_initial": float(self.stored[0]),
            "stored_final": float(self.stored[-1]),
            "delta_stored": float(self.delta_stored[-1]),
            "e_in": float(self.e_in[-1]),
            "e_out": float(self.e_out[-1]),
            "e_loss": float(self.e_loss[-1]),
            "e_net_final": float(self.e_net[-1]),
            "energy_scale": self.scale,
            "max_relative_residual": self.max_relative_residual,
            "final_relative_e_net": self.final_relative_e_net,
            "tolerance": self.tolerance,
            "status": "PASS" if self.passed() else "FAIL",
        }

    def report(self) -> str:
        s = self.summary()
        rows = [
            ("Ledger", s["ledger"]),
            ("Stored energy (initial)", f"{s['stored_initial']:.9e} J/m^2"),
            ("Stored energy (final)", f"{s['stored_final']:.9e} J/m^2"),
            ("Input (pump + sources)", f"{s['e_in']:.9e} J/m^2"),
            ("Output (ports)", f"{s['e_out']:.9e} J/m^2"),
            ("Loss (materials)", f"{s['e_loss']:.9e} J/m^2"),
            ("E_net (final)", f"{s['e_net_final']:.3e} J/m^2"),
            ("max |E_net| / scale", f"{s['max_relative_residual']:.3e}"),
            ("Status", f"{s['status']} (tol {s['tolerance']:.0e})"),
        ]
        width = max(len(k) for k, _ in rows)
        return "ENERGY AUDIT\n" + "\n".join(f"  {k:<{width}}  {v}" for k, v in rows)
