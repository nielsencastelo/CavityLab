r"""Batched native 1D FDTD: B independent simulations advanced together (CPU or GPU).

Same scheme and the same discrete-exact energy ledger as :class:`~cavitylab.solvers.fdtd1d.FDTD1D`
(see that module for the equations). All arrays carry a leading batch axis
(B, N+1) / (B, N). Every member shares the grid, dx, dt, eps_static profile and
modulation profile, and has its own modulation depth, angular frequency, phase and
conductivity. This is the shape of parameter sweeps, Floquet maps, ensembles and
optimizer populations, and is where the GPU pays off.

Backend: ``xp="numpy"`` (default), ``xp="torch"`` or ``xp="cupy"`` (NVIDIA GPU, float64).
Results come back as NumPy arrays.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from cavitylab.core.constants import C0, EPS0, MU0


class _Backend:
    """Minimal array-backend shim: numpy (CPU), cupy or torch (CUDA GPU); float64 or float32."""

    def __init__(self, name: str, dtype: str = "float64") -> None:
        self.name = name
        self.dtype = dtype
        if name == "numpy":
            self.mod = np
        elif name == "cupy":
            import cupy

            self.mod = cupy
        elif name == "torch":
            import torch

            if not torch.cuda.is_available():
                raise RuntimeError("torch backend requested but CUDA is not available")
            self.mod = torch
            self.device = torch.device("cuda")
        else:
            raise ValueError(f"unknown array backend {name!r}")

    def zeros(self, shape):
        if self.name == "torch":
            return self.mod.zeros(shape, dtype=getattr(self.mod, self.dtype), device=self.device)
        return self.mod.zeros(shape, dtype=self.dtype)

    def asarray(self, a):
        if self.name == "torch":
            a = np.array(a)  # writable copy (torch refuses read-only buffers)
            dtype = getattr(self.mod, self.dtype) if a.dtype.kind == "f" else None
            return self.mod.as_tensor(a, dtype=dtype, device=self.device)
        a = self.mod.asarray(a)
        return a.astype(self.dtype) if a.dtype.kind == "f" else a

    def host(self, a) -> np.ndarray:
        if isinstance(a, np.ndarray):
            return a
        if self.name == "torch":
            return a.detach().cpu().numpy()
        if self.name == "cupy":
            return a.get()
        return np.asarray(a)


def gpu_available() -> bool:
    """True if a CUDA GPU is usable through torch (preferred) or cupy."""
    try:
        import torch

        if torch.cuda.is_available():
            return True
    except Exception:  # noqa: BLE001 - any import/driver failure means no GPU
        pass
    try:
        import cupy

        return cupy.cuda.runtime.getDeviceCount() > 0
    except Exception:  # noqa: BLE001
        return False


@dataclass
class BatchResult:
    t: np.ndarray  # (n_rec,)
    W: np.ndarray  # (n_rec, B) exact discrete stored energy
    pump: np.ndarray  # (n_rec, B) cumulative pump work
    loss: np.ndarray  # (n_rec, B) cumulative dissipation (outside port)
    port: np.ndarray  # (n_rec, B) cumulative dissipation in port
    naive_W: np.ndarray
    naive_pump: np.ndarray
    naive_loss: np.ndarray
    naive_port: np.ndarray
    modes: dict[int, np.ndarray] = field(default_factory=dict)  # (n_rec, B)
    modes_h: dict[int, np.ndarray] = field(default_factory=dict)
    captures: dict[str, np.ndarray] = field(default_factory=dict)  # per-member snapshots
    meta: dict = field(default_factory=dict)

    def e_net(self, naive: bool = False) -> np.ndarray:
        if naive:
            return self.naive_port + self.naive_loss + (self.naive_W - self.naive_W[0]) - self.naive_pump
        return self.port + self.loss + (self.W - self.W[0]) - self.pump

    def scale(self) -> np.ndarray:
        return np.maximum.reduce([np.abs(self.W).max(0), np.abs(self.pump).max(0),
                                  self.loss.max(0), self.port.max(0)])


class FDTD1DBatch:
    def __init__(
        self,
        length: float,
        n_cells: int,
        batch: int,
        courant: float = 0.5,
        eps_r: float | np.ndarray = 1.0,
        sigma: float | np.ndarray = 0.0,
        port_mask: np.ndarray | None = None,
        depth: float | np.ndarray = 0.0,
        angular_frequency: float | np.ndarray = 0.0,
        phase: float | np.ndarray = 0.0,
        profile: np.ndarray | None = None,
        cavity_length: float | None = None,
        xp: str = "numpy",
        dtype: str = "float64",
    ) -> None:
        self.xp_name = xp
        self.dtype = dtype
        self.bk = _Backend(xp, dtype)
        self.length = float(length)
        self.n_cells = int(n_cells)
        self.batch = int(batch)
        self.dx = self.length / self.n_cells
        n_nodes = self.n_cells + 1
        self.x_e = np.arange(n_nodes) * self.dx
        self.x_h = (np.arange(self.n_cells) + 0.5) * self.dx
        self.cavity_length = self.length if cavity_length is None else float(cavity_length)

        b = self.batch
        eps_s = EPS0 * np.broadcast_to(np.asarray(eps_r, float), (n_nodes,))
        self._eps_s = eps_s
        sig = np.asarray(sigma, float)
        sig = np.broadcast_to(sig[:, None] if sig.ndim == 1 else sig, (b, n_nodes))
        self._sigma = sig
        self._port = np.zeros(n_nodes, bool) if port_mask is None else np.asarray(port_mask, bool)
        self._depth = np.broadcast_to(np.asarray(depth, float), (b,)).copy()
        self._omega = np.broadcast_to(np.asarray(angular_frequency, float), (b,)).copy()
        self._phase = np.broadcast_to(np.asarray(phase, float), (b,)).copy()
        self._profile = np.ones(n_nodes) if profile is None else np.asarray(profile, float)

        min_eps = eps_s.min() * (1.0 - np.abs(self._depth).max() * np.abs(self._profile).max())
        if min_eps <= 0:
            raise ValueError("modulation depth makes the permittivity non-positive")
        self.c_max = C0 * math.sqrt(EPS0 / min_eps)
        self.courant = float(courant)
        self.dt = self.courant * self.dx / self.c_max

        self.E = self.bk.zeros((b, n_nodes))
        self.D = self.bk.zeros((b, n_nodes))
        self.H = self.bk.zeros((b, self.n_cells))
        self.step_index = 0

    @property
    def time(self) -> float:
        return self.step_index * self.dt

    def _eps(self, t: float):
        f = self._depth * np.sin(self._omega * t + self._phase)  # (B,) on host, cheap
        return self._eps_s_dev[None, :] * (1.0 + self.bk.asarray(f)[:, None] * self._profile_dev[None, :])

    def set_fields(self, E: np.ndarray, H_prev: np.ndarray) -> None:
        """E (B, N+1) or (N+1,), H_prev (B, N) or (N,) at t = step_index*dt and -dt/2."""
        bk = self.bk
        self._eps_s_dev = bk.asarray(self._eps_s)
        self._profile_dev = bk.asarray(self._profile)
        E = np.broadcast_to(np.asarray(E, float), (self.batch, self.n_cells + 1)).copy()
        E[:, 0] = E[:, -1] = 0.0
        self.E = bk.asarray(E)
        self.H = bk.asarray(np.broadcast_to(np.asarray(H_prev, float), (self.batch, self.n_cells)).copy())
        self.D = self._eps(self.time) * self.E

    def run(
        self,
        n_steps: int,
        record_every: int = 1,
        modes: Sequence[int] = (),
        capture_steps: np.ndarray | None = None,
        capture_modes: Sequence[int] = (),
    ) -> BatchResult:
        """Advance n_steps. ``capture_steps`` (B,) records the modal (E, H) projections of
        ``capture_modes`` for each member at its own step index (relative to this call)."""
        bk = self.bk
        dx, dt = self.dx, self.dt
        coef_h = dt / (MU0 * dx)
        coef_d = dt / dx
        sig = bk.asarray(self._sigma[:, 1:-1])
        port = bk.asarray(self._port[1:-1])
        notport = ~port
        has_loss = bool(np.any(self._sigma != 0.0))
        has_mod = bool(np.any(self._depth != 0.0))

        n_rec = n_steps // record_every + 1
        keys = ("W", "pump", "loss", "port", "naive_W", "naive_pump", "naive_loss", "naive_port")
        rec = {k: bk.zeros((n_rec, self.batch)) for k in keys}
        t_rec = np.zeros(n_rec)
        cum = {k: bk.zeros(self.batch) for k in ("pump", "loss", "port", "naive_pump", "naive_loss", "naive_port")}

        lc = self.cavity_length
        in_cav = self.x_e <= lc + 1e-12 * lc
        basis = {m: bk.asarray((2.0 / lc) * dx * np.sin(m * math.pi * self.x_e / lc) * in_cav) for m in modes}
        basis_h = {m: bk.asarray((2.0 / lc) * dx * np.cos(m * math.pi * self.x_h / lc) * (self.x_h <= lc))
                   for m in modes}
        mode_rec = {m: bk.zeros((n_rec, self.batch)) for m in modes}
        mode_h_rec = {m: bk.zeros((n_rec, self.batch)) for m in modes}

        cap_basis = {m: bk.asarray((2.0 / lc) * dx * np.sin(m * math.pi * self.x_e / lc) * in_cav)
                     for m in capture_modes}
        cap_basis_h = {m: bk.asarray((2.0 / lc) * dx * np.cos(m * math.pi * self.x_h / lc) * (self.x_h <= lc))
                       for m in capture_modes}
        captures = {}
        if capture_steps is not None:
            capture_steps = np.asarray(capture_steps, int)
            for m in capture_modes:
                captures[f"a{m}"] = bk.zeros(self.batch)
                captures[f"b{m}"] = bk.zeros(self.batch)
            captures["t"] = np.zeros(self.batch)
            captures["f"] = np.zeros(self.batch)  # eps factor (1 + delta sin) at capture time
            cap_set = {int(s) for s in capture_steps}

        E, D, H = self.E, self.D, self.H
        eps_cur = self._eps(self.time)
        r = 0
        for step in range(n_steps + 1):
            H_new = H + coef_h * (E[:, 1:] - E[:, :-1])
            if step % record_every == 0:
                h_avg = 0.5 * (H_new + H)
                t_rec[r] = self.time
                rec["W"][r] = 0.5 * dx * ((D * E).sum(1) + MU0 * (H_new * H).sum(1))
                rec["naive_W"][r] = 0.5 * dx * ((eps_cur * E * E).sum(1) + MU0 * (h_avg * h_avg).sum(1))
                for k in ("pump", "loss", "port", "naive_pump", "naive_loss", "naive_port"):
                    rec[k][r] = cum[k]
                for m in modes:
                    mode_rec[m][r] = E @ basis[m]
                    mode_h_rec[m][r] = h_avg @ basis_h[m]
                r += 1
            if capture_steps is not None and step in cap_set:
                sel = np.nonzero(capture_steps == step)[0]
                sel_d = bk.asarray(sel)
                h_avg = 0.5 * (H_new + H)
                for m in capture_modes:
                    captures[f"a{m}"][sel_d] = E[sel_d] @ cap_basis[m]
                    captures[f"b{m}"][sel_d] = h_avg[sel_d] @ cap_basis_h[m]
                captures["t"][sel] = self.time
                captures["f"][sel] = 1.0 + self._depth[sel] * np.sin(self._omega[sel] * self.time + self._phase[sel])
            if step == n_steps:
                break
            H = H_new
            eps_next = self._eps(self.time + dt) if has_mod else eps_cur
            rhs = D[:, 1:-1] + coef_d * (H[:, 1:] - H[:, :-1])
            e_old = E[:, 1:-1]
            if has_loss:
                rhs = rhs - 0.5 * dt * sig * e_old
                d_new = rhs / (1.0 + 0.5 * dt * sig / eps_next[:, 1:-1])
            else:
                d_new = rhs
            e_new = d_new / eps_next[:, 1:-1]
            if has_mod:
                cum["pump"] += 0.5 * dx * (D[:, 1:-1] * d_new * (1.0 / eps_next[:, 1:-1] - 1.0 / eps_cur[:, 1:-1])).sum(1)
                cum["naive_pump"] += -0.5 * dx * (e_old * e_old * (eps_next[:, 1:-1] - eps_cur[:, 1:-1])).sum(1)
            if has_loss:
                e_bar = 0.5 * (e_new + e_old)
                diss = dt * dx * sig * e_bar * e_bar
                diss_n = dt * dx * sig * e_old * e_old
                cum["loss"] += (diss * notport).sum(1)
                cum["port"] += (diss * port).sum(1)
                cum["naive_loss"] += (diss_n * notport).sum(1)
                cum["naive_port"] += (diss_n * port).sum(1)
            D[:, 1:-1] = d_new
            E[:, 1:-1] = e_new
            eps_cur = eps_next
            self.step_index += 1

        self.E, self.D, self.H = E, D, H

        host = bk.host
        return BatchResult(
            t=t_rec,
            **{k: host(v) for k, v in rec.items()},
            modes={m: host(v) for m, v in mode_rec.items()},
            modes_h={m: host(v) for m, v in mode_h_rec.items()},
            captures={k: host(v) for k, v in captures.items()},
            meta={"solver": "cavitylab.native_fdtd1d_batch", "backend": self.xp_name, "dx": dx, "dt": dt,
                  "n_cells": self.n_cells, "batch": self.batch, "courant": self.courant, "n_steps": n_steps,
                  "precision": self.dtype},
        )
