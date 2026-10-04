r"""Native 1D FDTD (Yee) solver with a discrete-exact energy ledger.

Polarization Ez/Hy on [0, L_domain] with PEC walls (Ez = 0 at both ends).
Ez/Dz live on integer nodes x_i = i*dx (i = 0..N), Hy on half nodes x_{i+1/2}.

Update (D-formulation, so that a time-varying permittivity is handled exactly):

    H^{n+1/2} = H^{n-1/2} + dt/(mu dx) (E^n_{i+1} - E^n_i)
    D^{n+1}   = D^n + dt/dx (H^{n+1/2}_{i+1/2} - H^{n+1/2}_{i-1/2}) - dt sigma Ebar - dt J^{n+1/2}
    E^{n+1}   = D^{n+1} / eps^{n+1},     Ebar = (E^n + E^{n+1}) / 2

The conduction term is semi-implicit (Crank-Nicolson in sigma).

Discrete energy (exactly conserved by the lossless, unmodulated, source-free scheme):

    W^n = dx/2 * sum_i D^n_i E^n_i + dx/2 * mu * sum_j H^{n+1/2}_j H^{n-1/2}_j

Exact per-step balance (see docs/theory/discrete_energy_identity.md):

    W^{n+1} - W^n = dx * sum_i [ 1/2 D^n D^{n+1} (1/eps^{n+1} - 1/eps^n)   (pump work)
                                - dt sigma Ebar^2                         (dissipation)
                                - dt J Ebar ]                             (source work)

Dissipation at nodes flagged by ``port_mask`` is booked as output (E_out); the rest
as loss (E_loss). Alongside the exact ledger the solver also books a *naive*
ledger (left-point rules, time-averaged H) to quantify how textbook estimators
create spurious energy residuals.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np

from cavitylab.core.constants import C0, EPS0, MU0


@dataclass
class Modulation:
    """Permittivity modulation eps(x, t) = eps_static(x) * (1 + depth * profile(x) * sin(Omega t + phase))."""

    depth: float
    angular_frequency: float
    phase: float = 0.0
    profile: np.ndarray | None = None  # node weights (N+1,), None means the whole domain
    ramp_time: float = 0.0  # optional smooth turn-on (sin^2 envelope) [s]

    def envelope(self, t: float) -> float:
        if self.ramp_time <= 0.0 or t >= self.ramp_time:
            return 1.0
        if t <= 0.0:
            return 0.0
        return math.sin(0.5 * math.pi * t / self.ramp_time) ** 2


@dataclass
class SoftSource:
    """Soft current-density source J_z(t) [A/m^2] at a single E node."""

    node: int
    waveform: Callable[[float], float]


@dataclass
class FDTDResult:
    """Time series recorded by :meth:`FDTD1D.run`. Cumulative quantities start at 0."""

    t: np.ndarray
    W: np.ndarray  # exact discrete stored energy [J/m^2]
    W_E: np.ndarray
    W_H: np.ndarray
    pump: np.ndarray  # cumulative work done by the modulation on the field
    loss: np.ndarray  # cumulative dissipation outside the port
    port: np.ndarray  # cumulative dissipation inside the port (output)
    source: np.ndarray  # cumulative work done by sources on the field
    naive_W: np.ndarray
    naive_pump: np.ndarray
    naive_loss: np.ndarray
    naive_port: np.ndarray
    naive_source: np.ndarray
    modes: dict[int, np.ndarray] = field(default_factory=dict)  # E projections a_m(t)
    meta: dict = field(default_factory=dict)


class FDTD1D:
    """1D Yee solver on a PEC-terminated domain with exact energy bookkeeping."""

    def __init__(
        self,
        length: float,
        n_cells: int,
        courant: float = 0.5,
        eps_r: float | np.ndarray = 1.0,
        sigma: float | np.ndarray = 0.0,
        port_mask: np.ndarray | None = None,
        modulation: Modulation | None = None,
        sources: Sequence[SoftSource] = (),
        cavity_length: float | None = None,
    ) -> None:
        if n_cells < 2:
            raise ValueError("n_cells must be >= 2")
        self.length = float(length)
        self.n_cells = int(n_cells)
        self.dx = self.length / self.n_cells
        self.x_e = np.arange(self.n_cells + 1) * self.dx
        self.x_h = (np.arange(self.n_cells) + 0.5) * self.dx
        n_nodes = self.n_cells + 1

        self.eps_static = EPS0 * np.broadcast_to(np.asarray(eps_r, float), (n_nodes,)).copy()
        self.sigma = np.broadcast_to(np.asarray(sigma, float), (n_nodes,)).copy()
        self.port_mask = (
            np.zeros(n_nodes, bool) if port_mask is None else np.asarray(port_mask, bool).copy()
        )
        self.modulation = modulation
        if modulation is not None and modulation.profile is None:
            modulation.profile = np.ones(n_nodes)
        self.sources = list(sources)
        self.cavity_length = self.length if cavity_length is None else float(cavity_length)

        min_eps = self.eps_static.min()
        if modulation is not None:
            min_eps *= 1.0 - abs(modulation.depth) * float(np.max(np.abs(modulation.profile)))
            if min_eps <= 0:
                raise ValueError("modulation depth makes the permittivity non-positive")
        self.c_max = C0 * math.sqrt(EPS0 / min_eps)
        self.courant = float(courant)
        self.dt = self.courant * self.dx / self.c_max

        self.step_index = 0
        self.E = np.zeros(n_nodes)
        self.H = np.zeros(self.n_cells)  # H^{n-1/2}
        self.D = np.zeros(n_nodes)

    # ------------------------------------------------------------------ state
    @property
    def time(self) -> float:
        return self.step_index * self.dt

    def eps_at(self, t: float) -> np.ndarray:
        m = self.modulation
        if m is None:
            return self.eps_static
        factor = m.depth * m.envelope(t) * math.sin(m.angular_frequency * t + m.phase)
        return self.eps_static * (1.0 + factor * m.profile)

    def set_fields(self, E: np.ndarray, H_prev: np.ndarray) -> None:
        """Set E^n at the current step and H^{n-1/2} (half a step earlier)."""
        E = np.asarray(E, float).copy()
        E[0] = E[-1] = 0.0
        self.E = E
        self.H = np.asarray(H_prev, float).copy()
        self.D = self.eps_at(self.time) * self.E

    # -------------------------------------------------------------------- run
    def run(self, n_steps: int, record_every: int = 1, modes: Sequence[int] = ()) -> FDTDResult:
        dx, dt = self.dx, self.dt
        coef_h = dt / (MU0 * dx)
        coef_d = dt / dx
        sig = self.sigma[1:-1]
        port = self.port_mask[1:-1]
        has_loss = bool(np.any(sig != 0.0))

        n_rec = n_steps // record_every + 1
        rec = {k: np.zeros(n_rec) for k in (
            "t", "W", "W_E", "W_H", "naive_W")}
        cum = np.zeros(4)  # pump, loss, port, source (exact)
        cum_naive = np.zeros(4)
        cum_rec = np.zeros((n_rec, 4))
        cum_naive_rec = np.zeros((n_rec, 4))

        lc = self.cavity_length
        in_cav = self.x_e <= lc + 1e-12 * lc
        basis = {
            m: (2.0 / lc) * dx * np.sin(m * math.pi * self.x_e / lc) * in_cav for m in modes
        }
        mode_rec = {m: np.zeros(n_rec) for m in modes}

        src_nodes = np.array([s.node for s in self.sources], int)
        E, D, H = self.E, self.D, self.H
        eps_cur = self.eps_at(self.time)
        r = 0
        for step in range(n_steps + 1):
            H_new = H + coef_h * (E[1:] - E[:-1])
            if step % record_every == 0:
                w_e = 0.5 * dx * np.dot(D, E)
                w_h = 0.5 * dx * MU0 * np.dot(H_new, H)
                h_avg = 0.5 * (H_new + H)
                rec["t"][r] = self.time
                rec["W_E"][r] = w_e
                rec["W_H"][r] = w_h
                rec["W"][r] = w_e + w_h
                rec["naive_W"][r] = 0.5 * dx * (np.dot(eps_cur * E, E) + MU0 * np.dot(h_avg, h_avg))
                cum_rec[r] = cum
                cum_naive_rec[r] = cum_naive
                for m, b in basis.items():
                    mode_rec[m][r] = np.dot(b, E)
                r += 1
            if step == n_steps:
                break  # keep H at n-1/2 so that the state stays consistent
            H = H_new

            t_next = self.time + dt
            eps_next = self.eps_at(t_next)
            J = np.zeros_like(E)
            if self.sources:
                t_half = self.time + 0.5 * dt
                for s in self.sources:
                    J[s.node] += s.waveform(t_half)

            rhs = D[1:-1] + coef_d * (H[1:] - H[:-1]) - dt * J[1:-1]
            if has_loss:
                rhs -= 0.5 * dt * sig * E[1:-1]
                D_new_i = rhs / (1.0 + 0.5 * dt * sig / eps_next[1:-1])
            else:
                D_new_i = rhs
            E_new_i = D_new_i / eps_next[1:-1]
            E_old_i = E[1:-1]
            e_bar = 0.5 * (E_new_i + E_old_i)

            # exact ledger
            if self.modulation is not None:
                cum[0] += 0.5 * dx * np.dot(
                    D[1:-1] * D_new_i, 1.0 / eps_next[1:-1] - 1.0 / eps_cur[1:-1]
                )
                cum_naive[0] += -0.5 * dx * np.dot(E_old_i**2, eps_next[1:-1] - eps_cur[1:-1])
            if has_loss:
                diss = dt * dx * sig * e_bar**2
                diss_naive = dt * dx * sig * E_old_i**2
                cum[1] += diss[~port].sum()
                cum[2] += diss[port].sum()
                cum_naive[1] += diss_naive[~port].sum()
                cum_naive[2] += diss_naive[port].sum()
            if self.sources:
                Js = J[src_nodes]
                cum[3] += -dt * dx * np.dot(Js, 0.5 * (E_new_i[src_nodes - 1] + E_old_i[src_nodes - 1]))
                cum_naive[3] += -dt * dx * np.dot(Js, E_old_i[src_nodes - 1])

            D[1:-1] = D_new_i
            E[1:-1] = E_new_i
            eps_cur = eps_next
            self.step_index += 1

        self.E, self.D, self.H = E, D, H
        return FDTDResult(
            t=rec["t"],
            W=rec["W"],
            W_E=rec["W_E"],
            W_H=rec["W_H"],
            pump=cum_rec[:, 0],
            loss=cum_rec[:, 1],
            port=cum_rec[:, 2],
            source=cum_rec[:, 3],
            naive_W=rec["naive_W"],
            naive_pump=cum_naive_rec[:, 0],
            naive_loss=cum_naive_rec[:, 1],
            naive_port=cum_naive_rec[:, 2],
            naive_source=cum_naive_rec[:, 3],
            modes=mode_rec,
            meta={
                "solver": "cavitylab.native_fdtd1d",
                "dx": dx,
                "dt": dt,
                "n_cells": self.n_cells,
                "courant": self.courant,
                "n_steps": n_steps,
                "record_every": record_every,
                "precision": str(E.dtype),
            },
        )
