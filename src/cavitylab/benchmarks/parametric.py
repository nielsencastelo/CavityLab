"""Floquet analysis of the native FDTD solver for a uniformly modulated cavity.

The monodromy matrix of mode m is measured directly from two FDTD runs started in
the two field quadratures (E ~ sin kx, H ~ 0) and (E ~ 0, H ~ cos kx). The modal
state is (p, q) = (f a_m, eta_s b_m), with a_m, b_m the projections of E and of the
time-averaged H. Both Floquet exponents can then be compared with the Hill ODE
(:func:`cavitylab.benchmarks.hill.floquet_growth_rate`), independently of which
branch a particular initial condition happens to excite.
"""

from __future__ import annotations

import math

import numpy as np

from cavitylab.core.constants import EPS0, MU0
from cavitylab.benchmarks.analytic1d import discrete_mode_fields
from cavitylab.solvers.fdtd1d import FDTD1D, Modulation


def fdtd_monodromy(
    length: float,
    n_cells: int,
    delta: float,
    nu: float,
    quality: float = math.inf,
    pump_periods: int = 20,
    courant: float = 0.5,
    mode: int = 1,
    phase: float = 0.0,
) -> tuple[np.ndarray, float]:
    """Measured transfer matrix of mode ``mode`` acting on (p, q), and the elapsed tau = omega_m t.

    ``nu`` is Omega / omega_m with omega_m the *continuum* mode frequency.
    """
    eta = math.sqrt(MU0 / EPS0)
    k = mode * math.pi / length
    w_m = k / math.sqrt(MU0 * EPS0)
    big_omega = nu * w_m
    sigma = 0.0 if math.isinf(quality) else w_m * EPS0 / quality

    states0, states1 = [], []
    t_end = None
    for quadrature in (0, 1):
        sim = FDTD1D(length, n_cells, courant, sigma=sigma,
                     modulation=Modulation(delta, big_omega, phase=phase))
        if quadrature == 0:
            E, H, _ = discrete_mode_fields(sim, mode)
        else:
            E = np.zeros_like(sim.x_e)
            H = (1.0 / eta) * np.cos(k * sim.x_h)
        sim.set_fields(E, H)
        if t_end is None:
            # integer number of steps closest to an integer number of pump periods
            n_steps = round(pump_periods * 2 * math.pi / big_omega / sim.dt)
            t_end = n_steps * sim.dt
        res = sim.run(n_steps, record_every=n_steps, modes=[mode])
        f = 1.0 + delta * np.sin(big_omega * res.t + phase)
        p = f * res.modes[mode]
        q = eta * res.modes_h[mode]
        states0.append([p[0], q[0]])
        states1.append([p[-1], q[-1]])

    y0 = np.array(states0).T
    y1 = np.array(states1).T
    return y1 @ np.linalg.inv(y0), w_m * t_end


def fdtd_floquet_exponents(*args, **kwargs) -> tuple[float, float]:
    """The two Floquet amplitude growth rates per unit tau, sorted descending.

    Same arguments as :func:`fdtd_monodromy`.
    """
    monodromy, tau_end = fdtd_monodromy(*args, **kwargs)
    mult = np.sort(np.abs(np.linalg.eigvals(monodromy)))[::-1]
    return float(math.log(mult[0]) / tau_end), float(math.log(mult[1]) / tau_end)
