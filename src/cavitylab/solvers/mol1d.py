r"""Independent reference solver: method of lines, 4th-order staggered FD + adaptive DOP853.

Used to validate the native Yee/leapfrog FDTD solver with a *different* spatial
discretization (4th-order staggered stencil with PEC image ghosts) and a different
time integrator (adaptive 8th-order Runge-Kutta). Normalized units: c = eps0 = mu0 = 1,
lengths in units of L, time in units of L/c.

    dD/dt = dH/dx                 (nodes, D = 0 at the PEC walls)
    dH/dt = d(D / eps(x, t))/dx   (half nodes)

Stored energy W = 1/2 sum D E dx + 1/2 sum H^2 dx (all fields at the same time).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.integrate import solve_ivp


def _d_half_from_nodes(e: np.ndarray, dx: float) -> np.ndarray:
    """4th-order derivative at half nodes from node values e[0..N] with e[0] = e[N] = 0 (odd images)."""
    ext = np.concatenate(([-e[1]], e, [-e[-2]]))  # ghosts e[-1] = -e[1], e[N+1] = -e[N-1]
    # half node j+1/2 uses e[j-1], e[j], e[j+1], e[j+2]  -> ext indices j, j+1, j+2, j+3
    return (27.0 * (ext[2:-1] - ext[1:-2]) - (ext[3:] - ext[:-3])) / (24.0 * dx)


def _d_nodes_from_half(h: np.ndarray, dx: float) -> np.ndarray:
    """4th-order derivative at interior nodes 1..N-1 from half-node values (even images at the walls)."""
    ext = np.concatenate(([h[0]], h, [h[-1]]))  # ghosts h[-1/2-1] = h[1/2], mirror at walls
    # node i (1..N-1) uses h[i-3/2], h[i-1/2], h[i+1/2], h[i+3/2] -> ext indices i-1, i, i+1, i+2
    n = len(h)
    i = np.arange(1, n)
    return (27.0 * (ext[i + 1] - ext[i]) - (ext[i + 2] - ext[i - 1])) / (24.0 * dx)


def run_mol(
    n_cells: int,
    eps_of: Callable[[np.ndarray, float], np.ndarray],
    e0: np.ndarray,
    h0: np.ndarray,
    t_eval: np.ndarray,
    rtol: float = 1e-10,
    atol: float = 1e-13,
) -> dict:
    """Integrate from t = 0 with E(x,0) = e0 (nodes) and H(x,0) = h0 (half nodes).

    ``eps_of(x, t)`` returns the relative permittivity at positions x and time t.
    Returns dict with t, W (stored energy) and final fields.
    """
    dx = 1.0 / n_cells
    x_e = np.arange(n_cells + 1) * dx
    x_i = x_e[1:-1]
    d0 = eps_of(x_i, 0.0) * e0[1:-1]

    def rhs(t, y):
        d = y[: n_cells - 1]
        h = y[n_cells - 1:]
        e = np.zeros(n_cells + 1)
        e[1:-1] = d / eps_of(x_i, t)
        return np.concatenate((_d_nodes_from_half(h, dx), _d_half_from_nodes(e, dx)))

    sol = solve_ivp(rhs, (0.0, float(t_eval[-1])), np.concatenate((d0, h0)), t_eval=t_eval,
                    method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    d = sol.y[: n_cells - 1]
    h = sol.y[n_cells - 1:]
    eps_t = np.stack([eps_of(x_i, t) for t in sol.t], axis=1)
    W = 0.5 * dx * (np.sum(d * d / eps_t, axis=0) + np.sum(h * h, axis=0))
    return {"t": sol.t, "W": W, "d_final": d[:, -1], "h_final": h[:, -1], "nfev": sol.nfev}
