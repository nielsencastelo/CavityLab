r"""Single-mode reference for a uniformly modulated, lossy 1D cavity (benchmark B-004).

With eps(t) = eps_s (1 + delta sin(Omega t + phi)) filling the whole PEC cavity, the
spatial modes decouple exactly. In dimensionless form (tau = omega_0 t,
p = d/eps_s, q = eta_s h, f = 1 + delta sin(nu tau + phi), nu = Omega/omega_0):

    dp/dtau = -q - p / (Q f)
    dq/dtau =  p / f

with normalized stored energy  w = p^2/f + q^2  (W = eps_s L/4 * w per unit area)
and normalized pump power  dw_pump/dtau = -(p/f)^2 * df/dtau.
This is a damped Hill equation; parametric resonance occurs near nu = 2/m.
For small delta at nu = 2 the amplitude growth rate is delta/4 - 1/(2Q),
so the threshold is delta_th = 2/Q.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp


def _rhs(tau, y, delta, nu, phi, inv_q):
    p, q, w_pump, w_loss = y
    f = 1.0 + delta * math.sin(nu * tau + phi)
    dfdtau = delta * nu * math.cos(nu * tau + phi)
    e = p / f
    return [-q - inv_q * e, e, -(e**2) * dfdtau, 2.0 * inv_q * e**2]


def integrate_mode(
    tau_end: float,
    delta: float,
    nu: float,
    quality: float = math.inf,
    phi: float = 0.0,
    p0: float = 1.0,
    q0: float = 0.0,
    n_eval: int = 2001,
    rtol: float = 1e-11,
):
    """Integrate the normalized mode ODE; returns dict with tau, p, q, w, w_pump, w_loss."""
    inv_q = 0.0 if math.isinf(quality) else 1.0 / quality
    tau = np.linspace(0.0, tau_end, n_eval)
    sol = solve_ivp(
        _rhs, (0.0, tau_end), [p0, q0, 0.0, 0.0], t_eval=tau, method="DOP853",
        rtol=rtol, atol=1e-14, args=(delta, nu, phi, inv_q),
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    p, q, w_pump, w_loss = sol.y
    f = 1.0 + delta * np.sin(nu * tau + phi)
    return {"tau": tau, "p": p, "q": q, "w": p**2 / f + q**2, "w_pump": w_pump, "w_loss": w_loss}


def floquet_growth_rate(delta: float, nu: float, quality: float = math.inf, phi: float = 0.0) -> float:
    """Largest Floquet amplitude growth rate per unit tau (negative = decaying)."""
    inv_q = 0.0 if math.isinf(quality) else 1.0 / quality
    period = 2.0 * math.pi / nu

    def rhs2(tau, y):
        p, q = y
        f = 1.0 + delta * math.sin(nu * tau + phi)
        e = p / f
        return [-q - inv_q * e, e]

    cols = []
    for y0 in ([1.0, 0.0], [0.0, 1.0]):
        sol = solve_ivp(rhs2, (0.0, period), y0, method="DOP853", rtol=1e-12, atol=1e-14)
        cols.append(sol.y[:, -1])
    monodromy = np.array(cols).T
    mult = np.max(np.abs(np.linalg.eigvals(monodromy)))
    return math.log(mult) / period
