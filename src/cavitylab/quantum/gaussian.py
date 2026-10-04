r"""Exact Gaussian-state dynamics of one parametrically modulated cavity mode.

Units: hbar = 1, tau = omega_0 t, vacuum quadrature variance 1/2. The mode is the
quantized version of the Hill system used for the classical cavity
(:mod:`cavitylab.benchmarks.hill`):

    H(tau) = p^2 / (2 f(tau)) + q^2 / 2,     f = 1 + delta sin(nu tau + phi)

(time-dependent "mass" f = eps(t)/eps_s). Loss is standard amplitude damping into a
thermal bath (Lindblad, rate kappa = 1/Q per tau, occupation n_th). For quadratic
Hamiltonians and linear Lindblad terms Gaussian states stay Gaussian, and the
first and second moments obey closed linear ODEs, so these results are exact.

    dV/dtau = A V + V A^T + kappa (n_th + 1/2) I,   A = [[-kappa/2, 1/f], [-1, -kappa/2]]

Energy ledger (exact, integrated alongside the moments):
    dE/dtau = P_pump + P_bath,   P_pump = -(V_pp + d_p^2) f' / (2 f^2)
E includes the zero-point energy 1/2. A T = 0 bath takes no energy from the vacuum
state (P_bath = 0 there): only excitations above the vacuum are dissipated.

Thermodynamics, using the bare (f = 1) mode Hamiltonian H0 = a^dag a as reference:
    n_tot      = (V_qq + V_pp - 1)/2 + |d|^2/2      (photons above vacuum)
    nu         = sqrt(det V)                         (symplectic eigenvalue)
    E_passive  = nu - 1/2                            (thermal state with the same entropy)
    ergotropy  = n_tot - E_passive
    S          = (nu + 1/2) ln(nu + 1/2) - (nu - 1/2) ln(nu - 1/2)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass
class ParametricModeResult:
    tau: np.ndarray
    V: np.ndarray  # (n, 2, 2) covariance in (q, p)
    d: np.ndarray  # (n, 2) means
    energy: np.ndarray  # <H(tau)> including zero-point energy
    pump_work: np.ndarray  # cumulative work done by the modulation
    bath_energy_out: np.ndarray  # cumulative energy delivered to the bath
    photons: np.ndarray
    ergotropy: np.ndarray
    passive_energy: np.ndarray
    entropy: np.ndarray

    @property
    def ledger_residual(self) -> np.ndarray:
        """E(tau) - E(0) - (W_pump - E_bath): zero up to integration error."""
        return self.energy - self.energy[0] - (self.pump_work - self.bath_energy_out)


def gaussian_thermo(V: np.ndarray, d: np.ndarray):
    """(photons, ergotropy, passive energy, entropy) of single-mode Gaussian states."""
    if V.ndim == 2:
        V = V[None]
    d = np.atleast_2d(d)
    n_tot = 0.5 * (V[:, 0, 0] + V[:, 1, 1] - 1.0) + 0.5 * (d[:, 0] ** 2 + d[:, 1] ** 2)
    nu = np.sqrt(np.clip(np.linalg.det(V), 0.25, None))
    passive = nu - 0.5
    with np.errstate(divide="ignore", invalid="ignore"):
        s = (nu + 0.5) * np.log(nu + 0.5) - np.where(nu > 0.5, (nu - 0.5) * np.log(np.maximum(nu - 0.5, 1e-300)), 0.0)
    return n_tot, n_tot - passive, passive, s


def _rhs(tau, y, delta, nu, phi, kappa, n_th):
    vqq, vqp, vpp, dq, dp, w_pump, e_bath = y
    f = 1.0 + delta * math.sin(nu * tau + phi)
    fdot = delta * nu * math.cos(nu * tau + phi)
    a11, a12, a21, a22 = -0.5 * kappa, 1.0 / f, -1.0, -0.5 * kappa
    noise = kappa * (n_th + 0.5)
    # V' = A V + V A^T + noise I
    dvqq = 2 * (a11 * vqq + a12 * vqp) + noise
    dvqp = a11 * vqp + a12 * vpp + a21 * vqq + a22 * vqp
    dvpp = 2 * (a21 * vqp + a22 * vpp) + noise
    ddq = a11 * dq + a12 * dp
    ddp = a21 * dq + a22 * dp
    p_pump = -(vpp + dp**2) * fdot / (2 * f * f)
    # energy change caused by the dissipator alone
    p_bath = 0.5 * ((-kappa * vpp + noise) / f + (-kappa * vqq + noise)) - 0.5 * kappa * (dp**2 / f + dq**2)
    return [dvqq, dvqp, dvpp, ddq, ddp, p_pump, -p_bath]


def simulate_parametric_mode(
    tau_end: float,
    delta: float,
    nu: float = 2.0,
    quality: float = math.inf,
    n_th: float = 0.0,
    phi: float = 0.0,
    V0: np.ndarray | None = None,
    d0: tuple[float, float] = (0.0, 0.0),
    n_eval: int = 2001,
    rtol: float = 1e-10,
) -> ParametricModeResult:
    """Integrate first/second moments and the energy ledger. Default start: vacuum."""
    kappa = 0.0 if math.isinf(quality) else 1.0 / quality
    V0 = 0.5 * np.eye(2) if V0 is None else np.asarray(V0, float)
    y0 = [V0[0, 0], V0[0, 1], V0[1, 1], d0[0], d0[1], 0.0, 0.0]
    tau = np.linspace(0.0, tau_end, n_eval)
    sol = solve_ivp(_rhs, (0.0, tau_end), y0, t_eval=tau, method="DOP853", rtol=rtol, atol=1e-13,
                    args=(delta, nu, phi, kappa, n_th))
    if not sol.success:
        raise RuntimeError(sol.message)
    vqq, vqp, vpp, dq, dp, w_pump, e_bath = sol.y
    V = np.stack([np.stack([vqq, vqp], -1), np.stack([vqp, vpp], -1)], -2)
    d = np.stack([dq, dp], -1)
    f = 1.0 + delta * np.sin(nu * tau + phi)
    energy = 0.5 * (vpp / f + vqq) + 0.5 * (dp**2 / f + dq**2)
    photons, erg, passive, s = gaussian_thermo(V, d)
    return ParametricModeResult(tau, V, d, energy, w_pump, e_bath, photons, erg, passive, s)


def photons_from_monodromy(M: np.ndarray, V0: np.ndarray | None = None) -> float:
    """Photons above vacuum after a linear symplectic map M acting on (q, p): V -> M V M^T.

    Lets a *classical* solver's measured transfer matrix (e.g. FDTD) predict the
    quantum photon number created from vacuum, which is exact for linear, lossless media.
    """
    V0 = 0.5 * np.eye(2) if V0 is None else V0
    V = M @ V0 @ M.T
    return 0.5 * (V[0, 0] + V[1, 1] - 1.0)
