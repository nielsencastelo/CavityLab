"""Analytical references for the 1D PEC cavity (benchmarks B-001..B-003)."""

from __future__ import annotations

import math

import numpy as np

from cavitylab.core.constants import EPS0, MU0


def yee_dispersion_omega(k: float, dx: float, dt: float, wave_speed: float) -> float:
    """Angular frequency of a Yee eigenmode: sin(w dt/2) = (v dt/dx) sin(k dx/2)."""
    s = wave_speed * dt / dx * math.sin(0.5 * k * dx)
    if abs(s) > 1.0:
        raise ValueError("unstable: Courant condition violated for this wavenumber")
    return 2.0 / dt * math.asin(s)


def discrete_mode_fields(sim, n: int, amplitude: float = 1.0, eps_r: float = 1.0):
    """Exact Yee eigenmode n of a homogeneous PEC cavity filling the whole domain.

    Returns (E^0 at nodes, H^{-1/2} at half nodes, discrete angular frequency).
    E_i^m = A sin(k x_i) cos(w m dt),  H^{m+1/2} = (A/eta) cos(k x_{i+1/2}) sin(w (m+1/2) dt).
    """
    k = n * math.pi / sim.length
    v = 1.0 / math.sqrt(MU0 * EPS0 * eps_r)
    eta = math.sqrt(MU0 / (EPS0 * eps_r))
    w = yee_dispersion_omega(k, sim.dx, sim.dt, v)
    E = amplitude * np.sin(k * sim.x_e)
    H = (amplitude / eta) * np.cos(k * sim.x_h) * math.sin(-0.5 * w * sim.dt)
    return E, H, w


def estimate_frequency(samples: np.ndarray, dt: float) -> float:
    """Frequency of a single-harmonic sampled signal from a_{n+1} + a_{n-1} = 2 cos(w dt) a_n.

    Least-squares estimate; exact (to round-off) for a pure undamped harmonic.
    """
    a = np.asarray(samples, float)
    num = np.dot(a[1:-1], a[2:] + a[:-2])
    den = 2.0 * np.dot(a[1:-1], a[1:-1])
    return math.acos(np.clip(num / den, -1.0, 1.0)) / dt


def odd_periodic_extension(f, length: float):
    """F(y): odd, 2L-periodic extension of f defined on [0, L] (PEC image method)."""

    def F(y):
        y = np.mod(np.asarray(y, float), 2.0 * length)
        return np.where(y <= length, f(np.minimum(y, length)), -f(np.clip(2.0 * length - y, 0.0, length)))

    return F


def gaussian(center: float, width: float, amplitude: float = 1.0):
    return lambda x: amplitude * np.exp(-0.5 * ((np.asarray(x) - center) / width) ** 2)


def pulse_fields(f, length: float, x, t: float, eps_r: float = 1.0):
    """Exact d'Alembert solution in a PEC cavity for E(x,0) = f(x), H(x,0) = 0.

    E(x,t) = [F(x-vt) + F(x+vt)]/2,  H(x,t) = [F(x+vt) - F(x-vt)]/(2 eta).
    """
    v = 1.0 / math.sqrt(MU0 * EPS0 * eps_r)
    eta = math.sqrt(MU0 / (EPS0 * eps_r))
    F = odd_periodic_extension(f, length)
    E = 0.5 * (F(x - v * t) + F(x + v * t))
    H = (F(x + v * t) - F(x - v * t)) / (2.0 * eta)
    return E, H
