"""Analytical benchmarks: Yee dispersion, convergence order, Hill/Floquet threshold."""

import math

import numpy as np
import pytest

from cavitylab import EPS0, FDTD1D, Modulation
from cavitylab.benchmarks.analytic1d import (
    discrete_mode_fields,
    estimate_frequency,
    gaussian,
    pulse_fields,
)
from cavitylab.benchmarks.hill import floquet_growth_rate, integrate_mode
from cavitylab.geometry import Cavity1D

L = 0.1


def test_discrete_eigenmode_frequency_matches_yee_dispersion():
    sim = FDTD1D(L, 100, 0.6)
    E, H, w = discrete_mode_fields(sim, 5)
    sim.set_fields(E, H)
    res = sim.run(4000, modes=[5])
    assert estimate_frequency(res.modes[5], sim.dt) == pytest.approx(w, rel=1e-11)
    assert w < Cavity1D(L).mode_angular_frequency(5)  # Yee is dispersive-slow


def _pulse_error(n_cells, courant=0.5, t_end_rt=1.37):
    sim = FDTD1D(L, n_cells, courant)
    f = gaussian(0.4 * L, L / 20)
    sim.set_fields(pulse_fields(f, L, sim.x_e, 0.0)[0], pulse_fields(f, L, sim.x_h, -0.5 * sim.dt)[1])
    n_steps = round(t_end_rt * 2 * L / 3e8 / sim.dt)
    sim.run(n_steps, record_every=n_steps)
    exact = pulse_fields(f, L, sim.x_e, sim.time)[0]
    return math.sqrt(np.mean((sim.E - exact) ** 2))


def test_second_order_convergence():
    errs = [_pulse_error(n) for n in (100, 200, 400)]
    orders = [math.log2(errs[i] / errs[i + 1]) for i in range(2)]
    assert all(1.9 < p < 2.1 for p in orders), orders


def test_mathieu_threshold_two_over_q():
    q = 200.0
    assert floquet_growth_rate(1.5 * 2 / q, 2.0, q) > 0
    assert floquet_growth_rate(0.5 * 2 / q, 2.0, q) < 0


def test_hill_ode_energy_balance():
    sol = integrate_mode(200.0, delta=0.05, nu=2.0, quality=100.0)
    residual = (sol["w"] - sol["w"][0]) - (sol["w_pump"] - sol["w_loss"])
    assert np.max(np.abs(residual)) < 1e-8 * np.max(sol["w"])


def test_fdtd_uniform_modulation_tracks_hill_ode():
    q, delta = 100.0, 0.06
    n = 400
    w0 = Cavity1D(L).mode_angular_frequency(1)
    sigma = w0 * EPS0 / q
    sim = FDTD1D(L, n, 0.5, sigma=sigma, modulation=Modulation(delta, 2 * w0))
    E, H, _ = discrete_mode_fields(sim, 1)
    sim.set_fields(E, H)
    periods = 30
    res = sim.run(round(periods * 2 * math.pi / w0 / sim.dt), record_every=50)
    sol = integrate_mode(w0 * res.t[-1], delta, 2.0, q, n_eval=len(res.t))
    ratio_fdtd = res.W[-1] / res.W[0]
    ratio_ode = sol["w"][-1] / sol["w"][0]
    assert ratio_fdtd == pytest.approx(ratio_ode, rel=2e-2)
