"""Gaussian quantum parametric mode: DCE growth, ledger, work audit, classical bridge."""

import math

import numpy as np

from cavitylab.benchmarks.parametric import fdtd_monodromy
from cavitylab.quantum import gaussian_thermo, photons_from_monodromy, simulate_parametric_mode


def test_vacuum_is_stationary_without_modulation():
    r = simulate_parametric_mode(100.0, delta=0.0, quality=50.0, n_eval=11)
    assert np.allclose(r.photons, 0.0, atol=1e-10)
    assert np.allclose(r.bath_energy_out, 0.0, atol=1e-10)  # zero-point energy does not flow to a T=0 bath


def test_resonant_dce_matches_rwa_and_stays_pure():
    d = 0.02
    r = simulate_parametric_mode(400.0, d, n_eval=5)
    assert abs(r.photons[-1] / math.sinh(d * 400 / 4) ** 2 - 1) < 5e-3
    assert r.ergotropy[-1] / r.photons[-1] > 0.9999
    assert np.abs(r.ledger_residual).max() < 1e-8


def test_lossy_ledger_and_ergotropy_bound():
    # stroboscopic grid tau = k*pi/2 where f = 1, so H(tau) equals the reference H0;
    # there ergotropy <= photons = W_pump - E_bath <= W_pump is a theorem
    r = simulate_parametric_mode(250 * math.pi, 0.03, quality=150.0, n_eval=501)
    assert np.abs(r.ledger_residual).max() < 1e-8 * r.pump_work.max()
    assert np.all(r.ergotropy <= r.pump_work + 1e-9)
    assert r.passive_energy[-1] > 0


def test_thermal_state_has_zero_ergotropy():
    n, erg, passive, s = gaussian_thermo(1.5 * np.eye(2), np.zeros(2))
    assert np.isclose(n[0], 1.0) and np.isclose(erg[0], 0.0) and np.isclose(passive[0], 1.0)
    assert s[0] > 0


def test_classical_fdtd_transfer_matrix_predicts_quantum_photons():
    M, tau = fdtd_monodromy(0.1, 200, 0.02, 2.0, pump_periods=100)
    n_fdtd = photons_from_monodromy(M[::-1, ::-1])
    n_q = simulate_parametric_mode(tau, 0.02, n_eval=3).photons[-1]
    assert abs(n_fdtd / n_q - 1) < 1e-3
    assert abs(np.linalg.det(M) - 1) < 1e-10
