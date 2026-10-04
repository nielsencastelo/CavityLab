"""The batched solver must reproduce the single-run solver and its ledger exactly."""

import math

import numpy as np
import pytest

from cavitylab import EPS0, FDTD1D, Modulation
from cavitylab.benchmarks.analytic1d import discrete_mode_fields
from cavitylab.benchmarks.parametric import batch_floquet_exponents
from cavitylab.benchmarks.hill import floquet_growth_rate
from cavitylab.solvers.fdtd1d_batch import FDTD1DBatch, gpu_available

L = 0.1
W1 = math.pi * 299792458.0 / L


def _reference(n=200, steps=3000):
    sig = W1 * EPS0 / 100
    s = FDTD1D(L, n, 0.5, sigma=sig, modulation=Modulation(0.1, 2 * W1, phase=0.3))
    E, H, _ = discrete_mode_fields(s, 1)
    s.set_fields(E, H)
    return s, s.run(steps, record_every=100, modes=[1]), E, H, sig


def _batch(xp, E, H, sig, n=200, steps=3000):
    b = FDTD1DBatch(L, n, 3, 0.5, sigma=[sig, sig, 0.0], depth=[0.1, 0.05, 0.1],
                    angular_frequency=[2 * W1, 2 * W1, 1.5 * W1], phase=[0.3, 0.0, 0.0], xp=xp)
    b.set_fields(E, H)
    return b, b.run(steps, record_every=100, modes=[1])


def test_batch_matches_single_run_numpy():
    s, r, E, H, sig = _reference()
    b, rb = _batch("numpy", E, H, sig)
    assert b.dt == s.dt
    assert np.abs(rb.W[:, 0] - r.W).max() < 1e-14 * r.W[0]
    assert np.abs(rb.pump[:, 0] - r.pump).max() < 1e-14 * r.W[0]
    assert np.abs(rb.naive_W[:, 0] - r.naive_W).max() < 1e-14 * r.W[0]
    assert np.abs(rb.modes[1][:, 0] - r.modes[1]).max() < 1e-12
    assert np.all(np.abs(rb.e_net()).max(0) / rb.scale() < 1e-12)


@pytest.mark.skipif(not gpu_available(), reason="no CUDA GPU")
def test_batch_gpu_matches_cpu():
    _, r, E, H, sig = _reference()
    _, rc = _batch("numpy", E, H, sig)
    xp = "torch"
    try:
        _, rg = _batch(xp, E, H, sig)
    except (ImportError, RuntimeError):
        pytest.skip("torch CUDA backend unavailable")
    assert np.abs(rg.W - rc.W).max() < 1e-12 * rc.W[0, 0]
    assert np.all(np.abs(rg.e_net()).max(0) / rg.scale() < 1e-12)


def test_batch_floquet_matches_theory():
    mu = batch_floquet_exponents(L, 128, np.array([0.02, 0.04]), np.array([2.0, 2.0]), 200.0, pump_periods=40)
    for i, d in enumerate((0.02, 0.04)):
        assert abs(mu[i, 0] - floquet_growth_rate(d, 2.0, 200.0)) < 2e-5
