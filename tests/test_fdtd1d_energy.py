"""Energy-ledger invariants of the native 1D FDTD solver."""

import math

import numpy as np
import pytest

from cavitylab import EnergyAudit, FDTD1D, Modulation, SoftSource
from cavitylab.benchmarks.analytic1d import discrete_mode_fields, gaussian, pulse_fields

L = 0.1


def pulse_sim(n_cells=200, **kw):
    sim = FDTD1D(L, n_cells, **kw)
    f = gaussian(0.4 * L, L / 30)
    E0, _ = pulse_fields(f, L, sim.x_e, 0.0)
    _, Hm = pulse_fields(f, L, sim.x_h, -0.5 * sim.dt)
    sim.set_fields(E0, Hm)
    return sim


def test_lossless_closed_cavity_conserves_discrete_energy():
    sim = pulse_sim()
    res = sim.run(20_000)
    audit = EnergyAudit.from_fdtd(res, tolerance=1e-12)
    assert audit.passed(), audit.report()
    assert np.ptp(res.W) / res.W[0] < 1e-12


def test_naive_energy_is_not_conserved():
    sim = pulse_sim()
    res = sim.run(5_000)
    assert np.ptp(res.naive_W) / res.naive_W[0] > 1e-6


@pytest.mark.parametrize("courant", [0.3, 0.7, 0.99])
def test_balance_closes_with_loss_source_modulation_and_port(courant):
    n = 300
    sim0 = FDTD1D(L, n, courant)
    port = sim0.x_e > 0.8 * L
    sigma = np.where(port, 2.0 * ((sim0.x_e - 0.8 * L) / (0.2 * L)) ** 2, 0.02)
    profile = (sim0.x_e > 0.6 * L) & (sim0.x_e < 0.7 * L)
    w = math.pi * 3e8 / L
    sim = FDTD1D(
        L, n, courant, eps_r=2.0, sigma=sigma, port_mask=port,
        modulation=Modulation(0.1, 1.7 * w, profile=profile.astype(float)),
        sources=[SoftSource(40, lambda t: 1e-2 * math.sin(w * t) * math.exp(-((t - 2e-9) / 1e-9) ** 2))],
    )
    E, H, _ = discrete_mode_fields(sim, 2, eps_r=2.0)
    sim.set_fields(E, H)
    res = sim.run(6_000)
    audit = EnergyAudit.from_fdtd(res, tolerance=1e-12)
    assert audit.passed(), audit.report()
    assert res.port[-1] > 0 and res.loss[-1] > 0 and res.source[-1] != 0 and res.pump[-1] != 0


def test_run_can_be_resumed_without_breaking_the_ledger():
    sim = pulse_sim(sigma=0.01)
    a = sim.run(1000)
    b = sim.run(1000)
    assert b.W[0] == pytest.approx(a.W[-1], rel=1e-14)
    assert b.W[0] - b.W[-1] == pytest.approx(b.loss[-1], rel=1e-10)
