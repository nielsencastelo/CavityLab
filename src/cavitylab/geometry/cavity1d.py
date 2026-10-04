"""1D cavity between perfect electric conductor (PEC) walls at x = 0 and x = L."""

from __future__ import annotations

import math
from dataclasses import dataclass

from cavitylab.core.constants import C0, EPS0


@dataclass(frozen=True)
class Cavity1D:
    """Homogeneous 1D cavity, Ez/Hy polarization, PEC walls.

    Attributes:
        length: cavity length L [m].
        eps_r: relative permittivity of the filling (lossless part).
        sigma: electric conductivity of the filling [S/m].
    """

    length: float
    eps_r: float = 1.0
    sigma: float = 0.0

    @property
    def wave_speed(self) -> float:
        return C0 / math.sqrt(self.eps_r)

    def wavenumber(self, n: int) -> float:
        return n * math.pi / self.length

    def mode_angular_frequency(self, n: int) -> float:
        """Analytical lossless eigenfrequency omega_n = n pi c / (L sqrt(eps_r))."""
        return self.wavenumber(n) * self.wave_speed

    def energy_decay_rate(self) -> float:
        """Energy decay rate gamma = sigma / eps for uniform conductivity [1/s]."""
        return self.sigma / (self.eps_r * EPS0)

    def quality_factor(self, n: int) -> float:
        """Q_n = omega_n / gamma (math.inf when lossless)."""
        gamma = self.energy_decay_rate()
        return math.inf if gamma == 0.0 else self.mode_angular_frequency(n) / gamma
