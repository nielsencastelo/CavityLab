"""CavityLab: reproducible, energy-audited simulation of electromagnetic and quantum cavities.

The core is physics-neutral. It does not assume net energy generation from the vacuum;
it provides infrastructure to simulate, reproduce, audit and falsify hypotheses.
"""

from cavitylab.core.constants import C0, EPS0, ETA0, MU0
from cavitylab.energy.audit import EnergyAudit
from cavitylab.geometry.cavity1d import Cavity1D
from cavitylab.solvers.fdtd1d import FDTD1D, Modulation, SoftSource

__version__ = "0.1.0.dev0"

__all__ = [
    "C0",
    "EPS0",
    "ETA0",
    "MU0",
    "Cavity1D",
    "EnergyAudit",
    "FDTD1D",
    "Modulation",
    "SoftSource",
    "__version__",
]
