"""
This file contains the definition of the light source used in CDSAXS.
"""
# For now, we only consider P-polarized planewave.

from dataclasses import dataclass
import numpy as np

@dataclass
class SourceConfig:
    thetaarray: list[float]
    lambda0: float
    #E0: float # This is amplitude. For P-polarized planewave, incident electric field Ei=E0*(1,0,0).

    @property
    def k0(self) -> float:
        return 2 * np.pi / self.lambda0
