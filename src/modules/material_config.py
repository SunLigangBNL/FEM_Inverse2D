"""
This file contains all the material parameters to set up the CDSAXS problem.
"""
# coordinate system: sample frame.

from dataclasses import dataclass
import numpy as np

@dataclass
class MaterialConfig:
    layerepsrarray: list[complex] # ordered along +x direction.
    scatterepsrmat: list[list[complex]] # ordered along (+z, +x) direction.

    def __post_init__(self):
        # Ensure there is at least one layer and (optionally) at least one scatterer
        assert len(self.layerepsrarray) >= 1, "Must have at least one background layer."

    @property
    def layerrefindexarray(self) -> list[complex]:
        return [np.sqrt(eps) for eps in self.layerepsrarray]

    @property
    def scatrefindexarray(self) -> list[list[complex]]:
        return [
            [np.sqrt(eps) for eps in row]
            for row in self.scatterepsrmat
        ]
