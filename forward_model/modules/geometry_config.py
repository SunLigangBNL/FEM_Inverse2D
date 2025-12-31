"""
This file contains the layout of the background medium and the scatterers in CDSAXS.
"""
# Coordinate system: sample frame.
# Unit: SI.
# (1) z-coordinates of N layered medium which are stacked along +z direction.
# (2) Scatterers are embeded in layer L.
# (3) x and z coordinates of the vertices of those scatterers are stored.

import numpy as np
from dataclasses import dataclass

@dataclass
class Scatterer:
    vertexmat: np.ndarray  # shape (N, 2), each row is [x, z]
    def xarray(self) -> np.ndarray:
        return self.vertexmat[:, 0]
    def zarray(self) -> np.ndarray:
        return self.vertexmat[:, 1]

@dataclass
class GeometryConfig:
    layerzarray: list[float]          # z-coordinates of all layer interfaces.
    scatterlayer: list[int]             # index of the layers that contain scatterers.
    scattercoordmat: list[list[Scatterer]]    # all Scatterer objects in that single layer

    @property
    def layerindexarray(self) -> list[int]:
        return list(range(len(self.layerzarray)+1))

        