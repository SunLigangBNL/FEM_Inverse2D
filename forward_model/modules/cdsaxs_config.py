"""
This file is the main CDSAXS configration file, it calls the source/geometry/material sub-configration files.
"""
# coordinate system: sample frame.
# unit: SI.

import numpy as np
from dataclasses import dataclass
from .source_config   import SourceConfig
from .geometry_config import GeometryConfig
from .material_config import MaterialConfig

@dataclass
class CDSAXSConfig:
    source: SourceConfig
    geometry: GeometryConfig
    material: MaterialConfig

    # only the two extra CDSAXS‑specific params for now:
    pitch: float # unit: m.
    dimqxbranch: int
    qxpeakarray: np.ndarray

    @property
    def deltaqx(self) -> float:
        return 2*np.pi/self.pitch

    @property
    def centerindex(self) -> int:
        return self.dimqxbranch
        
