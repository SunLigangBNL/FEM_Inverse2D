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
    obindexmask: np.ndarray # indeices of those peaks (which will be included in the output intmat) in the default qxrefarray.
    flipswitch: int=0 # 1 means geometry is flipped downward (X-ray hitting the substrate first).
    geometrymode: int=0
    # Remark of geometry model:
    # 0: direct polygon model defined in config.
    # 1: symmetric single trapezoid model with rounding.
    # 2: asymmetric single trapezoid model with rounding.
    # 3: asymmetric stacked trapezoids model with rounding.
    decayswitch: int=0 # apply DW factors or not.
    NLswitch: int=0 # normalization switch to set I(Qx=deltaqx,Qz=0)=1.

    @property
    def deltaqx(self) -> float:
        return 2*np.pi/self.pitch

    @property
    def centerindex(self) -> int:
        return self.dimqxbranch
        
