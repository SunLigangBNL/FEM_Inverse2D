#! /usr/bin/env python3
"""
Constructs a series of polygons from input data.
"""
import numpy as np

def polygons(keys):
    """
    Construct the material scatterers and the computational domain.
    """
    cd = keys["cd"]
    h = keys["h"]
    swa = keys["swa"]
    r_top = keys["r_top"]
    r_bot = keys["r_bot"]

    # offsets
    o1 = keys["height_offset_substrate"]
    o2 = keys["height_offset_air"]

    y1=(-1)*o1
    y2=0
    y3=h
    y4=h+o2

    # computational domain
    comp_domain = [
        -0.5*keys["pitch"], y1,
        0.5*keys["pitch"], y1,
        0.5*keys["pitch"], y4,
        -0.5*keys["pitch"], y4,
    ]

    # scatterer domain.
    y_ox_0 = 0
    y_ox_h = h    
    
    dx=0.5*h*np.tan(np.deg2rad(90-swa)) # always positive.
    x1 = -0.5*keys["pitch"]
    x2 = -0.5*cd-dx
    x3 = -0.5*cd+dx
    x4 = 0.5*cd-dx
    x5 = 0.5*cd+dx
    x6 = 0.5*keys["pitch"]

    scatter_poly = [
    	x1, y1,
        x6, y1,
        x6, y2,
        x5, y2,
        x4, y3,
        x3, y3,
        x2, y2,
        x1, y2
    ]

    # substrate domain.
    substrate_poly=[
    	x1, y1,
        x6, y1,
        x6, y2,
        x1, y2
    ]

    return comp_domain, scatter_poly, substrate_poly 

