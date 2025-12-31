#! /usr/bin/env python3
"""
Constructs a series of polygons from input data.
"""
import numpy as np

def polygons(keys):
    """
    Construct the material scatterers and the computational domain.
    """
    pitch=keys["pitch"]
    cd = keys["cd"]
    h = keys["h"]
    swaleft = keys["swaleft"]
    swaright = keys["swaright"]
    r_topleft = keys["r_topleft"]
    r_topright = keys["r_topright"]
    r_botleft = keys["r_botleft"]
    r_botright = keys["r_botright"]

    # offsets
    o1 = keys["height_offset_bot"]
    o2 = keys["height_offset_top"]

    y1=(-1)*o1
    y2=0
    y3=h
    y4=h+o2

    # computational domain
    comp_domain = [
        -0.5*pitch, y1,
        0.5*pitch, y1,
        0.5*pitch, y4,
        -0.5*pitch, y4,
    ]

    # scatterer domain.
    y_ox_0 = 0
    y_ox_h = h

    # coordinates of the primary trapezoid.
    shift1=h/(2*np.tan(np.radians(swaleft)))
    shift2=h/(2*np.tan(np.radians(swaright)))

    x1=-pitch/2   
    x2=-cd/2-shift1
    x3=-cd/2+shift1
    x4=cd/2-shift2
    x5=cd/2+shift2
    x6=pitch/2

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

