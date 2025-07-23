#! /usr/bin/env python3

# some initial values
layout_geometry_keys = {
    "pitch": 80,
    "cd": 30, 
    "h": 30, 
    "swa": 71.56505117707799,
    "r_top": 8.0,
    "r_bot": 5.0
}

layout_keys = {
    "msl": 1, # maximum side length.
    "height_offset_substrate": 1, # always positive.
    "height_offset_air": 1, # always positive.
    "rad_upper_n": 5, # must be integer.
    "rad_lower_n": 5, # must be integer.
}

sources_keys = {
    "theta": 0,
    "phi": 0,
    "lambda0": 1e-9
}

materials_keys = {
    "epsilonr_sub": 0.999435,   # epsr for substrate.
    "epsilonr_scat": 9.99435131564694e-1+3.278431066866658e-5*1j  # epsr for scatterer. Note the Python format here!!
}

project_keys = {
    "numerical_aperture": 1
}


keys = dict()

keys.update(layout_geometry_keys)
keys.update(layout_keys)
keys.update(sources_keys)
keys.update(materials_keys)
keys.update(project_keys)
