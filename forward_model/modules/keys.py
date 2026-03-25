#! /usr/bin/env python3

# some initial values
layout_geometry_keys = {
    "pitch": 80,
    "cd": 30, 
    "h": 30, 
    "swal": 88,
    "swar": 89,
    "rtl": 2.0,
    "rtr": 2.0,
    "rbl": 2.0,
    "rbr": 2.0,
    "t": 1.0,
    "thickness_Cr":1.0,
}

layout_keys = {
    "msl": 1, # maximum side length.
    "height_offset_bot": 1, # always positive.
    "height_offset_top": 1, # always positive.
    "rad_upper_n": 5, # must be integer.
    "rad_lower_n": 5, # must be integer.
}

sources_keys = {
    "theta": 0,
    "phi": 0,
    "lambda0": 1e-9
}

materials_keys = {
    "epsilonr_sub": 9.99435131564694e-1,   # epsr for substrate.
    "epsilonr_scat_re": 9.99435131564694e-1,  # real part of epsr for scatterer.
    "epsilonr_scat_im": 3.278431066866658e-5,  # imaginary part of epsr for scatterer.
    "epsilonr_thinlayer_re": 9.99435131564694e-1,  # real part of epsr for thin layer.
    "epsilonr_thinlayer_im": 3.278431066866658e-5  # imaginary part of epsr for thin layer.
}

project_keys = {
    "numerical_aperture": 1,
    "dwfacx": 0,
    "dwfacz": 0,
}

keys = dict()

keys.update(layout_geometry_keys)
keys.update(layout_keys)
keys.update(sources_keys)
keys.update(materials_keys)
keys.update(project_keys)
