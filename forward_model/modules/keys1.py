#! /usr/bin/env python3


# some initial values
layout_geometry_keys = {
    "pitch": 84, #80
    "cd": 42, #40
    "h": 32, #30
    "swa": 75, #71.5651
    #"t": 5,
    "r_top": 9,
    "r_bot": 6
}

layout_keys = {
    "msl": 1,
    "height_offset_substrate": 10, # always positive.
    "height_offset_air": 40,
    "rad_upper_n": 5,
    "rad_lower_n": 5,
}

project_keys = {
    "info_level": 3,
    "cur_deg": 1,
    "fem_deg": 4,
    "precision": 1e-3,
}

sources_keys = {
    "theta": 45,
    "phi": 0,
    #"vacuum_wavelength": 266e-9
}

materials_keys = {
    "epsilonr_1": 0.999435,   # substrate: Si.
    "epsilonr_2": 0.99964+0.0000295399*1j,  # line: PMMA. Note the format here!!
    "epsilonr_3": 1,                     # Air.
}

keys = dict()

keys.update(layout_geometry_keys)
keys.update(layout_keys)
keys.update(sources_keys)
keys.update(materials_keys)
keys.update(project_keys)
