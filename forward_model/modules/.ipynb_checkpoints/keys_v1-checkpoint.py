#! /usr/bin/env python3


# some initial values
layout_geometry_keys = {
    "pitch": 80,
    "cd": 40,
    "h": 30,
    "swa": 86,
    #"t": 0,
    "r_top": 10,
    "r_bot": 5
}

layout_keys = {
    "msl": 10,
    "height_offset_substrate": 10, # no negative sign here.
    "height_offset_air": 40,
    "rad_upper_n": 8,
    "rad_lower_n": 8,
}

project_keys = {
    "info_level": 3,
    "cur_deg": 1,
    "fem_deg": 3,
    "precision": 1e-4,
}

sources_keys = {
    "theta": 45,
    "phi": 45
    #,
    #"vacuum_wavelength": 266e-9
}

materials_keys = {
    "epsilonr_1": 0.999435,   # substrate: Si
    "epsilonr_2": (0.99964, 0.0000295399),  # line: PMMA.
    "epsilonr_3": 1,  # Air
}


keys = dict()

keys.update(layout_geometry_keys)
keys.update(layout_keys)
keys.update(sources_keys)
keys.update(materials_keys)
keys.update(project_keys)
