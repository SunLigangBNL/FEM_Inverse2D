#! /usr/bin/env python3

# refkeys={'cd': 40, 'h': 100, 'swaleft': 88, 'swaright': 88, 'r_topleft': 5.0,
#           'r_topright':5.0, 'r_botleft':5.0, 'r_botright':5.0, 'pitch': 100, 
#          'epsilonr_scat_re': np.real(epsr_W), 'epsilonr_scat_im':np.imag(epsr_W),
#         'lambda0': lambda0, 'msl': msl, 'dwfacx': dwfacxref, 'dwfacz': dwfaczref}


# some initial values
layout_geometry_keys = {
    "pitch": 80,
    "cd": 30, 
    "h": 30, 
    "swaleft": 88,
    "swaright": 89,
    "r_topleft": 5.0,
    "r_topleft": 5.0,
    "r_botleft": 10.0,
    "r_botleft": 10.0,
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
    "epsilonr_scat_im": 3.278431066866658e-5  # real part of epsr for scatterer.
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
