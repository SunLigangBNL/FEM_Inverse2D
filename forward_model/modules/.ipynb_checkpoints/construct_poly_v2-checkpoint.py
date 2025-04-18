#! /usr/bin/env python3
"""Constructs a series of polygons from input data.

"""
import numpy as np


def polygons(keys):
    """Construct the material line and substrate and the computational domain.

    :param keys: 

    """
    cd = keys["cd"]
    h = keys["h"]
    swa = keys["swa"]
    #t = keys["t"]
    r_top = keys["r_top"]
    r_bot = keys["r_bot"]

    # offsets
    o1 = keys["height_offset_substrate"] # 10
    o2 = keys["height_offset_air"] # 30

    # computational domain:
    comp_dom = [
        -.5 * keys["pitch"], -o1,
        .5 * keys["pitch"], -o1,
        .5 * keys["pitch"], o2,
        -.5 * keys["pitch"], o2,
    ]

    # Line domain:
    y_ox_0 = 0
    y_ox_h = h

    yc = 0.5 * (y_ox_h - y_ox_0)
    hc1 = y_ox_0 - yc
    hc2 = y_ox_h - yc
    dx1 = hc1 * np.tan(np.deg2rad(90 - swa))
    dx2 = hc2 * np.tan(np.deg2rad(90 - swa))

    x1 = -.5 * keys["pitch"]
    x2 = -.5 * cd - dx2
    x3 = -.5 * cd - dx1
    x4 = .5 * cd + dx1
    x5 = .5 * cd + dx2
    x6 = .5 * keys["pitch"]

    line_poly = [
        x6, -o1,
        x6, y_ox_0,
        x5, y_ox_0,
        x4, y_ox_h,
        x3, y_ox_h,
        x2, y_ox_0,
        x1, y_ox_0,
        x1, -o1
    ]

    substrate_poly=[
        x6, -o1,
        x6, y_ox_0,
        x1, y_ox_0,
        x1, -o1
    ]

    # # regular material layer
    # y_0 = y_ox_0 - t
    # y_h = y_ox_h - t

    # angle_d = 0.5 * (180 - swa)
    # dx = t / np.tan(np.deg2rad(angle_d))

    # x1 = x1
    # x2 = x2 + dx
    # x3 = x3 + dx
    # x4 = x4 - dx
    # x5 = x5 - dx
    # x6 = x6

    # material_poly = [
    #     x6, -o1,
    #     x6, y_0,
    #     x5, y_0,
    #     x4, y_h,
    #     x3, y_h,
    #     x2, y_0,
    #     x1, y_0,
    #     x1, -o1
    # ]

    return substrate_poly, line_poly, comp_dom


if __name__ == "__main__":
    keys = dict()
    keys["cd"] = 40
    keys["h"] = 30
    keys["swa"] = 86
    #keys["t"] = 4.9393
    keys["r_top"] = 10
    keys["r_bot"] = 5
    keys["pitch"] = 80
    keys["height_offset_substrate"] = 10
    keys["height_offset_air"] = 40

    m, o, c = polygons(keys)
    print("mat: {}".format(m))
    print("ox: {}".format(o))
