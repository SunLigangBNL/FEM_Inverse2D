#! /usr/bin/env python3

import scipy
from scipy import io
import numpy as np
import pathlib


def load_data():
    """Load the data from the matlab file."""
    datapath = (
        pathlib.Path(__file__).absolute().parent / "../data/measdata.mat"
    )
    data = scipy.io.loadmat(datapath)
    phis = [0, 0, 90, 90]
    pols = ["P", "S", "P", "S"]

    out_data = list()

    for ii, phi, pol in zip(range(len(phis)), phis, pols):
        thetas = data["data"][0][ii][1:, 0]
        thetas = np.round(thetas)
        intensities = data["data"][0][ii][1:, 1]

        data_dict = dict()
        data_dict["phi"] = phi
        data_dict["polarization"] = pol
        data_dict["thetas"] = thetas
        data_dict["intensities"] = intensities

        out_data.append(data_dict)

    return out_data
