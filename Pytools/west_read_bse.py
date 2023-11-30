#!/usr/bin/python3

#
# Copyright (C) 2015-2023 M. Govoni
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
# This file is part of WEST.
#

import h5py
import numpy as np


def bse_read_qb_param(*args, **kwargs):
    #
    fileName = args[0]
    #
    data = {}
    #
    with h5py.File(fileName, "r") as f:
        data["nwfcs"] = f["wfcs"].attrs.get("nwfcs")
        data["nx"] = f["wfcs"].attrs.get("nx")
        data["ny"] = f["wfcs"].attrs.get("ny")
        data["nz"] = f["wfcs"].attrs.get("nz")
    # return dictionary
    return data


def bse_read_qb_wfc(*args, **kwargs):
    #
    fileName = args[0]
    iwfc = int(args[1])
    #
    with h5py.File(fileName, "r") as f:
        wfc = f[f"wfcs/wfc{iwfc}"][:]
    # return numpy ndarray
    return wfc


def test():
    #
    fileName = "qb_wfc.1"
    nwfcs = 4
    nx = 30
    ny = 30
    nz = 30
    data = np.ones(nx * ny * nz, dtype="float64")
    #
    with h5py.File(fileName, "w") as f:
        # metadata
        wfcs = f.create_group("wfcs")
        wfcs.attrs.create("nwfcs", nwfcs)
        wfcs.attrs.create("nx", nx)
        wfcs.attrs.create("ny", ny)
        wfcs.attrs.create("nz", nz)
        # data
        for iwfc in range(nwfcs):
            wfcs.create_dataset(f"wfc{iwfc+1}", data=data)
    #
    print(bse_read_qb_param(fileName))
    print(bse_read_qb_wfc(fileName, 1))


if __name__ == "__main__":
    # execute only if run as a script
    test()
