import numpy as np
import json


def read_dipole_from_json(fileName):
    """
    Read transition dipole moments from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    dip = {}
    for key in raw_["output"]["D"]["K000001"]["dipole"]:
        dip[key] = np.array(raw_["output"]["D"]["K000001"]["dipole"][key], dtype=float)

    return dip


def test_dipole():
    """
    Test transition dipole moments.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_dip = read_dipole_from_json("./test008/ref/westpp.json")
    test_dip = read_dipole_from_json("./test008/test.westpp.save/westpp.json")

    for key in ref_dip:
        np.testing.assert_almost_equal(
            np.abs(ref_dip[key]),
            np.abs(test_dip[key]),
            decimal=-np.log10(float(parameters["tolerance"]["westpp"])),
        )
