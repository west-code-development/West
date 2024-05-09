import numpy as np
import json


def read_localization_from_json(fileName):
    """
    Read localization factor from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    return np.array(raw_["output"]["L"]["K000001"]["local_factor"], dtype=float)


def read_ipr_from_json(fileName):
    """
    Read inverse participation ratio (IPR) from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    return np.array(raw_["output"]["L"]["K000001"]["ipr"], dtype=float)


def test_localization():
    """
    Test localization factor.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_loc = read_localization_from_json("./test008/ref/westpp.json")
    test_loc = read_localization_from_json("./test008/test.westpp.save/westpp.json")

    np.testing.assert_almost_equal(
        ref_loc,
        test_loc,
        decimal=-np.log10(float(parameters["tolerance"]["localization"])),
    )


def test_ipr():
    """
    Test IPR.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_ipr = read_ipr_from_json("./test008/ref/westpp.json")
    test_ipr = read_ipr_from_json("./test008/test.westpp.save/westpp.json")

    np.testing.assert_almost_equal(
        ref_ipr,
        test_ipr,
        decimal=-np.log10(float(parameters["tolerance"]["localization"])),
    )
