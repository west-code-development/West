import numpy as np
import json


def read_trans_matrix_from_json(fileName):
    """
    Read unitary transformation matrix from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    return np.array(raw_["output"]["B"]["K000001"]["trans_matrix"], dtype=float)


def read_wannier_center_from_json(fileName):
    """
    Read Wannier centers from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    return np.array(raw_["output"]["B"]["K000001"]["wan_center"], dtype=float)


def test_trans_matrix():
    """
    Test unitary transformation matrix.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_trans = read_trans_matrix_from_json("./test015/ref/westpp.json")
    test_trans = read_trans_matrix_from_json("./test015/test.westpp.save/westpp.json")

    np.testing.assert_almost_equal(
        np.abs(ref_trans),
        np.abs(test_trans),
        decimal=-np.log10(float(parameters["tolerance"]["localization"])),
    )


def test_wannier_center():
    """
    Test Wannier centers.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_wanc = read_wannier_center_from_json("./test015/ref/westpp.json")
    test_wanc = read_wannier_center_from_json("./test015/test.westpp.save/westpp.json")

    np.testing.assert_almost_equal(
        ref_wanc,
        test_wanc,
        decimal=-np.log10(float(parameters["tolerance"]["localization"])),
    )
