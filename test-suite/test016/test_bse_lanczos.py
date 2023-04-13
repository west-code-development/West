import numpy as np
import json


def read_beta_from_json(fileName):
    """
    Read beta values from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    return np.array(raw_["output"]["lanczos"]["K000001"]["XX"]["beta"], dtype=float)


def read_zeta_from_json(fileName):
    """
    Read zeta values from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    return np.array(raw_["output"]["lanczos"]["K000001"]["XX"]["zeta"], dtype=float)


def test_beta():
    """
    Test beta values.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_beta = read_beta_from_json("./test016/ref/wbse.json")
    test_beta = read_beta_from_json("./test016/test.wbse.save/wbse.json")

    np.testing.assert_almost_equal(
        ref_beta,
        test_beta,
        decimal=-np.log10(float(parameters["tolerance"]["bse"])),
    )


def test_zeta():
    """
    Test zeta values.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_zeta = read_zeta_from_json("./test016/ref/wbse.json")
    test_zeta = read_zeta_from_json("./test016/test.wbse.save/wbse.json")

    np.testing.assert_almost_equal(
        ref_zeta,
        test_zeta,
        decimal=-np.log10(float(parameters["tolerance"]["bse"])),
    )
