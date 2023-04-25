import numpy as np
import json


def read_wbse_eigenvalues(fileName):
    """
    Read eigenvalues from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    return np.array(raw_["exec"]["davitr"][-1]["ev"], dtype=float)


def test_wbse_eigenvalues():
    """
    Test eigenvalues from JSON file.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_bse_eig = read_wbse_eigenvalues("./test019/ref/wbse.json")
    test_bse_eig = read_wbse_eigenvalues("./test019/test.wbse.save/wbse.json")

    np.testing.assert_almost_equal(
        ref_bse_eig,
        test_bse_eig,
        decimal=-np.log10(float(parameters["tolerance"]["bse"])),
    )
