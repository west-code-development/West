import numpy as np
import json


def read_proj_matrix_from_json(fileName):
    """
    Read projection matrix from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    proj = {}
    nexc = raw_["input"]["westpp_control"]["westpp_n_liouville_to_use"]

    for iexc in range(nexc):
        label = f"E{(iexc+1):06d}"
        proj[label] = np.array(
            raw_["output"][label]["K000001"]["projection"]["vals"], dtype=float
        )

    return proj


def test_proj_matrix():
    """
    Test projection matrix.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_proj = read_proj_matrix_from_json("./test027/ref/westpp.json")
    test_proj = read_proj_matrix_from_json("./test027/test.westpp.save/westpp.json")

    for iexc in ref_proj:
        np.testing.assert_almost_equal(
            np.abs(ref_proj[iexc]),
            np.abs(test_proj[iexc]),
            decimal=-np.log10(float(parameters["tolerance"]["localization"])),
        )
