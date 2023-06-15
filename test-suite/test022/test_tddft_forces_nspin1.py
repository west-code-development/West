import numpy as np
import json


def read_wbse_forces(fileName):
    """
    Read forces from JSON file.
    """

    with open(fileName, "r") as f:
        raw_ = json.load(f)

    forces_drhox1 = np.array(raw_["exec"]["forces"]["forces_drhox1"], dtype=float)
    forces_drhox2 = np.array(raw_["exec"]["forces"]["forces_drhox2"], dtype=float)
    forces_drhoz = np.array(raw_["exec"]["forces"]["forces_drhoz"], dtype=float)
    forces_total = np.array(raw_["exec"]["forces"]["forces_total"], dtype=float)
    forces_collect = np.array(raw_["exec"]["forces"]["forces_corrected"], dtype=float)

    forces = np.array([forces_drhox1, forces_drhox2, forces_drhoz, forces_total, forces_collect])

    return forces


def test_wbse_forces():
    """
    Test eigenvalues from JSON file.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    ref_forces = read_wbse_forces("./test022/ref/wbse.json")
    test_forces = read_wbse_forces("./test022/test.wbse.save/wbse.json")

    np.testing.assert_almost_equal(
        ref_forces,
        test_forces,
        decimal=-np.log10(float(parameters["tolerance"]["forces"])),
    )
