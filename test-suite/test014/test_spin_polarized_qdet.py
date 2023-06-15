import numpy as np
import json


def read_parameters_from_JSON(filename):
    """
    Read basic calculation parameters from JSON file.
    """

    with open(filename, "r") as f:
        raw_ = json.load(f)

    indexmap = np.array(raw_["output"]["Q"]["indexmap"], dtype=int)

    npair = len(indexmap)
    nspin = int(raw_["system"]["electron"]["nspin"])
    bands = np.array(raw_["input"]["wfreq_control"]["qp_bands"], dtype=int)

    return nspin, npair, bands, indexmap


def read_one_body_terms_from_JSON(filename):
    """
    Read one-body terms from JSON file.
    """

    with open(filename, "r") as f:
        raw_ = json.load(f)

    # read parameters from JSON file
    indexmap = np.array(raw_["output"]["Q"]["indexmap"], dtype=int)
    npair = len(indexmap)
    nspin = int(raw_["system"]["electron"]["nspin"])
    bands = np.array(raw_["input"]["wfreq_control"]["qp_bands"], dtype=int)

    # allocate one-body terms in basis of pairs of KS states
    h1e_pair = np.zeros((nspin, npair))

    # read one-body terms from file
    for ispin in range(nspin):
        string1 = "K" + format(ispin + 1, "06d")
        h1e_pair[ispin, :] = np.array(raw_["qdet"]["h1e"][string1], dtype=float)

    return h1e_pair


def read_two_body_terms_from_JSON(filename):
    """
    Read two-body terms from JSON file.
    """

    with open(filename, "r") as f:
        raw_ = json.load(f)

    # read parameters from JSON file
    indexmap = np.array(raw_["output"]["Q"]["indexmap"], dtype=int)
    npair = len(indexmap)
    nspin = int(raw_["system"]["electron"]["nspin"])
    bands = np.array(raw_["input"]["wfreq_control"]["qp_bands"], dtype=int)

    # allocate one- and two-body terms in basis of pairs of KS states
    eri_pair = np.zeros((nspin, nspin, npair, npair))

    # read two-body terms from file
    for ispin1 in range(nspin):
        string1 = "K" + format(ispin1 + 1, "06d")
        for ispin2 in range(nspin):
            string2 = "K" + format(ispin2 + 1, "06d")

            for ipair in range(npair):
                string3 = "pair" + format(ipair + 1, "06d")
                eri_pair[ispin1, ispin2, ipair, :] = np.array(
                    raw_["qdet"]["eri_w"][string1][string2][string3], dtype=float
                )

    return eri_pair


def test_qdet_parameters():
    """
    Test indexmap and other parameters in JSON file/
    """

    ref_data = read_parameters_from_JSON("./test014/ref/wfreq.json")
    test_data = read_parameters_from_JSON("./test014/test.wfreq.save/wfreq.json")

    # test npair
    assert ref_data[0] == test_data[0]
    # test nspin
    assert ref_data[1] == test_data[1]
    # test bands
    np.testing.assert_array_equal(ref_data[2], test_data[2])
    # test indexmap
    np.testing.assert_array_equal(ref_data[3], test_data[3])


def test_qdet_one_body_terms():
    """
    Test one-body terms of spin-polarized QDET calculation.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    h1e_ref = read_one_body_terms_from_JSON("./test014/ref/wfreq.json")
    h1e_test = read_one_body_terms_from_JSON("./test014/test.wfreq.save/wfreq.json")

    np.testing.assert_almost_equal(
        np.abs(h1e_ref),
        np.abs(h1e_test),
        decimal=-np.log10(float(parameters["tolerance"]["pdep_eigenvalue"])),
    )


def test_qdet_two_body_terms():
    """
    Test two-body terms of spin-polarized QDET calculation.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    eri_ref = read_two_body_terms_from_JSON("./test014/ref/wfreq.json")
    eri_test = read_two_body_terms_from_JSON("./test014/test.wfreq.save/wfreq.json")

    np.testing.assert_almost_equal(
        np.abs(eri_ref),
        np.abs(eri_test),
        decimal=-np.log10(float(parameters["tolerance"]["pdep_eigenvalue"])),
    )
