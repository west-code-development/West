import numpy as np
import json


def read_indexmap_from_json(filename):
    """
    Read indexmap between pairs (i,j) and the Kohn-Sham states i and j from JSON
    files.
    """
    with open(filename, "r") as f:
        raw_ = json.load(f)

    return np.array(raw_["output"]["Q"]["indexmap"], dtype=int)


def read_eri_from_json(keyword, filename):
    """
    Read different kind of two-body terms (screened, partially screened, and
    unscreened) from JSON file.
    """
    with open(filename, "r") as f:
        raw_ = json.load(f)

    indexmap = np.array(raw_["output"]["Q"]["indexmap"], dtype=int)
    n_pairs = len(indexmap)

    eri = np.zeros((n_pairs, n_pairs))

    for i in range(n_pairs):
        string = "pair" + format(i + 1, "06d")
        eri[i, :] = np.array(
            raw_["qdet"][keyword]["K000001"]["K000001"][string], dtype=float
        )

    return eri


def test_indexmap():
    """
    Test length and content of indexmap
    """

    ref_map = read_indexmap_from_json("./test013/ref/wfreq.json")
    test_map = read_indexmap_from_json("./test013/test.wfreq.save/wfreq.json")

    np.testing.assert_array_equal(ref_map, test_map)


def test_eri():
    """
    Test whether matrix elements of the screened Coulomb potential are correct.
    """
    # get parameters from JSON file
    with open("./parameters.json", "r") as f:
        parameters = json.load(f)

    for keyword in ["eri_w", "eri_w_full", "eri_vc"]:

        ref_eri = read_eri_from_json(keyword, "./test013/ref/wfreq.json")
        test_eri = read_eri_from_json(keyword, "./test013/test.wfreq.save/wfreq.json")

        np.testing.assert_almost_equal(
            ref_eri,
            test_eri,
            decimal=np.log10(float(parameters["tolerance"]["pdep_eigenvalue"])),
        )
