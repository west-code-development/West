#!/usr/bin/python3

#
# Copyright (C) 2015-2025 M. Govoni
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
# This file is part of WEST.
#

import sys
from os import path, remove
import yaml
import json

rytoev = 13.6056980659  # Ry to eV conversion factor

#########################
# STATIC DEFAULT VALUES #
#########################

default = {}
# input_west
default["input_west"] = {}
default["input_west"]["qe_prefix"] = "pwscf"
default["input_west"]["west_prefix"] = "west"
default["input_west"]["outdir"] = "./"
# wstat_control
default["wstat_control"] = {}
default["wstat_control"]["wstat_calculation"] = "S"
default["wstat_control"]["n_pdep_eigen"] = 1  # dynamically set to the number of electrons
default["wstat_control"]["n_pdep_times"] = 4
default["wstat_control"]["n_pdep_maxiter"] = 100
default["wstat_control"]["n_dfpt_maxiter"] = 250
default["wstat_control"]["n_pdep_read_from_file"] = 0
default["wstat_control"]["n_steps_write_restart"] = 1
default["wstat_control"]["trev_pdep"] = 1.0e-3
default["wstat_control"]["trev_pdep_rel"] = 1.0e-1
default["wstat_control"]["tr2_dfpt"] = 1.0e-16
default["wstat_control"]["l_kinetic_only"] = False
default["wstat_control"]["l_minimize_exx_if_active"] = False
default["wstat_control"]["n_exx_lowrank"] = 0  # dynamically set to the number of bands
default["wstat_control"]["l_use_ecutrho"] = False
default["wstat_control"]["qlist"] = [1]  # dynamically set to the actual number of q
# wfreq_control
default["wfreq_control"] = {}
default["wfreq_control"]["wfreq_calculation"] = "XWGQ"
default["wfreq_control"]["n_pdep_eigen_to_use"] = 1  # dynamically set to the number of electrons
default["wfreq_control"]["qp_bandrange"] = [1, 2]
default["wfreq_control"]["qp_bands"] = [0]
default["wfreq_control"]["macropol_calculation"] = "C"
default["wfreq_control"]["n_lanczos"] = 30
default["wfreq_control"]["n_imfreq"] = 128
default["wfreq_control"]["n_refreq"] = 272
default["wfreq_control"]["ecut_imfreq"] = 25.0  # dynamically set to ecutrho
default["wfreq_control"]["ecut_refreq"] = 2.0
default["wfreq_control"]["wfreq_eta"] = 0.05 / rytoev
default["wfreq_control"]["n_secant_maxiter"] = 21
default["wfreq_control"]["trev_secant"] = 0.05 / rytoev
default["wfreq_control"]["l_enable_lanczos"] = True
default["wfreq_control"]["l_qdet_verbose"] = False
default["wfreq_control"]["l_enable_off_diagonal"] = False
default["wfreq_control"]["ecut_spectralf"] = [-2.0, 1.0]
default["wfreq_control"]["n_spectralf"] = 204
# westpp_control
default["westpp_control"] = {}
default["westpp_control"]["westpp_calculation"] = "R"
default["westpp_control"]["westpp_range"] = [1, 2]
default["westpp_control"]["westpp_format"] = "C"
default["westpp_control"]["westpp_sign"] = False
default["westpp_control"]["westpp_n_pdep_eigen_to_use"] = 1
default["westpp_control"]["westpp_r0"] = [0.0, 0.0, 0.0]
default["westpp_control"]["westpp_nr"] = 100
default["westpp_control"]["westpp_rmax"] = 1.0
default["westpp_control"]["westpp_epsinfty"] = 1.0
default["westpp_control"]["westpp_box"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
default["westpp_control"]["westpp_n_liouville_to_use"] = 1
default["westpp_control"]["westpp_l_spin_flip"] = False
default["westpp_control"]["westpp_l_compute_tdm"] = False
default["westpp_control"]["westpp_wannier_tr_rel"] = 1.0e-6
default["westpp_control"]["westpp_l_dipole_realspace"] = False
# server_control
default["server_control"] = {}
default["server_control"]["document"] = "{}"
# wbse_init_control
default["wbse_init_control"] = {}
default["wbse_init_control"]["wbse_init_calculation"] = "S"
default["wbse_init_control"]["solver"] = "BSE"
default["wbse_init_control"]["bse_method"] = "PDEP"
default["wbse_init_control"]["n_pdep_eigen_to_use"] = 1  # dynamically set to the number of electrons
default["wbse_init_control"]["localization"] = "N"
default["wbse_init_control"]["wannier_tr_rel"] = 1.0e-6
default["wbse_init_control"]["wfc_from_qbox"] = "qb_wfc"
default["wbse_init_control"]["bisection_info"] = "bis_info"
default["wbse_init_control"]["chi_kernel"] = "CHI"
default["wbse_init_control"]["overlap_thr"] = 0.0
default["wbse_init_control"]["spin_channel"] = 0
default["wbse_init_control"]["n_trunc_bands"] = 0
# wbse_control
default["wbse_control"] = {}
default["wbse_control"]["wbse_calculation"] = "D"
default["wbse_control"]["qp_correction"] = ""
default["wbse_control"]["scissor_ope"] = 0.0
default["wbse_control"]["n_liouville_eigen"] = 1
default["wbse_control"]["n_liouville_times"] = 4
default["wbse_control"]["n_liouville_maxiter"] = 100
default["wbse_control"]["n_liouville_read_from_file"] = 0
default["wbse_control"]["trev_liouville"] = 0.001
default["wbse_control"]["trev_liouville_rel"] = 0.1
default["wbse_control"]["n_lanczos"] = 1000
default["wbse_control"]["n_steps_write_restart"] = 100
default["wbse_control"]["wbse_ipol"] = "XX"
default["wbse_control"]["l_dipole_realspace"] = False
default["wbse_control"]["wbse_epsinfty"] = 1.0
default["wbse_control"]["spin_excitation"] = "S"
default["wbse_control"]["l_preconditioning"] = True
default["wbse_control"]["l_pre_shift"] = False
default["wbse_control"]["l_spin_flip"] = False
default["wbse_control"]["l_spin_flip_kernel"] = False
default["wbse_control"]["l_spin_flip_alda0"] = False
default["wbse_control"]["l_print_spin_flip_kernel"] = False
default["wbse_control"]["spin_flip_cut"] = 1.0e-3
default["wbse_control"]["l_forces"] = False
default["wbse_control"]["forces_state"] = 1
default["wbse_control"]["forces_zeq_cg_tr"] = 1.0e-10
default["wbse_control"]["forces_zeq_n_cg_maxiter"] = 500
default["wbse_control"]["ddvxc_fd_coeff"] = 0.01
default["wbse_control"]["forces_inexact_krylov"] = 0
default["wbse_control"]["forces_inexact_krylov_tr"] = 1.0e-16
default["wbse_control"]["l_minimize_exx_if_active"] = False
default["wbse_control"]["n_exx_lowrank"] = 0  # dynamically set to the number of bands

############################
# DYNAMICAL DEFAULT VALUES #
############################


def update_default_values(key, kwargs):
    assert key in default.keys()
    #
    if key == "wstat_control":
        #
        assert ("nq") in kwargs.keys()
        nq = kwargs["nq"]
        default[key]["qlist"] = [i + 1 for i in range(nq)]
        #
        assert ("nelec") in kwargs.keys()
        nelec = kwargs["nelec"]
        default[key]["n_pdep_eigen"] = int(nelec)
        #
        assert ("nbnd") in kwargs.keys()
        nbnd = kwargs["nbnd"]
        default[key]["n_exx_lowrank"] = int(nbnd)
    #
    if key == "wfreq_control":
        #
        assert ("nelec") in kwargs.keys()
        nelec = kwargs["nelec"]
        default[key]["n_pdep_eigen_to_use"] = int(nelec)
        #
        assert ("ecutrho") in kwargs.keys()
        ecutrho = kwargs["ecutrho"]
        default[key]["ecut_imfreq"] = ecutrho
    #
    if key == "wbse_init_control":
        #
        assert ("nelec") in kwargs.keys()
        nelec = kwargs["nelec"]
        default[key]["n_pdep_eigen_to_use"] = int(nelec)
    #
    if key == "wbse_control":
        #
        assert ("nbnd") in kwargs.keys()
        nbnd = kwargs["nbnd"]
        default[key]["n_exx_lowrank"] = int(nbnd)


################
# OPEN & PARSE #
################


def open_and_parse_file(fileName="west.in"):
    """Opens a file and parses it using the YAML syntax.

    :param fileName: name of the file
    :type fileName: ``string``
    :return: parsed data
    :rtype: ``dict``

    """
    try:
        with open(fileName, "r") as stream:
            try:
                data = yaml.load(stream, Loader=yaml.SafeLoader)
            except:
                print(f"Cannot parse file: {fileName}")
    except:
        print(f"Cannot open file: {fileName}")
    #
    # Stop if input is not read successfully
    #
    assert isinstance(data, dict)
    #
    if "server_control" in data.keys():
        if "document" in data["server_control"].keys():
            jsonText = json.dumps(data["server_control"]["document"])
            data["server_control"]["document"] = jsonText
    #
    return data


##############
# CHECK DICT #
##############


def check_dict(parsed_data={}, default_data={}):
    """Check data: returns a dictionary with the same keys of default_data. If keys are matching, values of default_data are replaced with those of parsed_data.

    :param parsed_data: parsed data
    :type parsed_data: ``dict``
    :param default_data: default data
    :type default_data: ``dict``
    :return: checked data
    :rtype: ``dict``

    """
    #
    data = {}
    #
    for key in default_data.keys():
        if key in parsed_data.keys():
            data[key] = parsed_data[key]
        else:
            data[key] = default_data[key]
    #
    return data


############
# QP BANDS #
############


def parse_qp_bands(data_in, keyword, kwargs):
    """Parse qp_bands.

    :param data_in: input data
    :type data_in: ``dict``
    :param keyword: keyword
    :type keyword: ``string``
    :param kwargs: kwargs dictionary
    :type kwargs: ``dict``
    :return: updated data
    :rtype: ``dict``
    """
    #
    if keyword == "wfreq_control":
        #
        assert ("nspin") in kwargs.keys()
        nspin = kwargs["nspin"]
        qp_bandrange = data_in["qp_bandrange"]
        qp_bands = data_in["qp_bands"]
        qp_bands_new = []
        #
        if qp_bands == [0]:
            qp_bands_new = [list(range(qp_bandrange[0], qp_bandrange[1] + 1))] * nspin
        else:
            if isinstance(qp_bands[0], list):
                if nspin == 1:
                    qp_bands_new = qp_bands
                elif nspin == 2:
                    if len(qp_bands) == 1:
                        qp_bands_new = qp_bands * nspin
                    else:
                        assert isinstance(qp_bands[1], list)
                        assert len(qp_bands[0]) == len(qp_bands[1])
                        qp_bands_new = qp_bands
            else:
                qp_bands_new = [qp_bands] * nspin
        #
        data_in["qp_bands"] = qp_bands_new
    #
    return data_in


###########
# SUPPORT #
###########


def print_bar(prefix="", nmarks=92):
    """Prints bar.

    :param prefix: prefix
    :type prefix: ``string``
    :param nmarks: number of marks
    :type nmarks: ``int``
    """
    #
    s = prefix
    for i in range(nmarks):
        s += "-"
    print(s)


#########
# PRINT #
#########


def print_dict(title="input_west", data={}):
    """Prints data.

    :param title: title
    :type title: ``string``
    :param data: data to print
    :type data: ``dict``

    """
    #
    nmarks = 92
    nspaces = 5
    s = ""
    for i in range(nspaces):
        s += " "
    #
    print("")
    print_bar(s, nmarks)
    print(s + "I/O Summary : " + str(title))
    print_bar(s, nmarks)
    for key in data.keys():
        print(s + key, ":", data[key])
    print_bar(s, nmarks)
    sys.stdout.flush()


def print_json(fileName="west.json", title="input_west", data={}):
    """Prints data.

    :param fileName: name of the file
    :type title: ``string``
    :param title: title
    :type title: ``string``
    :param data: data to print
    :type data: ``dict``

    """
    #
    try:
        with open(fileName, "r") as f:
            j = json.load(f)
    except:
        j = {}
    #
    if not "input" in j.keys():
        j["input"] = {}
    j["input"][title] = data
    #
    with open(fileName, "w") as f:
        json.dump(j, f, indent=2)


#############
# INTERFACE #
#############


def read_keyword_from_file(*args, **kwargs):
    """Read keyword from file.

    :return: read data
    :rtype: ``dict``

    """
    #
    fileName = args[0]
    keyword = args[1]
    verbose = args[2]
    logfile = args[3]
    #
    # Assign static & dynamical defaults
    #
    default_data = default[keyword]
    update_default_values(keyword, kwargs)
    #
    # Read input file
    #
    input_data = open_and_parse_file(fileName)
    parsed_data = {}
    if keyword in input_data.keys():
        parsed_data = input_data[keyword]
    #
    # Compare defaults and input variables
    #
    data = check_dict(parsed_data, default_data)
    #
    # Parse qp_bandrange and qp_bands
    #
    data = parse_qp_bands(data, keyword, kwargs)
    #
    # Print
    #
    if verbose:
        print_dict(keyword, data)
        print_json(logfile, keyword, data)
    #
    return data


########
# TEST #
########


def test():
    #
    fileName = "west.in"
    jsonName = "west.json"
    #
    with open(fileName, "w") as file:
        file.write(
            """
input_west :
   qe_prefix : molecule
   west_prefix : molecule
   outdir : "./"
wstat_control :
   wstat_calculation : R # this is a comment
   unknown_key : value # this line will be read but not passed
server_control :
   document : {}
"""
        )
    #
    read_keyword_from_file(fileName, "input_west", True, jsonName)
    read_keyword_from_file(fileName, "wstat_control", True, jsonName, nq=20, nelec=10, nbnd=30)
    read_keyword_from_file(fileName, "wfreq_control", True, jsonName, nelec=10, ecutrho=30.0, nspin=2)
    read_keyword_from_file(fileName, "westpp_control", True, jsonName)
    read_keyword_from_file(fileName, "server_control", True, jsonName)
    read_keyword_from_file(fileName, "wbse_init_control", True, jsonName, nelec=10)
    read_keyword_from_file(fileName, "wbse_control", True, jsonName, nbnd=30)
    #
    remove(fileName)
    remove(jsonName)


if __name__ == "__main__":
    # execute only if run as a script
    test()
