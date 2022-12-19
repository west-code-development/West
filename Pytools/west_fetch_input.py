#!/usr/bin/python3

#
# Copyright (C) 2015-2022 M. Govoni
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

rytoev = 13.6056980659 # Ry to eV conversion factor

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
default["wstat_control"]["n_pdep_eigen"] = 1 # dynamically set to the number of electrons
default["wstat_control"]["n_pdep_times"] = 4
default["wstat_control"]["n_pdep_maxiter"] = 100
default["wstat_control"]["n_dfpt_maxiter"] = 250
default["wstat_control"]["n_pdep_read_from_file"] = 0
default["wstat_control"]["n_steps_write_restart"] = 1
default["wstat_control"]["trev_pdep"] = 1.e-3
default["wstat_control"]["trev_pdep_rel"] = 1.e-1
default["wstat_control"]["tr2_dfpt"] = 1.e-12
default["wstat_control"]["l_kinetic_only"] = False
default["wstat_control"]["l_minimize_exx_if_active"] = False
default["wstat_control"]["l_use_ecutrho"] = False
default["wstat_control"]["qlist"] = [ 1 ] # dynamically set to the actual number of q
# wfreq_control
default["wfreq_control"] = {}
default["wfreq_control"]["wfreq_calculation"] = "XWGQ"
default["wfreq_control"]["n_pdep_eigen_to_use"] = 1 # dynamically set to the number of electrons
default["wfreq_control"]["qp_bandrange"] = [1, 2]
default["wfreq_control"]["qp_bands"] = [0]
default["wfreq_control"]["macropol_calculation"] = "N"
default["wfreq_control"]["n_lanczos"] = 30
default["wfreq_control"]["n_imfreq"] = 128
default["wfreq_control"]["n_refreq"] = 272
default["wfreq_control"]["ecut_imfreq"] = 25. # dynamically set to ecutrho
default["wfreq_control"]["ecut_refreq"] = 2.
default["wfreq_control"]["wfreq_eta"] = 0.05 / rytoev
default["wfreq_control"]["n_secant_maxiter"] = 21
default["wfreq_control"]["trev_secant"] = 0.05 / rytoev
default["wfreq_control"]["l_enable_lanczos"] = True
default["wfreq_control"]["l_qdet_verbose"] = False
default["wfreq_control"]["l_enable_off_diagonal"] = False
default["wfreq_control"]["o_restart_time"] = 0.
default["wfreq_control"]["ecut_spectralf"] = [-2., 1.]
default["wfreq_control"]["n_spectralf"] = 204
# westpp_control
default["westpp_control"] = {}
default["westpp_control"]["westpp_calculation"] = "R"
default["westpp_control"]["westpp_range"] = [1, 2]
default["westpp_control"]["westpp_format"] = "C"
default["westpp_control"]["westpp_sign"] = False
default["westpp_control"]["westpp_n_pdep_eigen_to_use"] = 1
default["westpp_control"]["westpp_r0"] = [0., 0., 0.]
default["westpp_control"]["westpp_nr"] = 100
default["westpp_control"]["westpp_rmax"] = 1.
default["westpp_control"]["westpp_epsinfty"] = 1.
default["westpp_control"]["westpp_box"] = [0., 0., 0., 0., 0., 0.]
default["westpp_control"]["westpp_n_liouville_to_use"] = 1
# server_control
default["server_control"] = {}
default["server_control"]["document"] = "{}"
# wbse_init_control
default["wbse_init_control"] = {}
default["wbse_init_control"]["wbse_init_calculation"] = "S"
default["wbse_init_control"]["localization"] = "N"
default["wbse_init_control"]["wfc_from_qbox"] = "qb_wfc"
default["wbse_init_control"]["bisection_info"] = "bis_info"
default["wbse_init_control"]["chi_kernel"] = "CHI"
default["wbse_init_control"]["overlap_thr"] = 0.
default["wbse_init_control"]["spin_channel"] = 0
# wbse_control
default["wbse_control"] = {}
default["wbse_control"]["wbse_calculation"] = "D"
default["wbse_control"]["solver"] = "BSE"
default["wbse_control"]["qp_correction"] = ""
default["wbse_control"]["scissor_ope"] = 0.
default["wbse_control"]["n_liouville_eigen"] = 1
default["wbse_control"]["n_liouville_times"] = 4
default["wbse_control"]["n_liouville_maxiter"] = 100
default["wbse_control"]["n_liouville_read_from_file"] = 0
default["wbse_control"]["trev_liouville"] = 0.001
default["wbse_control"]["trev_liouville_rel"] = 0.1
default["wbse_control"]["n_lanczos"] = 1000
default["wbse_control"]["n_steps_write_restart"] = 100
default["wbse_control"]["wbse_ipol"] = "XX"
default["wbse_control"]["macropol_calculation"] = "N"
default["wbse_control"]["wbse_epsinfty"] = 1.
default["wbse_control"]["spin_excitation"] = "S"
default["wbse_control"]["l_preconditioning"] = False
default["wbse_control"]["l_reduce_io"] = True

############################
# DYNAMICAL DEFAULT VALUES #
############################

def update_default_values(key,kwargs) :
    assert key in default.keys()
    #
    if key == "wstat_control" :
       #
       assert("nq") in kwargs.keys()
       nq = kwargs["nq"]
       default[key]["qlist"] = [ i+1 for i in range(nq) ]
       #
       assert("nelec") in kwargs.keys()
       nelec = kwargs["nelec"]
       default[key]["n_pdep_eigen"] = int(nelec)
    #
    if key == "wfreq_control" :
       #
       assert("nelec") in kwargs.keys()
       nelec = kwargs["nelec"]
       default[key]["n_pdep_eigen_to_use"] = int(nelec)
       #
       assert("ecutrho") in kwargs.keys()
       ecutrho = kwargs["ecutrho"]
       default[key]["ecut_imfreq"] = ecutrho

################
# OPEN & PARSE #
################

def open_and_parse_file(fileName="west.in") :
    """Opens a file and parses it using the YAML sintax

    :param fileName: name of the file
    :type fileName: ``string``
    :return: parsed data
    :rtype: ``dict``

    """
    data = {}
    try :
       with open(fileName, "r") as stream:
           try:
              data = yaml.load(stream,Loader=yaml.SafeLoader)
           except:
              print("Cannot parse file")
    except :
       print("Cannot open file : ",fileName)
    #
    if "server_control" in data.keys() :
       if "document" in data["server_control"].keys() :
          jsonText = json.dumps(data["server_control"]["document"])
          data["server_control"]["document"] = jsonText
    #
    return data

##############
# CHECK DICT #
##############

def check_dict(parsed_data={}, default_data={}) :
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
    for key in default_data.keys() :
        if key in parsed_data.keys() :
           data[key] = parsed_data[key]
        else :
           data[key] = default_data[key]
    #
    return data

###########
# SUPPORT #
###########

def print_bar(prefix="",nmarks=92) :
    """Prints bar.

    :param prefix: prefix
    :type prefix: ``string``
    :param nmarks: number of marks
    :type nmarks: ``int``
    """
    #
    s = prefix
    for i in range(nmarks) :
       s+="-"
    print(s)

#########
# PRINT #
#########

def print_dict(title="input_west", data={}) :
    """Prints data.

    :param title: title
    :type title: ``string``
    :param data: data to print
    :type default_data: ``dict``

    """
    #
    nmarks = 92
    nspaces = 5
    s = ""
    for i in range(nspaces) :
       s+=" "
    #
    print("")
    print_bar(s,nmarks)
    print(s+"I/O Summary : "+str(title))
    print_bar(s,nmarks)
    for key in data.keys() :
       print(s+key,":",data[key])
    print_bar(s,nmarks)
    sys.stdout.flush()

#############
# INTERFACE #
#############

def read_keyword_from_file(*args, **kwargs):
    """Read keyword from file

    :return: read data
    :rtype: ``dict``

    """
    #
    fileName = args[0]
    keyword = args[1]
    verbose = args[2]
    #
    # Assign static & dynamical defaults
    #
    default_data = default[keyword]
    update_default_values(keyword,kwargs)
    #
    # Read input file
    #
    input_data = open_and_parse_file(fileName)
    parsed_data = {}
    if keyword in input_data.keys() :
       parsed_data = input_data[keyword]
    #
    # Compare defaults and input variables
    #
    data = check_dict( parsed_data, default_data )
    #
    # Print
    #
    if (verbose) :
       print_dict(keyword, data)
    #
    return data

########
# TEST #
########

def test() :
    #
    fileName = "west.in"
    #
    with open(fileName, "w") as file :
       file.write("""
input_west :
   qe_prefix : molecule
   west_prefix : molecule
   outdir : "./"
wstat_control :
   wstat_calculation : R # this is a comment
   unknown_key : value # this line will be read but not passed
server_control :
   document : {}
""")
    #
    read_keyword_from_file(fileName,"input_west",True)
    read_keyword_from_file(fileName,"wstat_control",True,nq=20,nelec=10)
    read_keyword_from_file(fileName,"wfreq_control",True,nelec=10,ecutrho=30.)
    read_keyword_from_file(fileName,"westpp_control",True)
    read_keyword_from_file(fileName,"server_control",True)
    read_keyword_from_file(fileName,"wbse_init_control",True)
    read_keyword_from_file(fileName,"wbse_control",True)
    #
    remove(fileName)

if __name__ == "__main__":
    # execute only if run as a script
    test()

