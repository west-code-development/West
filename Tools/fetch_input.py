#!/usr/bin/python3

# File: fetch_input.py
# Test: python3 fetch_input.py

from __future__ import print_function
import sys
import os
import yaml
import json

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
default["wstat_control"]["n_pdep_eigen"] = 1
default["wstat_control"]["n_pdep_times"] = 4
default["wstat_control"]["n_pdep_maxiter"] = 100
default["wstat_control"]["n_dfpt_maxiter"] = 250
default["wstat_control"]["n_pdep_read_from_file"] = 0
default["wstat_control"]["trev_pdep"] = 1.e-3
default["wstat_control"]["trev_pdep_rel"] = 1.e-1
default["wstat_control"]["tr2_dfpt"] = 1.e-12
default["wstat_control"]["l_kinetic_only"] = False
default["wstat_control"]["l_minimize_exx_if_active"] = False
default["wstat_control"]["l_use_ecutrho"] = False
default["wstat_control"]["qlist"] = [ 1 ]

############################
# DYNAMICAL DEFAULT VALUES #
############################

def update_default_values(key,kwargs) :
    assert key in default.keys()
    #
    if key == "wstat_control" :
       assert("nq") in kwargs.keys()
       nq = kwargs["nq"] 
       default[key]["qlist"] = [ i+1 for i in range(nq) ]

def open_and_parse_file(fileName="west.in") :
    """Opens a file and parses it using the YAML sintax 

    :param fileName: name of the file
    :type fileName: ``string``
    :return: parsed data
    :rtype: ``dict``

    """
    data = {}
    try : 
       with open(fileName, 'r') as stream:
           try:
              data = yaml.load(stream,Loader=yaml.SafeLoader)
           except:
              print("Cannot parse file")
    except : 
       print("Cannot open file : ",fileName)
    #
    return data

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
    print_bar(s,nmarks)
    print(s+"I/O Summary : "+str(title))
    print_bar(s,nmarks)
    for key in data.keys() :
       print(s+key,"=",data[key])
    print_bar(s,nmarks)
    sys.stdout.flush()

def read_keyword_from_file(*args, **kwargs):
    """Read keyword from file  

    :return: read data
    :rtype: ``dict``

    """
    #
    fileName = args[0]
    keyword = args[1] 
    #
    parsed_data = open_and_parse_file(fileName)[keyword]
    default_data = default[keyword]
    #
    update_default_values(keyword,kwargs) 
    #
    data = check_dict( parsed_data, default_data )
    #
    print_dict(keyword, data)
    #
    return data

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
   wstat_calculation : R
""")
    #
    read_keyword_from_file(fileName,"input_west")
    read_keyword_from_file(fileName,"wstat_control",nq=20)
    #
    os.remove(fileName)

if __name__ == "__main__":
    # execute only if run as a script
    test()

