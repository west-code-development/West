#!/usr/bin/python

# File: fetch_input.py

from __future__ import print_function
import sys
import os
import yaml
import json

def open_and_parse_file(fileName) :
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
    return data

def check_data(keyword, parsed_data, default) : 
    #
    data = {}
    #
    if keyword in parsed_data.keys() : 
       for key in default.keys() : 
          if key in parsed_data[keyword].keys() : 
             data[key] = parsed_data[keyword][key]
          else : 
             data[key] = default[key]
    else : 
       print("Missing keyword in input file : ",keyword)
    return data

def read_input_west(*args, **kwargs):
    #
    fileName = args[0]
    #
    default = {}
    default["qe_prefix"] = "pwscf"
    default["west_prefix"] = "west"
    default["outdir"] = "./"
    #
    parsed_data = open_and_parse_file(fileName)
    #
    checked_data = check_data( "input_west", parsed_data, default )
    print( checked_data )
    return checked_data

def read_wstat_control(*args, **kwargs):
    #
    fileName = args[0]
    nq = args[1]
    #
    default = {}
    default["wstat_calculation"] = "S"
    default["n_pdep_eigen"] = 1
    default["n_pdep_times"] = 4
    default["n_pdep_maxiter"] = 100
    default["n_dfpt_maxiter"] = 250
    default["n_pdep_read_from_file"] = 0
    default["trev_pdep"] = 1.e-3
    default["trev_pdep_rel"] = 1.e-1
    default["tr2_dfpt"] = 1.e-12
    default["l_kinetic_only"] = False
    default["l_minimize_exx_if_active"] = False
    default["l_use_ecutrho"] = False
    default["qlist"] = [ i+1 for i in range(nq) ]
    #
    parsed_data = open_and_parse_file(fileName)
    #
    checked_data = check_data( "wstat_control", parsed_data, default )
    print( checked_data )
    return checked_data
