#!/usr/bin/python3

from __future__ import print_function
from os import mkdir 
import sys

#############
# INTERFACE #
#############

def my_mkdir(*args, **kwargs):
    #
    path = args[0]
    #
    try:
        mkdir(path)
    except OSError:
        #print (f"Creation of the directory {path} failed")
        pass
    else:
        #print (f"Successfully created the directory {path} ")
        pass
    sys.stdout.flush()

########
# TEST #
########

def test() :
    #
    dirname = "./wstat.save"
    #
    my_mkdir(dirname)

if __name__ == "__main__":
    # execute only if run as a script
    test()

