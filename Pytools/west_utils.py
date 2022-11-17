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

