#!/usr/bin/python3

#
# Copyright (C) 2015-2024 M. Govoni
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
# This file is part of WEST.
#

from time import sleep, perf_counter as pc
from os import path, remove
from abc import ABC, abstractmethod
import json

##############
# SUPERCLASS #
##############


class ClientServer(ABC):
    #
    def __init__(self, client_lockfile, maxsec=21600, sleepsec=10, document={}):
        #
        self.client_lockfile = client_lockfile
        self.maxsec = maxsec
        self.sleepsec = sleepsec
        self.document = document
        super().__init__()

    #
    @abstractmethod
    def before_sleep(self):  # subclass needs to implement this method
        pass

    #
    @abstractmethod
    def after_sleep(self):  # subclass needs to implement this method
        pass

    #
    @abstractmethod
    def awake_condition(self):  # subclass needs to implement this method
        pass

    #
    def start(self):
        #
        # ====================================
        self.before_sleep()
        # ====================================
        #
        awake = 1  # I am awake if this is zero
        t0 = pc()
        while pc() - t0 <= self.maxsec:
            #
            # =================================
            self.awake_condition()
            # =================================
            #
            exists = path.exists(self.client_lockfile)
            if not exists:
                awake = 0
                break
            else:
                sleep(self.sleepsec)
        #
        # ====================================
        self.after_sleep()
        # ====================================
        return awake


###############
# SERVERCLASS #
###############


class QboxServer(ClientServer):
    #
    def before_sleep(self):
        #
        # Determine the name of the server file
        #
        client_image = self.client_lockfile.split(".")[1]
        self.server_inputfile = f"qb.{client_image}.in"  # we assume that server_number = client_image
        #
        # Prepare reponse commands
        #
        command_suffix = ""
        perturbation_list = []
        #
        if "response" in self.document.keys():
            response = self.document["response"]
            if "approximation" in response:
                if response["approximation"] == "RPA":
                    command_suffix += "-RPA "
                if response["approximation"] == "IPA":
                    command_suffix += "-IPA "
            if "amplitude" in response:
                command_suffix += f'-amplitude {response["amplitude"]} '
            if "nitscf" in response:
                command_suffix += f'{response["nitscf"]} '
            else:
                command_suffix += "100 "
            if "nite" in response:
                command_suffix += f'{response["nite"]}'
            #
            # Read list of perturbation files from lockfile
            #
            with open(self.client_lockfile, "r") as f:
                for cnt, line in enumerate(f):
                    perturbation_list.append(line.replace("\n", ""))
        #
        # Create the INPUT file for the server
        #
        with open(self.server_inputfile, "w") as f:
            #
            # First attempt to write SCRIPT
            #
            if "script" in self.document.keys():
                for line in self.document["script"]:
                    f.write(line + "\n")
            #
            # Then attempt to write RESPONSE commands (many)
            #
            if "response" in self.document.keys():
                for pert in perturbation_list:
                    f.write(f"response -vext {pert} " + command_suffix + "\n")
        #
        # Awake server, by removing its lockfile
        #
        if path.exists(self.server_inputfile + ".lock"):
            remove(self.server_inputfile + ".lock")

    #
    def awake_condition(self):
        #
        # If server gets to sleeps, awake the client
        #
        if path.exists(self.server_inputfile + ".lock"):
            remove(self.client_lockfile)

    #
    def after_sleep(self):
        pass


#############
# INTERFACE #
#############


def sleep_and_wait(*args, **kwargs):
    #
    client_lockfile = args[0]  # name of client lockfile
    maxsec = 12 * 60 * 60  # 12 hours, Max sleep time (in s)
    sleepsec = 1  # 1 second, Sleep interval (in s)
    document = {}
    #
    # change defaults
    #
    if "maxsec" in kwargs.keys():
        maxsec = kwargs["maxsec"]
    if "sleepsec" in kwargs.keys():
        sleepsec = kwargs["sleepsec"]
    if "document" in kwargs.keys():
        document = json.loads(kwargs["document"])
    #
    #  consider_only allows to define a document with selected keys
    #
    if "consider_only" in kwargs.keys():
        consider_only_list = json.loads(kwargs["consider_only"])
        temp_document = {}
        for key in consider_only_list:
            if key in document.keys():
                temp_document[key] = document[key]
            else:
                print(f"Could not consider missing key: {key}")
        # replace document
        document = temp_document
    #
    server = QboxServer(client_lockfile, maxsec, sleepsec, document)
    return_int = server.start()
    #
    return return_int


########
# TEST #
########


def test():
    with open("I.1.lock", "w") as f:
        f.write("I.1_P.1.xml")
    sleep_and_wait(
        "I.1.lock",
        maxsec=60,
        sleepsec=2,
        document='{"response": {"amplitude": 0, "nitscf": 20, "nite": 0}}',
        consider_only='["response"]',
    )


if __name__ == "__main__":
    # execute only if run as a script
    test()
