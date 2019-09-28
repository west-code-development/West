#!/usr/bin/python3

from time import sleep, perf_counter as pc
from os import path, remove

def server_input_file(client_lockfile) : 
    #
    # Determine the name of the server file given the name of the lockfile
    #
    client_image = client_lockfile.split(".")[1] 
    return f"qb.{client_image}.in" # we assume that server_number = client_image


def create_input_file_for_server(client_lockfile) :
    #
    # List of perturbation files 
    #
    perturbation_list = []
    with open(client_lockfile,"r") as f:
       for cnt, line in enumerate(f):
          perturbation_list.append(line)
    #
    # Create the input file for the server 
    # 
    with open(server_input_file(client_lockfile),"w") as f: 
       f.write("load gs.xml\n")
       for pert in perturbation_list : 
           f.write(f"response -vext filename {pert}")

 
def before_sleep(client_lockfile) :
    #
    create_input_file_for_server(client_lockfile)
    #
    # Awake server, by removing its lockfile 
    #
    if(path.exists(server_input_file(client_lockfile)+".lock")) :
       remove(server_input_file(client_lockfile)+".lock")


def awake_condition(client_lockfile) : 
    #
    # If server gets to sleeps, awake the client 
    #
    if( path.exists(server_input_file(client_lockfile)+".lock")) :  
       remove(client_lockfile)

def after_sleep() : 
    #
    # HACK here
    #
    pass
 


def sleep_and_wait_for_lock_to_be_removed(*args, **kwargs):
    #
    # Name of lockfile
    #
    client_lockfile = args[0]
    #
    # Max sleep time (in s)
    # 
    maxsec = 12 * 60 * 60 # 12 hours 
    if "maxsec" in kwargs.keys() : 
       maxsec = kwargs["maxsec"]
    #
    # Sleep interval (in s) 
    #
    sleepsec = 1 # 1 second
    if "sleepsec" in kwargs.keys() : 
       sleepsec = kwargs["sleepsec"]
    #
    # ====================================
    before_sleep(client_lockfile)
    # ====================================
    #
    awake = 1 # I am awake if this is zero
    t0 = pc()
    while (pc()-t0 <= maxsec) :
       #
       # =================================
       awake_condition(client_lockfile)
       # =================================
       #
       exists = path.exists(client_lockfile)
       if (not exists) : 
          awake = 0 
          break
       else : 
          sleep(sleepsec)
    #
    # ====================================
    after_sleep()
    # ====================================
    #
    return awake
    

def test() :
    with open("I.1.lock","w") as f :
       f.write(" ")
    sleep_and_wait_for_lock_to_be_removed("I.1.lock",maxsec=60,sleepsec=2)

if __name__ == "__main__":
    # execute only if run as a script
    test()

