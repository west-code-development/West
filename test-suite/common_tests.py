import numpy as np 
import json
import os
from xml.etree import ElementTree as ET


def test_files_exist_and_job_done(list_of_files):
    """ 
    Check that every file in the list exists and contains the string JOB DONE
    """
    for file in list_of_files : 
        assert os.path.isfile(file), f"{file} is missing"
        with open(file,'r') as f:
           assert 'JOB DONE' in f.read(), f"{file} does not have the string JOB DONE"


def read_total_energy(fileName):
    """
    Reads the PW total energy from file.
    """

    xmlData = ET.parse(fileName)
    total_energy = np.array([(xmlData.find('output').find('total_energy').find('etot').text)],dtype='f8')

    return total_energy


def read_and_test_total_energies(fileA,fileB,tol):
    """
    Reads and tests the total energies
    """
    test_energy = read_total_energy(fileA)
    ref_energy = read_total_energy(fileB) 
    assert np.allclose(test_energy,ref_energy,rtol=0,atol=tol), "Total energies changed"


def read_wfreq_energies(fileName):
     """ reads the energies in wfreq
     """
     #
     with open(fileName, 'r') as f:
         data = json.load(f)

     en = {}
     for ik in range(1,data['system']['electron']['nkstot']+1):
         kindex = f"K{ik:06d}"
         en[ik] = {} 
         for key in data["output"]["Q"][kindex] :
            if 'im' in data["output"]["Q"][kindex][key] and 're' in data["output"]["Q"][kindex][key]: 
               en[ik][key] = np.array(data["output"]["Q"][kindex][key]["re"],dtype='c16')
               en[ik][key] += 1j * np.array(data["output"]["Q"][kindex][key]["im"],dtype='c16')
            else: 
               en[ik][key] = np.array(data["output"]["Q"][kindex][key],dtype='f8') 
     return en 



def read_and_test_wfreq_energies(fileA,fileB,tol):
    """
    Reads and tests single-particle energies
    """
    test_en = read_wfreq_energies(fileA)
    ref_en = read_wfreq_energies(fileB)
    for ik in test_en : 
       for key in test_en[ik] :
           assert np.allclose(test_en[ik][key],ref_en[ik][key],rtol=0,atol=tol), f"Single-particle energies changed, ik {ik}, field {key}"
