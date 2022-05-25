import numpy as np
import json
import os
import pytest
from xml.etree import ElementTree as ET


def check_files_exist_and_job_done(list_of_files):
    """
    Check that every file in the list exists and contains the string JOB DONE
    """

    for file in list_of_files:
        assert os.path.isfile(file),f'{file} is missing'
        with open(file,'r') as f:
            assert 'JOB DONE' in f.read(),f'{file} does not have the string JOB DONE'


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
    maxDiff = np.amax(np.abs(np.subtract(test_energy,ref_energy)))
    print(f'\npwscf max diff: {maxDiff}')

    assert np.allclose(test_energy,ref_energy,rtol=0,atol=tol),'Total energies changed'


def read_wstat_eigenvalues(fileName):
    """
    Reads the eigenvalues in wstat
    """

    with open(fileName,'r') as f:
        data = json.load(f)

    eig = {}
    if 'Q1' in data['exec']['davitr']:
        for iq in data['exec']['davitr']:
            eig[iq] = np.array(data['exec']['davitr'][iq][-1]['ev'],dtype='f8')
    else:
        eig['Q1'] = np.array(data['exec']['davitr'][-1]['ev'],dtype='f8')

    return eig


def read_and_test_wstat_eigenvalues(fileA,fileB,tol):
    """
    Reads and tests PDEP eigenvalues
    """

    test_eig = read_wstat_eigenvalues(fileA)
    ref_eig = read_wstat_eigenvalues(fileB)

    maxDiff = 0.0
    for iq in test_eig:
        maxDiff = max(maxDiff,np.amax(np.abs(np.subtract(test_eig[iq],ref_eig[iq]))))
    print(f'wstat max diff: {maxDiff}')

    for iq in test_eig:
        assert np.allclose(test_eig[iq],ref_eig[iq],rtol=0,atol=tol),f'PDEP eigenvalues changed, iq {iq}'


def read_wfreq_energies(fileName):
    """
    Reads the energies in wfreq
    """

    with open(fileName,'r') as f:
        data = json.load(f)

    en = {}
    for ik in range(1,data['system']['electron']['nkstot']+1):
        kindex = f'K{ik:06d}'
        en[ik] = {}
        for key in data['output']['Q'][kindex]:
            if 'im' in data['output']['Q'][kindex][key] and 're' in data['output']['Q'][kindex][key]:
                en[ik][key] = np.array(data['output']['Q'][kindex][key]['re'],dtype='c16')
                en[ik][key] += 1j * np.array(data['output']['Q'][kindex][key]['im'],dtype='c16')
            else:
                en[ik][key] = np.array(data['output']['Q'][kindex][key],dtype='f8')

    return en


def read_and_test_wfreq_energies(fileA,fileB,tol):
    """
    Reads and tests single-particle energies
    """

    test_en = read_wfreq_energies(fileA)
    ref_en = read_wfreq_energies(fileB)

    maxDiff = 0.0
    for ik in test_en:
        for key in test_en[ik]:
            maxDiff = max(maxDiff,np.amax(np.abs(np.subtract(test_en[ik][key],ref_en[ik][key]))))
    print(f'wfreq max diff: {maxDiff}')

    for ik in test_en:
        for key in test_en[ik]:
            assert np.allclose(test_en[ik][key],ref_en[ik][key],rtol=0,atol=tol),f'Single-particle energies changed, ik {ik}, field {key}'


@pytest.mark.parametrize('testcase',['test001','test002','test003','test004','test005','test006','test007'])
def test_west(testcase):

    with open('test_parameters.json','r') as f:
        parameters = json.load(f)

    check_files_exist_and_job_done([testcase+'/pw.out',testcase+'/wstat.out',testcase+'/wfreq.out'])

    read_and_test_total_energies(testcase+'/test.save/data-file-schema.xml',testcase+'/ref/data-file-schema.xml',float(parameters['tolerance']['total_energy']))

    read_and_test_wstat_eigenvalues(testcase+'/test.wstat.save/wstat.json',testcase+'/ref/wstat.json',float(parameters['tolerance']['pdep_eigenvalue']))

    read_and_test_wfreq_energies(testcase+'/test.wfreq.save/wfreq.json',testcase+'/ref/wfreq.json',float(parameters['tolerance']['singleparticle_energy']))
