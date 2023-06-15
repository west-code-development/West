import numpy as np
import json
import pytest
from xml.etree import ElementTree as ET

######################
# REUSABLE FUNCTIONS #
######################


def read_total_energy(fileName):
    """
    Reads the PW total energy from file
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

    maxDiff = np.amax(np.abs(test_energy-ref_energy))
    print(f'Total energy (pwscf) max diff: {maxDiff}')

    assert np.allclose(test_energy,ref_energy,rtol=0,atol=tol),'Total energies changed'


def read_wstat_eigenvalues(fileName):
    """
    Reads the eigenvalues in wstat
    """

    with open(fileName,'r') as f:
        data = json.load(f)

    pdep_eig = {}
    if 'Q1' in data['exec']['davitr']:
        for iq in data['exec']['davitr']:
            pdep_eig[iq] = np.array(data['exec']['davitr'][iq][-1]['ev'],dtype='f8')
    else:
        pdep_eig['Q1'] = np.array(data['exec']['davitr'][-1]['ev'],dtype='f8')

    return pdep_eig


def read_and_test_wstat_eigenvalues(fileA,fileB,tol):
    """
    Reads and tests PDEP eigenvalues
    """

    test_pdep_eig = read_wstat_eigenvalues(fileA)
    ref_pdep_eig = read_wstat_eigenvalues(fileB)

    maxDiff = 0.0
    for iq in test_pdep_eig:
        maxDiff = max(maxDiff,np.amax(np.abs(test_pdep_eig[iq]-ref_pdep_eig[iq])))
    print(f'PDEP eigenvalues (wstat) max diff: {maxDiff}')

    for iq in test_pdep_eig:
        assert np.allclose(test_pdep_eig[iq],ref_pdep_eig[iq],rtol=0,atol=tol),f'PDEP eigenvalues changed, iq {iq}'


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
                if key != 'occupation':
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
            maxDiff = max(maxDiff,np.amax(np.abs(np.abs(test_en[ik][key])-np.abs(ref_en[ik][key]))))
    print(f'Single-particle energy (wfreq) max diff: {maxDiff}')

    for ik in test_en:
        for key in test_en[ik]:
            assert np.allclose(np.abs(test_en[ik][key]),np.abs(ref_en[ik][key]),rtol=0,atol=tol),f'Single-particle energies changed, ik {ik}, field {key}'


#########
# TESTS #
#########


@pytest.mark.parametrize('testdir',['test001','test002','test003','test004','test005','test006','test007','test008','test009','test010','test011','test012','test013','test014','test015','test016','test017','test018','test019','test020','test021','test022','test023','test024','test025', 'test026'])
def test_totalEnergy(testdir):
    with open('parameters.json','r') as f:
        parameters = json.load(f)
    read_and_test_total_energies(testdir+'/test.save/data-file-schema.xml',testdir+'/ref/pw.xml',float(parameters['tolerance']['total_energy']))


@pytest.mark.parametrize('testdir',['test001','test002','test003','test004','test005','test006','test007','test009','test010','test011','test012','test013','test014','test016','test017','test020'])
def test_pdepEigen(testdir):
    with open('parameters.json','r') as f:
        parameters = json.load(f)
    read_and_test_wstat_eigenvalues(testdir+'/test.wstat.save/wstat.json',testdir+'/ref/wstat.json',float(parameters['tolerance']['pdep_eigenvalue']))


@pytest.mark.parametrize('testdir',['test001','test002','test003','test004','test006','test007','test009','test010','test011','test012','test013','test014','test020'])
def test_singleparticleEnergy(testdir):
    with open('parameters.json','r') as f:
        parameters = json.load(f)
    read_and_test_wfreq_energies(testdir+'/test.wfreq.save/wfreq.json',testdir+'/ref/wfreq.json',float(parameters['tolerance']['singleparticle_energy']))
