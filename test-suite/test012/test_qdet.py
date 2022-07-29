import numpy as np 
import json 
import os 
import pytest

def read_h1e_from_json(filename):
    """
    Read QDET one-body terms from JSON file.
    """
    with open(filename, 'r') as f:
        raw_ = json.load(f)
    
    return np.array(raw_['qdet']['h1e']['K000001'], dtype=float)

def test_h1e():     
    """     
    Test whether matrix elements of the screened Coulomb potential are correct.     
    """     
    # get parameters from JSON file     
    with open('./parameters.json', 'r') as f:         
        parameters = json.load(f)     
    
    ref_h1e = read_h1e_from_json('./test012/ref/wfreq.json')     
    test_h1e = read_h1e_from_json('./test012/test.wfreq.save/wfreq.json')     
    
    np.testing.assert_almost_equal(ref_h1e, test_h1e,
            decimal=np.log10(float(parameters['tolerance']['pdep_eigenvalue'])))
