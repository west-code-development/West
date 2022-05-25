import numpy as np 
 
exec(open("../common_tests.py").read())


if __name__ == "__main__":
    #
    # Get parameters
    #
    with open('../test_parameters.json', 'r') as f:
       parameters = json.load(f)
    #
    # Test files 
    #
    test_files_exist_and_job_done(['pw.out','wstat.out','wfreq.out'])
    #
    # Test total energies
    #
    read_and_test_total_energies('test.save/data-file-schema.xml','ref/data-file-schema.xml',float(parameters['tolerance']['total_energy']))
    #
    # Test PDEP
    #
    # TDB
    #
    # Test Wfreq
    #
    read_and_test_wfreq_energies('ref/wfreq.json','ref/wfreq.json',float(parameters['tolerance']['singleparticle_energy']))
