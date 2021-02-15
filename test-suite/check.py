import json
import os
import shutil
import sys

activeTests = [1,2,3,4,5]
tolerance = 0.0001

def check_pw(prefix):
    # Check "JOB DONE"
    fileName = prefix+'/pw.out'

    if not os.path.isfile(fileName):
        print(prefix+': pw.out not found')
        return 1

    ok = False
    test = 0.0
    with open(fileName,'r') as f:
        for line in f:
            if 'JOB DONE' in line:
                ok = True

            if '!    total energy' in line:
                test = float(line.split()[4])

    if not ok:
        print(prefix+': pwscf failed')
        return 1

    # Remove stderr (may exceed artifact size limit)
    fileName = prefix+'/pw.err'

    if os.path.isfile(fileName):
        os.remove(fileName)

    # Load reference
    fileName = prefix+'/ref/pw.ref'

    with open(fileName,'r') as f:
        for line in f:
            if '!    total energy' in line:
                ref = float(line.split()[4])

    # Total energy error
    absDiff = abs(ref-test)

    if absDiff > tolerance:
        print(prefix+': pwscf bad (diff '+"{:.2e}".format(absDiff)+')')
        return 1
    else:
        print(prefix+': pwscf good (diff '+"{:.2e}".format(absDiff)+')')
        return 0

def check_wstat(prefix):
    # Check "JOB DONE"
    fileName = prefix+'/wstat.out'

    if not os.path.isfile(fileName):
        print(prefix+': wstat.out not found')
        return 1

    ok = False
    with open(fileName,'r') as f:
        for line in f:
            if 'JOB DONE' in line:
                ok = True
                break

    if not ok:
        print(prefix+': wstat failed')
        return 1

    # Remove stderr (may exceed artifact size limit)
    fileName = prefix+'/wstat.err'

    if os.path.isfile(fileName):
        os.remove(fileName)

    # Load reference
    fileName = prefix+'/ref/summary.json'

    with open(fileName,'r') as f:
        jsonData = json.load(f)

    # Load test
    fileName = prefix+'/test.wstat.save/summary.json'

    if not os.path.isfile(fileName):
        print(prefix+': wstat output not found')
        return 1

    ok = True
    maxDiff = 0.0

    with open(fileName,'r') as f:
        jsonData2 = json.load(f)

    for ii in range(len(jsonData['dielectric_matrix']['pdep'])):
        ref = jsonData['dielectric_matrix']['pdep'][ii]['eigenval']
        test = jsonData2['dielectric_matrix']['pdep'][ii]['eigenval']

        for jj in range(len(ref)):
            absDiff = abs(ref[jj]-test[jj])
            maxDiff = max(maxDiff,absDiff)

            if absDiff > tolerance:
                ok = False

    if ok:
        print(prefix+': wstat good (max diff '+"{:.2e}".format(maxDiff)+')')
        return 0
    else:
        print(prefix+': wstat bad (max diff '+"{:.2e}".format(maxDiff)+')')
        return 1

def check_wfreq(prefix):
    # Check "JOB DONE"
    fileName = prefix+'/wfreq.out'

    if not os.path.isfile(fileName):
        print(prefix+': wfreq.out not found')
        return 1

    ok = False
    with open(fileName,'r') as f:
        for line in f:
            if 'JOB DONE' in line:
                ok = True
                break

    if not ok:
        print(prefix+': wfreq failed')
        return 1

    # Remove stderr (may exceed artifact size limit)
    fileName = prefix+'/wfreq.err'

    if os.path.isfile(fileName):
        os.remove(fileName)

    # Load reference
    fileName = prefix+'/ref/wfreq.json'

    with open(fileName,'r') as f:
        jsonData = json.load(f)

    # Load test
    fileName = prefix+'/test.wfreq.save/wfreq.json'

    if not os.path.isfile(fileName):
        print(prefix+': wfreq output not found')
        return 1

    ok = True
    maxDiff = 0.0

    with open(fileName,'r') as f:
        jsonData2 = json.load(f)

    for ii in range(jsonData['system']['electron']['nkstot']):
        ref = jsonData['output']['Q']['K'+str(ii+1).zfill(6)]['sigmax']
        test = jsonData2['output']['Q']['K'+str(ii+1).zfill(6)]['sigmax']

        for jj in range(len(ref)):
            absDiff = abs(ref[jj]-test[jj])
            maxDiff = max(maxDiff,absDiff)

            if absDiff > tolerance:
                ok = False

        ref = jsonData['output']['Q']['K'+str(ii+1).zfill(6)]['sigmac_eqpSec']['re']
        test = jsonData2['output']['Q']['K'+str(ii+1).zfill(6)]['sigmac_eqpSec']['re']

        for jj in range(len(ref)):
            absDiff = abs(ref[jj]-test[jj])
            maxDiff = max(maxDiff,absDiff)

            if absDiff > tolerance:
                ok = False

        ref = jsonData['output']['Q']['K'+str(ii+1).zfill(6)]['sigmac_eqpSec']['im']
        test = jsonData2['output']['Q']['K'+str(ii+1).zfill(6)]['sigmac_eqpSec']['im']

        for jj in range(len(ref)):
            absDiff = abs(ref[jj]-test[jj])
            maxDiff = max(maxDiff,absDiff)

            if absDiff > tolerance:
                ok = False

    if ok:
        print(prefix+': wfreq good (max diff '+"{:.2e}".format(maxDiff)+')')
        return 0
    else:
        print(prefix+': wfreq bad (max diff '+"{:.2e}".format(maxDiff)+')')
        return 1

if __name__ == "__main__":
    ok = True

    for ii in activeTests:
        print()

        prefix = 'test'+str(ii).zfill(3)

        err = check_pw(prefix)
        if err != 0:
            ok = False

        err = check_wstat(prefix)
        if err != 0:
            ok = False

        err = check_wfreq(prefix)
        if err != 0:
            ok = False

    if not ok:
        sys.exit(1)
