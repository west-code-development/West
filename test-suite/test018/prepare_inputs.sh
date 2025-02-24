#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/H_ONCV_PBE-1.0.upf
${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Si_ONCV_PBE-1.1.upf

cat > pw.in << EOF
&control
calculation  = 'scf'
restart_mode = 'from_scratch'
pseudo_dir   = './'
outdir       = './'
prefix       = 'test'
/
&system
ibrav           = 1
celldm(1)       = 20
nat             = 5
ntyp            = 2
ecutwfc         = 25.0
nbnd            = 30
assume_isolated = 'mp'
input_dft       = 'pbe0'
/
&electrons
diago_full_acc = .true.
/
ATOMIC_SPECIES
Si 28.0855  Si_ONCV_PBE-1.1.upf
H  1.00794   H_ONCV_PBE-1.0.upf
ATOMIC_POSITIONS bohr
Si      10.000000   10.000000  10.000000
H       11.614581   11.614581  11.614581
H        8.385418    8.385418  11.614581
H        8.385418   11.614581   8.385418
H       11.614581    8.385418   8.385418
K_POINTS gamma
EOF


cat > wbse_init.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

wbse_init_control:
  wbse_init_calculation: S
  solver: TDDFT
EOF


cat > wbse.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

wbse_init_control:
  wbse_init_calculation: S
  solver: TDDFT

wbse_control:
  wbse_calculation: L
  l_dipole_realspace: True
  n_lanczos: 200
EOF
