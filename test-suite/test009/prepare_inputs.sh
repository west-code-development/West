#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/O_ONCV_PBE-1.0.upf

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
nat             = 2
ntyp            = 1
ecutwfc         = 25.0
nbnd            = 12
assume_isolated = 'mp'
occupations     = 'from_input'
/
&electrons
diago_full_acc = .true.
/
ATOMIC_SPECIES
O 15.9994 O_ONCV_PBE-1.0.upf
ATOMIC_POSITIONS angstrom
O    6.25001033  6.25001033  5.64063432
O    6.25001033  6.25001033  6.84863432
K_POINTS gamma
OCCUPATIONS
2.00  2.00  2.00  2.00  2.00  1.00  1.00  0.00  0.00  0.00
0.00  0.00
EOF


cat > wstat.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

wstat_control:
  wstat_calculation: S
  n_pdep_eigen: 50
EOF


cat > wfreq.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

wstat_control:
  wstat_calculation: S
  n_pdep_eigen: 50

wfreq_control:
  wfreq_calculation: XWGQ
  macropol_calculation: N
  n_pdep_eigen_to_use: 50
  qp_bandrange: [1,7]
  n_refreq: 300
  ecut_refreq: 2.0
EOF
