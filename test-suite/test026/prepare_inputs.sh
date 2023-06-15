#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/O_ONCV_PBE-1.0.upf

cat > pw.in << EOF
&control
calculation  = 'scf'
restart_mode = 'from_scratch'
pseudo_dir   = './'
outdir       = './'
prefix       = 'test'
wf_collect   = .true.
/
&system
ibrav             = 1
celldm(1)         = 20
nat               = 2
ntyp              = 1
nspin             = 2
ecutwfc           = 25
nbnd              = 16
tot_magnetization = 2.
input_dft         = 'LDA'
/
&electrons
diago_full_acc = .true.
/
ATOMIC_SPECIES
O 16.00  O_ONCV_PBE-1.0.upf
ATOMIC_POSITIONS crystal
O        0.460000000   0.500000000   0.500000000
O        0.540000000   0.500000000   0.500000000
K_POINTS gamma
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
  wbse_calculation: D
  n_liouville_eigen: 4
  trev_liouville: 0.00000001
  trev_liouville_rel: 0.000001
  l_pre_shift: True
  l_spin_flip: True
  l_spin_flip_kernel: True
  l_spin_flip_alda0: False
  l_forces: True
  forces_state: 4
EOF
