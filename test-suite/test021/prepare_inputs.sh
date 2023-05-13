#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/O_ONCV_PBE-1.0.upf

cat > pw.in << EOF
&control
    calculation='scf',
    prefix='test'
    pseudo_dir='./'
    outdir = './',
    verbosity='high',
    wf_collect=.TRUE.
/
&system
    ibrav=0,
    nat=2,
    ntyp=1,
    ecutwfc=25,
    nbnd=16
    nspin=2
    tot_magnetization=2
    input_dft='LDA'
/
&electrons
    diago_full_acc=.true.,
    conv_thr=1.0D-8
/
ATOMIC_SPECIES
O 16.00  O_ONCV_PBE-1.0.upf
K_POINTS gamma
CELL_PARAMETERS {angstrom}
15 0 0
0 15 0
0 0 15
ATOMIC_POSITIONS {crystal}
O        0.480000000   0.500000000   0.500000000
O        0.560000000   0.500000000   0.500000000
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
  scissor_ope: 0.0
  n_liouville_eigen: 10
  n_liouville_times: 4
  n_liouville_maxiter: 1000
  n_liouville_read_from_file: 0
  trev_liouville: 0.00000001
  trev_liouville_rel: 0.000001
  wbse_epsinfty: 1.0
  spin_excitation: S
  l_preconditioning: True
  l_pre_shift: True
  l_spin_flip: True
  l_spin_flip_kernel: True
  l_spin_flip_alda0: False
  l_print_spin_flip_kernel: False
  spin_flip_cut1: 1000
  l_reduce_io: True
EOF
