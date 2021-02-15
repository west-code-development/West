#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/H_ONCV_PBE-1.0.upf
${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/C_ONCV_PBE-1.0.upf

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
nat               = 5
ntyp              = 2
nspin             = 2
ecutwfc           = 25.0
nbnd              = 10
assume_isolated   = 'mp'
tot_magnetization = 0.
/
&electrons
diago_full_acc = .true.
/
ATOMIC_SPECIES
C  12.0107   C_ONCV_PBE-1.0.upf
H  1.00794   H_ONCV_PBE-1.0.upf
ATOMIC_POSITIONS angstrom
C  0.0000  0.0000  0.0000
H  0.6276 -0.6275  0.6276
H -0.6276  0.6276  0.6276
H -0.6276 -0.6276 -0.6276
H  0.6276  0.6276 -0.6276
K_POINTS {gamma}
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
  n_pdep_eigen_to_use: 50
  qp_bandrange: [1,5]
  n_refreq: 300
  ecut_refreq: 2.0
EOF
