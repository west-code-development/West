#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/C_ONCV_PBE-1.0.upf
${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/N_ONCV_PBE-1.0.upf

cat > pw.in << EOF
&control
calculation  = 'scf'
restart_mode = 'from_scratch'
pseudo_dir   = './'
outdir       = './'
prefix       = 'test'
/
&system
ibrav             = 0
nat               = 15
ntyp              = 2
ecutwfc           = 25
nosym             = .true.
tot_charge        = -1
nspin             = 2
nbnd              = 40
tot_magnetization = 2
input_dft         = 'pbe0'
exx_fraction      = 0.18
/
&electrons
diago_full_acc = .true.
/
ATOMIC_SPECIES
C  12.0107  C_ONCV_PBE-1.0.upf
N  14.0067  N_ONCV_PBE-1.0.upf
ATOMIC_POSITIONS crystal
C            -0.0011453699       -0.0011377611        0.0067006607
C            -0.0002744852       -0.0002672037        0.5008176648
C            -0.0011044829        0.4955610337        0.0066821802
C             0.0450632327        0.4851167525        0.4846789951
C             0.4955573223       -0.0010991005        0.0066828147
C             0.4851269633        0.0450567053        0.4846639674
N             0.4994339805        0.4994639176        0.0018025833
C             0.4850798144        0.4850821394        0.4846967685
C             0.1230453760        0.1230396916        0.1309068249
C             0.1272210852        0.1272104167        0.6248827206
C             0.1356116446        0.6167545386        0.1309034755
C             0.1272035150        0.6206680852        0.6249000996
C             0.6167608167        0.1356077552        0.1308943700
C             0.6206777059        0.1271982848        0.6248928293
C             0.6167428815        0.6167447445        0.1308940454
K_POINTS gamma
CELL_PARAMETERS angstrom
-0.000000 3.566790 3.566790
 3.566790 0.000000 3.566790
 3.566790 3.566790 0.000000
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
  wbse_calculation: D
  n_liouville_eigen: 2
  n_liouville_times: 20
  trev_liouville: 0.00000001
  trev_liouville_rel: 0.000001
  l_pre_shift: False
  l_forces: True
  forces_state: 1
EOF
