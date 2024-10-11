#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Mg_ONCV_PBE-1.2.upf
${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/O_ONCV_PBE-1.2.upf

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
ibrav   = 0
nat     = 15
ntyp    = 2
ecutwfc = 20.0
nbnd    = 80
/
&electrons
/
ATOMIC_SPECIES
Mg 24.3050  Mg_ONCV_PBE-1.2.upf
O  15.9994   O_ONCV_PBE-1.2.upf
ATOMIC_POSITIONS crystal
Mg       0.000000000   0.000000000   0.000000000
Mg      -0.002103461  -0.002103461   0.502103461
Mg      -0.002103461   0.502103461  -0.002103461
Mg      -0.002103461   0.502103461   0.502103461
Mg       0.502103461  -0.002103461  -0.002103461
Mg       0.502103461  -0.002103461   0.502103461
Mg       0.502103461   0.502103461  -0.002103461
Mg       0.500000000   0.500000000   0.500000000
O        0.250000000   0.250000000   0.750000000
O        0.250000000   0.750000000   0.250000000
O        0.250000000   0.750000000   0.750000000
O        0.750000000   0.250000000   0.250000000
O        0.750000000   0.250000000   0.750000000
O        0.750000000   0.750000000   0.250000000
O        0.750000000   0.750000000   0.750000000
K_POINTS {gamma}
CELL_PARAMETERS angstrom
0.000000 4.249236 4.249236
4.249236 0.000000 4.249236
4.249236 4.249236 0.000000
EOF


cat > westpp.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

westpp_control:
  westpp_calculation: LBD
  westpp_range: [61,80]
  westpp_box: [4.015, 12.045, 4.015, 12.045, 0.00, 8.03]
EOF
