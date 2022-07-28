#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/H_ONCV_PBE-1.0.upf
${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Si_ONCV_PBE-1.1.upf

cat > pw.in << EOF
&CONTROL
  calculation = 'scf',
  pseudo_dir = './'
  wf_collect = .TRUE.
/
&SYSTEM
  ibrav = 0,
  nat = 15,
  ntyp = 2,
  input_dft = 'PBE'
  ecutwfc = 80
  nspin = 1
  nbnd = 80
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Mg  24.3050 Mg_ONCV_PBE-1.2.upf
  O  15.9994 O_ONCV_PBE-1.2.upf
ATOMIC_POSITIONS crystal
Mg      -0.000000000   0.000000000  -0.000000000
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
K_POINTS gamma
CELL_PARAMETERS angstrom
  -0.000000 4.249236 4.249236
  4.249236 0.000000 4.249236
  4.249236 4.249236 0.000000
EOF


cat > wstat.in << EOF
wstat_control:
  wstat_calculation: S
  n_pdep_eigen: 61
  trev_pdep: 0.01
EOF


cat > wfreq.in << EOF
wstat_control: 
  wstat_calculation: S
  l_minimize_exx_if_active: true
  n_pdep_eigen: 61
  trev_pdep: 0.01

wfreq_control:
  qp_bands: [61, 63, 64, 65]
  wfreq_calculation: XWGQH
  n_pdep_eigen_to_use: 61
  l_enable_lanczos: true
  l_enable_off_diagonal: true
  macropol_calculation: C
EOF
