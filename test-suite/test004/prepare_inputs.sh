#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Si_ONCV_PBE-1.1.upf

cat > pw.in << EOF
&control
calculation  = 'scf'
restart_mode = 'from_scratch'
pseudo_dir   = './'
outdir       = './'
prefix       = 'test'
wf_collect   = .TRUE.
/
&SYSTEM
ibrav       = 2,
a           = 5.43,
nat         = 2,
ntyp        = 1,
ecutwfc     = 30.0,
nbnd        = 10
noinv       = .true.
nosym       = .true.
/
&ELECTRONS
diago_full_acc = .true.
conv_thr    = 1.D-12
/
ATOMIC_SPECIES
 Si 28.085   Si_ONCV_PBE-1.1.upf
ATOMIC_POSITIONS (crystal)
 Si  0.0000  0.0000  0.0000
 Si  0.2500  0.2500  0.2500
K_POINTS (automatic)
1 1 2   0 0 0
EOF


cat > wstat.in << EOF
{
  "input_west": {
    "qe_prefix": "test",
    "west_prefix": "test",
    "outdir": "./"
  },
  "wstat_control": {
    "wstat_calculation": "S",
    "n_pdep_eigen": 10
  }
}
EOF


cat > wfreq.in << EOF 
{
  "input_west": {
    "qe_prefix": "test",
    "west_prefix": "test",
    "outdir": "./"
  },
  "wstat_control": {
    "wstat_calculation": "S",
    "n_pdep_eigen": 10
  },
  "wfreq_control": {
    "wfreq_calculation": "XWGQ",
    "n_pdep_eigen_to_use": 10,
    "qp_bandrange": [1,5],
    "n_refreq": 300,
    "ecut_refreq": 2.0, 
    "macropol_calculation" : "C"
  }
}
EOF
