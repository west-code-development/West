#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Mg_ONCV_PBE-1.2.upf
${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/O_ONCV_PBE-1.2.upf

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
{                                                                               
  "input_west": {},                                                             
  "wstat_control": {                                                            
    "wstat_calculation": "S",                                                   
    "l_minimize_exx_if_active": true,                                           
    "n_pdep_eigen": 61,                                                         
    "trev_pdep": 0.01                                                           
  }                                                                             
}      
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
  wfreq_calculation: XWGQH
  l_enable_off_diagonal: true
  n_pdep_eigen_to_use: 50
  qp_bands: [8,9,10]
  n_refreq: 300
  ecut_refreq: 2.0
EOF
