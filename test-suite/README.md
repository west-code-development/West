# WEST test-suite

## List of tests

1. SiH4 molecule, GW, Gamma only
2. CH4 molecule, GW, Gamma only, JSON input
3. CH4 molecule, spin-polarized GW, Gamma only
4. Si bulk, GW, k-mesh 1x1x2
5. SiH4 molecule, westpp output PDEP, Gamma only
6. SiH4 molecule, GW hybrid ACE, Gamma only
7. Si bulk, GW hybrid ACE, k-mesh 1x1x2
8. MgO, westpp localization factor, Gamma only
9. O2 molecule, GW fractional occupation, Gamma only
10. SiH4 molecule, GW `l_off_diagonal`, Gamma only
11. SiH4 molecule, GW `qp_bands`, Gamma only
12. SiH4 molecule, QDET, Gamma only
13. SiH4 molecule, QDET verbosity, Gamma only
14. NV- diamond, spin-polarized QDET, Gamma only
15. SiH4 molecule, westpp Wannier localization, Gamma only
16. SiH4 molecule, BSE Lanczos, Gamma only
17. SiH4 molecule, BSE Davidson, Gamma only
18. SiH4 molecule, TDDFT (PBE0) Lanczos, Gamma only
19. SiH4 molecule, GW hybrid no ACE, Gamma only
20. Formaldehyde molecule, TDDFT (PBE) forces, Gamma only
21. NV- diamond spin-polarized TDDFT (PBE) forces, Gamma only
22. Formaldehyde molecule, TDDFT (PBE0) forces, Gamma only
23. NV- diamond spin-polarized TDDFT (DDH) forces, Gamma only
24. O2 molecule, spin-flip TDDFT (LDA) forces, Gamma only
25. O2 molecule, spin-flip TDDFT (PBE) forces, Gamma only
26. O2 molecule, spin-flip TDDFT (PBE0) forces, Gamma only
27. SiH4 molecule, QDET `n_pdep_eigen_off_diagonal`, Gamma only

## Executing tests

In order to manually run the test, simply run `make` or `make all` in the test-suite.

## Adding new tests

Additional tests can be added with the following procedure:

1) Create a new test directory with:
  - `prepare_inputs.sh`, to generate all inputs (including links to download pseudopotential files)
  - `ref` subdirectory, that contains reference data
  - `Makefile`, that contains instructions about the execution of the codes

We recommend that you copy and adjust an existing test script.

2) Add the new subdirectory name to the list in the `test-suite/Makefile`.

3) Customize the file `main_test.py`.

4) Add a short description of the new test to the list above.
