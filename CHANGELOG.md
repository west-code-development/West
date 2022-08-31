Change Log
==========

v5.2.0 (2022/xx/xx)
-------------------

- Release of Quantum Defect Embedding Theory (QDET, N. Sheng et al., JCTC 18 3512 (2022)), enabling the use of WEST to study strongly correlated states of defects in solids.
- Improved the performance of the k-point case of `wfreq`.
- Bug fix. Fixed the k-point case of `wfreq` when `qp_bands` is used or when `qp_bandrange` does not start from 1.
- Updated CI/CD. Added tests to cover the new functionalities.
- Updated documentation. Added more tutorials.

v5.1.0 (2022/08/19)
-------------------

- Full release of GPU-accelerated full-frequency GW implementation (V. Yu and M. Govoni, JCTC 18 4690 (2022)), e.g., enabling the use of WEST at OLCF/Summit, NERSC/Perlmutter, ALCF/Polaris. GPU acceleration of `wstat` and `wfreq` is restricted to NVIDIA GPUs.
- Introduced infrastructure changes to prepare for QDET release. Enabled calculation with fractional occupation, added the possibility to specify a set of bands, instead of a range in `wfreq`, enabled the computation of the off-diagonal matrix elements in G0W0.
- Code updated for compatibility with Quantum ESPRESSO 7.1. QE must be compiled without CMake.
- Added the calculation of inverse participation ratio and localization factor in `westpp`.
- Bug fix. Fixed calculation of dipole matrix elements in `westpp`.
- Updated CI/CD to use the pytest framework.
- Updated documentation. Updated build instructions for ALCF/Theta, NERSC/Cori, and UChicago/RCC/Midway3. Added build instructions for NERSC/Perlmutter, OLCF/Summit, ANL/LCRC/Bebop, ANL/LCRC/Swing, and NVIDIA DGX A100. Added more tutorials.

v5.0.0 (2022/05/13)
-------------------

- Updated WEST to be compatible with QE 7.0 (before it was compatible with QE 6.1.0). QE must be compiled without CMake.
- Activated pool parallelization in `wstat` and `wfreq` to distribute spin polarization (for systems with nspin=2).
- Activated and optimized band group parallelization in `wstat` and `wfreq` (was disabled in QE 6.1.0).
- Enabled the calculation of dipole matrix elements in `westpp`.
- Updated build system. Automatically refresh the `make.depend` files when running `make conf`.
- Added interface to the ELPA eigensolver.
- Added support to read the HDF5 output of pwscf.
- Removed large direct access I/O, resolving issues with the NVIDIA nvfortran compiler.
- Removed obsolete or non-standard Fortran code, reducing the number of warnings during compilation.
- Updated library dependency to Json-Fortran 8.3.0.
- Bug fix. Print an error message when the code fails to read the restart files.
- Bug fix. Fixed the `XwgQ` restart mode of `wfreq`, i.e., computing Q from previously completed W and G.
- Bug fix. Avoid overflow in `IO_kernel/wfreqio.f90` in large-scale runs.
- Bug fix. Check that the mandatory logicals `nosym` and `noinv` are set to true for systems with k-points.
- To avoid overwriting JSON files when files with the same name already exist, a suffix is appended to the name of the new file.
- Updated CI/CD. Adapted nightly tests to cover OpenMP and ScaLAPACK.
- Updated CI/CD. Added tests of hybrid functionals.
- Updated documentation.

v4.3.0 (2021/05/26)
-------------------

- Introduced new data layout. Parallelization over bands allows to distribute data in a more flexible way in `wstat`. This feature also helps reduce memory per image. Band parallelization is enabled by specifying `-nb xxx` from the command line.
- Introduced checkpointing in `wstat`. With the new keyword `n_steps_write_restart` one can control how often the code produces restarts (default value is `1`).
- Improved I/O in `wfreq`. The number of I/O operations is reduced in `solve_wfreq`, and in the gamma case of `solve_gfreq`.
- Updated library dependency to Json-Fortran 8.2.1, resolving compilation issues with PGI 19.10.
- Updated the initialization of Forpy, such that in case of a module import error, the code has a better chance to print a clear error message.
- Updated build. Now the code builds with the NVIDIA/PGI Fortran compiler.
- Bug fix. Fixed an undefined variable in `Tools/set_npwq.f90`. Added `IMPLICIT NONE to all program units to let the compiler catch such errors in the future.
- Bug fix. Avoid overflow in `Wstat/wstat_memory_report.f90` in large-scale runs.
- Bug fix. Reset permissions to all source files in Westpp. All files appeared to be executable before.
- Bug fix. Reset `make.depend` files, added support to `make -j`.
- Updated CI/CD. Use updated Docker images.
- Updated CI/CD. Check numerical results and fail the CI if results don't match (see `check.py`). DFT checks error in total energy, WSTAT checks maximum error in PDEP eigenvalues, WFREQ checks maximum error in QP energies.
- Updated CI/CD. Added tests of images and OpenMP threads to nightly tests.
- Updated documentation. Fixed doc build with sphinx 3.5.0+.
- Updated documentation. Updated build instructions for ALCF/Theta, NERSC/Cori, UChicago/RCC/Midway3, macOS.
- Updated documentation. Added more tutorials.
- Updated documentation. Updated manual.

v4.2.1 (2020/10/19)
-------------------

- Added support for python 3.8 (--embed)
- Updated scripts for RCC-Midway and MacOSX
- Solved bugs in reporting conf layer in Makefile
- Updated documentation

v4.2.0 (2020/07/03)
-------------------

- Introduced automatic installation of missing python packages
- Introduced the conf layer in the Makefile to ease installation
- Updated manual

v4.1.0 (2019/10/18)
-------------------

- Improved usability of client/server mode with `server_control`
- Reduced execution time of wfreq (W) for solids
- Added build instructions for RCC-Midway and MacOSX
- Updated build instructions for ALCF-Theta
- Updated manual

v4.0.0 (2019/09/30)
-------------------

- Added client/server mode
- Added coupling to Qbox code (http://qboxcode.org)
- Added python3 interface
- Simplified the input format (now accepting both JSON and YAML formats)
- Expanded documentation

v3.1.1 (2018/09/19)
-------------------

- Python suite for pre- and post process WEST calculations

v3.1.0 (2018/06/30)
-------------------

- Introduction of k-points sampling
- Porting to Intel KNL
- Migration of all developments to a private GitLab server, and master branch mirrored to GitHub
- Test suite and continuous integration in place to automatically test the integrity of the code at every addition
- Expansion of documentation, and streamlining of its generation using markup language (Sphinx)

v3.0.0 (2017/07/16)
-------------------

- Restructuring of I/O in JSON (JavaScript Object Notation) format, thus enabling seamless integration with pre-/postprocessing tools and compatibility with Jupyter electronic notebooks

v2.0.0 (2016/10/19)
-------------------

- Implementation of spin-orbit coupling
- Implementation of novel hybrid functionals derived from GW self-energy

v1.1.0 (2016/05/20)
-------------------

- Addition of postprocessing routines, forming the seed for WESTpy

v1.0.3 (2016/01/08)
-------------------

- Efficiency improvements and bug fixes to the contour deformation technique

v1.0.2 (2015/09/23)
-------------------

- Efficiency improvements and bug fixes to the PDEP algorithm

v1.0.1 (2015/06/20)
-------------------

- Initial beta release of WEST: GW without empty states
