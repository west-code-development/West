Change Log
==========

v5.0.0 (2021/05/XX)
-------------------

- New variables
   - Restart info was written out in every Davidson iteration. Now this can be changed by setting the `n_steps_write_restart` keyword. Default is `1`, i.e., same behavior as before.

- Bug fixes 
   - Update the initialization for Forpy, such that in case of a module import error, the code has a better chance to print a clear error message.
   - Fix an undefined variable in `Tools/set_npwq.f90`. Add `IMPLICIT NONE to all program units to let the compiler catch such errors in the future.
   - Avoid overflow in `Wstat/wstat_memory_report.f90` in large-scale runs.
   - Fix non-standard Fortran codes that do not work with NVIDIA/PGI Fortran compiler.
   - Fix file mode for source files in Westpp. All files were executable before.
   - Work around an input parsing issue encountered with pgfortran and nvfortran.

- CI
   - Use updated Docker images
   - Check numerical results and fail the CI if results don't match. (see check.py). DFT checks error in total energy, WSTAT checks maximum error in PDEP eigenvalues, WFREQ checks maximum error in QP energies
   - Add tests of images and OpenMP threads to nightly test

- Documentation
   - Fix doc build with sphinx 3.5.0+. (see also https://github.com/sphinx-doc/sphinx/issues/8885)
   - Update build instructions for ALCF/Theta, NERSC/Cori, UChicago/RCC/Midway3, macOS 
   - Updated documentation

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
