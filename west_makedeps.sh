#!/bin/sh
# compute dependencies for West

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
TOPDIR=`pwd`

if test $# = 0
then
# this is the list of all directories for which we want to find dependencies
# upon include files *.h or *.fh or modules. Note that libraries that are
# externally maintained should not go into this list

    dirs=" West/Coulomb_kernel West/DFPT_kernel West/FFT_kernel \
           West/Hamiltonian_kernel West/IO_kernel West/Modules \
           West/Para_kernel West/Tools West/Westpp West/Wfreq West/Wstat "

elif
    test $1 = "-addson"
then
    echo "The script for adding new dependencies is running"
    echo "Usage: $0 -addson DIR DEPENDENCY_DIRS"
    echo "$0 assumes that the new dependencies are in $TOPDIR/../"
    dirs=$2
    shift
    shift
    add_deps=$*
    echo "dependencies in $add_deps will be searched for $dirs"
else
    dirs=$*
fi


for dir in $dirs; do

    # the following command removes a trailing slash
    DIR=`echo ${dir%/}`

    # the following would also work
    #DIR=`echo $dir | sed "s,/$,,"`

    # set inter-directory dependencies - only directories containing
    # modules that are used, or files that are included, by routines
    # in directory DIR should be listed in DEPENDS
    # (directory DIR itself should not be listed in DEPENDS)
    LEVEL1=..
    LEVEL2=../..
    # default
    DEPENDS="$LEVEL1/include"
    # for convenience, used later
    DEPEND1="$LEVEL1/include $LEVEL1/iotk/src $LEVEL1/FFTXlib $LEVEL1/LAXlib"
    DEPEND2="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/FFTXlib $LEVEL2/LAXlib \
             $LEVEL2/Modules"
    case $DIR in
        West/Coulomb_kernel )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/DFPT_kernel $LEVEL1/FFT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/IO_kernel $LEVEL1/Modules $LEVEL1/Para_kernel $LEVEL1/Tools $LEVEL1/Westpp $LEVEL1/Wfreq $LEVEL1/Wstat" ;;
        West/DFPT_kernel )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/FFT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/IO_kernel $LEVEL1/Modules $LEVEL1/Para_kernel $LEVEL1/Tools $LEVEL1/Westpp $LEVEL1/Wfreq $LEVEL1/Wstat" ;;
        West/FFT_kernel )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/DFPT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/IO_kernel $LEVEL1/Modules $LEVEL1/Para_kernel $LEVEL1/Tools $LEVEL1/Westpp $LEVEL1/Wfreq $LEVEL1/Wstat" ;;
        West/Hamiltonian_kernel )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/DFPT_kernel $LEVEL1/FFT_kernel $LEVEL1/IO_kernel $LEVEL1/Modules $LEVEL1/Para_kernel $LEVEL1/Tools $LEVEL1/Westpp $LEVEL1/Wfreq $LEVEL1/Wstat" ;;
        West/IO_kernel )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/DFPT_kernel $LEVEL1/FFT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/Modules $LEVEL1/Para_kernel $LEVEL1/Tools $LEVEL1/Westpp $LEVEL1/Wfreq $LEVEL1/Wstat" ;;
        West/Modules )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/DFPT_kernel $LEVEL1/FFT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/IO_kernel $LEVEL1/Para_kernel $LEVEL1/Tools $LEVEL1/Westpp $LEVEL1/Wfreq $LEVEL1/Wstat" ;;
        West/Para_kernel )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/DFPT_kernel $LEVEL1/FFT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/IO_kernel $LEVEL1/Modules $LEVEL1/Tools $LEVEL1/Westpp $LEVEL1/Wfreq $LEVEL1/Wstat" ;;
        West/Tools )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/DFPT_kernel $LEVEL1/FFT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/IO_kernel $LEVEL1/Modules $LEVEL1/Para_kernel $LEVEL1/Westpp $LEVEL1/Wfreq $LEVEL1/Wstat" ;;
        West/Westpp )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/DFPT_kernel $LEVEL1/FFT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/IO_kernel $LEVEL1/Modules $LEVEL1/Para_kernel $LEVEL1/Tools $LEVEL1/Wfreq $LEVEL1/Wstat" ;;
        West/Wfreq )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/DFPT_kernel $LEVEL1/FFT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/IO_kernel $LEVEL1/Modules $LEVEL1/Para_kernel $LEVEL1/Tools $LEVEL1/Westpp $LEVEL1/Wstat" ;;
        West/Wstat )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL1/Libraries/Base64 $LEVEL1/Libraries/Forpy $LEVEL1/Libraries/Json $LEVEL1/Libraries/Json_test $LEVEL1/Coulomb_kernel $LEVEL1/DFPT_kernel $LEVEL1/FFT_kernel $LEVEL1/Hamiltonian_kernel $LEVEL1/IO_kernel $LEVEL1/Modules $LEVEL1/Para_kernel $LEVEL1/Tools $LEVEL1/Westpp $LEVEL1/Wfreq" ;;
    *)
# if addson needs a make.depend file
        DEPENDS="$DEPENDS $add_deps"

    esac

    # list of all system modules
    sysdeps="iso_c_binding iso_fortran_env f90_unix_io f90_unix_env \
             f90_unix_proc ifcore ifport"

    # list of all external library modules or include files
    libdeps="mpi omp_lib hdf5 mkl_dfti mkl_dfti.f90 fftw3.f03 fftw3.f \
             xc_version.h xc_f03_lib_m elpa elpa1 \
             mbd w90_io fox_dom fox_wxml m_common_io \
             device_fbuff_m device_memcpy_m device_auxfunc_m"

    # list of all cuda-related modules
    cudadeps="cublas cublas_v2 cudafor curand cufft flops_tracker cusolverdn \
              zhegvdx_gpu dsyevd_gpu dsygvdx_gpu eigsolve_vars \
              nvtx_inters"

    # generate dependencies file (only for directories that are present)

    if test -d $TOPDIR/../$DIR
    then
        cd $TOPDIR/../$DIR

        $TOPDIR/../install/moduledep.sh $DEPENDS > make.depend
        $TOPDIR/../install/includedep.sh $DEPENDS >> make.depend

        # remove unwanted dependency upon system and library modules
        for no_dep in $sysdeps $libdeps $cudadeps; do
            echo "/@$no_dep@/d" >> removedeps.tmp
        done
        sed -f removedeps.tmp make.depend > tmp; mv tmp make.depend
        /bin/rm removedeps.tmp

        # check for missing dependencies
        if grep @ make.depend
        then
            notfound=1
            echo WARNING: dependencies not found in directory $DIR
        else
           echo directory $DIR : ok
        fi
    else
        echo directory $DIR : not present in $TOPDIR
    fi
done
if test "$notfound" = ""
then
    echo all dependencies updated successfully
fi
