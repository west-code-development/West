PROGRAM test

USE json_module, ONLY : json_file

INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

TYPE(json_file) :: json
INTEGER :: i
INTEGER :: iiarg, nargs
INTEGER :: numsp
LOGICAL :: found
CHARACTER(LEN=512) :: input_file
CHARACTER(LEN=:),ALLOCATABLE :: cval
REAL(DP) :: rval 
INTEGER :: ival 
INTEGER,ALLOCATABLE :: ivec(:)
REAL(DP),ALLOCATABLE :: rvec(:)
LOGICAL :: lval 
!
INTEGER :: iunit
INTEGER :: i0, i0_
!
CHARACTER(LEN=512),PARAMETER :: ifile="input_file.json"
CHARACTER(LEN=512),PARAMETER :: ofile="output_file.json"
!
! INIT INPUT FILE
!
OPEN( NEWUNIT=iunit, FILE=TRIM(ifile) )
WRITE( iunit, * ) '{ "test" : { "i0" : 100} }' 
CLOSE( iunit )
!
! INIT OUTPUT FILE
!
OPEN( NEWUNIT=iunit, FILE=TRIM(ofile) )
WRITE( iunit, * ) '{}' 
CLOSE( iunit )
!
! READ
!
CALL json%initialize()
CALL json%load_file(filename = TRIM(ifile))
CALL json%get('test.i0', i0, found)
i0_ = 0
IF( found ) i0_ = i0
CALL json%destroy()
!
! WRITE
!
CALL json%initialize()
CALL json%load_file(filename = TRIM(ofile))
CALL json%add('test.i0', i0_)
!
OPEN( NEWUNIT=iunit, FILE=TRIM(ofile) )
CALL json%print_file( iunit )
CLOSE( iunit )
CALL json%destroy()

END PROGRAM
