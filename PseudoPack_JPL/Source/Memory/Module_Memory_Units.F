#include "pseudopack.h"


!***************************************************************************
! MODULE      :: Memory_Units
!                 
! PURPOSE     :: Unit of the system memory
!
! VARIABLE    :: (see below)
!
! SUBROUTINE  :: NONE
!
! ARGUMENTS   :: NONE

! USAGE       :: N/A
! 
! MODULE USED :: NONE
!
! AUTHOR      :: Wai Sun Don, PhD
! DATE        :: 3-14-1999
! VERSION     :: 2.0
!***************************************************************************


MODULE Memory_Units

implicit NONE

integer , parameter :: KBytes = 1024
integer , parameter :: MBytes = 1024*KBytes

#if defined (CRAY_C90)
REALTYPE, parameter :: MUnits =   ONE/MBytes
#else
# if defined (SINGLE_PRECISION)
REALTYPE, parameter :: MUnits =  FOUR/MBytes
# else
REALTYPE, parameter :: MUnits = EIGHT/MBytes
# endif
#endif

#if defined (CRAY_C90)
character(LEN=6), parameter :: MUnits_Name = 'MWords'
#else
character(LEN=6), parameter :: MUnits_Name = 'MBytes'
#endif

END MODULE Memory_Units
