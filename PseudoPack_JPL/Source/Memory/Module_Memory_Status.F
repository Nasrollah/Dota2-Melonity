#include "pseudopack.h"


!***************************************************************************
! MODULE      :: Memory_Allocation_Status
!                 
! PURPOSE     :: Memory Allocation/DeAllocation Error Check
!
! VARIABLE    :: M_Error  (Integer 
!                         Error flag from the allocation/deallocation call 
!                         with STAT=M_Error
!
!                (see below)
!
! SUBROUTINE  :: Memory_Error_Check
!
! ARGUMENTS   :: Case   (Integer
!
!                  0     Allocation error check
!                  1   DeAllocation error check
!
! USAGE       :: Standard
!
! MODULE USED :: PS_IO_Unit
!                Memory_Units
!
! AUTHOR      :: Wai Sun Don, PhD
! DATE        :: 3-14-1999
! VERSION     :: 2.0
!***************************************************************************


MODULE Memory_Allocation_Status

#if defined (DEBUG)
  USE Memory_Units

  USE PS_IO_Unit  , ONLY: lid6, lid99
  USE Processor   , ONLY: I_Am

implicit NONE

integer           :: Domain_Number    = -1
integer           :: M_Error          = 0
integer           :: Memory_Requested = 0
character(LEN=80) :: Subroutine_Name  = ' '
character(LEN=80) :: Subroutine_Call  = ' '

PRIVATE

PUBLIC  :: Domain_Number
PUBLIC  :: M_Error
PUBLIC  :: Memory_Requested
PUBLIC  :: Subroutine_Name
PUBLIC  :: Subroutine_Call

PUBLIC  :: Memory_Error_Check

CONTAINS

  Subroutine Memory_Error_Check (Case)

  integer           :: Case 
  character(LEN=12) :: Cases

  if (M_Error == 0) then
    Domain_Number    = -1
    Memory_Requested = 0
    Subroutine_Call  = ' '
    Subroutine_Name  = ' '

    RETURN
  endif

  if (Case == 0) Cases = 'ALLOCATION'
  if (Case == 1) Cases = 'DEALLOCATION'

  write (lid6 ,1000) TRIM(Cases), TRIM(Subroutine_Name),             &
                     Domain_Number, I_Am,                            &
                     CEILING(Memory_Requested*MUnits), MUnits_Name, M_Error
  write (lid99,1000) TRIM(Cases), TRIM(Subroutine_Name),             &
                     Domain_Number, I_Am,                            &
                     CEILING(Memory_Requested*MUnits), MUnits_Name, M_Error

  CLOSE (lid99)

  STOP

1000 format (1x/1x,72('*')/                           &
             1x,20x,'<Memory ',(a),' Error>'//        &
             1x,'   Subroutine Name  : ',(a)/         &
             1x,'   Domain Number    : ',i10/         &
             1x,'   Processor Number : ',i10/         &
             1x,'   Memory Requested : ',i10,1x,a6/   &
             1x,'   Error Number     : ',i10/         &
             1x,72('*')/)

  END Subroutine Memory_Error_Check
#else

integer           :: M_Error          = 0

PRIVATE

PUBLIC  :: M_Error
#endif

END MODULE Memory_Allocation_Status
