#include "pseudopack.h"
 

!***************************************************************************
! MODULE      :: Counters
!
! PURPOSE     :: Integer array of Counter 
!
! VARIABLE    :: Counter   (GLOBAL Integer Array of size (1:32)
!                          Initialized to 0
!
! SUBROUTINE  :: NONE
!
! ARGUMENTS   :: NONE
!
! AUTHOR      :: Wai Sun Don, PhD
! DATE        :: 3-14-1999
! VERSION     :: 2.0
!***************************************************************************


MODULE Counters

implicit NONE

integer, parameter :: Number_Of_Counter = 128

integer, dimension(Number_Of_Counter) :: PS_Counter = 0
logical, dimension(Number_Of_Counter) :: PS_Switch  = .FALSE.

END MODULE Counters
