#include "pseudopack.h"
 

!***************************************************************************
! MODULE      :: PS_IO_Unit
!
! PURPOSE     :: Assign logical IO unit for files
!
! VARIABLE    :: (See below)
!
! SUBROUTINE  :: NONE
!
! ARGUMENTS   :: NONE
!
! AUTHOR      :: Wai Sun Don, PhD
! DATE        :: 3-14-1999
! VERSION     :: 2.0
!***************************************************************************


MODULE PS_IO_Unit

implicit NONE

logical :: IO_Error = .FALSE. ! Flag indicating Input/Ouput Errors

  ! File Unit Number used for 

integer :: lid5  = 5     ! Standard Screen Input 
integer :: lid6  = 6     ! Standard Screen Output

integer :: lid98  = 98   ! Library call Internal Use 
integer :: lid99  = 99   ! Library call Error
integer :: lid100 = 100  ! Library call Internal Use 
integer :: lid101 = 101  ! Library call Internal Use 

END MODULE PS_IO_Unit
