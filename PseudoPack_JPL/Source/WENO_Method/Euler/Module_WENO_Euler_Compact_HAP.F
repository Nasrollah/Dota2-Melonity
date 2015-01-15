#include "pseudopack.h"


Module WENO_Euler_Compact_HAP

  USE WENO_Coefficients
  USE WENO_Option

implicit NONE

INTERFACE PS_WENO_Euler_Compact_HAP
  MODULE PROCEDURE PS_WENO_Euler_Compact_HAP_REAL
  MODULE PROCEDURE PS_WENO_Euler_Compact_HAP_INT
END INTERFACE

PRIVATE
PUBLIC  :: PS_WENO_Euler_Compact_HAP

CONTAINS

  Subroutine PS_WENO_Euler_Compact_HAP_REAL &
               (Order, N0, N5, i0, i1, F, F_Half, Power, WENO_Epsilon)

#include "WENO_Euler_Declaration.i"

  REALTYPE                   :: Power

  integer                    :: j, k, r

                             Epsilon_WENO = PS_WENO_Option%Epsilon
  if (PRESENT(WENO_Epsilon)) Epsilon_WENO = WENO_Epsilon

#include "WENO_Euler_Compact_HAP.i"

  END Subroutine PS_WENO_Euler_Compact_HAP_REAL 
!
! -----------------------------------------------------------------------------
!
  Subroutine PS_WENO_Euler_Compact_HAP_INT  &
               (Order, N0, N5, i0, i1, F, F_Half, Power, WENO_Epsilon)

#include "WENO_Euler_Declaration.i"

  integer                    :: Power

  integer                    :: j, k, r

                             Epsilon_WENO = PS_WENO_Option%Epsilon
  if (PRESENT(WENO_Epsilon)) Epsilon_WENO = WENO_Epsilon

#include "WENO_Euler_Compact_HAP.i"

  END Subroutine PS_WENO_Euler_Compact_HAP_INT  

END Module WENO_Euler_Compact_HAP
