#include "pseudopack.h"


Module WENO_Scalar_Short_Zico

  USE WENO_Coefficients
  USE WENO_Option

implicit NONE

INTERFACE PS_WENO_Scalar_Short_Zico
  MODULE PROCEDURE PS_WENO_Scalar_Short_Zico_REAL
  MODULE PROCEDURE PS_WENO_Scalar_Short_Zico_INT
END INTERFACE

PRIVATE
PUBLIC  :: PS_WENO_Scalar_Short_Zico

CONTAINS

  Subroutine PS_WENO_Scalar_Short_Zico_REAL &
               (Order, N0, N5, i0, i1, F, F_Half, Power, WENO_Epsilon)

#include "WENO_Scalar_Declaration.i"

  REALTYPE                   :: Power

                             Epsilon_WENO = PS_WENO_Option%Epsilon
  if (PRESENT(WENO_Epsilon)) Epsilon_WENO = WENO_Epsilon

#include "WENO_Scalar_Short_Zico.i"

  END Subroutine PS_WENO_Scalar_Short_Zico_REAL
!
! -----------------------------------------------------------------------------
!
  Subroutine PS_WENO_Scalar_Short_Zico_INT  &
               (Order, N0, N5, i0, i1, F, F_Half, Power, WENO_Epsilon)

#include "WENO_Scalar_Declaration.i"

  integer                    :: Power

                             Epsilon_WENO = PS_WENO_Option%Epsilon
  if (PRESENT(WENO_Epsilon)) Epsilon_WENO = WENO_Epsilon

#include "WENO_Scalar_Short_Zico.i"

  END Subroutine PS_WENO_Scalar_Short_Zico_INT

END Module WENO_Scalar_Short_Zico

