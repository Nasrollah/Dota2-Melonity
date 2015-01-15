#include "pseudopack.h"


MODULE Adjust_D_F_Epsilon

IMPLICIT NONE

PRIVATE
PUBLIC :: PS_Adjust_D_F_Epsilon

CONTAINS
!
! =======================================================================
!  
  Subroutine PS_Adjust_D_F_Epsilon (Adjust_D_Flux, Navier_Stokes, Epsilon, &
                                    D_F_Epsilon, dt, dx, dy, dz)

  logical            :: Navier_Stokes
  integer            :: Adjust_D_Flux, Power
  REALTYPE           :: Epsilon, D_F_Epsilon
  REALTYPE           :: dt, dx, d_min
  REALTYPE, OPTIONAL :: dy, dz

  Power = 1 ; if (Navier_Stokes) Power = 2

  SELECT CASE (Adjust_D_Flux)
    CASE (0)
      D_F_Epsilon = ZERO

    CASE (1)
      D_F_Epsilon = Epsilon

    CASE (2)
      D_F_Epsilon = Epsilon/dt

    CASE (3)
                       d_min =     dx
      if (PRESENT(dy)) d_min = MIN(dy,d_min)
      if (PRESENT(dz)) d_min = MIN(dz,d_min)

      D_F_Epsilon = Epsilon/(dt*d_min**Power)

  END SELECT

  END Subroutine PS_Adjust_D_F_Epsilon 

END MODULE Adjust_D_F_Epsilon
