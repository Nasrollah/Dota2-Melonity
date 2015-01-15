#include "pseudopack.h"

!----------------------------------------------------------------------------
!         f_Map_x    is the mapping function         (x, x0, x1, alpha, beta)
! Inverse_f_Map_x    is the inverse mapping function (y, x0, x1, alpha, beta)
!
!       D_f_Map_x_1  is the first  derivative        (x, x0, x1, alpha, beta)
!       D_f_Map_x_2  is the second derivative        (x, x0, x1, alpha, beta)
!       D_f_Map_x_3  is the third  derivative        (x, x0, x1, alpha, beta)
!       D_f_Map_x_4  is the fourth derivative        (x, x0, x1, alpha, beta)
!
! where  x is the computational grid point
!        y is the physical grid point
!       x0 is the Leftmost  boundary point of the domain
!       x1 is the Rightmost boundary point of the domain
!
! To define a mapping, only f_Map_x, Inverse_f_Map_x and as many derivatives of
!  f_Map_x as needed are required.
!
! The derivatives of f_Map_x, D_f_Map_n are based on the Computational Grid x,
!  that is, Chebyshev, Legendre and Fourier nodes before any mapping applied.
!----------------------------------------------------------------------------

MODULE User_Defined_Map_Function

implicit NONE

PRIVATE

PUBLIC :: f_Map_x, Inverse_f_Map_x, D_f_Map_x

CONTAINS

  REALTYPE FUNCTION f_Map_x (x, x0, x1, alpha, beta) 

  REALTYPE, INTENT(IN) :: x, x0, x1, alpha, beta

  f_Map_x = ASIN(alpha*x)/ASIN(alpha)

  END FUNCTION f_Map_x 

!******************************************************
  REALTYPE FUNCTION Inverse_f_Map_x (y, x0, x1, alpha, beta) 

  REALTYPE, INTENT(IN) :: y, x0, x1, alpha, beta

  Inverse_f_Map_x = SIN(ASIN(alpha)*y)/alpha

  END FUNCTION Inverse_f_Map_x

!******************************************************
  Subroutine D_f_Map_x (x, x0, x1, alpha, beta,                             &
                        D_f_Map_x_1, D_f_Map_x_2, D_f_Map_x_3, D_f_Map_x_4, &
                        Map) 

  integer , INTENT(INOUT) :: Map
  REALTYPE, INTENT(IN)    :: x, x0, x1, alpha, beta
  REALTYPE, INTENT(OUT)   :: D_f_Map_x_1, D_f_Map_x_2, D_f_Map_x_3, D_f_Map_x_4 

  REALTYPE                :: kappa, gamma

  Map = -5

  kappa = ASIN(alpha) ; gamma = -kappa*kappa

  D_f_Map_x_1 = (kappa/alpha)*SQRT(ONE-(alpha*x)**2)
  D_f_Map_x_2 = gamma*x
  D_f_Map_x_3 = gamma*D_f_Map_x_1 
  D_f_Map_x_4 = gamma*D_f_Map_x_2 

  END Subroutine D_f_Map_x 

END MODULE User_Defined_Map_Function
