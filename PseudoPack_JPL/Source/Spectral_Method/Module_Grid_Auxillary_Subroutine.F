#include "pseudopack.h"


MODULE Grid_Auxillary_Subroutine

implicit NONE

INTERFACE Map_Computational_2_Physical
  MODULE PROCEDURE Map_Computational_2_Physical_S
  MODULE PROCEDURE Map_Computational_2_Physical_V
END INTERFACE

INTERFACE Map_Physical_2_Computational
  MODULE PROCEDURE Map_Physical_2_Computational_S
  MODULE PROCEDURE Map_Physical_2_Computational_V
END INTERFACE

PRIVATE

PUBLIC  :: Map_Computational_2_Physical
PUBLIC  :: Map_Physical_2_Computational
PUBLIC  :: Setup_Metric_A
PUBLIC  :: Apply_Linear_Map

CONTAINS

! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Computational_2_Physical
! DESCRIPTION  : Map the Chebyshev grid points to the physical domain.
!        Mapped [c_Left,c_Right] -> [x_Left,x_Right]
! ----------------------------------------------------------------------
  Subroutine Map_Computational_2_Physical_V (N, c_Left,c_Right, c, &
                                                x_Left,x_Right, x)

  integer  :: N
  REALTYPE :: c_Left, c_Right, x_Left, x_Right, Slope

  REALTYPE, dimension(N) :: x, c

  Slope = (x_Right-x_Left)/(c_Right-c_Left)

      x = x_Right+Slope*(c-c_Right)

  END Subroutine Map_Computational_2_Physical_V

! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Physical_2_Computational
! DESCRIPTION  : Map the Chebyshev grid points to the physical domain.
!        Mapped [x_Left,x_Right] -> [c_Left,c_Right]
! ----------------------------------------------------------------------
  Subroutine Map_Physical_2_Computational_V (N, x_Left,x_Right, x, &
                                                c_Left,c_Right, c)

  integer  :: N
  REALTYPE :: c_Left, c_Right, x_Left, x_Right, Slope

  REALTYPE, dimension(N) :: x, c

  Slope  = (c_Right-c_Left)/(x_Right-x_Left)

  c = c_Right+Slope*(x-x_Right)

  END Subroutine Map_Physical_2_Computational_V

! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Computational_2_Physical
! DESCRIPTION  : Map the Chebyshev grid points to the physical domain.
!        Mapped [c_Left,c_Right] -> [x_Left,x_Right]
! ----------------------------------------------------------------------
  Subroutine Map_Computational_2_Physical_S (N, c_Left,c_Right, c, &
                                                x_Left,x_Right, x)

  integer  :: N
  REALTYPE :: c_Left, c_Right, x_Left, x_Right, Slope

  REALTYPE :: x, c

  Slope = (x_Right-x_Left)/(c_Right-c_Left)

      x = x_Right+Slope*(c-c_Right)

  END Subroutine Map_Computational_2_Physical_S

! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Physical_2_Computational
! DESCRIPTION  : Map the Chebyshev grid points to the physical domain.
!        Mapped [x_Left,x_Right] -> [c_Left,c_Right]
! ----------------------------------------------------------------------
  Subroutine Map_Physical_2_Computational_S (N, x_Left,x_Right, x, &
                                                c_Left,c_Right, c)

  integer  :: N
  REALTYPE :: c_Left, c_Right, x_Left, x_Right, Slope

  REALTYPE :: x, c

  Slope  = (c_Right-c_Left)/(x_Right-x_Left)

  c = c_Right+Slope*(x-x_Right)

  END Subroutine Map_Physical_2_Computational_S
! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Metric_A
! DESCRIPTION  : Setup the Transformation Metric A for differentiation
!                up to order 4.
!
! D_Map_D_x_k is an array of size N0+1 which stored the k_th derivative
!                of the mapping function  xi=xi(x)
! ----------------------------------------------------------------------
  Subroutine Setup_Metric_A (N0, Number_Of_Derivative,       &
                                 D_Map_D_x_1, D_Map_D_x_2,   &
                                 D_Map_D_x_3, D_Map_D_x_4, A)

  integer  :: N0, Number_Of_Derivative

  REALTYPE, dimension(0:N0,*) :: A
  REALTYPE, dimension(0:N0)   :: D_Map_D_x_1, D_Map_D_x_2,  &
                                 D_Map_D_x_3, D_Map_D_x_4

    A(:, 1) =       D_Map_D_x_1

  if (Number_Of_Derivative >= 2) then
    A(:, 2) =       D_Map_D_x_2
    A(:, 3) =       D_Map_D_x_1**2
  endif

  if (Number_Of_Derivative >= 3) then
    A(:, 4) =       D_Map_D_x_3
    A(:, 5) = THREE*D_Map_D_x_2*D_Map_D_x_1
    A(:, 6) =       D_Map_D_x_1**3
  endif

  if (Number_Of_Derivative >= 4) then
    A(:, 7) =       D_Map_D_x_4
    A(:, 8) =  FOUR*D_Map_D_x_3*D_Map_D_x_1+THREE*D_Map_D_x_2**2
    A(:, 9) =   SIX*D_Map_D_x_2*D_Map_D_x_1**2
    A(:,10) =       D_Map_D_x_1**4
  endif

  END Subroutine Setup_Metric_A 

! ----------------------------------------------------------------------
! FUNCTION NAME: Apply_Linear_Map
! DESCRIPTION  : Apply transformation metric of linear map Slope to Metric A
!
! ----------------------------------------------------------------------
  Subroutine Apply_Linear_Map (N0, Number_Of_Derivative, Slope, A)

  integer  :: N0, Number_Of_Derivative
  REALTYPE :: Slope, S1, S2, S3, S4

  REALTYPE, dimension(0:N0,*) :: A

    S1 = Slope

    A(:, 1) = A(:, 1)*S1

  if (Number_Of_Derivative >= 2) then
    S2 = S1*S1

    A(:, 2) = A(:, 2)*S2
    A(:, 3) = A(:, 3)*S2
  endif

  if (Number_Of_Derivative >= 3) then
    S3 = S2*S1

    A(:, 4) = A(:, 4)*S3
    A(:, 5) = A(:, 5)*S3
    A(:, 6) = A(:, 6)*S3
  endif

  if (Number_Of_Derivative >= 4) then
    S4 = S3*S1

    A(:, 7) = A(:, 7)*S4
    A(:, 8) = A(:, 8)*S4
    A(:, 9) = A(:, 9)*S4
    A(:,10) = A(:,10)*S4
  endif

  END Subroutine Apply_Linear_Map 

END MODULE Grid_Auxillary_Subroutine
