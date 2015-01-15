#include "pseudopack.h"


MODULE Central_Derivative

  USE Central_Differentiation

IMPLICIT NONE

INTERFACE PS_Central_Derivative
  MODULE PROCEDURE PS_Central_Derivative_1D
  MODULE PROCEDURE PS_Central_Derivative_2D
  MODULE PROCEDURE PS_Central_Derivative_3D

  MODULE PROCEDURE PS_Central_Derivative_0_1D
  MODULE PROCEDURE PS_Central_Derivative_0_2D
  MODULE PROCEDURE PS_Central_Derivative_0_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Central_Derivative

CONTAINS

!
!===================================================================
!
  Subroutine PS_Central_Derivative_1D (Index, Order, Degree, Nx, dx, &
                                       N0, N5, N2, N3,               &
                                       NV, Q, D_Q                    )

  integer  :: Index, Order, Degree, Nx, NV
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,NV) :: Q, D_Q

  integer  :: n

  do n = 1,NV
    call PS_Central_Diff (Index,                                   &
                          Order, Degree, Nx, N0, N5, N2, N3, dx,   &
                          Q(:,n), D_Q(:,n))
  enddo

  END Subroutine PS_Central_Derivative_1D 
!
!===================================================================
!
  Subroutine PS_Central_Derivative_2D (Index, Order, Degree, Nx, dx, &
                                       N0, N5, N2, N3,               &
                                       M0, M5, M2, M3,               &
                                       NV, Q, D_Q                    )

  integer  :: Index, Order, Degree, Nx, NV
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, D_Q

  integer  :: n

  do n = 1,NV
    call PS_Central_Diff (Index,                                   &
                          Order, Degree, Nx, N0, N5, N2, N3, dx,   &
                          Order, Degree, Nx, M0, M5, M2, M3, dx,   &
                          Q(:,:,n), D_Q(:,:,n))
  enddo

  END Subroutine PS_Central_Derivative_2D 
!
!===================================================================
!
  Subroutine PS_Central_Derivative_3D (Index, Order, Degree, Nx, dx, &
                                       N0, N5, N2, N3,               &
                                       M0, M5, M2, M3,               &
                                       K0, K5, K2, K3,               &
                                       NV, Q, D_Q                    )

  integer  :: Index, Order, Degree, Nx, NV
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Q

  integer  :: n

  do n = 1,NV
    call PS_Central_Diff (Index,                                   &
                          Order, Degree, Nx, N0, N5, N2, N3, dx,   &
                          Order, Degree, Nx, M0, M5, M2, M3, dx,   &
                          Order, Degree, Nx, K0, K5, K2, K3, dx,   &
                          Q(:,:,:,n), D_Q(:,:,:,n))
  enddo

  END Subroutine PS_Central_Derivative_3D 



!
!===================================================================
!
  Subroutine PS_Central_Derivative_0_1D (Index, Order, Degree, Nx, dx, &
                                         N0, N5, N2, N3,               &
                                         Q, D_Q                        )

  integer  :: Index, Order, Degree, Nx
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5) :: Q, D_Q

    call PS_Central_Diff (Index,                                   &
                          Order, Degree, Nx, N0, N5, N2, N3, dx,   &
                          Q(:), D_Q(:))

  END Subroutine PS_Central_Derivative_0_1D 
!
!===================================================================
!
  Subroutine PS_Central_Derivative_0_2D (Index, Order, Degree, Nx, dx, &
                                         N0, N5, N2, N3,               &
                                         M0, M5, M2, M3,               &
                                         Q, D_Q                        )

  integer  :: Index, Order, Degree, Nx
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,M0:M5) :: Q, D_Q

    call PS_Central_Diff (Index,                                   &
                          Order, Degree, Nx, N0, N5, N2, N3, dx,   &
                          Order, Degree, Nx, M0, M5, M2, M3, dx,   &
                          Q(:,:), D_Q(:,:))

  END Subroutine PS_Central_Derivative_0_2D 
!
!===================================================================
!
  Subroutine PS_Central_Derivative_0_3D (Index, Order, Degree, Nx, dx, &
                                         N0, N5, N2, N3,               &
                                         M0, M5, M2, M3,               &
                                         K0, K5, K2, K3,               &
                                         Q, D_Q                        )

  integer  :: Index, Order, Degree, Nx
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Q, D_Q

    call PS_Central_Diff (Index,                                   &
                          Order, Degree, Nx, N0, N5, N2, N3, dx,   &
                          Order, Degree, Nx, M0, M5, M2, M3, dx,   &
                          Order, Degree, Nx, K0, K5, K2, K3, dx,   &
                          Q(:,:,:), D_Q(:,:,:))

  END Subroutine PS_Central_Derivative_0_3D 

END MODULE Central_Derivative
