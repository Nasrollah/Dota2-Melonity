#include "pseudopack.h"

MODULE WENO_Derivative

  USE WENO_Differentiation

implicit NONE

INTERFACE PS_WENO_Derivative
  MODULE PROCEDURE PS_WENO_Derivative_1D
  MODULE PROCEDURE PS_WENO_Derivative_2D
  MODULE PROCEDURE PS_WENO_Derivative_3D

  MODULE PROCEDURE PS_WENO_Derivative_0_1D
  MODULE PROCEDURE PS_WENO_Derivative_0_2D
  MODULE PROCEDURE PS_WENO_Derivative_0_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_WENO_Derivative

CONTAINS

!
!===================================================================
!
  Subroutine PS_WENO_Derivative_1D (Index, Order, dx,           &
                                      N0, N5, N2, N3,           &
                                    NV, Q, D_Q                  )

  integer  :: Index, Order, NV
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,NV) :: Q, D_Q

  integer  :: n

  do n = 1,NV
    call PS_Diff_WENO   (Index, Order, N0, N5, N2, N3, dx,             &
                                Q(:,n), D_Q(:,n))
  enddo

  END Subroutine PS_WENO_Derivative_1D
!
!===================================================================
!
  Subroutine PS_WENO_Derivative_2D (Index, Order, dx,           &
                                      N0, N5, N2, N3,           &
                                      M0, M5, M2, M3,           &
                                    NV, Q, D_Q                  )

  integer  :: NV
  integer  :: Index, Order
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, D_Q

  integer  :: n

  do n = 1,NV
    call PS_Diff_WENO   (Index, Order, N0, N5, N2, N3, dx,             &
                                Order, M0, M5, M2, M3, dx,             &
                                Q(:,:,n), D_Q(:,:,n))
  enddo

  END Subroutine PS_WENO_Derivative_2D
!
!===================================================================
!
  Subroutine PS_WENO_Derivative_3D (Index, Order, dx,           &
                                      N0, N5, N2, N3,           &
                                      M0, M5, M2, M3,           &
                                      K0, K5, K2, K3,           &
                                    NV, Q, D_Q                  )

  integer  :: Index, Order, NV
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Q

  integer  :: n

  do n = 1,NV
    call PS_Diff_WENO   (Index, Order, N0, N5, N2, N3, dx,             &
                                Order, M0, M5, M2, M3, dx,             &
                                Order, K0, K5, K2, K3, dx,             &
                                Q(:,:,:,n), D_Q(:,:,:,n))
  enddo

  END Subroutine PS_WENO_Derivative_3D
!
!===================================================================
!
  Subroutine PS_WENO_Derivative_0_1D (Index, Order, dx,           &
                                        N0, N5, N2, N3,           &
                                      Q, D_Q                      )

  integer  :: Index, Order
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5) :: Q, D_Q

    call PS_Diff_WENO   (Index, Order, N0, N5, N2, N3, dx,             &
                                Q(:), D_Q(:))

  END Subroutine PS_WENO_Derivative_0_1D
!
!===================================================================
!
  Subroutine PS_WENO_Derivative_0_2D (Index, Order, dx,           &
                                        N0, N5, N2, N3,           &
                                        M0, M5, M2, M3,           &
                                      Q, D_Q                      )

  integer  :: NV
  integer  :: Index, Order
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,M0:M5) :: Q, D_Q

    call PS_Diff_WENO   (Index, Order, N0, N5, N2, N3, dx,             &
                                Order, M0, M5, M2, M3, dx,             &
                                Q(:,:), D_Q(:,:))

  END Subroutine PS_WENO_Derivative_0_2D
!
!===================================================================
!
  Subroutine PS_WENO_Derivative_0_3D (Index, Order, dx,           &
                                        N0, N5, N2, N3,           &
                                        M0, M5, M2, M3,           &
                                        K0, K5, K2, K3,           &
                                      Q, D_Q                      )

  integer  :: Index, Order
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Q, D_Q

    call PS_Diff_WENO   (Index, Order, N0, N5, N2, N3, dx,             &
                                Order, M0, M5, M2, M3, dx,             &
                                Order, K0, K5, K2, K3, dx,             &
                                Q(:,:,:), D_Q(:,:,:))

  END Subroutine PS_WENO_Derivative_0_3D

END MODULE WENO_Derivative
