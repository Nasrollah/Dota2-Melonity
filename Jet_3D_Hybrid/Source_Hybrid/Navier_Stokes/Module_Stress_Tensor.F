#include "pseudopack.h"


MODULE Stress_Tensor

  USE PseudoPack

  USE Navier_Stokes_Coefficients, ONLY: NS_N2, NS_N3, NS_VD_Order_x
  USE Navier_Stokes_Coefficients, ONLY: NS_M2, NS_M3, NS_VD_Order_y
  USE Navier_Stokes_Coefficients, ONLY: NS_K2, NS_K3, NS_VD_Order_z

  USE Viscous_Dissipation_Boundary_Condition

implicit NONE

INTERFACE PS_Stress_Tensor
  MODULE PROCEDURE Stress_Tensor_11
  MODULE PROCEDURE Stress_Tensor_21
  MODULE PROCEDURE Stress_Tensor_31
END INTERFACE

REALTYPE, parameter :: c21 = TWO
REALTYPE, parameter :: c23 = TWO/THREE

PRIVATE
PUBLIC :: PS_Stress_Tensor

CONTAINS
!
!===================================================================
!
  Subroutine Stress_Tensor_11 (N0, N5, N2, N3, CD_x,  &
                               U, Tau, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                      :: CD_x
  REALTYPE, dimension(N0:N5)                    :: U
  REALTYPE, dimension(N0:N5,1)                  :: Tau

    NS_N2 = N2 - NS_VD_Order_x/2 ; NS_N3 = N3 + NS_VD_Order_x/2

  call PS_CD_Diff (1, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      U, Tau(:,1))

  call PS_Stress_Tensor_BC    (N0, N5, N2, N3,         &
                               Tau, MPI_Comm_Type)

  END Subroutine Stress_Tensor_11
!
!===================================================================
!
  Subroutine Stress_Tensor_21 (N0, N5, N2, N3, CD_x,  &
                               M0, M5, M2, M3, CD_y,  &
                               U, V, Tau, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                      :: CD_x, CD_y
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V
  REALTYPE, dimension(N0:N5,M0:M5,3)            :: Tau

  REALTYPE, dimension(N0:N5,M0:M5)              :: Div_U

  Div_U = ZERO

    NS_N2 = N2 - NS_VD_Order_x/2 ; NS_N3 = N3 + NS_VD_Order_x/2
    NS_M2 = M2 - NS_VD_Order_y/2 ; NS_M3 = M3 + NS_VD_Order_y/2

  call PS_CD_Diff (1, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      U, Tau(:,:,1))             ! dUdx

  call PS_CD_Diff (2, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      V, Tau(:,:,2))             ! dVdy

  Div_U = Tau(:,:,1) + Tau(:,:,2)

  Tau(:,:,1) = c21*Tau(:,:,1) - c23*Div_U           ! Tau_xx
  Tau(:,:,2) = c21*Tau(:,:,2) - c23*Div_U           ! Tau_yy

! --------------------
  call PS_CD_Diff (2, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      U, Tau(:,:,3))             ! dUdy

  call PS_CD_Diff (1, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      V, Div_U)                  ! dVdx

  Tau(:,:,3) = Tau(:,:,3) + Div_U                   ! Tau_xy
! --------------------

  call PS_Stress_Tensor_BC    (N0, N5, N2, N3,         &
                               M0, M5, M2, M3,         &
                               Tau, MPI_Comm_Type)

  END Subroutine Stress_Tensor_21
!
!===================================================================
!
  Subroutine Stress_Tensor_31 (N0, N5, N2, N3, CD_x,  &
                               M0, M5, M2, M3, CD_y,  &
                               K0, K5, K2, K3, CD_z,  &
                               U, V, W, Tau, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                      :: CD_x, CD_y, CD_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,6)      :: Tau

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: Div_U

  Div_U = ZERO

    NS_N2 = N2 - NS_VD_Order_x/2 ; NS_N3 = N3 + NS_VD_Order_x/2
    NS_M2 = M2 - NS_VD_Order_y/2 ; NS_M3 = M3 + NS_VD_Order_y/2
    NS_K2 = K2 - NS_VD_Order_z/2 ; NS_K3 = K3 + NS_VD_Order_z/2

  call PS_CD_Diff (1, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      U, Tau(:,:,:,1))           ! dUdx

  call PS_CD_Diff (2, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      V, Tau(:,:,:,2))           ! dVdy

  call PS_CD_Diff (3, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      W, Tau(:,:,:,3))           ! dWdz

  Div_U = Tau(:,:,:,1) + Tau(:,:,:,2) + Tau(:,:,:,3)

  Tau(:,:,:,1) = c21*Tau(:,:,:,1) - c23*Div_U       ! Tau_xx
  Tau(:,:,:,2) = c21*Tau(:,:,:,2) - c23*Div_U       ! Tau_yy
  Tau(:,:,:,3) = c21*Tau(:,:,:,3) - c23*Div_U       ! Tau_zz

! --------------------
  call PS_CD_Diff (2, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      U, Tau(:,:,:,4))           ! dUdy

  call PS_CD_Diff (1, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      V, Div_U)                  ! dVdx

  Tau(:,:,:,4) = Tau(:,:,:,4) + Div_U               ! Tau_xy

! --------------------
  call PS_CD_Diff (3, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      U, Tau(:,:,:,5))           ! dUdz

  call PS_CD_Diff (1, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      W, Div_U)                  ! dWdx

  Tau(:,:,:,5) = Tau(:,:,:,5) + Div_U               ! Tau_xz

! --------------------
  call PS_CD_Diff (3, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      V, Tau(:,:,:,6))           ! dVdz

  call PS_CD_Diff (2, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      W, Div_U)                  ! dWdy

  Tau(:,:,:,6) = Tau(:,:,:,6) + Div_U               ! Tau_yz

  call PS_Stress_Tensor_BC    (N0, N5, N2, N3,         &
                               M0, M5, M2, M3,         &
                               K0, K5, K2, K3,         &
                               Tau, MPI_Comm_Type)

  END Subroutine Stress_Tensor_31

END MODULE Stress_Tensor
