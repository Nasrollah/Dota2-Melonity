! --------
#if 0
! --------
#include "pseudopack.h"


MODULE Thermal_Gradient_module

  USE PseudoPack

  USE Navier_Stokes_Coefficients, ONLY: NS_N2, NS_N3, NS_TC_Order_x
  USE Navier_Stokes_Coefficients, ONLY: NS_M2, NS_M3, NS_TC_Order_y
  USE Navier_Stokes_Coefficients, ONLY: NS_K2, NS_K3, NS_TC_Order_z

  USE Thermal_Gradient_Boundary_Condition

implicit NONE

INTERFACE PS_Thermal_Gradient
  MODULE PROCEDURE Thermal_Gradient_1D
  MODULE PROCEDURE Thermal_Gradient_2D
  MODULE PROCEDURE Thermal_Gradient_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Thermal_Gradient

CONTAINS
!
!===================================================================
!
  Subroutine Thermal_Gradient_1D (N0, N5, N2, N3, CD_x,  &
                                  kappa, T, Grad_T, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                    :: CD_x
  REALTYPE, dimension(N0:N5)                  :: kappa, T
  REALTYPE, dimension(N0:N5,1)                :: Grad_T

  Grad_T = ZERO

    NS_N2 = N2 - NS_TC_Order_x/2
    NS_N3 = N3 + NS_TC_Order_x/2

  call PS_CD_Diff (1, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      T, Grad_T(:,1))

  call PS_Thermal_Gradient_BC (N0, N5, N2, N3, &
                              Grad_T, MPI_Comm_Type)

  Grad_T(:,1) = kappa*Grad_T(:,1)

  END Subroutine Thermal_Gradient_1D
!
!===================================================================
!
  Subroutine Thermal_Gradient_2D (N0, N5, N2, N3, CD_x, &
                                  M0, M5, M2, M3, CD_y, &
                                  kappa, T, Grad_T, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                    :: CD_x, CD_y
  REALTYPE, dimension(N0:N5,M0:M5)            :: kappa, T
  REALTYPE, dimension(N0:N5,M0:M5,2)          :: Grad_T

  Grad_T = ZERO

    NS_N2 = N2 - NS_TC_Order_x/2
    NS_N3 = N3 + NS_TC_Order_x/2

    NS_M2 = M2 - NS_TC_Order_y/2
    NS_M3 = M3 + NS_TC_Order_y/2

  call PS_CD_Diff (1, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      T, Grad_T(:,:,1))

  call PS_CD_Diff (2, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      T, Grad_T(:,:,2))

  Grad_T(:,:,1) = kappa*Grad_T(:,:,1)
  Grad_T(:,:,2) = kappa*Grad_T(:,:,2)

  call PS_Thermal_Gradient_BC (N0, N5, N2, N3,         &
                               M0, M5, M2, M3,         &
                               Grad_T, MPI_Comm_Type)

  END Subroutine Thermal_Gradient_2D
!
!===================================================================
!
  Subroutine Thermal_Gradient_3D (N0, N5, N2, N3, CD_x,  &
                                  M0, M5, M2, M3, CD_y,  &
                                  K0, K5, K2, K3, CD_z,  &
                                  kappa, T, Grad_T, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                    :: CD_x, CD_y, CD_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)      :: kappa, T
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,3)    :: Grad_T

  Grad_T = ZERO

    NS_N2 = N2 - NS_TC_Order_x/2
    NS_N3 = N3 + NS_TC_Order_x/2

    NS_M2 = M2 - NS_TC_Order_y/2
    NS_M3 = M3 + NS_TC_Order_y/2

    NS_K2 = K2 - NS_TC_Order_z/2
    NS_K3 = K3 + NS_TC_Order_z/2

!-----------------------------
  call PS_CD_Diff (1, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      T, Grad_T(:,:,:,1))

  call PS_CD_Diff (2, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      T, Grad_T(:,:,:,2))

  call PS_CD_Diff (3, N0, N5, NS_N2, NS_N3, 1, CD_x,   &
                      M0, M5, NS_M2, NS_M3, 1, CD_y,   &
                      K0, K5, NS_K2, NS_K3, 1, CD_z,   &
                      T, Grad_T(:,:,:,3))

  call PS_Thermal_Gradient_BC (N0, N5, N2, N3,         &
                               M0, M5, M2, M3,         &
                               K0, K5, K2, K3,         &
                               Grad_T, MPI_Comm_Type)

  Grad_T(:,:,:,1) = kappa*Grad_T(:,:,:,1)
  Grad_T(:,:,:,2) = kappa*Grad_T(:,:,:,2)
  Grad_T(:,:,:,3) = kappa*Grad_T(:,:,:,3)

  END Subroutine Thermal_Gradient_3D

END MODULE Thermal_Gradient_module

! --------
#else
! --------

#include "pseudopack.h"


MODULE Thermal_Gradient_module

  USE PseudoPack

! USE Central_Difference_Operators

  USE Navier_Stokes_Coefficients, ONLY: NS_N2, NS_N3, NS_TC_Order_x
  USE Navier_Stokes_Coefficients, ONLY: NS_M2, NS_M3, NS_TC_Order_y
  USE Navier_Stokes_Coefficients, ONLY: NS_K2, NS_K3, NS_TC_Order_z

  USE Thermal_Gradient_Boundary_Condition

implicit NONE

INTERFACE PS_Thermal_Gradient
  MODULE PROCEDURE Thermal_Gradient_1D
  MODULE PROCEDURE Thermal_Gradient_2D
  MODULE PROCEDURE Thermal_Gradient_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Thermal_Gradient

CONTAINS
!
!===================================================================
!
  Subroutine Thermal_Gradient_1D (N0, N5, N2, N3, CD_x,  &
                                  kappa, T, Grad_T, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                    :: CD_x
  REALTYPE, dimension(N0:N5)                  :: kappa, T
  REALTYPE, dimension(N0:N5,1)                :: Grad_T

    NS_N2 = N2 - NS_TC_Order_x/2 ; NS_N3 = N3 + NS_TC_Order_x/2

!-----------------------------
  call PS_CD_Grad (N0, N5, NS_N2, NS_N3, CD_x,   &
                   T, Grad_T)

  call PS_Thermal_Gradient_BC (N0, N5, N2, N3,         &
                               Grad_T, MPI_Comm_Type)

  Grad_T(:,1) = kappa*Grad_T(:,1)

  END Subroutine Thermal_Gradient_1D
!
!===================================================================
!
  Subroutine Thermal_Gradient_2D (N0, N5, N2, N3, CD_x,  &
                                  M0, M5, M2, M3, CD_y,  &
                                  kappa, T, Grad_T, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                    :: CD_x, CD_y
  REALTYPE, dimension(N0:N5,M0:M5)            :: kappa, T
  REALTYPE, dimension(N0:N5,M0:M5,2)          :: Grad_T

    NS_N2 = N2 - NS_TC_Order_x/2 ; NS_N3 = N3 + NS_TC_Order_x/2
    NS_M2 = M2 - NS_TC_Order_y/2 ; NS_M3 = M3 + NS_TC_Order_y/2

!-----------------------------
  call PS_CD_Grad (N0, N5, NS_N2, NS_N3, CD_x,   &
                   M0, M5, NS_M2, NS_M3, CD_y,   &
                   T, Grad_T)

  call PS_Thermal_Gradient_BC (N0, N5, N2, N3,         &
                               M0, M5, M2, M3,         &
                               Grad_T, MPI_Comm_Type)

  Grad_T(:,:,1) = kappa*Grad_T(:,:,1)
  Grad_T(:,:,2) = kappa*Grad_T(:,:,2)

  END Subroutine Thermal_Gradient_2D
!
!===================================================================
!
  Subroutine Thermal_Gradient_3D (N0, N5, N2, N3, CD_x,  &
                                  M0, M5, M2, M3, CD_y,  &
                                  K0, K5, K2, K3, CD_z,  &
                                  kappa, T, Grad_T, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                    :: CD_x, CD_y, CD_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)      :: kappa, T
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,3)    :: Grad_T

    NS_N2 = N2 - NS_TC_Order_x/2 ; NS_N3 = N3 + NS_TC_Order_x/2
    NS_M2 = M2 - NS_TC_Order_y/2 ; NS_M3 = M3 + NS_TC_Order_y/2
    NS_K2 = K2 - NS_TC_Order_z/2 ; NS_K3 = K3 + NS_TC_Order_z/2

!-----------------------------
  call PS_CD_Grad (N0, N5, NS_N2, NS_N3, CD_x,   &
                   M0, M5, NS_M2, NS_M3, CD_y,   &
                   K0, K5, NS_K2, NS_K3, CD_z,   &
                   T, Grad_T)

  call PS_Thermal_Gradient_BC (N0, N5, N2, N3,         &
                               M0, M5, M2, M3,         &
                               K0, K5, K2, K3,         &
                               Grad_T, MPI_Comm_Type)

  Grad_T(:,:,:,1) = kappa*Grad_T(:,:,:,1)
  Grad_T(:,:,:,2) = kappa*Grad_T(:,:,:,2)
  Grad_T(:,:,:,3) = kappa*Grad_T(:,:,:,3)

  END Subroutine Thermal_Gradient_3D

END MODULE Thermal_Gradient_module

! --------
#endif
! --------
