#include "pseudopack.h"


MODULE Viscous_Dissipation_module

  USE PseudoPack
  USE Stress_Tensor

implicit NONE

INTERFACE PS_Viscous_Dissipation
  MODULE PROCEDURE Viscous_Term_1D
  MODULE PROCEDURE Viscous_Term_2D
  MODULE PROCEDURE Viscous_Term_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Viscous_Dissipation

CONTAINS
!
!===================================================================
!
  Subroutine Viscous_Term_1D (N0, N5, N2, N3, CD_x,   &
                              NV, mu, U, D_F, Grad_T, & 
                              MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: NV

  integer  , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                           :: CD_x
  REALTYPE, dimension(N0:N5)                         :: mu, U
  REALTYPE, dimension(N0:N5,NV)                      :: D_F
  REALTYPE, dimension(N0:N5,1),             OPTIONAL :: Grad_T

  REALTYPE, dimension(N0:N5,1)                :: Tau
  REALTYPE, dimension(N0:N5)                  :: Temp

  Tau = ZERO ; Temp = ZERO 

  call PS_Stress_Tensor (N0, N5, N2, N3, CD_x,   &
                         U, Tau, MPI_Comm_Type)

    Tau(:,1) = mu*Tau(:,1)

!-----------------------------
  call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                      Tau(:,1), D_F(:,2))

!-----------------------------
                       Tau(:,1) = U*Tau(:,1)                ! U*T_xx 
  if (PRESENT(Grad_T)) Tau(:,1) =   Tau(:,1) + Grad_T(:,1)  ! U*T_xx + Grad_T

  call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                      Tau(:,1), D_F(:,3))

  D_F(:,3) = D_F(:,3) + Temp

  END Subroutine Viscous_Term_1D
!
!===================================================================
!
  Subroutine Viscous_Term_2D (N0, N5, N2, N3, CD_x,       &
                              M0, M5, M2, M3, CD_y,       &
                              NV, mu, U, V, D_F, Grad_T,  &
                              MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV

  integer  , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                           :: CD_x, CD_y
  REALTYPE, dimension(N0:N5,M0:M5)                   :: mu, U, V
  REALTYPE, dimension(N0:N5,M0:M5,NV)                :: D_F
  REALTYPE, dimension(N0:N5,M0:M5,2),       OPTIONAL :: Grad_T

  REALTYPE, dimension(N0:N5,M0:M5,3)          :: Tau
  REALTYPE, dimension(N0:N5,M0:M5)            :: Temp

  integer  :: n

  Tau = ZERO ; Temp = ZERO 

  call PS_Stress_Tensor (N0, N5, N2, N3, CD_x,   &
                         M0, M5, M2, M3, CD_y,   &
                         U, V, Tau, MPI_Comm_Type)

  do n = 1,3
    Tau(:,:,n) = mu*Tau(:,:,n)
  enddo

!-----------------------------
  call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      Tau(:,:,1), D_F(:,:,2))

  call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      Tau(:,:,3), Temp)

  D_F(:,:,2) = D_F(:,:,2) + Temp

!-----------------------------
  call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      Tau(:,:,3), D_F(:,:,3))

  call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      Tau(:,:,2), Temp)

  D_F(:,:,3) = D_F(:,:,3) + Temp

!-----------------------------
  Tau(:,:,1) = U*Tau(:,:,1) + V*Tau(:,:,3)   ! U*T_xx + V*T_xy
  Tau(:,:,2) = U*Tau(:,:,3) + V*Tau(:,:,2)   ! U*T_xy + V*T_yy

  if (PRESENT(Grad_T)) Tau(:,:,1:2) = Tau(:,:,1:2) + Grad_T(:,:,1:2)

  call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      Tau(:,:,1), D_F(:,:,4))

  call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      Tau(:,:,2), Temp)

  D_F(:,:,4) = D_F(:,:,4) + Temp

  END Subroutine Viscous_Term_2D
!
!===================================================================
!
  Subroutine Viscous_Term_3D (N0, N5, N2, N3, CD_x,         &
                              M0, M5, M2, M3, CD_y,         &
                              K0, K5, K2, K3, CD_z,         &
                              NV, mu, U, V, W, D_F, Grad_T, &
                              MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV

  integer  , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                           :: CD_x, CD_y, CD_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)             :: mu, U, V, W
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)          :: D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,3), OPTIONAL :: Grad_T

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,6)    :: Tau
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)      :: Temp

  integer  :: n

  Tau = ZERO ; Temp = ZERO 

  call PS_Stress_Tensor (N0, N5, N2, N3, CD_x,   &
                         M0, M5, M2, M3, CD_y,   &
                         K0, K5, K2, K3, CD_z,   &
                         U, V, W, Tau, MPI_Comm_Type)

  do n = 1,6
    Tau(:,:,:,n) = mu*Tau(:,:,:,n)
  enddo


!-----------------------------
  call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,1), D_F(:,:,:,2))   ! d Txx /dx

  call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,4), Temp)           ! d Txy /dy

  D_F(:,:,:,2) = D_F(:,:,:,2) + Temp

  call PS_CD_Diff (3, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,5), Temp)           ! d Txz /dz

  D_F(:,:,:,2) = D_F(:,:,:,2) + Temp

!-----------------------------
  call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,4), D_F(:,:,:,3))   ! d Txy /dx

  call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,2), Temp)           ! d Tyy /dy

  D_F(:,:,:,3) = D_F(:,:,:,3) + Temp

  call PS_CD_Diff (3, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,6), Temp)           ! d Tyz /dz

  D_F(:,:,:,3) = D_F(:,:,:,3) + Temp

!-----------------------------
  call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,5), D_F(:,:,:,4))   ! d Tzx /dx

  call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,6), Temp)           ! d Tzy /dy

  D_F(:,:,:,4) = D_F(:,:,:,4) + Temp

  call PS_CD_Diff (3, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,3), Temp)           ! d Tzz /dz

  D_F(:,:,:,4) = D_F(:,:,:,4) + Temp

!-----------------------------
  Tau(:,:,:,1) = U*Tau(:,:,:,1) + V*Tau(:,:,:,4) + W*Tau(:,:,:,5) 
  Tau(:,:,:,2) = U*Tau(:,:,:,4) + V*Tau(:,:,:,2) + W*Tau(:,:,:,6)  
  Tau(:,:,:,3) = U*Tau(:,:,:,5) + V*Tau(:,:,:,6) + W*Tau(:,:,:,3)  

  if (PRESENT(Grad_T)) Tau(:,:,:,1:3) = Tau(:,:,:,1:3) + Grad_T(:,:,:,1:3)

  call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,1), D_F(:,:,:,5))

  call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,2), Temp)

  D_F(:,:,:,5) = D_F(:,:,:,5) + Temp

  call PS_CD_Diff (3, N0, N5, N2, N3, 1, CD_x,   &
                      M0, M5, M2, M3, 1, CD_y,   &
                      K0, K5, K2, K3, 1, CD_z,   &
                      Tau(:,:,:,2), Temp)

  D_F(:,:,:,5) = D_F(:,:,:,5) + Temp

  END Subroutine Viscous_Term_3D

END MODULE Viscous_Dissipation_module
