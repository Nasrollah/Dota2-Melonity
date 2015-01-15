#include "pseudopack.h"


MODULE Molecular_Diffusion_module

  USE PseudoPack

! USE Central_Difference_Operators

  USE Navier_Stokes_Coefficients, ONLY: NS_N2, NS_N3, NS_MD_Order_x
  USE Navier_Stokes_Coefficients, ONLY: NS_M2, NS_M3, NS_MD_Order_y
  USE Navier_Stokes_Coefficients, ONLY: NS_K2, NS_K3, NS_MD_Order_z

  USE Molecular_Diffusion_Boundary_Condition

implicit NONE

INTERFACE PS_Molecular_Diffusion
  MODULE PROCEDURE Molecular_Diffusion_1D
  MODULE PROCEDURE Molecular_Diffusion_2D
  MODULE PROCEDURE Molecular_Diffusion_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Molecular_Diffusion

CONTAINS
!
!===================================================================
!
  Subroutine Molecular_Diffusion_1D (N0, N5, N2, N3, CD_x, &
                                     NV, NF, NS, zeta, Q, D_F, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: NV, NF, NS

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                           :: CD_x
  REALTYPE, dimension(N0:N5)                         :: zeta
  REALTYPE, dimension(N0:N5,NV)                      :: Q, D_F
  REALTYPE, dimension(N0:N5)                         :: f
  REALTYPE, dimension(N0:N5,1)                       :: Grad_f

  integer  :: n, ms

    NS_N2 = N2 - NS_MD_Order_x/2 ; NS_N3 = N3 + NS_MD_Order_x/2

  do n = 1, NS-1
    ms = NF+n

    f = Q(:,ms)/Q(:,1)

    call PS_CD_Grad (N0, N5, NS_N2, NS_N3, CD_x,   &
                     f, Grad_f)

    call PS_Molecular_Diffusion_BC (N0, N5, N2, N3,         &
                                    Grad_f, MPI_Comm_Type)

    Grad_f(:,1) = zeta*Grad_f(:,1)*Q(:,1)

    call PS_CD_Div  (N0, N5, N2, N3, CD_x,   &
                     Grad_f, D_F(:,ms))
  enddo

  END Subroutine Molecular_Diffusion_1D
!
!===================================================================
!
  Subroutine Molecular_Diffusion_2D (N0, N5, N2, N3, CD_x, &
                                     M0, M5, M2, M3, CD_y, &
                                     NV, NF, NS, zeta, Q, D_F, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV, NF, NS

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                           :: CD_x, CD_y
  REALTYPE, dimension(N0:N5,M0:M5)                   :: zeta
  REALTYPE, dimension(N0:N5,M0:M5,NV)                :: Q, D_F
  REALTYPE, dimension(N0:N5,M0:M5)                   :: f
  REALTYPE, dimension(N0:N5,M0:M5,2)                 :: Grad_f

  integer  :: n, ms

    NS_N2 = N2 - NS_MD_Order_x/2 ; NS_N3 = N3 + NS_MD_Order_x/2
    NS_M2 = M2 - NS_MD_Order_y/2 ; NS_M3 = M3 + NS_MD_Order_y/2

  do n = 1, NS-1
    ms = NF+n

    f = Q(:,:,ms)/Q(:,:,1)

    call PS_CD_Grad (N0, N5, NS_N2, NS_N3, CD_x,   &
                     M0, M5, NS_M2, NS_M3, CD_y,   &
                     f, Grad_f)

    call PS_Molecular_Diffusion_BC (N0, N5, N2, N3,         &
                                    M0, M5, M2, M3,         &
                                    Grad_f, MPI_Comm_Type)

    Grad_f(:,:,1) = zeta*Grad_f(:,:,1)*Q(:,:,1)
    Grad_f(:,:,2) = zeta*Grad_f(:,:,2)*Q(:,:,1)

    call PS_CD_Div  (N0, N5, N2, N3, CD_x,   &
                     M0, M5, M2, M3, CD_y,   &
                     Grad_f, D_F(:,:,ms))
  enddo

  END Subroutine Molecular_Diffusion_2D
!
!===================================================================
!
  Subroutine Molecular_Diffusion_3D (N0, N5, N2, N3, CD_x, &
                                     M0, M5, M2, M3, CD_y, &
                                     K0, K5, K2, K3, CD_z, &
                                     NV, NF, NS, zeta, Q, D_F, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV, NF, NS

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)                           :: CD_x, CD_y, CD_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)             :: zeta
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)          :: Q, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)             :: f
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,3)           :: Grad_f

  integer  :: n, ms

    NS_N2 = N2 - NS_MD_Order_x/2 ; NS_N3 = N3 + NS_MD_Order_x/2
    NS_M2 = M2 - NS_MD_Order_y/2 ; NS_M3 = M3 + NS_MD_Order_y/2
    NS_K2 = K2 - NS_MD_Order_z/2 ; NS_K3 = K3 + NS_MD_Order_z/2

  do n = 1, NS-1
    ms = NF+n

    f = Q(:,:,:,ms)/Q(:,:,:,1)

    call PS_CD_Grad (N0, N5, NS_N2, NS_N3, CD_x,   &
                     M0, M5, NS_M2, NS_M3, CD_y,   &
                     K0, K5, NS_K2, NS_K3, CD_z,   &
                     f, Grad_f)

    call PS_Molecular_Diffusion_BC (N0, N5, N2, N3,         &
                                    M0, M5, M2, M3,         &
                                    K0, K5, K2, K3,         &
                                    Grad_f, MPI_Comm_Type)

    Grad_f(:,:,:,1) = zeta*Grad_f(:,:,:,1)*Q(:,:,:,1)
    Grad_f(:,:,:,2) = zeta*Grad_f(:,:,:,2)*Q(:,:,:,1)
    Grad_f(:,:,:,3) = zeta*Grad_f(:,:,:,3)*Q(:,:,:,1)

    call PS_CD_Div  (N0, N5, N2, N3, CD_x,   &
                     M0, M5, M2, M3, CD_y,   &
                     K0, K5, K2, K3, CD_z,   &
                     Grad_f, D_F(:,:,:,ms))
  enddo

  END Subroutine Molecular_Diffusion_3D

END MODULE Molecular_Diffusion_module

