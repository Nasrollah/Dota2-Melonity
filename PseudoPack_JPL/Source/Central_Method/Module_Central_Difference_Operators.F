#include "pseudopack.h"


MODULE Central_Difference_Operators

  USE Central_Difference_Function

implicit NONE

INTERFACE PS_CD_Grad
  MODULE PROCEDURE CD_Gradient_0D
  MODULE PROCEDURE CD_Gradient_1D
  MODULE PROCEDURE CD_Gradient_2D
  MODULE PROCEDURE CD_Gradient_3D
END INTERFACE

INTERFACE PS_CD_Div
  MODULE PROCEDURE CD_Divergence_0D
  MODULE PROCEDURE CD_Divergence_1D
  MODULE PROCEDURE CD_Divergence_2D
  MODULE PROCEDURE CD_Divergence_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_CD_Grad
PUBLIC :: PS_CD_Div

CONTAINS
!
!===================================================================
!
  Subroutine CD_Gradient_0D (N0, N5, N2, N3, CD_x, &
                             f, Grad_f)

  integer  :: N0, N5, N2, N3

  REALTYPE, dimension(:,:)                           :: CD_x
  REALTYPE, dimension(N0:N5)                         :: f
  REALTYPE, dimension(N0:N5)                         :: Grad_f

  Grad_f = ZERO 

    call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                        f, Grad_f)

  END Subroutine CD_Gradient_0D
!
!===================================================================
!
  Subroutine CD_Gradient_1D (N0, N5, N2, N3, CD_x, &
                             f, Grad_f)

  integer  :: N0, N5, N2, N3

  REALTYPE, dimension(:,:)                           :: CD_x
  REALTYPE, dimension(N0:N5)                         :: f
  REALTYPE, dimension(N0:N5,1)                       :: Grad_f

  Grad_f = ZERO 

    call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                        f, Grad_f(:,1))

  END Subroutine CD_Gradient_1D
!
!===================================================================
!
  Subroutine CD_Gradient_2D (N0, N5, N2, N3, CD_x, &
                             M0, M5, M2, M3, CD_y, &
                             f, Grad_f)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3

  REALTYPE, dimension(:,:)                           :: CD_x, CD_y
  REALTYPE, dimension(N0:N5,M0:M5)                   :: f
  REALTYPE, dimension(N0:N5,M0:M5,2)                 :: Grad_f

  Grad_f = ZERO 

    call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        f, Grad_f(:,:,1))

    call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        f, Grad_f(:,:,2))

  END Subroutine CD_Gradient_2D
!
!===================================================================
!
  Subroutine CD_Gradient_3D (N0, N5, N2, N3, CD_x, &
                             M0, M5, M2, M3, CD_y, &
                             K0, K5, K2, K3, CD_z, &
                             f, Grad_f)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  REALTYPE, dimension(:,:)                           :: CD_x, CD_y, CD_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)             :: f
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,3)           :: Grad_f

  Grad_f = ZERO 

    call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        K0, K5, K2, K3, 1, CD_z,   &
                        f, Grad_f(:,:,:,1))

    call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        K0, K5, K2, K3, 1, CD_z,   &
                        f, Grad_f(:,:,:,2))

    call PS_CD_Diff (3, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        K0, K5, K2, K3, 1, CD_z,   &
                        f, Grad_f(:,:,:,3))

  END Subroutine CD_Gradient_3D
!
!===================================================================
!
  Subroutine CD_Divergence_0D (N0, N5, N2, N3, CD_x, &
                               F, Div_F)

  integer  :: N0, N5, N2, N3

  REALTYPE, dimension(:,:)                           :: CD_x
  REALTYPE, dimension(N0:N5)                         :: F
  REALTYPE, dimension(N0:N5)                         :: Div_F

    Div_F = ZERO

    call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                        F, Div_F)

  END Subroutine CD_Divergence_0D
!
!===================================================================
!
  Subroutine CD_Divergence_1D (N0, N5, N2, N3, CD_x, &
                               F, Div_F)

  integer  :: N0, N5, N2, N3

  REALTYPE, dimension(:,:)                           :: CD_x
  REALTYPE, dimension(N0:N5,1)                       :: F
  REALTYPE, dimension(N0:N5)                         :: Div_F

    Div_F = ZERO

    call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                        F(:,1), Div_F)

  END Subroutine CD_Divergence_1D
!
!===================================================================
!
  Subroutine CD_Divergence_2D (N0, N5, N2, N3, CD_x, &
                               M0, M5, M2, M3, CD_y, &
                               F, Div_F)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3

  REALTYPE, dimension(:,:)                           :: CD_x, CD_y
  REALTYPE, dimension(N0:N5,M0:M5,2)                 :: F
  REALTYPE, dimension(N0:N5,M0:M5)                   :: Div_F
  REALTYPE, dimension(N0:N5,M0:M5)                   :: dF_1

    Div_F = ZERO ; dF_1 = ZERO 

    call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        F(:,:,1), dF_1)

    call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        F(:,:,2), Div_F)

    Div_F = Div_F + dF_1 

  END Subroutine CD_Divergence_2D
!
!===================================================================
!
  Subroutine CD_Divergence_3D (N0, N5, N2, N3, CD_x, &
                               M0, M5, M2, M3, CD_y, &
                               K0, K5, K2, K3, CD_z, &
                               F, Div_F)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  REALTYPE, dimension(:,:)                           :: CD_x, CD_y, CD_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,3)           :: F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)             :: Div_F

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)             :: dF_1, dF_2

    Div_F = ZERO ; dF_1 = ZERO ; dF_2 = ZERO 

    call PS_CD_Diff (1, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        K0, K5, K2, K3, 1, CD_z,   &
                        F(:,:,:,1), dF_1)

    call PS_CD_Diff (2, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        K0, K5, K2, K3, 1, CD_z,   &
                        F(:,:,:,2), dF_2)

    call PS_CD_Diff (3, N0, N5, N2, N3, 1, CD_x,   &
                        M0, M5, M2, M3, 1, CD_y,   &
                        K0, K5, K2, K3, 1, CD_z,   &
                        F(:,:,:,3), Div_F)

    Div_F = Div_F + dF_1 + dF_2

  END Subroutine CD_Divergence_3D

END MODULE Central_Difference_Operators
