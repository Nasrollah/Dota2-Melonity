#include "pseudopack.h"


MODULE Central_Euler_Flux

  USE Euler_Function
  USE Specific_Heat_Gamma
  USE Central_Difference_Function

IMPLICIT NONE

integer , parameter :: Order = 1

INTERFACE PS_CD_Euler_Flux
  MODULE PROCEDURE Euler_Flux_1D_0
  MODULE PROCEDURE Euler_Flux_2D_0
  MODULE PROCEDURE Euler_Flux_3D_0

  MODULE PROCEDURE Euler_Flux_1D_1
  MODULE PROCEDURE Euler_Flux_2D_1
  MODULE PROCEDURE Euler_Flux_3D_1

  MODULE PROCEDURE Euler_Flux_1D_2
  MODULE PROCEDURE Euler_Flux_2D_2
  MODULE PROCEDURE Euler_Flux_3D_2

  MODULE PROCEDURE Euler_Flux_1D_3
  MODULE PROCEDURE Euler_Flux_2D_3
  MODULE PROCEDURE Euler_Flux_3D_3

  MODULE PROCEDURE Euler_Flux_1D_a
  MODULE PROCEDURE Euler_Flux_2D_a
  MODULE PROCEDURE Euler_Flux_3D_a

  MODULE PROCEDURE Euler_Flux_1D_b
  MODULE PROCEDURE Euler_Flux_2D_b
  MODULE PROCEDURE Euler_Flux_3D_b

  MODULE PROCEDURE Euler_Flux_1D_c
  MODULE PROCEDURE Euler_Flux_2D_c
  MODULE PROCEDURE Euler_Flux_3D_c

  MODULE PROCEDURE Euler_Flux_1D_d
  MODULE PROCEDURE Euler_Flux_2D_d
  MODULE PROCEDURE Euler_Flux_3D_d
END INTERFACE

PRIVATE
PUBLIC :: PS_CD_Euler_Flux

CONTAINS
!
!==========================================================================
!
  Subroutine Euler_Flux_1D_0 (N0, N5, N2, N3, NV,              &
                              Q, U, P,                         &
                              CD_x, D_F, alpha, D_F_Epsilon, S)

  integer  :: N0, N5, N2, N3
  integer  :: NV

  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: U, P
  REALTYPE, dimension(N0:N5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x
  REALTYPE, dimension(1)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO 

  call PS_Euler_Flux (1, N0, N5, NV, Q, U, P, F, S)
  call PS_CD_Diff    (1, N0, N5, N2, N3, Order, CD_x,                 &
                         NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  END Subroutine Euler_Flux_1D_0
!
!==========================================================================
!
  Subroutine Euler_Flux_2D_0 (N0, N5, N2, N3,                  &
                              M0, M5, M2, M3, NV,              &
                              Q, U, V, P,                      &
                              CD_x, CD_y, D_F, alpha, D_F_Epsilon, S)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, P
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x, CD_y
  REALTYPE, dimension(2)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO ; D_G = ZERO

  call PS_Euler_Flux (1, N0, N5, M0, M5, NV, Q, U, V, P, F, S)
  call PS_CD_Diff    (1, N0, N5, N2, N3, Order, CD_x,                 &
                         M0, M5, M2, M3, Order, CD_y,                 &
                         NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  call PS_Euler_Flux (2, N0, N5, M0, M5, NV, Q, U, V, P, F, S)
  call PS_CD_Diff    (2, N0, N5, N2, N3, Order, CD_x,                 &
                         M0, M5, M2, M3, Order, CD_y,                 &
                         NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G

  END Subroutine Euler_Flux_2D_0 
!
!==========================================================================
!
  Subroutine Euler_Flux_3D_0 (N0, N5, N2, N3,                  &
                              M0, M5, M2, M3,                  &
                              K0, K5, K2, K3, NV,              &
                              Q, U, V, W, P,                   &
                              CD_x, CD_y, CD_z, D_F, alpha, D_F_Epsilon, S)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x, CD_y, CD_z
  REALTYPE, dimension(3)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO ; D_G = ZERO

  call PS_Euler_Flux (1, N0, N5, M0, M5, K0, K5, NV, Q, U, V, W, P, F, S)
  call PS_CD_Diff    (1, N0, N5, N2, N3, Order, CD_x,                 &
                         M0, M5, M2, M3, Order, CD_y,                 &
                         K0, K5, K2, K3, Order, CD_z,                 &
                         NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  call PS_Euler_Flux (2, N0, N5, M0, M5, K0, K5, NV, Q, U, V, W, P, F, S)
  call PS_CD_Diff    (2, N0, N5, N2, N3, Order, CD_x,                 &
                         M0, M5, M2, M3, Order, CD_y,                 &
                         K0, K5, K2, K3, Order, CD_z,                 &
                         NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G

  call PS_Euler_Flux (3, N0, N5, M0, M5, K0, K5, NV, Q, U, V, W, P, F, S)
  call PS_CD_Diff    (3, N0, N5, N2, N3, Order, CD_x,                 &
                         M0, M5, M2, M3, Order, CD_y,                 &
                         K0, K5, K2, K3, Order, CD_z,                 &
                         NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G 

  END Subroutine Euler_Flux_3D_0
!
!==========================================================================
!
  Subroutine Euler_Flux_1D_1 (N0, N5, N2, N3, NV,              &
                              Gamma, Q,                        &
                              CD_x, D_F, alpha, D_F_Epsilon, S)

  integer  :: N0, N5, N2, N3
  integer  :: NV

  REALTYPE :: Gamma

  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: U, P
  REALTYPE, dimension(N0:N5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x
  REALTYPE, dimension(1)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO 

  call PS_Euler_Primitive_Variable (   N0, N5, NV, Gamma, Q, U, P      )

  call PS_Euler_Flux               (1, N0, N5, NV,        Q, U, P, F, S)
  call PS_CD_Diff                  (1, N0, N5, N2, N3, Order, CD_x,    &
                                       NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  END Subroutine Euler_Flux_1D_1
!
!==========================================================================
!
  Subroutine Euler_Flux_2D_1 (N0, N5, N2, N3,                  &
                              M0, M5, M2, M3, NV,              &
                              Gamma, Q,                        &
                              CD_x, CD_y, D_F, alpha, D_F_Epsilon, S)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV

  REALTYPE :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, P
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x, CD_y
  REALTYPE, dimension(2)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO ; D_G = ZERO 

  call PS_Euler_Primitive_Variable (   N0, N5, M0, M5, NV,            &
                                    Gamma, Q, U, V, P)

  call PS_Euler_Flux               (1, N0, N5, M0, M5, NV,            &
                                           Q, U, V, P, F, S)
  call PS_CD_Diff                  (1, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  call PS_Euler_Flux               (2, N0, N5, M0, M5, NV,            &
                                           Q, U, V, P, F, S)
  call PS_CD_Diff                  (2, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G

  END Subroutine Euler_Flux_2D_1 
!
!==========================================================================
!
  Subroutine Euler_Flux_3D_1 (N0, N5, N2, N3,                  &
                              M0, M5, M2, M3,                  &
                              K0, K5, K2, K3, NV,              &
                              Gamma, Q,                        &
                              CD_x, CD_y, CD_z, D_F, alpha, D_F_Epsilon, S)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV

  REALTYPE :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x, CD_y, CD_z
  REALTYPE, dimension(3)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO ; D_G = ZERO 

  call PS_Euler_Primitive_Variable (   N0, N5, M0, M5, K0, K5, NV,    &
                                    Gamma, Q, U, V, W, P)

  call PS_Euler_Flux               (1, N0, N5, M0, M5, K0, K5, NV,    &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (1, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       K0, K5, K2, K3, Order, CD_z,   &
                                       NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  call PS_Euler_Flux               (2, N0, N5, M0, M5, K0, K5, NV,    &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (2, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       K0, K5, K2, K3, Order, CD_z,   &
                                       NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G

  call PS_Euler_Flux               (3, N0, N5, M0, M5, K0, K5, NV,    &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (3, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       K0, K5, K2, K3, Order, CD_z,   &
                                       NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G 

  END Subroutine Euler_Flux_3D_1
!
!==========================================================================
!
  Subroutine Euler_Flux_1D_2 (N0, N5, N2, N3, NV,              &
                              Gamma, Q,                        &
                              CD_x, D_F, alpha, D_F_Epsilon, S)

  integer  :: N0, N5, N2, N3
  integer  :: NV

  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: U, P, Gamma
  REALTYPE, dimension(N0:N5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x
  REALTYPE, dimension(1)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO 

  call PS_Euler_Primitive_Variable (   N0, N5, NV,                    &
                                    Gamma, Q, U, P)

  call PS_Euler_Flux               (1, N0, N5, NV,        Q, U, P, F, S)
  call PS_CD_Diff                  (1, N0, N5, N2, N3, Order, CD_x,   &
                                       NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  END Subroutine Euler_Flux_1D_2
!
!==========================================================================
!
  Subroutine Euler_Flux_2D_2 (N0, N5, N2, N3,                  &
                              M0, M5, M2, M3, NV,              &
                              Gamma, Q,                        &
                              CD_x, CD_y, D_F, alpha, D_F_Epsilon, S)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, P, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x, CD_y
  REALTYPE, dimension(2)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO ; D_G = ZERO 

  call PS_Euler_Primitive_Variable (   N0, N5, M0, M5, NV,            &
                                    Gamma, Q, U, V, P)

  call PS_Euler_Flux               (1, N0, N5, M0, M5, NV,            &
                                           Q, U, V, P, F, S)
  call PS_CD_Diff                  (1, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  call PS_Euler_Flux               (2, N0, N5, M0, M5, NV,            &
                                           Q, U, V, P, F, S)
  call PS_CD_Diff                  (2, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G

  END Subroutine Euler_Flux_2D_2 
!
!==========================================================================
!
  Subroutine Euler_Flux_3D_2 (N0, N5, N2, N3,                  &
                              M0, M5, M2, M3,                  &
                              K0, K5, K2, K3, NV,              &
                              Gamma, Q,                        &
                              CD_x, CD_y, CD_z, D_F, alpha, D_F_Epsilon, S)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x, CD_y, CD_z
  REALTYPE, dimension(3)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO ; D_G = ZERO 

  call PS_Euler_Primitive_Variable (   N0, N5, M0, M5, K0, K5, NV,    &
                                    Gamma, Q, U, V, W, P)

  call PS_Euler_Flux               (1, N0, N5, M0, M5, K0, K5, NV,    &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (1, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       K0, K5, K2, K3, Order, CD_z,   &
                                       NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  call PS_Euler_Flux               (2, N0, N5, M0, M5, K0, K5, NV,    &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (2, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       K0, K5, K2, K3, Order, CD_z,   &
                                       NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G

  call PS_Euler_Flux               (3, N0, N5, M0, M5, K0, K5, NV,    &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (3, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       K0, K5, K2, K3, Order, CD_z,   &
                                       NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G 

  END Subroutine Euler_Flux_3D_2
!
!==========================================================================
!
  Subroutine Euler_Flux_1D_3 (N0, N5, N2, N3, NV,               &
                              Gamma_Formulation, Gamma_Gas, Q,  &
                              CD_x, D_F, alpha, D_F_Epsilon, S, &
                              Specific_Heat)

  integer  :: N0, N5, N2, N3
  integer  :: NV

  integer                  :: Gamma_Formulation
  REALTYPE, dimension(:,:) :: Gamma_Gas

  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: U, P, Gamma
  REALTYPE, dimension(N0:N5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x
  REALTYPE, dimension(1)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  EXTERNAL    Specific_Heat
  OPTIONAL :: Specific_Heat

  D_F = ZERO 

  if (PRESENT(Specific_Heat)) then
    call    Specific_Heat          (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  else
    call PS_Specific_Heat_Gamma    (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  endif

  call PS_Euler_Primitive_Variable (   N0, N5, NV,                    &
                                    Gamma, Q, U, P)

  call PS_Euler_Flux               (1, N0, N5, NV,        Q, U, P, F, S)
  call PS_CD_Diff                  (1, N0, N5, N2, N3, Order, CD_x,   &
                                       NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  END Subroutine Euler_Flux_1D_3
!
!==========================================================================
!
  Subroutine Euler_Flux_2D_3 (N0, N5, N2, N3,                         &
                              M0, M5, M2, M3, NV,                     &
                              Gamma_Formulation, Gamma_Gas, Q,        &
                              CD_x, CD_y, D_F, alpha, D_F_Epsilon, S, &
                              Specific_Heat)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV

  integer                  :: Gamma_Formulation
  REALTYPE, dimension(:,:) :: Gamma_Gas

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, P, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x, CD_y
  REALTYPE, dimension(2)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  EXTERNAL    Specific_Heat
  OPTIONAL :: Specific_Heat

  D_F = ZERO ; D_G = ZERO 

  if (PRESENT(Specific_Heat)) then
    call    Specific_Heat          (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  else
    call PS_Specific_Heat_Gamma    (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  endif

  call PS_Euler_Primitive_Variable (   N0, N5, M0, M5, NV,            &
                                    Gamma, Q, U, V, P)

  call PS_Euler_Flux               (1, N0, N5, M0, M5, NV,            &
                                           Q, U, V, P, F, S)
  call PS_CD_Diff                  (1, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  call PS_Euler_Flux               (2, N0, N5, M0, M5, NV,            &
                                           Q, U, V, P, F, S)
  call PS_CD_Diff                  (2, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G

  END Subroutine Euler_Flux_2D_3
!
!==========================================================================
!
  Subroutine Euler_Flux_3D_3 (N0, N5, N2, N3,                               &
                              M0, M5, M2, M3,                               &
                              K0, K5, K2, K3, NV,                           &
                              Gamma_Formulation, Gamma_Gas, Q,              &
                              CD_x, CD_y, CD_z, D_F, alpha, D_F_Epsilon, S, &
                              Specific_Heat)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV

  integer                  :: Gamma_Formulation
  REALTYPE, dimension(:,:) :: Gamma_Gas

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F, D_G

  REALTYPE, dimension(:,:)            :: CD_x, CD_y, CD_z
  REALTYPE, dimension(3)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  EXTERNAL    Specific_Heat
  OPTIONAL :: Specific_Heat

  D_F = ZERO ; D_G = ZERO 

  if (PRESENT(Specific_Heat)) then
    call    Specific_Heat          (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  else
    call PS_Specific_Heat_Gamma    (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  endif

  call PS_Euler_Primitive_Variable (   N0, N5, M0, M5, K0, K5, NV,    &
                                    Gamma, Q, U, V, W, P)

  call PS_Euler_Flux               (1, N0, N5, M0, M5, K0, K5, NV,    &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (1, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       K0, K5, K2, K3, Order, CD_z,   &
                                       NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  call PS_Euler_Flux               (2, N0, N5, M0, M5, K0, K5, NV,    &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (2, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       K0, K5, K2, K3, Order, CD_z,   &
                                       NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G

  call PS_Euler_Flux               (3, N0, N5, M0, M5, K0, K5, NV,    &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (3, N0, N5, N2, N3, Order, CD_x,   &
                                       M0, M5, M2, M3, Order, CD_y,   &
                                       K0, K5, K2, K3, Order, CD_z,   &
                                       NV, F, D_G, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G) <  D_F_Epsilon) D_G = ZERO 

  D_F = D_F + D_G 

  END Subroutine Euler_Flux_3D_3
!
!==========================================================================
!
  Subroutine Euler_Flux_1D_a (Index, N0, N5, N2, N3, NV,              &
                              Q, U, P,                                &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: NV

  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: U, P
  REALTYPE, dimension(N0:N5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(1)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO 

  call PS_Euler_Flux (Index, N0, N5, NV,                                &
                             Q, U, P, F, S)
  call PS_CD_Diff    (Index, N0, N5, N2, N3, Order, CD,                 &
                             NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_1D_a
!
!==========================================================================
!
  Subroutine Euler_Flux_2D_a (Index, N0, N5, N2, N3,                  &
                                     M0, M5, M2, M3, NV,              &
                              Q, U, V, P,                             &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, P
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(2)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO 

  call PS_Euler_Flux (Index, N0, N5, M0, M5, NV,                        &
                             Q, U, V, P, F, S)
  call PS_CD_Diff    (Index, N0, N5, N2, N3, Order, CD,                 &
                             M0, M5, M2, M3, Order, CD,                 &
                             NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_2D_a 
!
!==========================================================================
!
  Subroutine Euler_Flux_3D_a (Index, N0, N5, N2, N3,                  &
                                     M0, M5, M2, M3,                  &
                                     K0, K5, K2, K3, NV,              &
                             Q, U, V, W, P,                           &
                             CD, D_F, alpha, D_F_Epsilon, Flux, S)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(3)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO 

  call PS_Euler_Flux (Index, N0, N5, M0, M5, K0, K5, NV,                &
                             Q, U, V, W, P, F, S)
  call PS_CD_Diff    (Index, N0, N5, N2, N3, Order, CD,                 &
                             M0, M5, M2, M3, Order, CD,                 &
                             K0, K5, K2, K3, Order, CD,                 &
                             NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_3D_a
!
!==========================================================================
!
  Subroutine Euler_Flux_1D_b (Index, N0, N5, N2, N3, NV,              &
                              Gamma, Q,                               &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: NV

  REALTYPE :: Gamma

  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: U, P
  REALTYPE, dimension(N0:N5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(1)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO 

  call PS_Euler_Primitive_Variable (       N0, N5, NV,                  &
                                    Gamma, Q, U, P)

  call PS_Euler_Flux               (Index, N0, N5, NV,                  &
                                           Q, U, P, F, S)
  call PS_CD_Diff                  (Index, N0, N5, N2, N3, Order, CD,   &
                                           NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_1D_b
!
!==========================================================================
!
  Subroutine Euler_Flux_2D_b (Index, N0, N5, N2, N3,                  &
                                     M0, M5, M2, M3, NV,              &
                              Gamma, Q,                               &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV

  REALTYPE :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, P
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(2)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO  

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, NV,          &
                                    Gamma, Q, U, V, P)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, NV,          &
                                           Q, U, V, P, F, S)
  call PS_CD_Diff                  (Index, N0, N5, N2, N3, Order, CD,   &
                                           M0, M5, M2, M3, Order, CD,   &
                                           NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_2D_b 
!
!==========================================================================
!
  Subroutine Euler_Flux_3D_b (Index, N0, N5, N2, N3,                  &
                                     M0, M5, M2, M3,                  &
                                     K0, K5, K2, K3, NV,              &
                              Gamma, Q,                               &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV

  REALTYPE :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(3)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO  

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, K0, K5, NV,  &
                                    Gamma, Q, U, V, W, P)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, K0, K5, NV,  &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (Index, N0, N5, N2, N3, Order, CD,   &
                                           M0, M5, M2, M3, Order, CD,   &
                                           K0, K5, K2, K3, Order, CD,   &
                                           NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_3D_b
!
!==========================================================================
!
  Subroutine Euler_Flux_1D_c (Index, N0, N5, N2, N3, NV,              &
                              Gamma, Q,                               &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: NV

  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: U, P, Gamma
  REALTYPE, dimension(N0:N5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(1)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO 

  call PS_Euler_Primitive_Variable (       N0, N5, NV,                  &
                                    Gamma, Q, U, P)

  call PS_Euler_Flux               (Index, N0, N5, NV,                  &
                                           Q, U, P, F, S)
  call PS_CD_Diff                  (Index, N0, N5, N2, N3, Order, CD,   &
                                           NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_1D_c
!
!==========================================================================
!
  Subroutine Euler_Flux_2D_c (Index, N0, N5, N2, N3,                  &
                                     M0, M5, M2, M3, NV,              &
                              Gamma, Q,                               &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, P, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(2)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO  

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, NV,          &
                                    Gamma, Q, U, V, P)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, NV,          &
                                           Q, U, V, P, F, S)
  call PS_CD_Diff                  (Index, N0, N5, N2, N3, Order, CD,   &
                                           M0, M5, M2, M3, Order, CD,   &
                                           NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_2D_c 
!
!==========================================================================
!
  Subroutine Euler_Flux_3D_c (Index, N0, N5, N2, N3,                  &
                                     M0, M5, M2, M3,                  &
                                     K0, K5, K2, K3, NV,              &
                              Gamma, Q,                               &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(3)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  D_F = ZERO  

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, K0, K5, NV,  &
                                    Gamma, Q, U, V, W, P)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, K0, K5, NV,  &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (Index, N0, N5, N2, N3, Order, CD,   &
                                           M0, M5, M2, M3, Order, CD,   &
                                           K0, K5, K2, K3, Order, CD,   &
                                           NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_3D_c
!
!==========================================================================
!
  Subroutine Euler_Flux_1D_d (Index, N0, N5, N2, N3, NV,              &
                              Gamma_Formulation, Gamma_Gas, Q,        &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S,   &
                              Specific_Heat)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: NV

  integer                  :: Gamma_Formulation
  REALTYPE, dimension(:,:) :: Gamma_Gas

  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: U, P, Gamma
  REALTYPE, dimension(N0:N5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(1)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  EXTERNAL    Specific_Heat
  OPTIONAL :: Specific_Heat

  D_F = ZERO 

  if (PRESENT(Specific_Heat)) then
    call    Specific_Heat          (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  else
    call PS_Specific_Heat_Gamma    (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  endif

  call PS_Euler_Primitive_Variable (       N0, N5, NV,                  &
                                    Gamma, Q, U, P)

  call PS_Euler_Flux               (Index, N0, N5, NV,                  &
                                           Q, U, P, F, S)
  call PS_CD_Diff                  (Index, N0, N5, N2, N3, Order, CD,   &
                                           NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_1D_d
!
!==========================================================================
!
  Subroutine Euler_Flux_2D_d (Index, N0, N5, N2, N3,                  &
                                     M0, M5, M2, M3, NV,              &
                              Gamma_Formulation, Gamma_Gas, Q,        &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S,   &
                              Specific_Heat)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV

  integer                  :: Gamma_Formulation
  REALTYPE, dimension(:,:) :: Gamma_Gas

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, P, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(2)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  EXTERNAL    Specific_Heat
  OPTIONAL :: Specific_Heat

  D_F = ZERO  

  if (PRESENT(Specific_Heat)) then
    call    Specific_Heat          (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  else
    call PS_Specific_Heat_Gamma    (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  endif

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, NV,          &
                                    Gamma, Q, U, V, P)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, NV,          &
                                           Q, U, V, P, F, S)
  call PS_CD_Diff                  (Index, N0, N5, N2, N3, Order, CD,   &
                                           M0, M5, M2, M3, Order, CD,   &
                                           NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_2D_d
!
!==========================================================================
!
  Subroutine Euler_Flux_3D_d (Index, N0, N5, N2, N3,                  &
                                     M0, M5, M2, M3,                  &
                                     K0, K5, K2, K3, NV,              &
                              Gamma_Formulation, Gamma_Gas, Q,        &
                              CD, D_F, alpha, D_F_Epsilon, Flux, S,   &
                              Specific_Heat)

  integer  :: Index
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV

  integer                  :: Gamma_Formulation
  REALTYPE, dimension(:,:) :: Gamma_Gas

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), OPTIONAL :: Flux

  REALTYPE, dimension(:,:)            :: CD
  REALTYPE, dimension(3)  , OPTIONAL  :: S

  REALTYPE, OPTIONAL :: D_F_Epsilon
  REALTYPE, OPTIONAL :: alpha

  EXTERNAL    Specific_Heat
  OPTIONAL :: Specific_Heat

  D_F = ZERO  

  if (PRESENT(Specific_Heat)) then
    call    Specific_Heat          (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  else
    call PS_Specific_Heat_Gamma    (Gamma_Formulation, Q, Gamma_Gas, Gamma)
  endif

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, K0, K5, NV,  &
                                    Gamma, Q, U, V, W, P)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, K0, K5, NV,  &
                                           Q, U, V, W, P, F, S)
  call PS_CD_Diff                  (Index, N0, N5, N2, N3, Order, CD,   &
                                           M0, M5, M2, M3, Order, CD,   &
                                           K0, K5, K2, K3, Order, CD,   &
                                           NV, F, D_F, alpha)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_F) <  D_F_Epsilon) D_F = ZERO 

  if (PRESENT(Flux)) Flux = F

  END Subroutine Euler_Flux_3D_d

END MODULE Central_Euler_Flux
