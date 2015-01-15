#include "pseudopack.h"


MODULE WENO_Euler_Flux_General

  USE Euler_Function
  USE WENO_Option
  USE WENO_Euler
  USE Backward_Differencing

IMPLICIT NONE

INTERFACE PS_WENO_Euler_Flux
  MODULE PROCEDURE WENO_Euler_Flux_1D_C_a
  MODULE PROCEDURE WENO_Euler_Flux_1D_P_a

  MODULE PROCEDURE WENO_Euler_Flux_2D_C_a
  MODULE PROCEDURE WENO_Euler_Flux_2D_P_a

  MODULE PROCEDURE WENO_Euler_Flux_3D_C_a
  MODULE PROCEDURE WENO_Euler_Flux_3D_P_a

  MODULE PROCEDURE WENO_Euler_Flux_1D_C_b
  MODULE PROCEDURE WENO_Euler_Flux_1D_P_b

  MODULE PROCEDURE WENO_Euler_Flux_2D_C_b
  MODULE PROCEDURE WENO_Euler_Flux_2D_P_b

  MODULE PROCEDURE WENO_Euler_Flux_3D_C_b
  MODULE PROCEDURE WENO_Euler_Flux_3D_P_b
END INTERFACE

PRIVATE
PUBLIC :: PS_WENO_Euler_Flux

CONTAINS
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_C_a (Index_x, Distributed_x, Order_x, dx, &
                                            N0, N5, N2, N3, NV,           &
                                     Gamma_Gas, Gamma, Q, D_Flux,         &
                                     Lambda, D_F_Epsilon,                 &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE, dimension(256,5) :: Gamma_Gas
  REALTYPE, dimension(N0:N5) :: Gamma

  REALTYPE, dimension(N0:N5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5)        :: U, P, H, C

  REALTYPE,                OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(1) , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(NV), OPTIONAL :: Lambda
  integer ,                OPTIONAL :: MPI_Comm_Type

  if (Index_x == 0) RETURN

  call PS_Euler_Primitive_Variable (       N0, N5, NV, Gamma, Q, U, P, H, C)

  call PS_Euler_Flux               (Index_x,                               &
                                           N0, N5, NV,        Q, U, P, F,  &
                                    V_Moving_Frame)

     call WENO_Euler_Flux_1D_P_a   (Index_x, Distributed_x, Order_x, dx,   &
                                           N0, N5, N2, N3, NV,             &
                                    Gamma_Gas, Gamma,                      &
                                    Q, U, P, H, C, F, D_Flux,              &
                                    Lambda, D_F_Epsilon,                   &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_1D_C_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_P_a (Index_x, Distributed_x, Order_x, dx,  &
                                            N0, N5, N2, N3, NV,            &
                                     Gamma_Gas, Gamma,                     &
                                     Q, U, P, H, C, F, D_Flux,             &
                                     Lambda, D_F_Epsilon,                  &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE, dimension(256,5) :: Gamma_Gas
  REALTYPE, dimension(N0:N5) :: Gamma

  REALTYPE, dimension(N0:N5,NV)     :: Q, F, D_Flux
  REALTYPE, dimension(N0:N5)        :: U, P, H, C

  REALTYPE,                OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(1) , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(NV), OPTIONAL :: Lambda
  integer ,                OPTIONAL :: MPI_Comm_Type

  if (Index_x == 0) RETURN

      call WENO_Euler_Flux_1D_x_a (Index_x, Distributed_x, Order_x, dx,    &
                                         N0, N5, N2, N3,  NV,              &
                                   Gamma_Gas, Gamma,                       &
                                   Q, U, P, H, C, F, D_Flux,               &
                                   Lambda, V_Moving_Frame, MPI_Comm_Type)

    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO

  END Subroutine WENO_Euler_Flux_1D_P_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_x_a  (Index, Distributed, Order, dx,    &
                                             N0, N5, N2, N3, NV,        &
                                      Gamma_Gas, Gamma,                 &
                                      Q, U, P, H, C, F, D_Flux,         &
                                      Lambda_Max_p,                     &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE, dimension(256,5) :: Gamma_Gas
  REALTYPE, dimension(N0:N5) :: Gamma

  REALTYPE, dimension(1),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5)        :: U, P, H, C, f1
  REALTYPE, dimension(N0:N5,NV,NV)  :: E_L, E_R

  REALTYPE, dimension(NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(NV)           :: Lambda_i
  REALTYPE, dimension(NV)           :: Lambda_Max

  integer  :: n, i, i_0, i_1
  REALTYPE, dimension(1), parameter :: Normal = (/ ONE /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  i_0 = N2-1 ; i_1 = N3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, NV,                      &
                                       U, C, Lambda_i, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  f1 = Q(:,4)/Q(:,1)

    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_i    (:)
    else
      Lambda_Max = Lambda_Max_p(:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0, i_1, NV,               &
                                   Normal, Gamma_Gas,                         &
                                   Q(:,1), U, H, f1, E_R, E_L                 ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0, i_1, NV,               &
                                   Q, F, Lambda_Max, E_L, E_R, F_Half         ) 

    call PS_Backward_Difference   (N0, N5, i_0, i_1, NV,                      &
                                   dx, F_Half, D_Flux, -ONE                   )

  END Subroutine WENO_Euler_Flux_1D_x_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_C_a (Index,                               &
                                     Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,             &
                                     Index_y, Distributed_y, Order_y, dy, &
                                              M0, M5, M2, M3,  NV,        &
                                     Gamma_Gas, Gamma, Q, D_Flux,         &
                                     Lambda, D_F_Epsilon,                 &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: Index_y, Order_y
  logical  :: Distributed_y
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: NV

  REALTYPE, dimension(256,5)       :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE,                  OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(2)   , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(:,:) , OPTIONAL :: Lambda
  integer ,                  OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, NV,             &
                                    Gamma, Q, U, V, P, H, C)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, NV,             &
                                           Q, U, V, P, F,                  &
                                           V_Moving_Frame)

  call WENO_Euler_Flux_2D_P_a      (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                    Index_y, Distributed_y, Order_y, dy, &
                                           M0, M5, M2, M3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, P, H, C, F, D_Flux,         &
                                    Lambda, D_F_Epsilon,                 &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_2D_C_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_P_a (Index_x, Distributed_x, Order_x, dx, &
                                            N0, N5, N2, N3,               &
                                     Index_y, Distributed_y, Order_y, dy, &
                                            M0, M5, M2, M3,  NV,          &
                                     Gamma_Gas, Gamma,                    &
                                     Q, U, V, P, H, C, F, D_Flux,         &
                                     Lambda, D_F_Epsilon,                 &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: Index_y, Order_y
  logical  :: Distributed_y
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: NV

  REALTYPE, dimension(256,5)       :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(2)  , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(:,:), OPTIONAL :: Lambda
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux, D_G
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F

  SELECT CASE (Index_x)
    CASE (1)
      call WENO_Euler_Flux_2D_x_a (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, P, H, C, F, D_Flux,              &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_2D_y_a (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, P, H, C, F, D_Flux,              &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)
  END SELECT 
  
  if (Index_x /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO
  else 
    D_Flux = ZERO
  endif

  SELECT CASE (Index_y)
    CASE (1)
      call WENO_Euler_Flux_2D_x_a (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, P, H, C, F, D_G,                 &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_2D_y_a (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, P, H, C, F, D_G,                 &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

  END SELECT

  if (Index_y /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO

    D_Flux = D_Flux + D_G
  endif

  END Subroutine WENO_Euler_Flux_2D_P_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_x_a (Index, Distributed, Order, dx,       &
                                            N0, N5, N2, N3,               &
                                            M0, M5, M2, M3,  NV,          &
                                     Gamma_Gas, Gamma,                    &
                                     Q, U, V, P, H, C, F, D_Flux,         &
                                     Lambda_Max_p,                        &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE, dimension(256,5)       :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE, dimension(2),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C, f1
  REALTYPE, dimension(N0:N5,NV,NV)        :: E_L, E_R

  REALTYPE, dimension(M0:M5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(M0:M5,NV)           :: Lambda_i
  REALTYPE, dimension(NV)                 :: Lambda_Max

  integer  :: n, j, i_0, i_1
  REALTYPE, dimension(2), parameter :: Normal = (/ ONE, ZERO /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  i_0 = N2-1 ; i_1 = N3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, NV,                 &
                                       U, V, C, Lambda_i, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  f1 = Q(:,:,5)/Q(:,:,1)

  do j = M2,M3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_i    (j,:)
    else
      Lambda_Max = Lambda_Max_p(j,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0, i_1, NV,               &
                                   Normal, Gamma_Gas,                         &
                                   Q(:,j,1), U(:,j), V(:,j), h(:,j), f1(:,j), &
                                   E_R, E_L ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0, i_1, NV,               &
                                   Q(:,j,:), F(:,j,:),                        &
                                   Lambda_Max, E_L, E_R, F_Half(:,j,:)        ) 

    call PS_Backward_Difference   (       N0, N5, i_0, i_1, NV,               &
                                   dx, F_Half(:,j,:), D_Flux(:,j,:), -ONE     )
  enddo

  END Subroutine WENO_Euler_Flux_2D_x_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_y_a  (Index, Distributed, Order, dy,       &
                                             N0, N5, N2, N3,               &
                                             M0, M5, M2, M3,  NV,          &
                                      Gamma_Gas, Gamma,                    &
                                      Q, U, V, P, H, C, F, D_Flux,         &
                                      Lambda_Max_p,                        &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV
  REALTYPE :: dy

  REALTYPE, dimension(256,5)       :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE, dimension(2),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C, f1
  REALTYPE, dimension(M0:M5,NV,NV)        :: E_L, E_R

  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(N0:N5,NV)           :: Lambda_j
  REALTYPE, dimension(NV)                 :: Lambda_Max

  integer  :: n, i, j_0, j_1
  REALTYPE, dimension(2), parameter :: Normal = (/ ZERO, ONE /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  j_0 = M2-1 ; j_1 = M3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, NV,                 &
                                       U, V, C, Lambda_j, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  f1 = Q(:,:,5)/Q(:,:,1)

  do i = N2,N3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_j    (i,:)
    else
      Lambda_Max = Lambda_Max_p(i,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       M0, M5, j_0, j_1, NV,               &
                                   Normal, Gamma_Gas,                         &
                                   Q(i,:,1), U(i,:), V(i,:), h(i,:), f1(i,:), &
                                   E_R, E_L ) 

    call PS_WENO_Euler            (Order, M0, M5, j_0, j_1, NV,               &
                                   Q(i,:,:), F(i,:,:),                        &
                                   Lambda_Max, E_L, E_R, F_Half(i,:,:)        ) 

    call PS_Backward_Difference   (       M0, M5, j_0, j_1, NV,               &
                                   dy, F_Half(i,:,:), D_Flux(i,:,:), -ONE     )
  enddo

  END Subroutine WENO_Euler_Flux_2D_y_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_C_a (Index,                               &
                                     Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,             &
                                     Index_y, Distributed_y, Order_y, dy, &
                                              M0, M5, M2, M3,             &
                                     Index_z, Distributed_z, Order_z, dz, &
                                              K0, K5, K2, K3,  NV,        &
                                     Gamma_Gas, Gamma, Q, D_Flux,         &
                                     Lambda, D_F_Epsilon,                 &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: Index_y, Order_y
  logical  :: Distributed_y
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: Index_z, Order_z
  logical  :: Distributed_z
  integer  :: K0, K5, K2, K3
  REALTYPE :: dz

  integer  :: NV

  REALTYPE, dimension(256,5)             :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE,                   OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(3)    , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(:,:,:), OPTIONAL :: Lambda
  integer ,                   OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, K0, K5, NV,      &
                                    Gamma, Q, U, V, W, P, H, C)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, K0, K5, NV,      &
                                           Q, U, V, W, P, F,                &
                                           V_Moving_Frame)

  call WENO_Euler_Flux_3D_P_a      (Index_x, Distributed_x, Order_x, dx,    &
                                           N0, N5, N2, N3,                  &
                                    Index_y, Distributed_y, Order_y, dy,    &
                                           M0, M5, M2, M3,                  &
                                    Index_z, Distributed_z, Order_z, dz,    &
                                           K0, K5, K2, K3,  NV,             &
                                    Gamma_Gas, Gamma,                       &
                                    Q, U, V, W, P, H, C, F, D_Flux,         &
                                    Lambda, D_F_Epsilon,                    &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_3D_C_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_P_a (Index_x, Distributed_x, Order_x, dx,    &
                                            N0, N5, N2, N3,                  &
                                     Index_y, Distributed_y, Order_y, dy,    &
                                            M0, M5, M2, M3,                  &
                                     Index_z, Distributed_z, Order_z, dz,    &
                                            K0, K5, K2, K3,  NV,             &
                                     Gamma_Gas, Gamma,                       &
                                     Q, U, V, W, P, H, C, F, D_Flux,         &
                                     Lambda, D_F_Epsilon,                    &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: Index_y, Order_y
  logical  :: Distributed_y
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: Index_z, Order_z
  logical  :: Distributed_z
  integer  :: K0, K5, K2, K3
  REALTYPE :: dz

  integer  :: NV

  REALTYPE, dimension(256,5)             :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE,                   OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(3)    , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(:,:,:), OPTIONAL :: Lambda
  integer ,                   OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux, D_G
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F

  SELECT CASE (Index_x)
    CASE (1)
      call WENO_Euler_Flux_3D_x_a (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, W, P, H, C, F, D_Flux,           &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_3D_y_a (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, W, P, H, C, F, D_Flux,           &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (3)
      call WENO_Euler_Flux_3D_z_a (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, W, P, H, C, F, D_Flux,           &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

  END SELECT 
  
  if (Index_x /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO
  else 
    D_Flux = ZERO
  endif

  SELECT CASE (Index_y)
    CASE (1)
      call WENO_Euler_Flux_3D_x_a (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, W, P, H, C, F, D_G,              &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_3D_y_a (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, W, P, H, C, F, D_G,              &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (3)
      call WENO_Euler_Flux_3D_z_a (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, W, P, H, C, F, D_G,              &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

  END SELECT

  if (Index_y /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO

    D_Flux = D_Flux + D_G
  endif

  SELECT CASE (Index_z)
    CASE (1)
      call WENO_Euler_Flux_3D_x_a (Index_z, Distributed_z, Order_z, dz,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, W, P, H, C, F, D_G,              &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_3D_y_a (Index_z, Distributed_z, Order_z, dz,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, W, P, H, C, F, D_G,              &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (3)
      call WENO_Euler_Flux_3D_z_a (Index_z, Distributed_z, Order_z, dz,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma_Gas, Gamma,                         &
                                   Q, U, V, W, P, H, C, F, D_G,              &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

  END SELECT

  if (Index_z /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO

    D_Flux = D_Flux + D_G
  endif

  END Subroutine WENO_Euler_Flux_3D_P_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_x_a  (Index, Distributed, Order, dx,          &
                                             N0, N5, N2, N3,                  &
                                             M0, M5, M2, M3,                  &
                                             K0, K5, K2, K3,  NV,             &
                                      Gamma_Gas, Gamma,                       &
                                      Q, U, V, W, P, H, C, F, D_Flux,         &
                                      Lambda_Max_p,                           &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE, dimension(256,5)             :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE, dimension(3),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C, f1
  REALTYPE, dimension(N0:N5,NV,NV)              :: E_L, E_R

  REALTYPE, dimension(M0:M5,K0:K5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(M0:M5,K0:K5,NV)           :: Lambda_i
  REALTYPE, dimension(NV)                       :: Lambda_Max

  integer  :: n, j, k, i_0, i_1
  REALTYPE, dimension(3), parameter :: Normal = (/ ONE, ZERO, ZERO /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  i_0 = N2-1 ; i_1 = N3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, K0, K5, NV,            &
                                       U, V, W, C, Lambda_i, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  f1 = Q(:,:,:,6)/Q(:,:,:,1)

  do k = K2,K3
  do j = M2,M3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_i    (j,k,:)
    else
      Lambda_Max = Lambda_Max_p(j,k,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0, i_1, NV,               &
                                   Normal, Gamma_Gas,                         &
                                   Q(:,j,k,1), U(:,j,k), V(:,j,k), W(:,j,k),  &
                                               h(:,j,k),f1(:,j,k),            &
                                   E_R, E_L ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0, i_1, NV,               &
                                   Q(:,j,k,:), F(:,j,k,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(:,j,k,:)      ) 

    call PS_Backward_Difference   (       N0, N5, i_0, i_1, NV,               &
                                   dx, F_Half(:,j,k,:), D_Flux(:,j,k,:), -ONE )
  enddo
  enddo

  END Subroutine WENO_Euler_Flux_3D_x_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_y_a  (Index, Distributed, Order, dx,          &
                                             N0, N5, N2, N3,                  &
                                             M0, M5, M2, M3,                  &
                                             K0, K5, K2, K3,  NV,             &
                                      Gamma_Gas, Gamma,                       &
                                      Q, U, V, W, P, H, C, F, D_Flux,         &
                                      Lambda_Max_p,                           &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE, dimension(256,5)             :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE, dimension(3),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C, f1
  REALTYPE, dimension(M0:M5,NV,NV)              :: E_L, E_R

  REALTYPE, dimension(N0:N5,K0:K5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(N0:N5,K0:K5,NV)           :: Lambda_j
  REALTYPE, dimension(NV)                       :: Lambda_Max

  integer  :: n, i, k, j_0, j_1
  REALTYPE, dimension(3), parameter :: Normal = (/ ZERO, ONE, ZERO /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  j_0 = M2-1 ; j_1 = M3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, K0, K5, NV,            &
                                       U, V, W, C, Lambda_j, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  f1 = Q(:,:,:,6)/Q(:,:,:,1)

  do k = K2,K3
  do i = N2,N3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_j    (i,k,:)
    else
      Lambda_Max = Lambda_Max_p(i,k,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       M0, M5, j_0, j_1, NV,               &
                                   Normal, Gamma_Gas,                         &
                                   Q(i,:,k,1), U(i,:,k), V(i,:,k), W(i,:,k),  &
                                               h(i,:,k),f1(i,:,k),            &
                                   E_R, E_L ) 

    call PS_WENO_Euler            (Order, M0, M5, j_0, j_1, NV,               &
                                   Q(i,:,k,:), F(i,:,k,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(i,:,k,:)      ) 

    call PS_Backward_Difference   (       M0, M5, j_0, j_1, NV,               &
                                   dx, F_Half(i,:,k,:), D_Flux(i,:,k,:), -ONE )
  enddo
  enddo

  END Subroutine WENO_Euler_Flux_3D_y_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_z_a  (Index, Distributed, Order, dx,          &
                                             N0, N5, N2, N3,                  &
                                             M0, M5, M2, M3,                  &
                                             K0, K5, K2, K3,  NV,             &
                                      Gamma_Gas, Gamma,                       &
                                      Q, U, V, W, P, H, C, F, D_Flux,         &
                                      Lambda_Max_p,                           &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE, dimension(256,5)             :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE, dimension(3),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C, f1
  REALTYPE, dimension(K0:K5,NV,NV)              :: E_L, E_R

  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Lambda_k
  REALTYPE, dimension(NV)                       :: Lambda_Max

  integer  :: n, i, j, k_0, k_1
  REALTYPE, dimension(3), parameter :: Normal = (/ ZERO, ZERO, ONE /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  k_0 = K2-1 ; k_1 = K3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, K0, K5, NV,            &
                                       U, V, W, C, Lambda_k, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  f1 = Q(:,:,:,6)/Q(:,:,:,1)

  do j = M2,M3
  do i = N2,N3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_k    (i,j,:)
    else
      Lambda_Max = Lambda_Max_p(i,j,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       K0, K5, k_0, k_1, NV,               &
                                   Normal, Gamma_Gas,                         &
                                   Q(i,j,:,1), U(i,j,:), V(i,j,:), W(i,j,:),  &
                                               h(i,j,:),f1(i,j,:),            &
                                   E_R, E_L ) 

    call PS_WENO_Euler            (Order, K0, K5, k_0, k_1, NV,               &
                                   Q(i,j,:,:), F(i,j,:,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(i,j,:,:)      ) 

    call PS_Backward_Difference   (       K0, K5, k_0, k_1, NV,               &
                                   dx, F_Half(i,j,:,:), D_Flux(i,j,:,:), -ONE )
  enddo
  enddo

  END Subroutine WENO_Euler_Flux_3D_z_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_C_b (Index_x, Distributed_x, Order_x, dx,     &
                                            N0, N5, N2, N3, NV,               &
                                     Gamma, Q, D_Flux,                        &
                                     Lambda, D_F_Epsilon,                     &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE                   :: Gamma

  REALTYPE, dimension(N0:N5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5)        :: U, P, H, C

  REALTYPE,                OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(1) , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(NV), OPTIONAL :: Lambda
  integer ,                OPTIONAL :: MPI_Comm_Type

  if (Index_x == 0) RETURN

  call PS_Euler_Primitive_Variable (       N0, N5, NV, Gamma, Q, U, P, H, C)

  call PS_Euler_Flux               (Index_x,                               &
                                           N0, N5, NV,        Q, U, P, F,  &
                                    V_Moving_Frame)

  call    WENO_Euler_Flux_1D_P_b   (Index_x, Distributed_x, Order_x, dx,   &
                                           N0, N5, N2, N3, NV,             &
                                    Gamma, Q, U, P, H, C, F, D_Flux,       &
                                    Lambda, D_F_Epsilon,                   &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_1D_C_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_P_b (Index_x, Distributed_x, Order_x, dx,     &
                                            N0, N5, N2, N3, NV,               &
                                     Gamma, Q, U, P, H, C, F, D_Flux,         &
                                     Lambda, D_F_Epsilon,                     &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE                   :: Gamma

  REALTYPE, dimension(N0:N5,NV)     :: Q, F, D_Flux
  REALTYPE, dimension(N0:N5)        :: U, P, H, C

  REALTYPE,                OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(1) , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(NV), OPTIONAL :: Lambda
  integer ,                OPTIONAL :: MPI_Comm_Type

  if (Index_x == 0) RETURN

      call WENO_Euler_Flux_1D_x_b (Index_x, Distributed_x, Order_x, dx,    &
                                         N0, N5, N2, N3,  NV,              &
                                   Gamma, Q, U, P, H, C, F, D_Flux,        &
                                   Lambda, V_Moving_Frame, MPI_Comm_Type)

    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO

  END Subroutine WENO_Euler_Flux_1D_P_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_x_b  (Index, Distributed, Order, dx,    &
                                             N0, N5, N2, N3, NV,        &
                                      Gamma, Q, U, P, H, C, F, D_Flux,  &
                                      Lambda_Max_p,                     &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE                   :: Gamma

  REALTYPE, dimension(1),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5)        :: U, P, H, C
  REALTYPE, dimension(N0:N5,NV,NV)  :: E_L, E_R

  REALTYPE, dimension(NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(NV)           :: Lambda_i
  REALTYPE, dimension(NV)           :: Lambda_Max

  integer  :: n, i, i_0, i_1
  REALTYPE, dimension(1), parameter :: Normal = (/ ONE /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  i_0 = N2-1 ; i_1 = N3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, NV,                      &
                                       U, C, Lambda_i, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_i    (:)
    else
      Lambda_Max = Lambda_Max_p(:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0, i_1, NV,               &
                                   Normal, Gamma, Q(:,1), U, H, E_R, E_L      ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0, i_1,                   &
                                   NV, Q, F, Lambda_Max, E_L, E_R, F_Half     ) 

    call PS_Backward_Difference   (       N0, N5, i_0, i_1, NV,               &
                                   dx, F_Half, D_Flux, -ONE                   )

  END Subroutine WENO_Euler_Flux_1D_x_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_C_b (Index,                               &
                                     Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,             &
                                     Index_y, Distributed_y, Order_y, dy, &
                                              M0, M5, M2, M3,  NV,        &
                                     Gamma, Q, D_Flux,                    &
                                     Lambda, D_F_Epsilon,                 &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: Index_y, Order_y
  logical  :: Distributed_y
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: NV

  REALTYPE                         :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(2)  , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(:,:), OPTIONAL :: Lambda
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, NV,             &
                                    Gamma, Q, U, V, P, H, C)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, NV,             &
                                           Q, U, V, P, F,                  &
                                           V_Moving_Frame)

  call WENO_Euler_Flux_2D_P_b      (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                    Index_y, Distributed_y, Order_y, dy, &
                                           M0, M5, M2, M3,  NV,          &
                                    Gamma, Q, U, V, P, H, C, F, D_Flux,  &
                                    Lambda, D_F_Epsilon,                 &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_2D_C_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_P_b (Index_x, Distributed_x, Order_x, dx, &
                                            N0, N5, N2, N3,               &
                                     Index_y, Distributed_y, Order_y, dy, &
                                            M0, M5, M2, M3,  NV,          &
                                     Gamma, Q, U, V, P, H, C, F, D_Flux,  &
                                     Lambda, D_F_Epsilon,                 &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: Index_y, Order_y
  logical  :: Distributed_y
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: NV

  REALTYPE                         :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(2)  , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(:,:), OPTIONAL :: Lambda
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux, D_G
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F

  SELECT CASE (Index_x)
    CASE (1)
      call WENO_Euler_Flux_2D_x_b (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,  NV,               &
                                   Gamma, Q, U, V, P, H, C, F, D_Flux,       &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_2D_y_b (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,  NV,               &
                                   Gamma, Q, U, V, P, H, C, F, D_Flux,       &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)
  END SELECT 
  
  if (Index_x /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO
  endif

  SELECT CASE (Index_y)
    CASE (1)
      call WENO_Euler_Flux_2D_x_b (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,  NV,               &
                                   Gamma, Q, U, V, P, H, C, F, D_G,          &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_2D_y_b (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,  NV,               &
                                   Gamma, Q, U, V, P, H, C, F, D_G,          &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

  END SELECT

  if (Index_y /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO

    D_Flux = D_Flux + D_G
  endif

  END Subroutine WENO_Euler_Flux_2D_P_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_x_b  (Index, Distributed, Order, dx,       &
                                             N0, N5, N2, N3,               &
                                             M0, M5, M2, M3,  NV,          &
                                      Gamma, Q, U, V, P, H, C, F, D_Flux,  &
                                      Lambda_Max_p,                        &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE                         :: Gamma

  REALTYPE, dimension(2),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C
  REALTYPE, dimension(N0:N5,NV,NV)        :: E_L, E_R

  REALTYPE, dimension(M0:M5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(M0:M5,NV)           :: Lambda_i
  REALTYPE, dimension(NV)                 :: Lambda_Max

  integer  :: n, j, i_0, i_1
  REALTYPE, dimension(2), parameter :: Normal = (/ ONE, ZERO /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  i_0 = N2-1 ; i_1 = N3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, NV,                 &
                                       U, V, C, Lambda_i, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  do j = M2,M3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_i    (j,:)
    else
      Lambda_Max = Lambda_Max_p(j,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0, i_1, NV,               &
                                   Normal, Gamma,                             &
                                   Q(:,j,1), U(:,j), V(:,j), h(:,j), E_R, E_L ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0, i_1, NV,               &
                                   Q(:,j,:), F(:,j,:),                        &
                                   Lambda_Max, E_L, E_R, F_Half(:,j,:)        ) 

    call PS_Backward_Difference   (       N0, N5, i_0, i_1, NV,               &
                                   dx, F_Half(:,j,:), D_Flux(:,j,:), -ONE     )
  enddo

  END Subroutine WENO_Euler_Flux_2D_x_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_y_b  (Index, Distributed, Order, dy,       &
                                             N0, N5, N2, N3,               &
                                             M0, M5, M2, M3,  NV,          &
                                      Gamma, Q, U, V, P, H, C, F, D_Flux,  &
                                      Lambda_Max_p,                        &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV
  REALTYPE :: dy

  REALTYPE                         :: Gamma

  REALTYPE, dimension(2),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C
  REALTYPE, dimension(M0:M5,NV,NV)        :: E_L, E_R

  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(N0:N5,NV)           :: Lambda_j
  REALTYPE, dimension(NV)                 :: Lambda_Max

  integer  :: n, i, j_0, j_1
  REALTYPE, dimension(2), parameter :: Normal = (/ ZERO, ONE /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  j_0 = M2-1 ; j_1 = M3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, NV,                 &
                                       U, V, C, Lambda_j, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  do i = N2,N3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_j    (i,:)
    else
      Lambda_Max = Lambda_Max_p(i,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       M0, M5, j_0, j_1, NV,               &
                                   Normal, Gamma,                             &
                                   Q(i,:,1), U(i,:), V(i,:), h(i,:), E_R, E_L ) 

    call PS_WENO_Euler            (Order, M0, M5, j_0, j_1, NV,               &
                                   Q(i,:,:), F(i,:,:),                        &
                                   Lambda_Max, E_L, E_R, F_Half(i,:,:)        ) 

    call PS_Backward_Difference   (       M0, M5, j_0, j_1, NV,               &
                                   dy, F_Half(i,:,:), D_Flux(i,:,:), -ONE     )
  enddo

  END Subroutine WENO_Euler_Flux_2D_y_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_C_b (Index,                               &
                                     Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,             &
                                     Index_y, Distributed_y, Order_y, dy, &
                                              M0, M5, M2, M3,             &
                                     Index_z, Distributed_z, Order_z, dz, &
                                              K0, K5, K2, K3,  NV,        &
                                     Gamma, Q, D_Flux,                    &
                                     Lambda, D_F_Epsilon,                 &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: Index_y, Order_y
  logical  :: Distributed_y
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: Index_z, Order_z
  logical  :: Distributed_z
  integer  :: K0, K5, K2, K3
  REALTYPE :: dz

  integer  :: NV

  REALTYPE                               :: Gamma

  REALTYPE,                   OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(3)    , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(:,:,:), OPTIONAL :: Lambda
  integer ,                   OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, K0, K5, NV,      &
                                    Gamma, Q, U, V, W, P, H, C)

  call PS_Euler_Flux               (Index, N0, N5, M0, M5, K0, K5, NV,      &
                                           Q, U, V, W, P, F,                &
                                           V_Moving_Frame)

  call WENO_Euler_Flux_3D_P_b      (Index_x, Distributed_x, Order_x, dx,    &
                                           N0, N5, N2, N3,                  &
                                    Index_y, Distributed_y, Order_y, dy,    &
                                           M0, M5, M2, M3,                  &
                                    Index_z, Distributed_z, Order_z, dz,    &
                                           K0, K5, K2, K3,  NV,             &
                                    Gamma, Q, U, V, W, P, H, C, F, D_Flux,  &
                                    Lambda, D_F_Epsilon,                    &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_3D_C_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_P_b (Index_x, Distributed_x, Order_x, dx,    &
                                            N0, N5, N2, N3,                  &
                                     Index_y, Distributed_y, Order_y, dy,    &
                                            M0, M5, M2, M3,                  &
                                     Index_z, Distributed_z, Order_z, dz,    &
                                            K0, K5, K2, K3,  NV,             &
                                     Gamma, Q, U, V, W, P, H, C, F, D_Flux,  &
                                     Lambda, D_F_Epsilon,                    &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: Index_y, Order_y
  logical  :: Distributed_y
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: Index_z, Order_z
  logical  :: Distributed_z
  integer  :: K0, K5, K2, K3
  REALTYPE :: dz

  integer  :: NV

  REALTYPE                               :: Gamma

  REALTYPE,                   OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(3)    , OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(:,:,:), OPTIONAL :: Lambda
  integer ,                   OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux, D_G
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F

  SELECT CASE (Index_x)
    CASE (1)
      call WENO_Euler_Flux_3D_x_b (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma, Q, U, V, W, P, H, C, F, D_Flux,    &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_3D_y_b (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma, Q, U, V, W, P, H, C, F, D_Flux,    &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (3)
      call WENO_Euler_Flux_3D_z_b (Index_x, Distributed_x, Order_x, dx,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma, Q, U, V, W, P, H, C, F, D_Flux,    &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

  END SELECT 
  
  if (Index_x /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO
  endif

  SELECT CASE (Index_y)
    CASE (1)
      call WENO_Euler_Flux_3D_x_b (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma, Q, U, V, W, P, H, C, F, D_G,       &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_3D_y_b (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma, Q, U, V, W, P, H, C, F, D_G,       &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (3)
      call WENO_Euler_Flux_3D_z_b (Index_y, Distributed_y, Order_y, dy,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma, Q, U, V, W, P, H, C, F, D_G,       &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

  END SELECT

  if (Index_y /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G)    <  D_F_Epsilon) D_G    = ZERO

    D_Flux = D_Flux + D_G
  endif

  SELECT CASE (Index_z)
    CASE (1)
      call WENO_Euler_Flux_3D_x_b (Index_z, Distributed_z, Order_z, dz,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma, Q, U, V, W, P, H, C, F, D_G,       &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (2)
      call WENO_Euler_Flux_3D_y_b (Index_z, Distributed_z, Order_z, dz,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma, Q, U, V, W, P, H, C, F, D_G,       &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

    CASE (3)
      call WENO_Euler_Flux_3D_z_b (Index_z, Distributed_z, Order_z, dz,      &
                                          N0, N5, N2, N3,                    &
                                          M0, M5, M2, M3,                    &
                                          K0, K5, K2, K3,  NV,               &
                                   Gamma, Q, U, V, W, P, H, C, F, D_G,       &
                                   Lambda, V_Moving_Frame,                   &
                                   MPI_Comm_Type)

  END SELECT

  if (Index_z /= 0) then
    if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_G)    <  D_F_Epsilon) D_G    = ZERO

    D_Flux = D_Flux + D_G
  endif

  END Subroutine WENO_Euler_Flux_3D_P_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_x_b  (Index, Distributed, Order, dx,          &
                                             N0, N5, N2, N3,                  &
                                             M0, M5, M2, M3,                  &
                                             K0, K5, K2, K3,  NV,             &
                                      Gamma, Q, U, V, W, P, H, C, F, D_Flux,  &
                                      Lambda_Max_p,                           &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE                               :: Gamma

  REALTYPE, dimension(3),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C
  REALTYPE, dimension(N0:N5,NV,NV)              :: E_L, E_R

  REALTYPE, dimension(M0:M5,K0:K5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(M0:M5,K0:K5,NV)           :: Lambda_i
  REALTYPE, dimension(NV)                       :: Lambda_Max

  integer  :: n, j, k, i_0, i_1
  REALTYPE, dimension(3), parameter :: Normal = (/ ONE, ZERO, ZERO /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  i_0 = N2-1 ; i_1 = N3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, K0, K5, NV,            &
                                       U, V, W, C, Lambda_i, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  do k = K2,K3
  do j = M2,M3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_i    (j,k,:)
    else
      Lambda_Max = Lambda_Max_p(j,k,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0, i_1, NV,               &
                                   Normal, Gamma,                             &
                                   Q(:,j,k,1), U(:,j,k), V(:,j,k), W(:,j,k),  &
                                               h(:,j,k),                      &
                                   E_R, E_L ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0, i_1, NV,               &
                                   Q(:,j,k,:), F(:,j,k,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(:,j,k,:)      ) 

    call PS_Backward_Difference   (       N0, N5, i_0, i_1, NV,               &
                                   dx, F_Half(:,j,k,:), D_Flux(:,j,k,:), -ONE )
  enddo
  enddo

  END Subroutine WENO_Euler_Flux_3D_x_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_y_b  (Index, Distributed, Order, dx,          &
                                             N0, N5, N2, N3,                  &
                                             M0, M5, M2, M3,                  &
                                             K0, K5, K2, K3,  NV,             &
                                      Gamma, Q, U, V, W, P, H, C, F, D_Flux,  &
                                      Lambda_Max_p,                           &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE                               :: Gamma

  REALTYPE, dimension(3),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C
  REALTYPE, dimension(M0:M5,NV,NV)              :: E_L, E_R

  REALTYPE, dimension(N0:N5,K0:K5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(N0:N5,K0:K5,NV)           :: Lambda_j
  REALTYPE, dimension(NV)                       :: Lambda_Max

  integer  :: n, i, k, j_0, j_1
  REALTYPE, dimension(3), parameter :: Normal = (/ ZERO, ONE, ZERO /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  j_0 = M2-1 ; j_1 = M3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, K0, K5, NV,            &
                                       U, V, W, C, Lambda_j, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  do k = K2,K3
  do i = N2,N3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_j    (i,k,:)
    else
      Lambda_Max = Lambda_Max_p(i,k,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (        M0, M5, j_0, j_1, NV,              &
                                   Normal, Gamma,                             &
                                   Q(i,:,k,1), U(i,:,k), V(i,:,k), W(i,:,k),  &
                                               h(i,:,k),                      &
                                   E_R, E_L ) 

    call PS_WENO_Euler            (Order, M0, M5, j_0, j_1, NV,               &
                                   Q(i,:,k,:), F(i,:,k,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(i,:,k,:)      ) 

    call PS_Backward_Difference   (       M0, M5, j_0, j_1, NV,               &
                                   dx, F_Half(i,:,k,:), D_Flux(i,:,k,:), -ONE )
  enddo
  enddo

  END Subroutine WENO_Euler_Flux_3D_y_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_z_b  (Index, Distributed, Order, dx,          &
                                             N0, N5, N2, N3,                  &
                                             M0, M5, M2, M3,                  &
                                             K0, K5, K2, K3,  NV,             &
                                      Gamma, Q, U, V, W, P, H, C, F, D_Flux,  &
                                      Lambda_Max_p,                           &
                                      V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index, Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer  :: NV
  REALTYPE :: dx

  REALTYPE                               :: Gamma

  REALTYPE, dimension(3),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C
  REALTYPE, dimension(K0:K5,NV,NV)              :: E_L, E_R

  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Lambda_k
  REALTYPE, dimension(NV)                       :: Lambda_Max

  integer  :: n, i, j, k_0, k_1
  REALTYPE, dimension(3), parameter :: Normal = (/ ZERO, ZERO, ONE /)

  REALTYPE :: Lambda_Factor, Sound_Speed_Min 

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  k_0 = K2-1 ; k_1 = K3

  if (.NOT. PRESENT(Lambda_Max_p)) then
    call PS_Maximum_Eigenvalue (Index, N0, N5, M0, M5, K0, K5, NV,            &
                                       U, V, W, C, Lambda_k, V_Moving_Frame,  &
                                MPI_Comm_Type)
  endif

  do j = M2,M3
  do i = N2,N3
    if (.NOT. PRESENT(Lambda_Max_p)) then
      Lambda_Max = Lambda_k    (i,j,:)
    else
      Lambda_Max = Lambda_Max_p(i,j,:)
    endif

    Lambda_Max = MAX(Sound_Speed_Min, Lambda_Max)*Lambda_Factor

    call PS_Euler_Roe_Eigenvector (       K0, K5, k_0, k_1, NV,               &
                                   Normal, Gamma,                             &
                                   Q(i,j,:,1), U(i,j,:), V(i,j,:), W(i,j,:),  &
                                               h(i,j,:),                      &
                                   E_R, E_L ) 

    call PS_WENO_Euler            (Order, K0, K5, k_0, k_1, NV,               &
                                   Q(i,j,:,:), F(i,j,:,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(i,j,:,:)      ) 

    call PS_Backward_Difference   (       K0, K5, k_0, k_1, NV,               &
                                   dx, F_Half(i,j,:,:), D_Flux(i,j,:,:), -ONE )
  enddo
  enddo

  END Subroutine WENO_Euler_Flux_3D_z_b

END MODULE WENO_Euler_Flux_General
