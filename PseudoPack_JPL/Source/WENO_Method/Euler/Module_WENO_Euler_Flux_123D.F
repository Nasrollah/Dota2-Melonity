#include "pseudopack.h"


MODULE WENO_Euler_Flux_123D

  USE Euler_Function
  USE WENO_Option
  USE WENO_Euler
  USE Backward_Differencing

IMPLICIT NONE

INTERFACE PS_WENO_Euler_Flux_1D_x
  MODULE PROCEDURE WENO_Euler_Flux_1D_x_a
  MODULE PROCEDURE WENO_Euler_Flux_1D_x_b
END INTERFACE 

INTERFACE PS_WENO_Euler_Flux_2D_x
  MODULE PROCEDURE WENO_Euler_Flux_2D_x_a
  MODULE PROCEDURE WENO_Euler_Flux_2D_x_b
END INTERFACE 

INTERFACE PS_WENO_Euler_Flux_2D_y
  MODULE PROCEDURE WENO_Euler_Flux_2D_y_a
  MODULE PROCEDURE WENO_Euler_Flux_2D_y_b
END INTERFACE 

INTERFACE PS_WENO_Euler_Flux_3D_x
  MODULE PROCEDURE WENO_Euler_Flux_3D_x_a
  MODULE PROCEDURE WENO_Euler_Flux_3D_x_b
END INTERFACE 

INTERFACE PS_WENO_Euler_Flux_3D_y
  MODULE PROCEDURE WENO_Euler_Flux_3D_y_a
  MODULE PROCEDURE WENO_Euler_Flux_3D_y_b
END INTERFACE 

INTERFACE PS_WENO_Euler_Flux_3D_z
  MODULE PROCEDURE WENO_Euler_Flux_3D_z_a
  MODULE PROCEDURE WENO_Euler_Flux_3D_z_b
END INTERFACE

PRIVATE
PUBLIC :: PS_WENO_Euler_Flux_1D_x

PUBLIC :: PS_WENO_Euler_Flux_2D_x
PUBLIC :: PS_WENO_Euler_Flux_2D_y

PUBLIC :: PS_WENO_Euler_Flux_3D_x
PUBLIC :: PS_WENO_Euler_Flux_3D_y
PUBLIC :: PS_WENO_Euler_Flux_3D_z

CONTAINS
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

  REALTYPE, dimension(:,:)   :: Gamma_Gas
  REALTYPE, dimension(N0:N5) :: Gamma

  REALTYPE, dimension(1),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5)        :: U, P, H, C
  REALTYPE, dimension(N0:N5)        :: f1
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

    if (NV == 4) f1 = Q(:,NV)/Q(:,1)

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0  , i_1, NV,             &
                                   Normal, Gamma_Gas,                         &
                                   Q(:,1), U, H, f1, E_L, E_R                 ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0  , i_1, NV,             &
                                   Q, F, Lambda_Max, E_L, E_R, F_Half         ) 

    call PS_Backward_Difference   (       N0, N5, i_0+1, i_1, NV,             &
                                   dx, F_Half, D_Flux, -ONE                   )

  END Subroutine WENO_Euler_Flux_1D_x_a
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

  REALTYPE, dimension(:,:)         :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE, dimension(2),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C
  REALTYPE, dimension(N0:N5)              :: f1
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

    if (NV == 5) f1 = Q(:,j,NV)/Q(:,j,1)

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0  , i_1, NV,             &
                                   Normal, Gamma_Gas,                         &
                                   Q(:,j,1), U(:,j), V(:,j), H(:,j), f1,      &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0  , i_1, NV,             &
                                   Q(:,j,:), F(:,j,:),                        &
                                   Lambda_Max, E_L, E_R, F_Half(:,j,:)        ) 

    call PS_Backward_Difference   (       N0, N5, i_0+1, i_1, NV,             &
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

  REALTYPE, dimension(:,:)         :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE, dimension(2),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5)        :: U, V, P, H, C
  REALTYPE, dimension(M0:M5)              :: f1
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

    if (NV == 5) f1 = Q(i,:,NV)/Q(i,:,1)

    call PS_Euler_Roe_Eigenvector (       M0, M5, j_0  , j_1, NV,             &
                                   Normal, Gamma_Gas,                         &
                                   Q(i,:,1), U(i,:), V(i,:), H(i,:), f1,      &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, M0, M5, j_0  , j_1, NV,             &
                                   Q(i,:,:), F(i,:,:),                        &
                                   Lambda_Max, E_L, E_R, F_Half(i,:,:)        ) 

    call PS_Backward_Difference   (       M0, M5, j_0+1, j_1, NV,             &
                                   dy, F_Half(i,:,:), D_Flux(i,:,:), -ONE     )
  enddo

  END Subroutine WENO_Euler_Flux_2D_y_a
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

  REALTYPE, dimension(:,:)               :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE, dimension(3),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C
  REALTYPE, dimension(N0:N5)                    :: f1
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

    if (NV == 6) f1 = Q(:,j,k,NV)/Q(:,j,k,1)

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0  , i_1, NV,             &
                                   Normal, Gamma_Gas,                         &
                                   Q(:,j,k,1), U(:,j,k), V(:,j,k), W(:,j,k),  &
                                               H(:,j,k), f1,                  &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0  , i_1, NV,             &
                                   Q(:,j,k,:), F(:,j,k,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(:,j,k,:)      ) 

    call PS_Backward_Difference   (       N0, N5, i_0+1, i_1, NV,             &
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

  REALTYPE, dimension(:,:)               :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE, dimension(3),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C
  REALTYPE, dimension(M0:M5)                    :: f1
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

    if (NV == 6) f1 = Q(i,:,k,NV)/Q(i,:,k,1)

    call PS_Euler_Roe_Eigenvector (       M0, M5, j_0  , j_1, NV,             &
                                   Normal, Gamma_Gas,                         &
                                   Q(i,:,k,1), U(i,:,k), V(i,:,k), W(i,:,k),  &
                                               H(i,:,k), f1,                  &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, M0, M5, j_0  , j_1, NV,             &
                                   Q(i,:,k,:), F(i,:,k,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(i,:,k,:)      ) 

    call PS_Backward_Difference   (       M0, M5, j_0+1, j_1, NV,             &
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

  REALTYPE, dimension(:,:)               :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE, dimension(3),  OPTIONAL :: V_Moving_Frame
  integer ,                OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)     :: F, F_Half
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)        :: U, V, W, P, H, C
  REALTYPE, dimension(K0:K5)                    :: f1
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

    if (NV == 6) f1 = Q(i,j,:,NV)/Q(i,j,:,1)

    call PS_Euler_Roe_Eigenvector (       K0, K5, k_0  , k_1, NV,             &
                                   Normal, Gamma_Gas,                         &
                                   Q(i,j,:,1), U(i,j,:), V(i,j,:), W(i,j,:),  &
                                               H(i,j,:), f1,                  &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, K0, K5, k_0  , k_1, NV,             &
                                   Q(i,j,:,:), F(i,j,:,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(i,j,:,:)      ) 

    call PS_Backward_Difference   (       K0, K5, k_0+1, k_1, NV,             &
                                   dx, F_Half(i,j,:,:), D_Flux(i,j,:,:), -ONE )
  enddo
  enddo

  END Subroutine WENO_Euler_Flux_3D_z_a
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
  REALTYPE, dimension(N0:N5)        :: f1
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

    if (NV == 4) f1 = Q(:,NV)/Q(:,1)

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0  , i_1, NV,             &
                                   Normal, Gamma, Q(:,1), U, H, f1, E_L, E_R  ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0  , i_1, NV,             &
                                   Q, F, Lambda_Max, E_L, E_R, F_Half         )

    call PS_Backward_Difference   (       N0, N5, i_0+1, i_1, NV,             &
                                   dx, F_Half, D_Flux, -ONE                   )

  END Subroutine WENO_Euler_Flux_1D_x_b
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
  REALTYPE, dimension(N0:N5)              :: f1
  REALTYPE, dimension(N0:N5,NV,NV)        :: E_L, E_R

  REALTYPE, dimension(M0:M5,NV), OPTIONAL :: Lambda_Max_p
  REALTYPE, dimension(M0:M5,NV)           :: Lambda_i
  REALTYPE, dimension(NV)                 :: Lambda_Max

  REALTYPE, dimension(2), parameter :: Normal = (/ ONE, ZERO /)
  integer                :: n, j, i_0, i_1

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

    if (NV == 5) f1 = Q(:,j,NV)/Q(:,j,1)

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0  , i_1, NV,             &
                                   Normal, Gamma,                             &
                                   Q(:,j,1), U(:,j), V(:,j), H(:,j), f1,      &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0  , i_1, NV,             &
                                   Q(:,j,:), F(:,j,:),                        &
                                   Lambda_Max, E_L, E_R, F_Half(:,j,:)        ) 

    call PS_Backward_Difference   (       N0, N5, i_0+1, i_1, NV,             &
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
  REALTYPE, dimension(M0:M5)              :: f1
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

    if (NV == 5) f1 = Q(i,:,NV)/Q(i,:,1)

    call PS_Euler_Roe_Eigenvector (       M0, M5, j_0  , j_1, NV,             &
                                   Normal, Gamma,                             &
                                   Q(i,:,1), U(i,:), V(i,:), H(i,:), f1,      &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, M0, M5, j_0  , j_1, NV,             &
                                   Q(i,:,:), F(i,:,:),                        &
                                   Lambda_Max, E_L, E_R, F_Half(i,:,:)        ) 

    call PS_Backward_Difference   (       M0, M5, j_0+1, j_1, NV,             &
                                   dy, F_Half(i,:,:), D_Flux(i,:,:), -ONE     )
  enddo

  END Subroutine WENO_Euler_Flux_2D_y_b
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
  REALTYPE, dimension(N0:N5)                    :: f1
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

    if (NV == 6) f1 = Q(:,j,k,NV)/Q(:,j,k,1)

    call PS_Euler_Roe_Eigenvector (       N0, N5, i_0  , i_1, NV,             &
                                   Normal, Gamma,                             &
                                   Q(:,j,k,1), U(:,j,k), V(:,j,k), W(:,j,k),  &
                                               H(:,j,k), f1,                  &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, N0, N5, i_0  , i_1, NV,             &
                                   Q(:,j,k,:), F(:,j,k,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(:,j,k,:)      ) 

    call PS_Backward_Difference   (       N0, N5, i_0+1, i_1, NV,             &
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
  REALTYPE, dimension(M0:M5)                    :: f1
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

    if (NV == 6) f1 = Q(i,:,k,NV)/Q(i,:,k,1)

    call PS_Euler_Roe_Eigenvector (       M0, M5, j_0  , j_1, NV,             &
                                   Normal, Gamma,                             &
                                   Q(i,:,k,1), U(i,:,k), V(i,:,k), W(i,:,k),  &
                                               H(i,:,k), f1,                  &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, M0, M5, j_0  , j_1, NV,             &
                                   Q(i,:,k,:), F(i,:,k,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(i,:,k,:)      ) 

    call PS_Backward_Difference   (       M0, M5, j_0+1, j_1, NV,             &
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
  REALTYPE, dimension(K0:K5)                    :: f1
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

    if (NV == 6) f1 = Q(i,j,:,NV)/Q(i,j,:,1)

    call PS_Euler_Roe_Eigenvector (       K0, K5, k_0  , k_1, NV,             &
                                   Normal, Gamma,                             &
                                   Q(i,j,:,1), U(i,j,:), V(i,j,:), W(i,j,:),  &
                                               H(i,j,:), f1,                  &
                                   E_L, E_R ) 

    call PS_WENO_Euler            (Order, K0, K5, k_0  , k_1, NV,             &
                                   Q(i,j,:,:), F(i,j,:,:),                    &
                                   Lambda_Max, E_L, E_R, F_Half(i,j,:,:)      ) 

    call PS_Backward_Difference   (       K0, K5, k_0+1, k_1, NV,             &
                                   dx, F_Half(i,j,:,:), D_Flux(i,j,:,:), -ONE )
  enddo
  enddo

  END Subroutine WENO_Euler_Flux_3D_z_b

END MODULE WENO_Euler_Flux_123D
