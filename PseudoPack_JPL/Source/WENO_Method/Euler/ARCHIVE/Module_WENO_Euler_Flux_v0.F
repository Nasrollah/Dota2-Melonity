#include "pseudopack.h"

#ifdef NO_SIZE_IN_DIMENSION
#define _SIZE(n)	n+2
#else
#define _SIZE(n)	SIZE(Q,DIM=n+1)
#endif

MODULE WENO_Euler_Flux

#if defined (PARALLEL_MPI)
  USE Processor
#endif
  USE Euler_Function
  USE WENO_Euler
  USE WENO_Option

implicit NONE

INTERFACE PS_WENO_Euler_Flux
  MODULE PROCEDURE PS_WENO_Euler_Flux_1D
  MODULE PROCEDURE PS_WENO_Euler_Flux_2D
  MODULE PROCEDURE PS_WENO_Euler_Flux_3D
END INTERFACE

REALTYPE :: Lambda_Factor   = ELEVEN/TEN
REALTYPE :: Sound_Speed_Min = TEN*EPSILON(ONE)

PRIVATE
PUBLIC :: PS_WENO_Euler_Flux

PUBLIC :: PS_WENO_Euler_Flux_1D
PUBLIC :: PS_WENO_Euler_Flux_1D_1

PUBLIC :: PS_WENO_Euler_Flux_2D
PUBLIC :: PS_WENO_Euler_Flux_2D_1
PUBLIC :: PS_WENO_Euler_Flux_2D_2

PUBLIC :: PS_WENO_Euler_Flux_3D
PUBLIC :: PS_WENO_Euler_Flux_3D_1
PUBLIC :: PS_WENO_Euler_Flux_3D_2
PUBLIC :: PS_WENO_Euler_Flux_3D_3

CONTAINS

!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_1D (Dimensions,                                &
                                    Index_x, Distributed_x, Order_x, dx,       &
                                             N0, N5, N2, N3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)

  integer  :: Dimensions 
  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx, Gamma

  REALTYPE, dimension(Dimensions  ), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions+2), OPTIONAL :: Lambda
  integer ,                          OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,Dimensions+2) :: Q, D_Flux

  call PS_WENO_Euler_Flux_1D_1 (Dimensions,                                &
                                Index_x, Distributed_x, Order_x, dx,       &
                                         N0, N5, N2, N3,                   &
                                Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                MPI_Comm_Type)

  END Subroutine PS_WENO_Euler_Flux_1D

!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_2D (Dimensions,                                &
                                    Index_x, Distributed_x, Order_x, dx,       &
                                             N0, N5, N2, N3,                   &
                                    Index_y, Distributed_y, Order_y, dy,       &
                                             M0, M5, M2, M3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)

  integer  :: Dimensions
  integer  :: Index_x, Order_x
  integer  :: Index_y, Order_y
  logical  :: Distributed_x
  logical  :: Distributed_y
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dx, dy, Gamma

  REALTYPE, dimension(Dimensions  ), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions+2), OPTIONAL :: Lambda
  integer ,                          OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,Dimensions+2) :: Q, D_Flux

  SELECT CASE (Index_x)
    CASE DEFAULT
      call PS_WENO_Euler_Flux_2D_1 (Dimensions,                                &
                                    Index_x, Distributed_x, Order_x, dx,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)
    CASE (2)
      call PS_WENO_Euler_Flux_2D_2 (Dimensions,                                &
                                    Index_x, Distributed_x, Order_x, dx,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)
  END SELECT

  SELECT CASE (Index_y)
    CASE (1)
      call PS_WENO_Euler_Flux_2D_1 (Dimensions,                                &
                                    Index_y, Distributed_y, Order_y, dy,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)
    CASE DEFAULT
      call PS_WENO_Euler_Flux_2D_2 (Dimensions,                                &
                                    Index_y, Distributed_y, Order_y, dy,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)
  END SELECT

  END Subroutine PS_WENO_Euler_Flux_2D
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_3D (Dimensions,                                &
                                    Index_x, Distributed_x, Order_x, dx,       &
                                             N0, N5, N2, N3,                   &
                                    Index_y, Distributed_y, Order_y, dy,       &
                                             M0, M5, M2, M3,                   &
                                    Index_z, Distributed_z, Order_z, dz,       &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)

  integer  :: Dimensions
  integer  :: Index_x, Order_x
  integer  :: Index_y, Order_y
  integer  :: Index_z, Order_z
  logical  :: Distributed_x
  logical  :: Distributed_y
  logical  :: Distributed_z
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dx, dy, dz, Gamma

  REALTYPE, dimension(Dimensions  ), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions+2), OPTIONAL :: Lambda
  integer ,                          OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Dimensions+2) :: Q, D_Flux

  SELECT CASE (Index_x)
    CASE DEFAULT
      call PS_WENO_Euler_Flux_3D_1 (Dimensions,                                &
                                    Index_x, Distributed_x, Order_x, dx,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D_2 (Dimensions,                                &
                                    Index_x, Distributed_x, Order_x, dx,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D_3 (Dimensions,                                &
                                    Index_x, Distributed_x, Order_x, dx,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)
  END SELECT

  SELECT CASE (Index_y)
    CASE (1)
      call PS_WENO_Euler_Flux_3D_1 (Dimensions,                                &
                                    Index_y, Distributed_y, Order_y, dy,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)

    CASE DEFAULT
      call PS_WENO_Euler_Flux_3D_2 (Dimensions,                                &
                                    Index_y, Distributed_y, Order_y, dy,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D_3 (Dimensions,                                &
                                    Index_y, Distributed_y, Order_y, dy,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)
  END SELECT

  SELECT CASE (Index_z)
    CASE (1)
      call PS_WENO_Euler_Flux_3D_1 (Dimensions,                                &
                                    Index_z, Distributed_z, Order_z, dz,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D_2 (Dimensions,                                &
                                    Index_z, Distributed_z, Order_z, dz,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)

    CASE DEFAULT
      call PS_WENO_Euler_Flux_3D_3 (Dimensions,                                &
                                    Index_z, Distributed_z, Order_z, dz,       &
                                             N0, N5, N2, N3,                   &
                                             M0, M5, M2, M3,                   &
                                             K0, K5, K2, K3,                   &
                                    Gamma, Q, D_Flux, V_Moving_Frame, Lambda,  &
                                    MPI_Comm_Type)
  END SELECT

  END Subroutine PS_WENO_Euler_Flux_3D 
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_1D_1 (Dimensions,                            &
                                      Index, Distributed, Order, dx,         &
                                      N0, N5, i_0, i_1,                      &
                                      Gamma, Q, D_Flux,                      &
                                      V_Moving_Frame, Lambda_o, MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  REALTYPE :: dx, Gamma

  REALTYPE, dimension(Dimensions  ), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions+2), OPTIONAL :: Lambda_o
  integer ,                          OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions)             :: S

  REALTYPE, dimension(N0:N5,Dimensions+2                       ) :: Q, D_Flux

  REALTYPE, dimension(N0:N5,_SIZE(Dimensions)                  ) :: F, F_Half
  REALTYPE, dimension(N0:N5,_SIZE(Dimensions),_SIZE(Dimensions)) :: E_L, E_R

  REALTYPE, dimension(N0:N5) :: U, P, C, r, h 

  REALTYPE, dimension(      _SIZE(Dimensions)                  ) :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(      _SIZE(Dimensions)                  ) :: Lambda_i
#endif

  integer  :: i
  REALTYPE :: s0, s1, t1, Um, hm, qm, cm, rcm, b1, b2, gm1

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  gm1 = Gamma - ONE

!$OMP PARALLEL DO
  do i = i_0-1, i_1
    E_R(i,1,1) = ONE ; E_R(i,1,2) = ONE  ; E_R(i,1,3) = ONE  
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_o)) then
    if (Distributed) call Eigenvalues
  endif
#endif

!$OMP PARALLEL DO 
  do i = N0, N5
    U(i) = Q(i,2)/Q(i,1)
    P(i) = gm1*(Q(i,3) - HALF*Q(i,1)*(U(i)**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C(i) = SQRT(ABS(Gamma*P(i)/Q(i,1)))
#else
                           C(i) = SQRT(ABS(Gamma*P(i)/Q(i,1)))
#endif

    r(i) = SQRT(ABS(Q(i,1)))
    h(i) = (Q(i,3)+P(i))/Q(i,1)

    F(i,1) =   (Q(i,1)         * U(i)       )
    F(i,2) =   (Q(i,2)         * U(i) + P(i))
    F(i,3) =   (Q(i,3) + P(i)) * U(i)
  enddo
!$OMP END PARALLEL DO

  if (S(1) /= ZERO) then
    do i = N0, N5
      F(i,1) = F(i,1) - S(1)*Q(i,1)
      F(i,2) = F(i,2) - S(1)*Q(i,2)
      F(i,3) = F(i,3) - S(1)*Q(i,3)
    enddo
  endif

  if (.NOT. PRESENT(Lambda_o)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)-C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)  )))
      Lambda(3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)+C)))
#if defined (PARALLEL_MPI)
    endif
#endif
  else
    Lambda = Lambda_o
  endif

  Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, hm, qm, cm) &
!$OMP             PRIVATE (rcm, b1, b2, t1)        
    do i = i_0-1, i_1
      s0 = r(i)/(r(i) + r(i+1)) ; s1 = ONE-s0

      Um = s0*U(i) + s1*U(i+1)
      hm = s0*h(i) + s1*h(i+1)

      qm = HALF*(Um*Um)
      cm = SQRT(ABS(gm1*(hm - qm)))

      E_R(i,2,1) = Um - cm
      E_R(i,2,2) = Um
      E_R(i,2,3) = Um + cm

      E_R(i,3,1) = hm - Um*cm
      E_R(i,3,2) = qm
      E_R(i,3,3) = hm + Um*cm

      rcm = ONE/cm
       b1 = gm1*rcm**2 ; b2 = qm*b1 ; t1 = Um*b1 

      E_L(i,1,1) =  HALF*( b2 + Um*rcm )
      E_L(i,1,2) = -HALF*( t1 +    rcm )
      E_L(i,1,3) =  HALF*  b1

      E_L(i,2,1) =   ONE - b2
      E_L(i,2,2) =         t1
      E_L(i,2,3) =       - b1

      E_L(i,3,1) =  HALF*( b2 - Um*rcm )
      E_L(i,3,2) = -HALF*( t1 -    rcm )
      E_L(i,3,3) =  HALF*  b1

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler (Order, N0, N5, i_0-1, i_1, _SIZE(Dimensions), &
                        Q(:,:), F, Lambda, E_L, E_R, F_Half) 

#if defined (PARALLEL_OPENMP)
!$OMP PARALLEL DO
    do i = i_0, i_1
      D_Flux(i,:) = (F_Half(i-1,:) - F_Half(i,:))/dx
    enddo
!$OMP END PARALLEL DO
#else
    D_Flux(i_0:i_1,:) = (F_Half(i_0-1:i_1-1,:) - F_Half(i_0:i_1,:))/dx
#endif

#if defined (PARALLEL_MPI)
  CONTAINS

    Subroutine Eigenvalues

    REALTYPE, dimension(_SIZE(Dimensions)) :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO
    do i = N0, N5
      U(i) = Q(i,2)/Q(i,1)
      P(i) = gm1*(Q(i,3) - HALF*Q(i,1)*(U(i)**2))

      C(i) = SQRT(ABS(Gamma*P(i)/Q(i,1)))
    enddo
!$OMP END PARALLEL DO

    Lambda_Local(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)-C)))
    Lambda_Local(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)  )))
    Lambda_Local(3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)+C)))

    if (N_Processor == 1) then
      Lambda_i = Lambda_Local
    else
                                  MPI_Communicator = MPI_Comm_World
      if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

      if (PRESENT(MPI_Comm_Type)) then
        call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

        ALLOCATE (Sub_Grid(N_Dims))

        Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

        call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                       MPI_Error_Status)

        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator_1D,     &
                            MPI_Error_Status)

        DEALLOCATE (Sub_Grid)

        call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
      else
        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator,        &
                            MPI_Error_Status)
      endif

    endif

    END Subroutine Eigenvalues
#endif

  END Subroutine PS_WENO_Euler_Flux_1D_1
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_2D_1 (Dimensions,                            &
                                      Index, Distributed, Order, dx,         &
                                      N0, N5, i_0, i_1,                      &
                                      M0, M5, j_0, j_1,                      &
                                      Gamma, Q, D_Flux,                      &
                                      V_Moving_Frame, Lambda_o, MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  REALTYPE :: dx, Gamma

  REALTYPE, dimension(Dimensions  ), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions+2), OPTIONAL :: Lambda_o
  integer ,                          OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions)             :: S

  REALTYPE, dimension(N0:N5,M0:M5,Dimensions+2                 ) :: Q, D_Flux

  REALTYPE, dimension(N0:N5,_SIZE(Dimensions)                  ) :: F, F_Half
  REALTYPE, dimension(N0:N5,_SIZE(Dimensions),_SIZE(Dimensions)) :: E_L, E_R

  REALTYPE, dimension(N0:N5) :: U, V, P, C, r, h 

  REALTYPE, dimension(         _SIZE(Dimensions)               ) :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(j_0:j_1, _SIZE(Dimensions)               ) :: Lambda_i
#endif

  integer  :: i, j
  REALTYPE :: s0, s1, t1, t2, Um, Vm, hm, qm, cm, rcm, b1, b2, gm1

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  gm1 = Gamma - ONE

!$OMP PARALLEL DO
  do i = i_0-1, i_1
    E_R(i,1,1) = ONE ; E_R(i,1,2) = ZERO ; E_R(i,1,3) = ONE  ; E_R(i,1,4) = ONE
                       E_R(i,2,2) = ZERO
                       E_R(i,3,2) = ONE

                       E_L(i,2,2) = ZERO ; E_L(i,2,3) = ONE  ; E_L(i,2,4) = ZERO
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_o)) then
    if (Distributed) call Eigenvalues
  endif
#endif

!--------------------- begin of outer loop in y-dir. -------------------

!$OMP PARALLEL DO PRIVATE (U, V, P, C, r, h, F, F_Half, Lambda)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)

  do j = j_0, j_1

    U = Q(:,j,2)/Q(:,j,1)
    V = Q(:,j,3)/Q(:,j,1)
    P = gm1*(Q(:,j,4) - HALF*Q(:,j,1)*(U**2+V**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma*P/Q(:,j,1)))
#else
                           C = SQRT(ABS(Gamma*P/Q(:,j,1)))
#endif

    r = SQRT(ABS(Q(:,j,1)))
    h = (Q(:,j,4)+P)/Q(:,j,1)

    F(:,1) =   (Q(:,j,1)      * U    )
    F(:,2) =   (Q(:,j,2)      * U + P)
    F(:,3) =   (Q(:,j,3)      * U    )
    F(:,4) =   (Q(:,j,4) + P) * U

    if (S(1) /= ZERO) then
      F(:,1) = F(:,1) - S(1)*Q(:,j,1)
      F(:,2) = F(:,2) - S(1)*Q(:,j,2)
      F(:,3) = F(:,3) - S(1)*Q(:,j,3)
      F(:,4) = F(:,4) - S(1)*Q(:,j,4)
    endif

  if (.NOT. PRESENT(Lambda_o)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(j,:)
    else
#endif

      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)-C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)  )))
      Lambda(3) = Lambda(2)
      Lambda(4) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)+C)))
#if defined (PARALLEL_MPI)
    endif
#endif
  else
    Lambda = Lambda_o
  endif

    Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, hm, qm, cm) &
!$OMP             PRIVATE (rcm, b1, b2, t1, t2)
    do i = i_0-1, i_1
      s0 = r(i)/(r(i) + r(i+1)) ; s1 = ONE-s0

      Um = s0*U(i) + s1*U(i+1)
      Vm = s0*V(i) + s1*V(i+1)
      hm = s0*h(i) + s1*h(i+1)

      qm = HALF*(Um*Um + Vm*Vm)
      cm = SQRT(ABS(gm1*(hm - qm)))

      E_R(i,2,1) = Um - cm
      E_R(i,2,3) = Um
      E_R(i,2,4) = Um + cm

      E_R(i,3,1) = Vm
      E_R(i,3,3) = Vm
      E_R(i,3,4) = Vm

      E_R(i,4,1) = hm - Um*cm
      E_R(i,4,2) = Vm
      E_R(i,4,3) = qm
      E_R(i,4,4) = hm + Um*cm

      rcm = ONE/cm
       b1 = gm1*rcm**2 ; b2 = qm * b1 ; t1 = Um * b1 ; t2 = Vm * b1 

      E_L(i,1,1) =  HALF*( b2 + Um*rcm )
      E_L(i,1,2) = -HALF*( t1 +    rcm )
      E_L(i,1,3) = -HALF*  t2
      E_L(i,1,4) =  HALF*  b1

      E_L(i,2,1) =       - Vm

      E_L(i,3,1) =   ONE - b2
      E_L(i,3,2) =         t1
      E_L(i,3,3) =         t2
      E_L(i,3,4) =       - b1

      E_L(i,4,1) =  HALF*( b2 - Um*rcm )
      E_L(i,4,2) = -HALF*( t1 -    rcm )
      E_L(i,4,3) = -HALF*  t2
      E_L(i,4,4) =  HALF*  b1

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler (Order, N0, N5, i_0-1, i_1, _SIZE(Dimensions), &
                        Q(:,j,:), F, Lambda, E_L, E_R, F_Half)

#if defined (PARALLEL_OPENMP)
!$OMP PARALLEL DO
    do i = i_0, i_1
      D_Flux(i,j,:) = (F_Half(i-1,:) - F_Half(i,:))/dx
    enddo
!$OMP END PARALLEL DO
#else
    D_Flux(i_0:i_1,j,:) = (F_Half(i_0-1:i_1-1,:) - F_Half(i_0:i_1,:))/dx
#endif

  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  CONTAINS

    Subroutine Eigenvalues

    REALTYPE, dimension(j_0:j_1,_SIZE(Dimensions)) :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO PRIVATE (U, V, P, C)
    do j = j_0, j_1

      U = Q(:,j,2)/Q(:,j,1)
      V = Q(:,j,3)/Q(:,j,1)
      P = gm1*(Q(:,j,4) - HALF*Q(:,j,1)*(U**2+V**2))

      C = SQRT(ABS(Gamma*P/Q(:,j,1)))

      Lambda_Local(j,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)-C)))
      Lambda_Local(j,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)  )))
      Lambda_Local(j,3) = Lambda_Local(j,2)
      Lambda_Local(j,4) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)+C)))

    enddo
!$OMP END PARALLEL DO

    if (N_Processor == 1) then
      Lambda_i = Lambda_Local
    else
                                  MPI_Communicator = MPI_Comm_World
      if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

      if (PRESENT(MPI_Comm_Type)) then
        call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

        ALLOCATE (Sub_Grid(N_Dims))

        Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

        call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                       MPI_Error_Status)

        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator_1D,     &
                            MPI_Error_Status)

        DEALLOCATE (Sub_Grid)

        call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
      else
        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator,        &
                            MPI_Error_Status)
      endif

    endif

    END Subroutine Eigenvalues
#endif

  END Subroutine PS_WENO_Euler_Flux_2D_1
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_2D_2 (Dimensions,                            &
                                      Index, Distributed, Order, dy,         &
                                      N0, N5, i_0, i_1,                      &
                                      M0, M5, j_0, j_1,                      &
                                      Gamma, Q, D_Flux,                      &
                                      V_Moving_Frame, Lambda_o, MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  REALTYPE :: dy, Gamma

  REALTYPE, dimension(Dimensions  ), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions+2), OPTIONAL :: Lambda_o
  integer ,                          OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions)             :: S

  REALTYPE, dimension(N0:N5,M0:M5,Dimensions+2                 ) :: Q, D_Flux

  REALTYPE, dimension(M0:M5,_SIZE(Dimensions)                  ) :: F, F_Half
  REALTYPE, dimension(M0:M5,_SIZE(Dimensions),_SIZE(Dimensions)) :: E_L, E_R

  REALTYPE, dimension(M0:M5) :: U, V, P, C, r, h 

  REALTYPE, dimension(         _SIZE(Dimensions)               ) :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(i_0:i_1, _SIZE(Dimensions)               ) :: Lambda_i
#endif

  integer  :: i, j
  REALTYPE :: s0, s1, t1, t2, Um, Vm, hm, qm, cm, rcm, b1, b2, gm1

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  gm1 = Gamma - ONE

!$OMP PARALLEL DO
  do j = j_0-1, j_1
    E_R(j,1,1) = ONE ; E_R(j,1,2) = ZERO ; E_R(j,1,3) = ONE  ; E_R(j,1,4) = ONE
                       E_R(j,2,2) = ONE
                       E_R(j,3,2) = ZERO

                       E_L(j,2,2) = ONE  ; E_L(j,2,3) = ZERO ; E_L(j,2,4) = ZERO
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_o)) then
    if (Distributed) call Eigenvalues
  endif
#endif

!--------------------- begin of outer loop in x-dir. -------------------

!$OMP PARALLEL DO PRIVATE (U, V, P, C, r, h, F, F_Half, Lambda)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)
  do i = i_0, i_1

    U = Q(i,:,2)/Q(i,:,1)
    V = Q(i,:,3)/Q(i,:,1)
    P = gm1*(Q(i,:,4) - HALF*Q(i,:,1)*(U**2+V**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma*P/Q(i,:,1)))
#else
                           C = SQRT(ABS(Gamma*P/Q(i,:,1)))
#endif

    r = SQRT(ABS(Q(i,:,1)))
    h = (Q(i,:,4)+P)/Q(i,:,1)

    F(:,1) =  (Q(i,:,1)      * V    )
    F(:,2) =  (Q(i,:,2)      * V    )
    F(:,3) =  (Q(i,:,3)      * V + P)
    F(:,4) =  (Q(i,:,4) + P) * V

    if (S(2) /= ZERO) then
      F(:,1) = F(:,1) - S(2)*Q(i,:,1)
      F(:,2) = F(:,2) - S(2)*Q(i,:,2)
      F(:,3) = F(:,3) - S(2)*Q(i,:,3)
      F(:,4) = F(:,4) - S(2)*Q(i,:,4)
    endif

  if (.NOT. PRESENT(Lambda_o)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(i,:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)-C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)  )))
      Lambda(3) = Lambda(2)
      Lambda(4) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)+C)))
#if defined (PARALLEL_MPI)
    endif
#endif
  else
    Lambda = Lambda_o
  endif

    Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, hm, qm, cm) &
!$OMP             PRIVATE (rcm, b1, b2, t1, t2)
    do j = j_0-1, j_1
      s0 = r(j)/(r(j) + r(j+1)) ; s1 = ONE-s0

      Um = s0*U(j) + s1*U(j+1)
      Vm = s0*V(j) + s1*V(j+1)
      hm = s0*h(j) + s1*h(j+1)

      qm = HALF*(Um*Um + Vm*Vm)
      cm = SQRT(ABS(gm1*(hm - qm)))

      E_R(j,2,1) = Um
      E_R(j,2,3) = Um
      E_R(j,2,4) = Um

      E_R(j,3,1) = Vm - cm
      E_R(j,3,3) = Vm
      E_R(j,3,4) = Vm + cm

      E_R(j,4,1) = hm - Vm*cm
      E_R(j,4,2) = Um
      E_R(j,4,3) = qm
      E_R(j,4,4) = hm + Vm*cm

      rcm = ONE/cm
       b1 = gm1*rcm**2 ; b2 = qm * b1 ; t1 = Um * b1 ; t2 = Vm * b1 

      E_L(j,1,1) =  HALF*( b2 + Vm*rcm )
      E_L(j,1,2) = -HALF*  t1
      E_L(j,1,3) = -HALF*( t2 +    rcm )
      E_L(j,1,4) =  HALF*  b1

      E_L(j,2,1) =       - Um

      E_L(j,3,1) =   ONE - b2
      E_L(j,3,2) =         t1
      E_L(j,3,3) =         t2
      E_L(j,3,4) =       - b1

      E_L(j,4,1) =  HALF*( b2 - Vm*rcm )
      E_L(j,4,2) = -HALF*  t1
      E_L(j,4,3) = -HALF*( t2 -    rcm )
      E_L(j,4,4) =  HALF*  b1

    enddo
!$OMP END PARALLEL DO
     call PS_WENO_Euler (Order, M0, M5, j_0-1, j_1, _SIZE(Dimensions), &
                         Q(i,:,:), F, Lambda, E_L, E_R, F_Half)

#if defined (PARALLEL_OPENMP)
!$OMP PARALLEL DO
    do j = j_0, j_1
      D_Flux(i,j,:) = D_Flux(i,j,:) + (F_Half(j-1,:) - F_Half(j,:))/dy
    enddo
!$OMP END PARALLEL DO
#else
    D_Flux(i,j_0:j_1,:) = D_Flux(i,j_0:j_1,:) + &
                         (F_Half(j_0-1:j_1-1,:) - F_Half(j_0:j_1,:))/dy
#endif

  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  CONTAINS

    Subroutine Eigenvalues 

    REALTYPE, dimension(i_0:i_1,_SIZE(Dimensions)) :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO PRIVATE (U, V, P, C) 
    do i = i_0, i_1

      U = Q(i,:,2)/Q(i,:,1)
      V = Q(i,:,3)/Q(i,:,1)
      P = gm1*(Q(i,:,4) - HALF*Q(i,:,1)*(U**2+V**2))

      C = SQRT(ABS(Gamma*P/Q(i,:,1)))

      Lambda_Local(i,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)-C)))
      Lambda_Local(i,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)  )))
      Lambda_Local(i,3) = Lambda_Local(i,2)
      Lambda_Local(i,4) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)+C)))

    enddo
!$OMP END PARALLEL DO

    if (N_Processor == 1) then
      Lambda_i = Lambda_Local
    else
                                  MPI_Communicator = MPI_Comm_World
      if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

      if (PRESENT(MPI_Comm_Type)) then
        call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

        ALLOCATE (Sub_Grid(N_Dims))

        Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

        call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                       MPI_Error_Status)

        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator_1D,     &
                            MPI_Error_Status)

        DEALLOCATE (Sub_Grid)

        call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
      else
        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator,        &
                            MPI_Error_Status)
      endif

    endif

    END Subroutine Eigenvalues 
#endif

  END Subroutine PS_WENO_Euler_Flux_2D_2 
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_3D_1 (Dimensions,                            &
                                      Index, Distributed, Order, dx,         &
                                      N0, N5, i_0, i_1,                      &
                                      M0, M5, j_0, j_1,                      &
                                      K0, K5, k_0, k_1,                      &
                                      Gamma, Q, D_Flux,                      &
                                      V_Moving_Frame, Lambda_o, MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1
  REALTYPE :: dx, Gamma

  REALTYPE, dimension(Dimensions  ), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions+2), OPTIONAL :: Lambda_o
  integer ,                          OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions)             :: S

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Dimensions+2           ) :: Q, D_Flux

  REALTYPE, dimension(N0:N5,_SIZE(Dimensions)                  ) :: F, F_Half
  REALTYPE, dimension(N0:N5,_SIZE(Dimensions),_SIZE(Dimensions)) :: E_L, E_R

  REALTYPE, dimension(N0:N5) :: U, V, W, P, C, r, h 

  REALTYPE, dimension(                 _SIZE(Dimensions)) :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(j_0:j_1,k_0:k_1, _SIZE(Dimensions)) :: Lambda_i
#endif

  integer  :: i, j, k
  REALTYPE :: s0, s1, t1, t2, t3, Um, Vm, Wm, hm, qm, cm, rcm, b1, b2, gm1

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  gm1 = Gamma - ONE

!$OMP PARALLEL DO
  do i = i_0-1, i_1
    E_R(i,1,1) = ONE ; E_R(i,1,2) = ZERO ; E_R(i,1,3) = ZERO 
                                           E_R(i,1,4) = ONE  ; E_R(i,1,5) = ONE
                       E_R(i,2,2) = ZERO ; E_R(i,2,3) = ZERO
                       E_R(i,3,2) = ONE  ; E_R(i,3,3) = ZERO
                       E_R(i,4,2) = ZERO ; E_R(i,4,3) = ONE

                       E_L(i,2,2) = ZERO ; E_L(i,2,3) = ONE  
                                           E_L(i,2,4) = ZERO ; E_L(i,2,5) = ZERO
                       E_L(i,3,2) = ZERO ; E_L(i,3,3) = ZERO 
                                           E_L(i,3,4) = ONE  ; E_L(i,3,5) = ZERO
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_o)) then
    if (Distributed) call Eigenvalues
  endif
#endif

!--------------------- begin of outer loop in y-dir and z-dir ------------

!$OMP PARALLEL DO  
  do k = k_0, k_1
!$OMP PARALLEL DO PRIVATE (U, V, W, P, C, r, h, F, F_Half, Lambda)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)                                
  do j = j_0, j_1

    U = Q(:,j,k,2)/Q(:,j,k,1)
    V = Q(:,j,k,3)/Q(:,j,k,1)
    W = Q(:,j,k,4)/Q(:,j,k,1)
    P = gm1*(Q(:,j,k,5) - HALF*Q(:,j,k,1)*(U**2+V**2+W**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma*P/Q(:,j,k,1)))
#else
                           C = SQRT(ABS(Gamma*P/Q(:,j,k,1)))
#endif

    r = SQRT(ABS(Q(:,j,k,1)))
    h = (Q(:,j,k,5)+P)/Q(:,j,k,1)

    F(:,1) =   (Q(:,j,k,1)      * U    )
    F(:,2) =   (Q(:,j,k,2)      * U + P)
    F(:,3) =   (Q(:,j,k,3)      * U    )
    F(:,4) =   (Q(:,j,k,4)      * U    )
    F(:,5) =   (Q(:,j,k,5) + P) * U

    if (S(1) /= ZERO) then
      F(:,1) = F(:,1) - S(1)*Q(:,j,k,1)
      F(:,2) = F(:,2) - S(1)*Q(:,j,k,2)
      F(:,3) = F(:,3) - S(1)*Q(:,j,k,3)
      F(:,4) = F(:,4) - S(1)*Q(:,j,k,4)
      F(:,5) = F(:,5) - S(1)*Q(:,j,k,5)
    endif

  if (.NOT. PRESENT(Lambda_o)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(j,k,:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)-C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)  )))
      Lambda(3) = Lambda(2)
      Lambda(4) = Lambda(2)
      Lambda(5) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)+C)))
#if defined (PARALLEL_MPI)
    endif
#endif
  else
    Lambda = Lambda_o
  endif

    Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, Wm, hm, qm, cm) &
!$OMP             PRIVATE (rcm, b1, b2, t1, t2, t3)
    do i = i_0-1, i_1
      s0 = r(i)/(r(i) + r(i+1)) ; s1 = ONE-s0

      Um = s0*U(i) + s1*U(i+1)
      Vm = s0*V(i) + s1*V(i+1)
      Wm = s0*W(i) + s1*W(i+1)
      hm = s0*h(i) + s1*h(i+1)

      qm = HALF*(Um*Um + Vm*Vm + Wm*Wm)
      cm = SQRT(ABS(gm1*(hm - qm)))

      E_R(i,2,1) = Um - cm
      E_R(i,2,4) = Um 
      E_R(i,2,5) = Um + cm

      E_R(i,3,1) = Vm
      E_R(i,3,4) = Vm
      E_R(i,3,5) = Vm

      E_R(i,4,1) = Wm
      E_R(i,4,4) = Wm
      E_R(i,4,5) = Wm

      E_R(i,5,1) = hm - Um*cm
      E_R(i,5,2) = Vm
      E_R(i,5,3) = Wm
      E_R(i,5,4) = qm
      E_R(i,5,5) = hm + Um*cm

      rcm = ONE/cm 
       b1 = gm1*rcm**2 ; b2 = qm*b1 ; t1 = Um*b1 ; t2 = Vm*b1 ; t3 = Wm*b1

      E_L(i,1,1) =  HALF*( b2 + Um*rcm )
      E_L(i,1,2) = -HALF*( t1 +    rcm )
      E_L(i,1,3) = -HALF*  t2
      E_L(i,1,4) = -HALF*  t3
      E_L(i,1,5) =  HALF*  b1

      E_L(i,2,1) =       - Vm

      E_L(i,3,1) =       - Wm

      E_L(i,4,1) =   ONE - b2
      E_L(i,4,2) =         t1
      E_L(i,4,3) =         t2
      E_L(i,4,4) =         t3
      E_L(i,4,5) =       - b1

      E_L(i,5,1) =  HALF*( b2 - Um*rcm )
      E_L(i,5,2) = -HALF*( t1 -    rcm )
      E_L(i,5,3) = -HALF*  t2
      E_L(i,5,4) = -HALF*  t3
      E_L(i,5,5) =  HALF*  b1

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler (Order, N0, N5, i_0-1, i_1, _SIZE(Dimensions), &
                        Q(:,j,k,:), F, Lambda, E_L, E_R, F_Half) 

#if defined (PARALLEL_OPENMP)
!$OMP PARALLEL DO
    do i = i_0, i_1
      D_Flux(i,j,k,:) = (F_Half(i-1,:) - F_Half(i,:))/dx
    enddo
!$OMP END PARALLEL DO
#else
    D_Flux(i_0:i_1,j,k,:) = (F_Half(i_0-1:i_1-1,:) - F_Half(i_0:i_1,:))/dx
#endif

  enddo
!$OMP END PARALLEL DO
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  CONTAINS

    Subroutine Eigenvalues

    REALTYPE, dimension(j_0:j_1,k_0:k_1,_SIZE(Dimensions)) :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO PRIVATE (U, V, W, P, C)
    do k = k_0, k_1
!$OMP PARALLEL DO
    do j = j_0, j_1

      U = Q(:,j,k,2)/Q(:,j,k,1)
      V = Q(:,j,k,3)/Q(:,j,k,1)
      W = Q(:,j,k,4)/Q(:,j,k,1)
      P = gm1*(Q(:,j,k,5) - HALF*Q(:,j,k,1)*(U**2+V**2+W**2))

      C = SQRT(ABS(Gamma*P/Q(:,j,k,1)))

      Lambda_Local(j,k,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)-C)))
      Lambda_Local(j,k,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)  )))
      Lambda_Local(j,k,3) = Lambda_Local(j,k,2)
      Lambda_Local(j,k,4) = Lambda_Local(j,k,2)
      Lambda_Local(j,k,5) = MAX(Sound_Speed_Min, MAXVAL(ABS(U-S(1)+C)))

    enddo
!$OMP END PARALLEL DO
    enddo
!$OMP END PARALLEL DO

    if (N_Processor == 1) then
      Lambda_i = Lambda_Local
    else
                                  MPI_Communicator = MPI_Comm_World
      if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

      if (PRESENT(MPI_Comm_Type)) then
        call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

        ALLOCATE (Sub_Grid(N_Dims))

        Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

        call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                       MPI_Error_Status)

        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator_1D,     &
                            MPI_Error_Status)

        DEALLOCATE (Sub_Grid)

        call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
      else
        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator,        &
                            MPI_Error_Status)
      endif

    endif

    END Subroutine Eigenvalues
#endif

  END Subroutine PS_WENO_Euler_Flux_3D_1
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_3D_2 (Dimensions,                            &
                                      Index, Distributed, Order, dy,         &
                                      N0, N5, i_0, i_1,                      &
                                      M0, M5, j_0, j_1,                      &
                                      K0, K5, k_0, k_1,                      &
                                      Gamma, Q, D_Flux,                      &
                                      V_Moving_Frame, Lambda_o, MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1
  REALTYPE :: dy, Gamma

  REALTYPE, dimension(Dimensions  ), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions+2), OPTIONAL :: Lambda_o
  integer ,                          OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions)             :: S

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Dimensions+2           ) :: Q, D_Flux

  REALTYPE, dimension(M0:M5,_SIZE(Dimensions)                  ) :: F, F_Half
  REALTYPE, dimension(M0:M5,_SIZE(Dimensions),_SIZE(Dimensions)) :: E_L, E_R

  REALTYPE, dimension(M0:M5) :: U, V, W, P, C, r, h 

  REALTYPE, dimension(                 _SIZE(Dimensions)) :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(i_0:i_1,k_0:k_1, _SIZE(Dimensions)) :: Lambda_i
#endif

  integer  :: i, j, k
  REALTYPE :: s0, s1, t1, t2, t3, Um, Vm, Wm, hm, qm, cm, rcm, b1, b2, gm1

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  gm1 = Gamma - ONE

!$OMP PARALLEL DO
  do j = j_0-1, j_1
    E_R(j,1,1) = ONE ; E_R(j,1,2) = ZERO ; E_R(j,1,3) = ZERO 
                                           E_R(j,1,4) = ONE  ; E_R(j,1,5) = ONE
                       E_R(j,2,2) = ONE  ; E_R(j,2,3) = ZERO
                       E_R(j,3,2) = ZERO ; E_R(j,3,3) = ZERO
                       E_R(j,4,2) = ZERO ; E_R(j,4,3) = ONE

                       E_L(j,2,2) = ONE  ; E_L(j,2,3) = ZERO 
                                           E_L(j,2,4) = ZERO ; E_L(j,2,5) = ZERO
                       E_L(j,3,2) = ZERO ; E_L(j,3,3) = ZERO 
                                           E_L(j,3,4) = ONE  ; E_L(j,3,5) = ZERO
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_o)) then
    if (Distributed) call Eigenvalues
  endif
#endif

!--------------------- begin of outer loop in x-dir and z-dir ------------

!$OMP PARALLEL DO PRIVATE (U, V, W, P, C, r, h, F, F_Half, Lambda)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)
  do k = k_0, k_1
!$OMP PARALLEL DO 
  do i = i_0, i_1

    U = Q(i,:,k,2)/Q(i,:,k,1)
    V = Q(i,:,k,3)/Q(i,:,k,1)
    W = Q(i,:,k,4)/Q(i,:,k,1)
    P = gm1*(Q(i,:,k,5) - HALF*Q(i,:,k,1)*(U**2+V**2+W**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma*P/Q(i,:,k,1)))
#else
                           C = SQRT(ABS(Gamma*P/Q(i,:,k,1)))
#endif

    r = SQRT(ABS(Q(i,:,k,1)))
    h = (Q(i,:,k,5)+P)/Q(i,:,k,1)

    F(:,1) =   (Q(i,:,k,1)      * V    )
    F(:,2) =   (Q(i,:,k,2)      * V    )
    F(:,3) =   (Q(i,:,k,3)      * V + P)
    F(:,4) =   (Q(i,:,k,4)      * V    )
    F(:,5) =   (Q(i,:,k,5) + P) * V

    if (S(2) /= ZERO) then
      F(:,1) = F(:,1) - S(2)*Q(i,:,k,1)
      F(:,2) = F(:,2) - S(2)*Q(i,:,k,2)
      F(:,3) = F(:,3) - S(2)*Q(i,:,k,3)
      F(:,4) = F(:,4) - S(2)*Q(i,:,k,4)
      F(:,5) = F(:,5) - S(2)*Q(i,:,k,5)
    endif

  if (.NOT. PRESENT(Lambda_o)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(i,k,:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)-C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)  )))
      Lambda(3) = Lambda(2)
      Lambda(4) = Lambda(2)
      Lambda(5) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)+C)))
#if defined (PARALLEL_MPI)
    endif
#endif
  else
    Lambda = Lambda_o
  endif

    Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, Wm, hm, qm, cm) &
!$OMP             PRIVATE (rcm, b1, b2, t1, t2, t3)
    do j = j_0-1, j_1
      s0 = r(j)/(r(j) + r(j+1)) ; s1 = ONE-s0

      Um = s0*U(j) + s1*U(j+1)
      Vm = s0*V(j) + s1*V(j+1)
      Wm = s0*W(j) + s1*W(j+1)
      hm = s0*h(j) + s1*h(j+1)

      qm = HALF*(Um*Um + Vm*Vm + Wm*Wm)
      cm = SQRT(ABS(gm1*(hm - qm)))

      E_R(j,2,1) = Um
      E_R(j,2,4) = Um
      E_R(j,2,5) = Um

      E_R(j,3,1) = Vm - cm
      E_R(j,3,4) = Vm
      E_R(j,3,5) = Vm + cm

      E_R(j,4,1) = Wm
      E_R(j,4,4) = Wm
      E_R(j,4,5) = Wm

      E_R(j,5,1) = hm - Vm*cm
      E_R(j,5,2) = Um
      E_R(j,5,3) = Wm
      E_R(j,5,4) = qm
      E_R(j,5,5) = hm + Vm*cm

      rcm = ONE/cm 
       b1 = gm1*rcm**2 ; b2 = qm*b1 ; t1 = Um*b1 ; t2 = Vm*b1 ; t3 = Wm*b1

      E_L(j,1,1) =  HALF*( b2 + Vm*rcm )
      E_L(j,1,2) = -HALF*  t1
      E_L(j,1,3) = -HALF*( t2 +    rcm )
      E_L(j,1,4) = -HALF*  t3
      E_L(j,1,5) =  HALF*  b1

      E_L(j,2,1) =       - Um

      E_L(j,3,1) =       - Wm

      E_L(j,4,1) =   ONE - b2
      E_L(j,4,2) =         t1
      E_L(j,4,3) =         t2
      E_L(j,4,4) =         t3
      E_L(j,4,5) =       - b1

      E_L(j,5,1) =  HALF*( b2 - Vm*rcm )
      E_L(j,5,2) = -HALF*  t1
      E_L(j,5,3) = -HALF*( t2 -    rcm )
      E_L(j,5,4) = -HALF*  t3
      E_L(j,5,5) =  HALF*  b1

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler (Order, M0, M5, j_0-1, j_1, _SIZE(Dimensions), &
                        Q(i,:,k,:), F, Lambda, E_L, E_R, F_Half) 

#if defined (PARALLEL_OPENMP)
!$OMP PARALLEL DO
    do j = j_0, j_1
      D_Flux(i,j,k,:) = D_Flux(i,j,k,:) + (F_Half(j-1,:) - F_Half(j,:))/dy
    enddo
!$OMP END PARALLEL DO
#else
    D_Flux(i,j_0:j_1,k,:) = D_Flux(i,j_0:j_1,k,:) + &
                            (F_Half(j_0-1:j_1-1,:) - F_Half(j_0:j_1,:))/dy
#endif

  enddo
!$OMP END PARALLEL DO
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  CONTAINS

    Subroutine Eigenvalues

    REALTYPE, dimension(i_0:i_1,k_0:k_1,_SIZE(Dimensions)) :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO PRIVATE (U, V, W, P, C)
    do k = k_0, k_1
!$OMP PARALLEL DO 
    do i = i_0, i_1

      U = Q(i,:,k,2)/Q(i,:,k,1)
      V = Q(i,:,k,3)/Q(i,:,k,1)
      W = Q(i,:,k,4)/Q(i,:,k,1)
      P = gm1*(Q(i,:,k,5) - HALF*Q(i,:,k,1)*(U**2+V**2+W**2))

      C = SQRT(ABS(Gamma*P/Q(i,:,k,1)))

      Lambda_Local(i,k,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)-C)))
      Lambda_Local(i,k,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)  )))
      Lambda_Local(i,k,3) = Lambda_Local(i,k,2)
      Lambda_Local(i,k,4) = Lambda_Local(i,k,2)
      Lambda_Local(i,k,5) = MAX(Sound_Speed_Min, MAXVAL(ABS(V-S(2)+C)))

    enddo
!$OMP END PARALLEL DO
    enddo
!$OMP END PARALLEL DO

    if (N_Processor == 1) then
      Lambda_i = Lambda_Local
    else
                                  MPI_Communicator = MPI_Comm_World
      if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

      if (PRESENT(MPI_Comm_Type)) then
        call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

        ALLOCATE (Sub_Grid(N_Dims))

        Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

        call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                       MPI_Error_Status)

        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator_1D,     &
                            MPI_Error_Status)

        DEALLOCATE (Sub_Grid)

        call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
      else
        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator,        &
                            MPI_Error_Status)
      endif

    endif

    END Subroutine Eigenvalues
#endif

  END Subroutine PS_WENO_Euler_Flux_3D_2 
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_3D_3 (Dimensions,                            &
                                      Index, Distributed, Order, dz,         &
                                      N0, N5, i_0, i_1,                      &
                                      M0, M5, j_0, j_1,                      &
                                      K0, K5, k_0, k_1,                      &
                                      Gamma, Q, D_Flux,                      &
                                      V_Moving_Frame, Lambda_o, MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1
  REALTYPE :: dz, Gamma

  REALTYPE, dimension(Dimensions  ), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions+2), OPTIONAL :: Lambda_o
  integer ,                          OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions)             :: S

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Dimensions+2           ) :: Q, D_Flux

  REALTYPE, dimension(K0:K5,_SIZE(Dimensions)                  ) :: F, F_Half
  REALTYPE, dimension(K0:K5,_SIZE(Dimensions),_SIZE(Dimensions)) :: E_L, E_R

  REALTYPE, dimension(K0:K5) :: U, V, W, P, C, r, h 

  REALTYPE, dimension(                 _SIZE(Dimensions)) :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(i_0:i_1,j_0:j_1, _SIZE(Dimensions)) :: Lambda_i
#endif

  integer  :: i, j, k
  REALTYPE :: s0, s1, t1, t2, t3, Um, Vm, Wm, hm, qm, cm, rcm, b1, b2, gm1

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  gm1 = Gamma - ONE

!$OMP PARALLEL DO
  do k = k_0-1, k_1
    E_R(k,1,1) = ONE ; E_R(k,1,2) = ZERO ; E_R(k,1,3) = ZERO 
                                           E_R(k,1,4) = ONE  ; E_R(k,1,5) = ONE
                       E_R(k,2,2) = ONE  ; E_R(k,2,3) = ZERO
                       E_R(k,3,2) = ZERO ; E_R(k,3,3) = ONE
                       E_R(k,4,2) = ZERO ; E_R(k,4,3) = ZERO

                       E_L(k,2,2) = ONE  ; E_L(k,2,3) = ZERO 
                                           E_L(k,2,4) = ZERO ; E_L(k,2,5) = ZERO
                       E_L(k,3,2) = ZERO ; E_L(k,3,3) = ONE  
                                           E_L(k,3,4) = ZERO ; E_L(k,3,5) = ZERO
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_o)) then
    if (Distributed) call Eigenvalues
  endif
#endif

!--------------------- begin of outer loop in y-dir and z-dir ------------

!$OMP PARALLEL DO PRIVATE (U, V, W, P, C, r, h, F, F_Half, Lambda)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)
  do j = j_0, j_1
!$OMP PARALLEL DO 
  do i = i_0, i_1

    U = Q(i,j,:,2)/Q(i,j,:,1)
    V = Q(i,j,:,3)/Q(i,j,:,1)
    W = Q(i,j,:,4)/Q(i,j,:,1)
    P = gm1*(Q(i,j,:,5) - HALF*Q(i,j,:,1)*(U**2+V**2+W**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma*P/Q(i,j,:,1)))
#else
                           C = SQRT(ABS(Gamma*P/Q(i,j,:,1)))
#endif

    r = SQRT(ABS(Q(i,j,:,1)))
    h = (Q(i,j,:,5)+P)/Q(i,j,:,1)

    F(:,1) =   (Q(i,j,:,1)      * W    )
    F(:,2) =   (Q(i,j,:,2)      * W    )
    F(:,3) =   (Q(i,j,:,3)      * W    )
    F(:,4) =   (Q(i,j,:,4)      * W + P)
    F(:,5) =   (Q(i,j,:,5) + P) * W

    if (S(3) /= ZERO) then
      F(:,1) = F(:,1) - S(3)*Q(i,j,:,1)
      F(:,2) = F(:,2) - S(3)*Q(i,j,:,2)
      F(:,3) = F(:,3) - S(3)*Q(i,j,:,3)
      F(:,4) = F(:,4) - S(3)*Q(i,j,:,4)
      F(:,5) = F(:,5) - S(3)*Q(i,j,:,5)
    endif

  if (.NOT. PRESENT(Lambda_o)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(i,j,:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(W-S(3)-C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(W-S(3)  )))
      Lambda(3) = Lambda(2)
      Lambda(4) = Lambda(2)
      Lambda(5) = MAX(Sound_Speed_Min, MAXVAL(ABS(W-S(3)+C)))
#if defined (PARALLEL_MPI)
    endif
#endif
  else
    Lambda = Lambda_o
  endif

    Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, Wm, hm, qm, cm) &
!$OMP             PRIVATE (rcm, b1, b2, t1, t2, t3)
    do k = k_0-1, k_1
      s0 = r(k)/(r(k) + r(k+1)) ; s1 = ONE-s0

      Um = s0*U(k) + s1*U(k+1)
      Vm = s0*V(k) + s1*V(k+1)
      Wm = s0*W(k) + s1*W(k+1)
      hm = s0*h(k) + s1*h(k+1)

      qm = HALF*(Um*Um + Vm*Vm + Wm*Wm)
      cm = SQRT(ABS(gm1*(hm - qm)))

      E_R(k,2,1) = Um
      E_R(k,2,4) = Um
      E_R(k,2,5) = Um

      E_R(k,3,1) = Vm
      E_R(k,3,4) = Vm
      E_R(k,3,5) = Vm

      E_R(k,4,1) = Wm - cm
      E_R(k,4,4) = Wm 
      E_R(k,4,5) = Wm + cm

      E_R(k,5,1) = hm - Wm*cm
      E_R(k,5,2) = Um
      E_R(k,5,3) = Vm
      E_R(k,5,4) = qm
      E_R(k,5,5) = hm + Wm*cm

      rcm = ONE/cm 
       b1 = gm1*rcm**2 ; b2 = qm*b1 ; t1 = Um*b1 ; t2 = Vm*b1 ; t3 = Wm*b1

      E_L(k,1,1) =  HALF*( b2 + Wm*rcm )
      E_L(k,1,2) = -HALF*  t1
      E_L(k,1,3) = -HALF*  t2
      E_L(k,1,4) = -HALF*( t3 +    rcm )
      E_L(k,1,5) =  HALF*  b1

      E_L(k,2,1) =       - Um

      E_L(k,3,1) =       - Vm

      E_L(k,4,1) =   ONE - b2
      E_L(k,4,2) =         t1
      E_L(k,4,3) =         t2
      E_L(k,4,4) =         t3
      E_L(k,4,5) =       - b1

      E_L(k,5,1) =  HALF*( b2 - Wm*rcm )
      E_L(k,5,2) = -HALF*  t1 
      E_L(k,5,3) = -HALF*  t2
      E_L(k,5,4) = -HALF*( t3 -    rcm )
      E_L(k,5,5) =  HALF*  b1

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler (Order, K0, K5, k_0-1, k_1, _SIZE(Dimensions), &
                        Q(i,j,:,:), F, Lambda, E_L, E_R, F_Half) 

#if defined (PARALLEL_OPENMP)
!$OMP PARALLEL DO
    do k = k_0, k_1
      D_Flux(i,j,k,:) = D_Flux(i,j,k,:) + (F_Half(k-1,:) - F_Half(k,:))/dz
    enddo
!$OMP END PARALLEL DO
#else
    D_Flux(i,j,k_0:k_1,:) = D_Flux(i,j,k_0:k_1,:) + &
                            (F_Half(k_0-1:k_1-1,:) - F_Half(k_0:k_1,:))/dz
#endif

  enddo
!$OMP END PARALLEL DO
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  CONTAINS

    Subroutine Eigenvalues

    REALTYPE, dimension(i_0:i_1,j_0:j_1,_SIZE(Dimensions)) :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO PRIVATE (U, V, W, P, C)
    do j = j_0, j_1
!$OMP PARALLEL DO
    do i = i_0, i_1

      U = Q(i,j,:,2)/Q(i,j,:,1)
      V = Q(i,j,:,3)/Q(i,j,:,1)
      W = Q(i,j,:,4)/Q(i,j,:,1)
      P = gm1*(Q(i,j,:,5) - HALF*Q(i,j,:,1)*(U**2+V**2+W**2))

      C = SQRT(ABS(Gamma*P/Q(i,j,:,1)))

      Lambda_Local(i,j,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(W-S(3)-C)))
      Lambda_Local(i,j,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(W-S(3)  )))
      Lambda_Local(i,j,3) = Lambda_Local(i,j,2)
      Lambda_Local(i,j,4) = Lambda_Local(i,j,2)
      Lambda_Local(i,j,5) = MAX(Sound_Speed_Min, MAXVAL(ABS(W-S(3)+C)))

    enddo
!$OMP END PARALLEL DO
    enddo
!$OMP END PARALLEL DO

    if (N_Processor == 1) then
      Lambda_i = Lambda_Local
    else
                                  MPI_Communicator = MPI_Comm_World
      if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

      if (PRESENT(MPI_Comm_Type)) then
        call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

        ALLOCATE (Sub_Grid(N_Dims))

        Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

        call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                       MPI_Error_Status)

        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator_1D,     &
                            MPI_Error_Status)

        DEALLOCATE (Sub_Grid)

        call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
      else
        call MPI_ALLREDUCE (Lambda_Local, Lambda_i, SIZE(Lambda_i),         &
                            MPI_REALTYPE, MPI_MAX, MPI_Communicator,        &
                            MPI_Error_Status)
      endif

    endif

    END Subroutine Eigenvalues
#endif

  END Subroutine PS_WENO_Euler_Flux_3D_3

END MODULE WENO_Euler_Flux
