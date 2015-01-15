#include "pseudopack.h"


MODULE WENO_Euler_Flux_TWO_Species

#if defined (PARALLEL_MPI)
   USE Processor
#endif
  USE Euler_Function
  USE WENO_Euler
  USE WENO_Option

implicit NONE

INTERFACE PS_WENO_Euler_Flux
  MODULE PROCEDURE PS_WENO_Euler_Flux_1D2S
  MODULE PROCEDURE PS_WENO_Euler_Flux_2D2S
  MODULE PROCEDURE PS_WENO_Euler_Flux_3D2S
END INTERFACE

REALTYPE :: Lambda_Factor   = ELEVEN/TEN
REALTYPE :: Sound_Speed_Min = TEN*EPSILON(ONE)

PRIVATE
PUBLIC :: PS_WENO_Euler_Flux

PUBLIC :: PS_WENO_Euler_Flux_1D2S
PUBLIC :: PS_WENO_Euler_Flux_1D2S_1

PUBLIC :: PS_WENO_Euler_Flux_2D2S
PUBLIC :: PS_WENO_Euler_Flux_2D2S_1
PUBLIC :: PS_WENO_Euler_Flux_2D2S_2

PUBLIC :: PS_WENO_Euler_Flux_3D2S
PUBLIC :: PS_WENO_Euler_Flux_3D2S_1
PUBLIC :: PS_WENO_Euler_Flux_3D2S_2
PUBLIC :: PS_WENO_Euler_Flux_3D2S_3

CONTAINS

!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_1D2S (Dimensions,                              &
                                      Index_x, Distributed_x, Order_x, dx,     &
                                               N0, N5, N2, N3, NV,             &
                                      Gamma_Formulation, Gamma_Gas,            &
                                      Q, D_Flux, V_Moving_Frame, Lambda,       &
                                      MPI_Comm_Type)

  integer  :: Dimensions

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  REALTYPE, dimension(N0:N5,NV) :: Q, D_Flux

  integer                       :: Gamma_Formulation
  REALTYPE, dimension(:,:)      :: Gamma_Gas

  REALTYPE, dimension(Dimensions), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(NV)        , OPTIONAL :: Lambda
  integer ,                        OPTIONAL :: MPI_Comm_Type

  SELECT CASE (Index_x)
    CASE (0)
      RETURN

    CASE DEFAULT
      call PS_WENO_Euler_Flux_1D2S_1 (Dimensions,                              &
                                      Index_x, Distributed_x, Order_x, dx,     &
                                               N0, N5, N2, N3, NV,             &
                                      Gamma_Formulation, Gamma_Gas,            &
                                      Q, D_Flux, V_Moving_Frame, Lambda,       &
                                      MPI_Comm_Type)
  END SELECT

  END Subroutine PS_WENO_Euler_Flux_1D2S
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_2D2S (Dimensions,                              &
                                      Index_x, Distributed_x, Order_x, dx,     &
                                               N0, N5, N2, N3,                 &
                                      Index_y, Distributed_y, Order_y, dy,     &
                                               M0, M5, M2, M3, NV,             &
                                      Gamma_Formulation, Gamma_Gas,            &
                                      Q, D_Flux, V_Moving_Frame, Lambda,       &
                                      MPI_Comm_Type)

  integer  :: Dimensions

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: Index_y, Order_y
  logical  :: Distributed_y
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, D_Flux

  integer                             :: Gamma_Formulation
  REALTYPE, dimension(:,:)            :: Gamma_Gas

  REALTYPE, dimension(Dimensions), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(NV)        , OPTIONAL :: Lambda
  integer ,                        OPTIONAL :: MPI_Comm_Type

  SELECT CASE (Index_x)
    CASE (0)
      RETURN

    CASE DEFAULT
      call PS_WENO_Euler_Flux_2D2S_1 (Dimensions,                             &
                                      Index_x, Distributed_x, Order_x, dx,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)
    CASE (2)
      call PS_WENO_Euler_Flux_2D2S_2 (Dimensions,                             &
                                      Index_x, Distributed_x, Order_x, dx,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)
  END SELECT

  SELECT CASE (Index_y)
    CASE (0)
      RETURN

    CASE (1)
      call PS_WENO_Euler_Flux_2D2S_1 (Dimensions,                             &
                                      Index_y, Distributed_y, Order_y, dy,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)

    CASE DEFAULT
      call PS_WENO_Euler_Flux_2D2S_2 (Dimensions,                             &
                                      Index_y, Distributed_y, Order_y, dy,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)
  END SELECT

  END Subroutine PS_WENO_Euler_Flux_2D2S
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_3D2S (Dimensions,                              &
                                      Index_x, Distributed_x, Order_x, dx,     &
                                               N0, N5, N2, N3,                 &
                                      Index_y, Distributed_y, Order_y, dy,     &
                                               M0, M5, M2, M3,                 &
                                      Index_z, Distributed_z, Order_z, dz,     &
                                               K0, K5, K2, K3, NV,             &
                                      Gamma_Formulation, Gamma_Gas,            &
                                      Q, D_Flux, V_Moving_Frame, Lambda,       &
                                      MPI_Comm_Type)

  integer  :: Dimensions

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

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Flux

  integer                                   :: Gamma_Formulation
  REALTYPE, dimension(:,:)                  :: Gamma_Gas

  REALTYPE, dimension(Dimensions), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(NV)        , OPTIONAL :: Lambda
  integer ,                        OPTIONAL :: MPI_Comm_Type

  SELECT CASE (Index_x)
    CASE (0)
      RETURN

    CASE DEFAULT
      call PS_WENO_Euler_Flux_3D2S_1 (Dimensions,                             &
                                      Index_x, Distributed_x, Order_x, dx,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3,                &
                                               K0, K5, K2, K3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D2S_2 (Dimensions,                             &
                                      Index_x, Distributed_x, Order_x, dx,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3,                &
                                               K0, K5, K2, K3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D2S_3 (Dimensions,                             &
                                      Index_x, Distributed_x, Order_x, dx,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3,                &
                                               K0, K5, K2, K3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)
  END SELECT

  SELECT CASE (Index_y)
    CASE (0)
      RETURN

    CASE (1)
      call PS_WENO_Euler_Flux_3D2S_1 (Dimensions,                             &
                                      Index_y, Distributed_y, Order_y, dy,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3,                &
                                               K0, K5, K2, K3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)

    CASE DEFAULT
      call PS_WENO_Euler_Flux_3D2S_2 (Dimensions,                             &
                                      Index_y, Distributed_y, Order_y, dy,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3,                &
                                               K0, K5, K2, K3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D2S_3 (Dimensions,                             &
                                      Index_y, Distributed_y, Order_y, dy,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3,                &
                                               K0, K5, K2, K3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)
  END SELECT

  SELECT CASE (Index_z)
    CASE (0)
      RETURN

    CASE (1)
      call PS_WENO_Euler_Flux_3D2S_1 (Dimensions,                             &
                                      Index_z, Distributed_z, Order_z, dz,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3,                &
                                               K0, K5, K2, K3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D2S_2 (Dimensions,                             &
                                      Index_z, Distributed_z, Order_z, dz,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3,                &
                                               K0, K5, K2, K3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)

    CASE DEFAULT
      call PS_WENO_Euler_Flux_3D2S_3 (Dimensions,                             &
                                      Index_z, Distributed_z, Order_z, dz,    &
                                               N0, N5, N2, N3,                &
                                               M0, M5, M2, M3,                &
                                               K0, K5, K2, K3, NV,            &
                                      Gamma_Formulation, Gamma_Gas,           &
                                      Q, D_Flux, V_Moving_Frame, Lambda,      &
                                      MPI_Comm_Type)
  END SELECT

  END Subroutine PS_WENO_Euler_Flux_3D2S 
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_1D2S_1 (Dimensions,                            &
                                        Index, Distributed, Order, dx,         &
                                        N0, N5, i_0, i_1, NV,                  &
                                        Gamma_Formulation, Gamma_Gas,          & 
                                        Q, D_Flux, V_Moving_Frame, Lambda_Max, &
                                        MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: NV
  REALTYPE :: dx

  REALTYPE :: Cp_1, Cv_1
  REALTYPE :: Cp_2, Cv_2

  REALTYPE :: N_x
  REALTYPE :: M_x, UdotN

  integer                          :: Gamma_Formulation
  REALTYPE, dimension(:,:)         :: Gamma_Gas
  REALTYPE, dimension(N0:N5)       :: Gamma

  REALTYPE, dimension(NV)        , OPTIONAL :: Lambda_Max
  integer ,                        OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions)           :: S

  REALTYPE, dimension(N0:N5,NV)    :: Q, D_Flux

  REALTYPE, dimension(N0:N5,NV)    :: F, F_Half
  REALTYPE, dimension(N0:N5,NV,NV) :: E_L, E_R

  REALTYPE, dimension(N0:N5) :: U, P, C, r, h 
  REALTYPE, dimension(N0:N5) :: f1, f2
  REALTYPE, dimension(N0:N5) :: U_dot_N

  REALTYPE, dimension(NV)    :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(NV) :: Lambda_i
#endif

  integer  :: i
  REALTYPE :: s0, s1, b1, b2, b3, b4, c1, u1, v1, w1
  REALTYPE :: Um, Vm, Wm, hm, f1m, f2m, gm, gm1, gmp, qm, cm, em, Xm, Ym

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  N_x = ONE

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_Max)) then
    if (Distributed) call Eigenvalues (Q)
  endif
#endif

!$OMP PARALLEL DO 
  do i = N0, N5
     U(i) = Q(i,2)/Q(i,1)
    f1(i) = Q(i,4)/Q(i,1)

    f1(i) = MIN(MAX(ZERO,f1(i)),ONE)
    f2(i) = ONE - f1(i)

    Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
    Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

    Gamma(i) = (f1(i)*Cp_1 + f2(i)*Cp_2)/(f1(i)*Cv_1 + f2(i)*Cv_2)

    P(i) = (Gamma(i)-ONE)*(Q(i,3) - HALF*Q(i,1)*(U(i)**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C(i) = SQRT(ABS(Gamma(i)*P(i)/Q(i,1)))
#else
                           C(i) = SQRT(ABS(Gamma(i)*P(i)/Q(i,1)))
#endif

    r(i) = SQRT(ABS(Q(i,1)))
    h(i) = (Q(i,3)+P(i))/Q(i,1)

    U_dot_N(i) = U(i)*N_x

    F(i,1) =   (Q(i,1)         * U(i)       )
    F(i,2) =   (Q(i,2)         * U(i) + P(i))
    F(i,3) =   (Q(i,3) + P(i)) * U(i)
    F(i,4) =   (Q(i,4)         * U(i)       )
  enddo
!$OMP END PARALLEL DO

  if (S(1) /= ZERO) then
    do i = N0, N5
      F(i,1) = F(i,1) - S(1)*Q(i,1)
      F(i,2) = F(i,2) - S(1)*Q(i,2)
      F(i,3) = F(i,3) - S(1)*Q(i,3)
      F(i,4) = F(i,4) - S(1)*Q(i,4)
    enddo
  endif

  if (.NOT. PRESENT(Lambda_Max)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)+C(i))))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)-C(i))))
      Lambda(3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)     )))
      Lambda(4) = Lambda(3)
#if defined (PARALLEL_MPI)
    endif
#endif
  else
    Lambda = Lambda_Max
  endif

  Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, Wm, hm)           &
!$OMP             PRIVATE (f1m, f2m, gm, gm1, gmp)           &
!$OMP             PRIVATE (qm, cm, em, Xm, Ym, UdotN)        &
!$OMP             PRIVATE (b1, b2, b3, b4, c1, u1, v1, w1)
    do i = i_0-1, i_1

      s0 = r(i)/(r(i) + r(i+1)) ; s1 = ONE-s0

       Um = s0* U(i) + s1* U(i+1)
       hm = s0* h(i) + s1* h(i+1)
      f1m = s0*f1(i) + s1*f1(i+1)

      f1m = MIN(MAX(ZERO,f1m),ONE)
      f2m = ONE - f1m

!     gm  = f1m*Gamma_Gas(1,1) + f2m*Gamma_Gas(2,1)
!     gmp = Gamma_Gas(1,1) - Gamma_Gas(2,1)

      gm  = (f1m*Cp_1 + f2m*Cp_2)/(f1m*Cv_1 + f2m*Cv_2)
      gmp = (Cp_1*Cv_2-Cp_2*Cv_1)/(f1m*Cv_1 + f2m*Cv_2)**2

      gm1 = gm-ONE

      qm = HALF*(Um*Um)
      cm = SQRT(ABS(gm1*(hm - qm)))
      em = cm**2/(gm*gm1)

      Xm  = gmp*em ; Ym = Xm/gm1

      UdotN = Um*N_x

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN 

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = ZERO

      E_R(i,2,1) =  Um     + cm*N_x
      E_R(i,2,2) =  Um     - cm*N_x
      E_R(i,2,3) = ONE
      E_R(i,2,4) = ZERO

      E_R(i,3,1) =  hm     + cm*UdotN
      E_R(i,3,2) =  hm     - cm*UdotN
      E_R(i,3,3) = -qm*M_x + Um         
      E_R(i,3,4) =         - Ym

      E_R(i,4,1) = f1m
      E_R(i,4,2) = f1m
      E_R(i,4,3) = f1m*M_x
      E_R(i,4,4) = ONE

      b1 = HALF*gm1/cm**2 ; b2 = qm*b1 ; b3 = Ym*b1 ; b4 = b2 - f1m*b3
      c1 = HALF/cm        ; u1 = Um*b1 

      E_L(i,1,1) =     - c1*UdotN  + b4
      E_L(i,1,2) =       c1*N_x    - u1
      E_L(i,1,3) =                   b1
      E_L(i,1,4) =                   b3

      E_L(i,2,1) =       c1*UdotN  + b4
      E_L(i,2,2) =     - c1*N_x    - u1
      E_L(i,2,3) =                   b1
      E_L(i,2,4) =                   b3

      E_L(i,3,1) =  UdotN*N_x-TWO*Um*b4
      E_L(i,3,2) =            TWO*Um*u1
      E_L(i,3,3) =           -TWO*Um*b1
      E_L(i,3,4) =           -TWO*Um*b3

      E_L(i,4,1) = -f1m
      E_L(i,4,2) =  ZERO
      E_L(i,4,3) =  ZERO
      E_L(i,4,4) =  ONE

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler   (Order, N0, N5, i_0-1, i_1, NV,            &
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

    Subroutine Eigenvalues (Q)

    REALTYPE, dimension(N0:N5,NV) :: Q
    REALTYPE, dimension(NV)       :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO
    do i = N0, N5
       U(i) = Q(i,2)/Q(i,1)
      f1(i) = Q(i,4)/Q(i,1)

      f1(i) = MIN(MAX(ZERO,f1(i)),ONE)
      f2(i) = ONE - f1(i)

      Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
      Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

      Gamma(i) = (f1(i)*Cp_1 + f2(i)*Cp_2)/(f1(i)*Cv_1 + f2(i)*Cv_2)

       P(i) = (Gamma(i)-ONE)*(Q(i,3) - HALF*Q(i,1)*(U(i)**2))

       C(i) = SQRT(ABS(Gamma(i)*P(i)/Q(i,1)))

       U_dot_N(i) = U(i)*N_x
    enddo
!$OMP END PARALLEL DO

    Lambda_Local(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)+C)))
    Lambda_Local(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)-C)))
    Lambda_Local(3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)  )))
    Lambda_Local(4) = Lambda_Local(3)

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

  END Subroutine PS_WENO_Euler_Flux_1D2S_1
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_2D2S_1 (Dimensions,                            &
                                        Index, Distributed, Order, dx,         &
                                        N0, N5, i_0, i_1,                      &
                                        M0, M5, j_0, j_1, NV,                  &
                                        Gamma_Formulation, Gamma_Gas,          &
                                        Q, D_Flux, V_Moving_Frame, Lambda_Max, &
                                        MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: NV
  REALTYPE :: dx

  REALTYPE :: Cp_1, Cv_1
  REALTYPE :: Cp_2, Cv_2

  REALTYPE :: N_x, N_y
  REALTYPE :: M_x, M_y, UdotN

  integer                          :: Gamma_Formulation
  REALTYPE, dimension(:,:)         :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE, dimension(NV)        , OPTIONAL :: Lambda_Max
  integer ,                        OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions)           :: S

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, D_Flux

  REALTYPE, dimension(N0:N5,NV)       :: F, F_Half
  REALTYPE, dimension(N0:N5,NV,NV)    :: E_L, E_R

  REALTYPE, dimension(N0:N5) :: U, V, P, C, r, h 
  REALTYPE, dimension(N0:N5) :: f1, f2
  REALTYPE, dimension(N0:N5) :: U_dot_N

  REALTYPE, dimension(NV)    :: Lambda

#if defined (PARALLEL_MPI)
  REALTYPE, dimension(j_0:j_1,NV) :: Lambda_i
#endif

  integer  :: i, j
  REALTYPE :: s0, s1, b1, b2, b3, b4, c1, u1, v1, w1
  REALTYPE :: Um, Vm, Wm, hm, f1m, f2m, gm, gm1, gmp, qm, cm, em, Xm, Ym

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  N_x = ONE  ; N_y = ZERO 

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_Max)) then
    if (Distributed) call Eigenvalues (Q)
  endif
#endif

!--------------------- begin of outer loop in y-dir. -------------------

!$OMP PARALLEL DO PRIVATE (U, V, P, C, r, h, F, F_Half, Lambda, MR_Flag_x)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)
  do j = j_0, j_1

    U  = Q(:,j,2)/Q(:,j,1)
    V  = Q(:,j,3)/Q(:,j,1)
    f1 = Q(:,j,5)/Q(:,j,1)

    f1 = MIN(MAX(ZERO,f1),ONE)
    f2 = ONE - f1

    Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
    Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

    Gamma(:,j) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

    P  = (Gamma(:,j)-ONE)*(Q(:,j,4) - HALF*Q(:,j,1)*(U**2+V**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma(:,j)*P/Q(:,j,1)))
#else
                           C = SQRT(ABS(Gamma(:,j)*P/Q(:,j,1)))
#endif

    r  = SQRT(ABS(Q(:,j,1)))
    h  = (Q(:,j,4)+P)/Q(:,j,1)

    F(:,1) =   (Q(:,j,1)      * U    )
    F(:,2) =   (Q(:,j,2)      * U + P)
    F(:,3) =   (Q(:,j,3)      * U    )
    F(:,4) =   (Q(:,j,4) + P) * U
    F(:,5) =   (Q(:,j,5)      * U    )

    if (S(1) /= ZERO) then
      F(:,1) = F(:,1) - S(1)*Q(:,j,1)
      F(:,2) = F(:,2) - S(1)*Q(:,j,2)
      F(:,3) = F(:,3) - S(1)*Q(:,j,3)
      F(:,4) = F(:,4) - S(1)*Q(:,j,4)
      F(:,5) = F(:,5) - S(1)*Q(:,j,5)
    endif

    U_dot_N = U*N_x + V*N_y

 if (.NOT. PRESENT(Lambda_Max)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(j,:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)+C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)-C)))
      Lambda(3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)  )))
      Lambda(4) = Lambda(3)
      Lambda(5) = Lambda(3)
#if defined (PARALLEL_MPI)
    endif
#endif
  else
    Lambda = Lambda_Max
  endif

    Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, Wm, hm)           &
!$OMP             PRIVATE (f1m, f2m, gm, gm1, gmp)           &
!$OMP             PRIVATE (qm, cm, em, Xm, Ym, UdotN)        &
!$OMP             PRIVATE (b1, b2, b3, b4, c1, u1, v1, w1)
    do i = i_0-1, i_1

      s0 = r(i)/(r(i) + r(i+1)) ; s1 = ONE-s0

       Um = s0* U(i) + s1* U(i+1)
       Vm = s0* V(i) + s1* V(i+1)
       hm = s0* h(i) + s1* h(i+1)
      f1m = s0*f1(i) + s1*f1(i+1)

      f1m = MIN(MAX(ZERO,f1m),ONE)
      f2m = ONE - f1m

!     gm  = f1m*Gamma_Gas(1,1) + f2m*Gamma_Gas(2,1)
!     gmp = Gamma_Gas(1,1) - Gamma_Gas(2,1)

      gm  = (f1m*Cp_1 + f2m*Cp_2)/(f1m*Cv_1 + f2m*Cv_2)
      gmp = (Cp_1*Cv_2-Cp_2*Cv_1)/(f1m*Cv_1 + f2m*Cv_2)**2

      gm1 = gm-ONE

      qm = HALF*(Um*Um + Vm*Vm)
      cm = SQRT(ABS(gm1*(hm - qm)))
      em = cm**2/(gm*gm1)

      Xm  = gmp*em ; Ym = Xm/gm1

      UdotN = Um*N_x + Vm*N_y

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = M_y
      E_R(i,1,5) = ZERO

      E_R(i,2,1) =  Um     + cm*N_x
      E_R(i,2,2) =  Um     - cm*N_x
      E_R(i,2,3) = ONE                  ! UdotN
      E_R(i,2,4) = ZERO
      E_R(i,2,5) = ZERO

      E_R(i,3,1) =  Vm     + cm*N_y
      E_R(i,3,2) =  Vm     - cm*N_y
      E_R(i,3,3) = ZERO
      E_R(i,3,4) = ONE                  ! UdotN
      E_R(i,3,5) = ZERO

      E_R(i,4,1) =  hm     + cm*UdotN
      E_R(i,4,2) =  hm     - cm*UdotN
      E_R(i,4,3) = -qm*M_x + Um         ! -qm*N_x + Um*UdotN 
      E_R(i,4,4) = -qm*M_y + Vm         ! -qm*N_y + Vm*UdotN 
      E_R(i,4,5) =         - Ym

      E_R(i,5,1) = f1m
      E_R(i,5,2) = f1m
      E_R(i,5,3) = f1m*M_x              ! f1m*N_x
      E_R(i,5,4) = f1m*M_y              ! f1m*N_y
      E_R(i,5,5) = ONE

      b1 = HALF*gm1/cm**2 ; b2 = qm*b1 ; b3 = Ym*b1 ; b4 = b2 - f1m*b3
      c1 = HALF/cm        ; u1 = Um*b1 ; v1 = Vm*b1

      E_L(i,1,1) =     - c1*UdotN  + b4
      E_L(i,1,2) =       c1*N_x    - u1
      E_L(i,1,3) =       c1*N_y    - v1
      E_L(i,1,4) =                   b1
      E_L(i,1,5) =                   b3

      E_L(i,2,1) =       c1*UdotN  + b4
      E_L(i,2,2) =     - c1*N_x    - u1
      E_L(i,2,3) =     - c1*N_y    - v1
      E_L(i,2,4) =                   b1
      E_L(i,2,5) =                   b3

      E_L(i,3,1) =  UdotN*N_x-TWO*Um*b4
      E_L(i,3,2) =            TWO*Um*u1 - N_x*N_x + ONE 
      E_L(i,3,3) =            TWO*Um*v1 - N_x*N_y
      E_L(i,3,4) =           -TWO*Um*b1
      E_L(i,3,5) =           -TWO*Um*b3

      E_L(i,4,1) =  UdotN*N_y-TWO*Vm*b4
      E_L(i,4,2) =            TWO*Vm*u1 - N_y*N_x
      E_L(i,4,3) =            TWO*Vm*v1 - N_y*N_y + ONE
      E_L(i,4,4) =           -TWO*Vm*b1
      E_L(i,4,5) =           -TWO*Vm*b3

      E_L(i,5,1) = -f1m
      E_L(i,5,2) =  ZERO
      E_L(i,5,3) =  ZERO
      E_L(i,5,4) =  ZERO
      E_L(i,5,5) =  ONE

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler   (Order, N0, N5, i_0-1, i_1, NV,            &
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

    Subroutine Eigenvalues (Q)

    REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
    REALTYPE, dimension(j_0:j_1,NV)     :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO PRIVATE (U, V, P, C, f1, f2, Cp_1, Cp_2, Cv_1, Cv_2) 
    do j = j_0, j_1

      U  = Q(:,j,2)/Q(:,j,1)
      V  = Q(:,j,3)/Q(:,j,1)
      f1 = Q(:,j,5)/Q(:,j,1)

      f1 = MIN(MAX(ZERO,f1),ONE)
      f2 = ONE - f1

      Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
      Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

      Gamma(:,j) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

      P = (Gamma(:,j)-ONE)*(Q(:,j,4) - HALF*Q(:,j,1)*(U**2+V**2))

      C = SQRT(ABS(Gamma(:,j)*P/Q(:,j,1)))

      U_dot_N = U*N_x + V*N_y

      Lambda_Local(j,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)+C)))
      Lambda_Local(j,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)-C)))
      Lambda_Local(j,3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)  )))
      Lambda_Local(j,4) = Lambda_Local(j,3)
      Lambda_Local(j,5) = Lambda_Local(j,3)

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

  END Subroutine PS_WENO_Euler_Flux_2D2S_1
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_2D2S_2 (Dimensions,                            &
                                        Index, Distributed, Order, dy,         &
                                        N0, N5, i_0, i_1,                      &
                                        M0, M5, j_0, j_1, NV,                  &
                                        Gamma_Formulation, Gamma_Gas,          &
                                        Q, D_Flux, V_Moving_Frame, Lambda_Max, &
                                        MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: NV
  REALTYPE :: dy

  REALTYPE :: Cp_1, Cv_1
  REALTYPE :: Cp_2, Cv_2

  REALTYPE :: N_x, N_y
  REALTYPE :: M_x, M_y, UdotN

  integer                          :: Gamma_Formulation
  REALTYPE, dimension(:,:)         :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE, dimension(NV)        , OPTIONAL :: Lambda_Max
  integer ,                        OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions)           :: S

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, D_Flux

  REALTYPE, dimension(M0:M5,NV)       :: F, F_Half
  REALTYPE, dimension(M0:M5,NV,NV)    :: E_L, E_R

  REALTYPE, dimension(M0:M5) :: U, V, P, C, r, h 
  REALTYPE, dimension(M0:M5) :: f1, f2
  REALTYPE, dimension(M0:M5) :: U_dot_N

  REALTYPE, dimension(NV)    :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(i_0:i_1,NV) :: Lambda_i
#endif

  integer  :: i, j
  REALTYPE :: s0, s1, b1, b2, b3, b4, c1, u1, v1, w1
  REALTYPE :: Um, Vm, Wm, hm, f1m, f2m, gm, gm1, gmp, qm, cm, em, Xm, Ym

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  N_x = ZERO ; N_y = ONE

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_Max)) then
    if (Distributed) call Eigenvalues (Q)
  endif
#endif

!--------------------- begin of outer loop in x-dir. -------------------

!$OMP PARALLEL DO PRIVATE (U, V, P, C, r, h, F, F_Half, Lambda, MR_Flag_y)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)
  do i = i_0, i_1

    U  = Q(i,:,2)/Q(i,:,1)
    V  = Q(i,:,3)/Q(i,:,1)
    f1 = Q(i,:,5)/Q(i,:,1)

    f1 = MIN(MAX(ZERO,f1),ONE)
    f2 = ONE - f1

    Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
    Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

    Gamma(i,:) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

    P  = (Gamma(i,:)-ONE)*(Q(i,:,4) - HALF*Q(i,:,1)*(U**2+V**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma(i,:)*P/Q(i,:,1)))
#else
                           C = SQRT(ABS(Gamma(i,:)*P/Q(i,:,1)))
#endif

    r  = SQRT(ABS(Q(i,:,1)))
    h  = (Q(i,:,4)+P)/Q(i,:,1)

    F(:,1) =  (Q(i,:,1)      * V    )
    F(:,2) =  (Q(i,:,2)      * V    )
    F(:,3) =  (Q(i,:,3)      * V + P)
    F(:,4) =  (Q(i,:,4) + P) * V
    F(:,5) =  (Q(i,:,5)      * V    )

    if (S(2) /= ZERO) then
      F(:,1) = F(:,1) - S(2)*Q(i,:,1)
      F(:,2) = F(:,2) - S(2)*Q(i,:,2)
      F(:,3) = F(:,3) - S(2)*Q(i,:,3)
      F(:,4) = F(:,4) - S(2)*Q(i,:,4)
      F(:,5) = F(:,5) - S(2)*Q(i,:,5)
    endif

    U_dot_N = U*N_x + V*N_y

  if (.NOT. PRESENT(Lambda_Max)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(i,:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)+C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)-C)))
      Lambda(3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)  )))
      Lambda(4) = Lambda(3)
      Lambda(5) = Lambda(3)
#if defined (PARALLEL_MPI)
    endif
#endif
  else
    Lambda = Lambda_Max
  endif

    Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, Wm, hm)           &
!$OMP             PRIVATE (f1m, f2m, gm, gm1, gmp)           &
!$OMP             PRIVATE (qm, cm, em, Xm, Ym, UdotN)        &
!$OMP             PRIVATE (b1, b2, b3, b4, c1, u1, v1, w1)
    do j = j_0-1, j_1

      s0 = r(j)/(r(j) + r(j+1)) ; s1 = ONE-s0

       Um = s0* U(j) + s1* U(j+1)
       Vm = s0* V(j) + s1* V(j+1)
       hm = s0* h(j) + s1* h(j+1)
      f1m = s0*f1(j) + s1*f1(j+1)

      f1m = MIN(MAX(ZERO,f1m),ONE)
      f2m = ONE - f1m

!     gm  = f1m*Gamma_Gas(1,1) + f2m*Gamma_Gas(2,1)
!     gmp = Gamma_Gas(1,1) - Gamma_Gas(2,1)

      gm  = (f1m*Cp_1 + f2m*Cp_2)/(f1m*Cv_1 + f2m*Cv_2)
      gmp = (Cp_1*Cv_2-Cp_2*Cv_1)/(f1m*Cv_1 + f2m*Cv_2)**2

      gm1 = gm-ONE

      qm = HALF*(Um*Um + Vm*Vm)
      cm = SQRT(ABS(gm1*(hm - qm)))
      em = cm**2/(gm*gm1)

      Xm  = gmp*em ; Ym = Xm/gm1

      UdotN = Um*N_x + Vm*N_y

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN

      E_R(j,1,1) = ONE
      E_R(j,1,2) = ONE
      E_R(j,1,3) = M_x
      E_R(j,1,4) = M_y
      E_R(j,1,5) = ZERO

      E_R(j,2,1) =  Um     + cm*N_x
      E_R(j,2,2) =  Um     - cm*N_x
      E_R(j,2,3) = ONE
      E_R(j,2,4) = ZERO
      E_R(j,2,5) = ZERO

      E_R(j,3,1) =  Vm     + cm*N_y
      E_R(j,3,2) =  Vm     - cm*N_y
      E_R(j,3,3) = ZERO
      E_R(j,3,4) = ONE
      E_R(j,3,5) = ZERO

      E_R(j,4,1) =  hm     + cm*UdotN
      E_R(j,4,2) =  hm     - cm*UdotN
      E_R(j,4,3) = -qm*M_x + Um
      E_R(j,4,4) = -qm*M_y + Vm
      E_R(j,4,5) =         - Ym

      E_R(j,5,1) = f1m
      E_R(j,5,2) = f1m
      E_R(j,5,3) = f1m*M_x
      E_R(j,5,4) = f1m*M_y
      E_R(j,5,5) = ONE

      b1 = HALF*gm1/cm**2 ; b2 = qm*b1 ; b3 = Ym*b1 ; b4 = b2 - f1m*b3
      c1 = HALF/cm        ; u1 = Um*b1 ; v1 = Vm*b1

      E_L(j,1,1) =     - c1*UdotN  + b4
      E_L(j,1,2) =       c1*N_x    - u1
      E_L(j,1,3) =       c1*N_y    - v1
      E_L(j,1,4) =                   b1
      E_L(j,1,5) =                   b3

      E_L(j,2,1) =       c1*UdotN  + b4
      E_L(j,2,2) =     - c1*N_x    - u1
      E_L(j,2,3) =     - c1*N_y    - v1
      E_L(j,2,4) =                   b1
      E_L(j,2,5) =                   b3

      E_L(j,3,1) =  UdotN*N_x-TWO*Um*b4
      E_L(j,3,2) =            TWO*Um*u1 - N_x*N_x + ONE
      E_L(j,3,3) =            TWO*Um*v1 - N_x*N_y
      E_L(j,3,4) =           -TWO*Um*b1
      E_L(j,3,5) =           -TWO*Um*b3

      E_L(j,4,1) =  UdotN*N_y-TWO*Vm*b4
      E_L(j,4,2) =            TWO*Vm*u1 - N_y*N_x
      E_L(j,4,3) =            TWO*Vm*v1 - N_y*N_y + ONE
      E_L(j,4,4) =           -TWO*Vm*b1
      E_L(j,4,5) =           -TWO*Vm*b3

      E_L(j,5,1) = -f1m
      E_L(j,5,2) =  ZERO
      E_L(j,5,3) =  ZERO
      E_L(j,5,4) =  ZERO
      E_L(j,5,5) =  ONE

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler   (Order, M0, M5, j_0-1, j_1, NV,            &
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

    Subroutine Eigenvalues (Q)

    REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
    REALTYPE, dimension(i_0:i_1,NV)     :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO PRIVATE (U, V, P, C, f1, f2, Cp_1, Cp_2, Cv_1, Cv_2) 
    do i = i_0, i_1

      U  = Q(i,:,2)/Q(i,:,1)
      V  = Q(i,:,3)/Q(i,:,1)
      f1 = Q(i,:,5)/Q(i,:,1)

      f1 = MIN(MAX(ZERO,f1),ONE)
      f2 = ONE - f1

      Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
      Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

      Gamma(i,:) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

      P = (Gamma(i,:)-ONE)*(Q(i,:,4) - HALF*Q(i,:,1)*(U**2+V**2))

      C = SQRT(ABS(Gamma(i,:)*P/Q(i,:,1)))

      U_dot_N = U*N_x + V*N_y

      Lambda_Local(i,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)+C)))
      Lambda_Local(i,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)-C)))
      Lambda_Local(i,3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)  )))
      Lambda_Local(i,4) = Lambda_Local(i,3)
      Lambda_Local(i,5) = Lambda_Local(i,3)

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

  END Subroutine PS_WENO_Euler_Flux_2D2S_2 
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_3D2S_1 (Dimensions,                            &
                                        Index, Distributed, Order, dx,         &
                                        N0, N5, i_0, i_1,                      &
                                        M0, M5, j_0, j_1,                      &
                                        K0, K5, k_0, k_1, NV,                  &
                                        Gamma_Formulation, Gamma_Gas,          &
                                        Q, D_Flux, V_Moving_Frame, Lambda_Max, &
                                        MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1
  integer  :: NV
  REALTYPE :: dx

  REALTYPE :: Cp_1, Cv_1
  REALTYPE :: Cp_2, Cv_2

  REALTYPE :: N_x, N_y, N_z
  REALTYPE :: M_x, M_y, M_z, UdotN

  integer                                :: Gamma_Formulation
  REALTYPE, dimension(:,:)               :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE, dimension(NV)        , OPTIONAL :: Lambda_Max
  integer ,                        OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions)           :: S

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Flux

  REALTYPE, dimension(N0:N5,NV)    :: F, F_Half
  REALTYPE, dimension(N0:N5,NV,NV) :: E_L, E_R

  REALTYPE, dimension(N0:N5) :: U, V, W, P, C, r, h 
  REALTYPE, dimension(N0:N5) :: f1, f2
  REALTYPE, dimension(N0:N5) :: U_dot_N

  REALTYPE, dimension(NV)    :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(j_0:j_1,k_0:k_1,NV) :: Lambda_i
#endif

  integer  :: i, j, k
  REALTYPE :: s0, s1, b1, b2, b3, b4, c1, u1, v1, w1
  REALTYPE :: Um, Vm, Wm, hm, f1m, f2m, gm, gm1, gmp, qm, cm, em, Xm, Ym

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  N_x = ONE  ; N_y = ZERO ; N_z = ZERO

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_Max)) then
    if (Distributed) call Eigenvalues (Q)
  endif
#endif

!--------------------- begin of outer loop in y-dir and z-dir ------------

!$OMP PARALLEL DO  
  do k = k_0, k_1
!$OMP PARALLEL DO PRIVATE (U, V, W, P, C, r, h, F, F_Half, Lambda, MR_Flag_x)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)                                
  do j = j_0, j_1

    U  = Q(:,j,k,2)/Q(:,j,k,1)
    V  = Q(:,j,k,3)/Q(:,j,k,1)
    W  = Q(:,j,k,4)/Q(:,j,k,1)
    f1 = Q(:,j,k,6)/Q(:,j,k,1)

    f1 = MIN(MAX(ZERO,f1),ONE)
    f2 = ONE - f1

    Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
    Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

    Gamma(:,j,k) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

    P = (Gamma(:,j,k)-ONE)*(Q(:,j,k,5) - HALF*Q(:,j,k,1)*(U**2+V**2+W**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma(:,j,k)*P/Q(:,j,k,1)))
#else
                           C = SQRT(ABS(Gamma(:,j,k)*P/Q(:,j,k,1)))
#endif

    r = SQRT(ABS(Q(:,j,k,1)))
    h = (Q(:,j,k,5)+P)/Q(:,j,k,1)

    F(:,1) =   (Q(:,j,k,1)      * U    )
    F(:,2) =   (Q(:,j,k,2)      * U + P)
    F(:,3) =   (Q(:,j,k,3)      * U    )
    F(:,4) =   (Q(:,j,k,4)      * U    )
    F(:,5) =   (Q(:,j,k,5) + P) * U
    F(:,6) =   (Q(:,j,k,6)      * U    )

    if (S(1) /= ZERO) then
      F(:,1) = F(:,1) - S(1)*Q(:,j,k,1)
      F(:,2) = F(:,2) - S(1)*Q(:,j,k,2)
      F(:,3) = F(:,3) - S(1)*Q(:,j,k,3)
      F(:,4) = F(:,4) - S(1)*Q(:,j,k,4)
      F(:,5) = F(:,5) - S(1)*Q(:,j,k,5)
      F(:,6) = F(:,6) - S(1)*Q(:,j,k,6)
    endif

    U_dot_N = U*N_x + V*N_y + W*N_z

  if (.NOT. PRESENT(Lambda_Max)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(j,k,:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)+C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)-C)))
      Lambda(3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)  )))
      Lambda(4) = Lambda(3)
      Lambda(5) = Lambda(3)
      Lambda(6) = Lambda(3)
#if defined (PARALLEL_MPI)
   endif
#endif
  else
    Lambda = Lambda_Max
  endif

  Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, Wm, hm)           &
!$OMP             PRIVATE (f1m, f2m, gm, gm1, gmp)           &
!$OMP             PRIVATE (qm, cm, em, Xm, Ym, UdotN)        &
!$OMP             PRIVATE (b1, b2, b3, b4, c1, u1, v1, w1)
    do i = i_0-1, i_1

      s0 = r(i)/(r(i) + r(i+1)) ; s1 = ONE-s0

       Um = s0* U(i) + s1* U(i+1)
       Vm = s0* V(i) + s1* V(i+1)
       Wm = s0* W(i) + s1* W(i+1)
       hm = s0* h(i) + s1* h(i+1)
      f1m = s0*f1(i) + s1*f1(i+1)

      f1m = MIN(MAX(ZERO,f1m),ONE)
      f2m = ONE - f1m

!     gm  = f1m*Gamma_Gas(1,1) + f2m*Gamma_Gas(2,1)
!     gmp = Gamma_Gas(1,1) - Gamma_Gas(2,1)

      gm  = (f1m*Cp_1 + f2m*Cp_2)/(f1m*Cv_1 + f2m*Cv_2)
      gmp = (Cp_1*Cv_2-Cp_2*Cv_1)/(f1m*Cv_1 + f2m*Cv_2)**2

      gm1 = gm-ONE

      qm = HALF*(Um*Um + Vm*Vm + Wm*Wm)
      cm = SQRT(ABS(gm1*(hm - qm)))
      em = cm**2/(gm*gm1)

      Xm  = gmp*em ; Ym = Xm/gm1

      UdotN = Um*N_x + Vm*N_y + Wm*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = M_y
      E_R(i,1,5) = M_z
      E_R(i,1,6) = ZERO

      E_R(i,2,1) =  Um     + cm*N_x
      E_R(i,2,2) =  Um     - cm*N_x
      E_R(i,2,3) = ONE
      E_R(i,2,4) = ZERO
      E_R(i,2,5) = ZERO
      E_R(i,2,6) = ZERO

      E_R(i,3,1) =  Vm     + cm*N_y
      E_R(i,3,2) =  Vm     - cm*N_y
      E_R(i,3,3) = ZERO
      E_R(i,3,4) = ONE
      E_R(i,3,5) = ZERO
      E_R(i,3,6) = ZERO

      E_R(i,4,1) =  Wm     + cm*N_z
      E_R(i,4,2) =  Wm     - cm*N_z
      E_R(i,4,3) = ZERO
      E_R(i,4,4) = ZERO
      E_R(i,4,5) = ONE
      E_R(i,4,6) = ZERO

      E_R(i,5,1) =  hm     + cm*UdotN
      E_R(i,5,2) =  hm     - cm*UdotN
      E_R(i,5,3) = -qm*M_x + Um
      E_R(i,5,4) = -qm*M_y + Vm
      E_R(i,5,5) = -qm*M_z + Wm
      E_R(i,5,6) =         - Ym

      E_R(i,6,1) = f1m
      E_R(i,6,2) = f1m
      E_R(i,6,3) = f1m*M_x
      E_R(i,6,4) = f1m*M_y
      E_R(i,6,5) = f1m*M_z
      E_R(i,6,6) = ONE

      b1 = HALF*gm1/cm**2 ; b2 = qm*b1 ; b3 = Ym*b1 ; b4 = b2 - f1m*b3
      c1 = HALF/cm        ; u1 = Um*b1 ; v1 = Vm*b1 ; w1 = Wm*b1

      E_L(i,1,1) =     - c1*UdotN  + b4
      E_L(i,1,2) =       c1*N_x    - u1
      E_L(i,1,3) =       c1*N_y    - v1
      E_L(i,1,4) =       c1*N_z    - w1
      E_L(i,1,5) =                   b1
      E_L(i,1,6) =                   b3

      E_L(i,2,1) =       c1*UdotN  + b4
      E_L(i,2,2) =     - c1*N_x    - u1
      E_L(i,2,3) =     - c1*N_y    - v1
      E_L(i,2,4) =     - c1*N_z    - w1
      E_L(i,2,5) =                   b1
      E_L(i,2,6) =                   b3

      E_L(i,3,1) =  UdotN*N_x-TWO*Um*b4
      E_L(i,3,2) =            TWO*Um*u1 - N_x*N_x + ONE
      E_L(i,3,3) =            TWO*Um*v1 - N_x*N_y
      E_L(i,3,4) =            TWO*Um*w1 - N_x*N_z
      E_L(i,3,5) =           -TWO*Um*b1
      E_L(i,3,6) =           -TWO*Um*b3

      E_L(i,4,1) =  UdotN*N_y-TWO*Vm*b4
      E_L(i,4,2) =            TWO*Vm*u1 - N_y*N_x
      E_L(i,4,3) =            TWO*Vm*v1 - N_y*N_y + ONE
      E_L(i,4,4) =            TWO*Vm*w1 - N_y*N_z
      E_L(i,4,5) =           -TWO*Vm*b1
      E_L(i,4,6) =           -TWO*Vm*b3

      E_L(i,5,1) =  UdotN*N_z-TWO*Wm*b4
      E_L(i,5,2) =            TWO*Wm*u1 - N_z*N_x
      E_L(i,5,3) =            TWO*Wm*v1 - N_z*N_y 
      E_L(i,5,4) =            TWO*Wm*w1 - N_z*N_z + ONE
      E_L(i,5,5) =           -TWO*Wm*b1
      E_L(i,5,6) =           -TWO*Wm*b3

      E_L(i,6,1) = -f1m
      E_L(i,6,2) =  ZERO
      E_L(i,6,3) =  ZERO
      E_L(i,6,4) =  ZERO
      E_L(i,6,5) =  ZERO
      E_L(i,6,6) =  ONE

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler   (Order, N0, N5, i_0-1, i_1, NV,            &
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

    Subroutine Eigenvalues (Q)

    REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
    REALTYPE, dimension(j_0:j_1,k_0:k_1,NV)   :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO PRIVATE (U, V, W, P, C, f1, f2, Cp_1, Cp_2, Cv_1, Cv_2)
    do k = k_0, k_1
!$OMP PARALLEL DO
    do j = j_0, j_1
      U  = Q(:,j,k,2)/Q(:,j,k,1)
      V  = Q(:,j,k,3)/Q(:,j,k,1)
      W  = Q(:,j,k,4)/Q(:,j,k,1)
      f1 = Q(:,j,k,6)/Q(:,j,k,1)

      f1 = MIN(MAX(ZERO,f1),ONE)
      f2 = ONE - f1

      Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
      Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

      Gamma(:,j,k) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

      P = (Gamma(:,j,k)-ONE)*(Q(:,j,k,5) - HALF*Q(:,j,k,1)*(U**2+V**2+W**2))

      C = SQRT(ABS(Gamma(:,j,k)*P/Q(:,j,k,1)))

      U_dot_N = U*N_x + V*N_y + W*N_z

      Lambda_Local(j,k,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)+C)))
      Lambda_Local(j,k,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)-C)))
      Lambda_Local(j,k,3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(1)  )))
      Lambda_Local(j,k,4) = Lambda_Local(j,k,3)
      Lambda_Local(j,k,5) = Lambda_Local(j,k,3)
      Lambda_Local(j,k,6) = Lambda_Local(j,k,3)

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

  END Subroutine PS_WENO_Euler_Flux_3D2S_1
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_3D2S_2 (Dimensions,                            &
                                        Index, Distributed, Order, dy,         &
                                        N0, N5, i_0, i_1,                      &
                                        M0, M5, j_0, j_1,                      &
                                        K0, K5, k_0, k_1, NV,                  &
                                        Gamma_Formulation, Gamma_Gas,          &
                                        Q, D_Flux, V_Moving_Frame, Lambda_Max, &
                                        MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1
  integer  :: NV
  REALTYPE :: dy

  REALTYPE :: Cp_1, Cv_1
  REALTYPE :: Cp_2, Cv_2

  REALTYPE :: N_x, N_y, N_z
  REALTYPE :: M_x, M_y, M_z, UdotN

  integer                          :: Gamma_Formulation
  REALTYPE, dimension(:,:)         :: Gamma_Gas

  REALTYPE, dimension(NV)        , OPTIONAL :: Lambda_Max
  integer ,                        OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions)           :: S

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: Gamma

  REALTYPE, dimension(M0:M5,NV)    :: F, F_Half
  REALTYPE, dimension(M0:M5,NV,NV) :: E_L, E_R

  REALTYPE, dimension(M0:M5) :: U, V, W, P, C, r, h
  REALTYPE, dimension(M0:M5) :: f1, f2
  REALTYPE, dimension(M0:M5) :: U_dot_N

  REALTYPE, dimension(NV)    :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(i_0:i_1,k_0:k_1,NV) :: Lambda_i
#endif

  integer  :: i, j, k
  REALTYPE :: s0, s1, b1, b2, b3, b4, c1, u1, v1, w1
  REALTYPE :: Um, Vm, Wm, hm, f1m, f2m, gm, gm1, gmp, qm, cm, em, Xm, Ym

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  N_x = ZERO ; N_y = ONE  ; N_z = ZERO

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_Max)) then
    if (Distributed) call Eigenvalues (Q)
  endif
#endif

!--------------------- begin of outer loop in x-dir and z-dir ------------

!$OMP PARALLEL DO PRIVATE (U, V, W, P, C, r, h, F, F_Half, Lambda, MR_Flag_y)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)
  do k = k_0, k_1
!$OMP PARALLEL DO 
  do i = i_0, i_1

    U  = Q(i,:,k,2)/Q(i,:,k,1)
    V  = Q(i,:,k,3)/Q(i,:,k,1)
    W  = Q(i,:,k,4)/Q(i,:,k,1)
    f1 = Q(i,:,k,6)/Q(i,:,k,1)

    f1 = MIN(MAX(ZERO,f1),ONE)
    f2 = ONE - f1

    Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
    Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

    Gamma(i,:,k) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

    P = (Gamma(i,:,k)-ONE)*(Q(i,:,k,5) - HALF*Q(i,:,k,1)*(U**2+V**2+W**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma(i,:,k)*P/Q(i,:,k,1)))
#else
                           C = SQRT(ABS(Gamma(i,:,k)*P/Q(i,:,k,1)))
#endif

    r = SQRT(ABS(Q(:,j,k,1)))
    r = SQRT(ABS(Q(i,:,k,1)))
    h = (Q(i,:,k,5)+P)/Q(i,:,k,1)

    F(:,1) =   (Q(i,:,k,1)      * V    )
    F(:,2) =   (Q(i,:,k,2)      * V    )
    F(:,3) =   (Q(i,:,k,3)      * V + P)
    F(:,4) =   (Q(i,:,k,4)      * V    )
    F(:,5) =   (Q(i,:,k,5) + P) * V
    F(:,6) =   (Q(i,:,k,6)      * V    )

    if (S(2) /= ZERO) then
      F(:,1) = F(:,1) - S(2)*Q(i,:,k,1)
      F(:,2) = F(:,2) - S(2)*Q(i,:,k,2)
      F(:,3) = F(:,3) - S(2)*Q(i,:,k,3)
      F(:,4) = F(:,4) - S(2)*Q(i,:,k,4)
      F(:,5) = F(:,5) - S(2)*Q(i,:,k,5)
      F(:,6) = F(:,6) - S(2)*Q(i,:,k,6)
    endif

    U_dot_N = U*N_x + V*N_y + W*N_z

  if (.NOT. PRESENT(Lambda_Max)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(i,k,:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)+C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)-C)))
      Lambda(3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)  )))
      Lambda(4) = Lambda(3)
      Lambda(5) = Lambda(3)
      Lambda(6) = Lambda(3)
#if defined (PARALLEL_MPI)
   endif
#endif
  else
    Lambda = Lambda_Max
  endif

    Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, Wm, hm)           &
!$OMP             PRIVATE (f1m, f2m, gm, gm1, gmp)           &
!$OMP             PRIVATE (qm, cm, em, Xm, Ym, UdotN)        &
!$OMP             PRIVATE (b1, b2, b3, b4, c1, u1, v1, w1)
    do j = j_0-1, j_1

      s0 = r(j)/(r(j) + r(j+1)) ; s1 = ONE-s0

       Um = s0* U(j) + s1* U(j+1)
       Vm = s0* V(j) + s1* V(j+1)
       Wm = s0* W(j) + s1* W(j+1)
       hm = s0* h(j) + s1* h(j+1)
      f1m = s0*f1(j) + s1*f1(j+1)

      f1m = MIN(MAX(ZERO,f1m),ONE)
      f2m = ONE - f1m

!     gm  = f1m*Gamma_Gas(1,1) + f2m*Gamma_Gas(2,1)
!     gmp = Gamma_Gas(1,1) - Gamma_Gas(2,1)

      gm  = (f1m*Cp_1 + f2m*Cp_2)/(f1m*Cv_1 + f2m*Cv_2)
      gmp = (Cp_1*Cv_2-Cp_2*Cv_1)/(f1m*Cv_1 + f2m*Cv_2)**2

      gm1 = gm-ONE

      qm = HALF*(Um*Um + Vm*Vm + Wm*Wm)
      cm = SQRT(ABS(gm1*(hm - qm)))
      em = cm**2/(gm*gm1)

      Xm  = gmp*em ; Ym = Xm/gm1

      UdotN = Um*N_x + Vm*N_y + Wm*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN

      E_R(j,1,1) = ONE
      E_R(j,1,2) = ONE
      E_R(j,1,3) = M_x
      E_R(j,1,4) = M_y
      E_R(j,1,5) = M_z
      E_R(j,1,6) = ZERO

      E_R(j,2,1) =  Um     + cm*N_x
      E_R(j,2,2) =  Um     - cm*N_x
      E_R(j,2,3) = ONE
      E_R(j,2,4) = ZERO
      E_R(j,2,5) = ZERO
      E_R(j,2,6) = ZERO

      E_R(j,3,1) =  Vm     + cm*N_y
      E_R(j,3,2) =  Vm     - cm*N_y
      E_R(j,3,3) = ZERO
      E_R(j,3,4) = ONE
      E_R(j,3,5) = ZERO
      E_R(j,3,6) = ZERO

      E_R(j,4,1) =  Wm     + cm*N_z
      E_R(j,4,2) =  Wm     - cm*N_z
      E_R(j,4,3) = ZERO
      E_R(j,4,4) = ZERO
      E_R(j,4,5) = ONE
      E_R(j,4,6) = ZERO

      E_R(j,5,1) =  hm     + cm*UdotN
      E_R(j,5,2) =  hm     - cm*UdotN
      E_R(j,5,3) = -qm*M_x + Um
      E_R(j,5,4) = -qm*M_y + Vm
      E_R(j,5,5) = -qm*M_z + Wm
      E_R(j,5,6) =         - Ym

      E_R(j,6,1) = f1m
      E_R(j,6,2) = f1m
      E_R(j,6,3) = f1m*M_x
      E_R(j,6,4) = f1m*M_y
      E_R(j,6,5) = f1m*M_z
      E_R(j,6,6) = ONE

      b1 = HALF*gm1/cm**2 ; b2 = qm*b1 ; b3 = Ym*b1 ; b4 = b2 - f1m*b3
      c1 = HALF/cm        ; u1 = Um*b1 ; v1 = Vm*b1 ; w1 = Wm*b1

      E_L(j,1,1) =     - c1*UdotN  + b4
      E_L(j,1,2) =       c1*N_x    - u1
      E_L(j,1,3) =       c1*N_y    - v1
      E_L(j,1,4) =       c1*N_z    - w1
      E_L(j,1,5) =                   b1
      E_L(j,1,6) =                   b3

      E_L(j,2,1) =       c1*UdotN  + b4
      E_L(j,2,2) =     - c1*N_x    - u1
      E_L(j,2,3) =     - c1*N_y    - v1
      E_L(j,2,4) =     - c1*N_z    - w1
      E_L(j,2,5) =                   b1
      E_L(j,2,6) =                   b3

      E_L(j,3,1) =  UdotN*N_x-TWO*Um*b4
      E_L(j,3,2) =            TWO*Um*u1 - N_x*N_x + ONE
      E_L(j,3,3) =            TWO*Um*v1 - N_x*N_y
      E_L(j,3,4) =            TWO*Um*w1 - N_x*N_z
      E_L(j,3,5) =           -TWO*Um*b1
      E_L(j,3,6) =           -TWO*Um*b3

      E_L(j,4,1) =  UdotN*N_y-TWO*Vm*b4
      E_L(j,4,2) =            TWO*Vm*u1 - N_y*N_x
      E_L(j,4,3) =            TWO*Vm*v1 - N_y*N_y + ONE
      E_L(j,4,4) =            TWO*Vm*w1 - N_y*N_z
      E_L(j,4,5) =           -TWO*Vm*b1
      E_L(j,4,6) =           -TWO*Vm*b3

      E_L(j,5,1) =  UdotN*N_z-TWO*Wm*b4
      E_L(j,5,2) =            TWO*Wm*u1 - N_z*N_x
      E_L(j,5,3) =            TWO*Wm*v1 - N_z*N_y 
      E_L(j,5,4) =            TWO*Wm*w1 - N_z*N_z + ONE
      E_L(j,5,5) =           -TWO*Wm*b1
      E_L(j,5,6) =           -TWO*Wm*b3

      E_L(j,6,1) = -f1m
      E_L(j,6,2) =  ZERO
      E_L(j,6,3) =  ZERO
      E_L(j,6,4) =  ZERO
      E_L(j,6,5) =  ZERO
      E_L(j,6,6) =  ONE

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler   (Order, M0, M5, j_0-1, j_1, NV,            &
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

    Subroutine Eigenvalues (Q)

    REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
    REALTYPE, dimension(i_0:i_1,k_0:k_1,NV)   :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

!$OMP PARALLEL DO PRIVATE (U, V, W, P, C, f1, f2, Cp_1, Cp_2, Cv_1, Cv_2)
    do k = k_0, k_1
!$OMP PARALLEL DO 
    do i = i_0, i_1

       U = Q(i,:,k,2)/Q(i,:,k,1)
       V = Q(i,:,k,3)/Q(i,:,k,1)
       W = Q(i,:,k,4)/Q(i,:,k,1)
      f1 = Q(i,:,k,6)/Q(i,:,k,1)

      f1 = MIN(MAX(ZERO,f1),ONE)
      f2 = ONE - f1

      Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
      Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

      Gamma(i,:,k) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

      P = (Gamma(i,:,k)-ONE)*(Q(i,:,k,5) - HALF*Q(i,:,k,1)*(U**2+V**2+W**2))

      C = SQRT(ABS(Gamma(i,:,k)*P/Q(i,:,k,1)))

      U_dot_N = U*N_x + V*N_y + W*N_z

      Lambda_Local(i,k,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)+C)))
      Lambda_Local(i,k,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)-C)))
      Lambda_Local(i,k,3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(2)  )))
      Lambda_Local(i,k,4) = Lambda_Local(i,k,3)
      Lambda_Local(i,k,5) = Lambda_Local(i,k,3)
      Lambda_Local(i,k,6) = Lambda_Local(i,k,3)

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

  END Subroutine PS_WENO_Euler_Flux_3D2S_2 
!
!===================================================================
!
  Subroutine PS_WENO_Euler_Flux_3D2S_3 (Dimensions,                            &
                                        Index, Distributed, Order, dz,         &
                                        N0, N5, i_0, i_1,                      &
                                        M0, M5, j_0, j_1,                      &
                                        K0, K5, k_0, k_1, NV,                  &
                                        Gamma_Formulation, Gamma_Gas,          &
                                        Q, D_Flux, V_Moving_Frame, Lambda_Max, &
                                        MPI_Comm_Type)

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1
  integer  :: NV
  REALTYPE :: dz

  REALTYPE :: Cp_1, Cv_1
  REALTYPE :: Cp_2, Cv_2

  REALTYPE :: N_x, N_y, N_z
  REALTYPE :: M_x, M_y, M_z, UdotN

  integer                                :: Gamma_Formulation
  REALTYPE, dimension(:,:)               :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE, dimension(NV)        , OPTIONAL :: Lambda_Max
  integer ,                        OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(Dimensions), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(Dimensions)           :: S

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Flux

  REALTYPE, dimension(K0:K5,NV)    :: F, F_Half
  REALTYPE, dimension(K0:K5,NV,NV) :: E_L, E_R

  REALTYPE, dimension(K0:K5) :: U, V, W, P, C, r, h
  REALTYPE, dimension(K0:K5) :: f1, f2
  REALTYPE, dimension(K0:K5) :: U_dot_N

  REALTYPE, dimension(NV)    :: Lambda
#if defined (PARALLEL_MPI)
  REALTYPE, dimension(i_0:i_1,j_0:j_1,NV) :: Lambda_i
#endif

  integer  :: i, j, k
  REALTYPE :: s0, s1, b1, b2, b3, b4, c1, u1, v1, w1
  REALTYPE :: Um, Vm, Wm, hm, f1m, f2m, gm, gm1, gmp, qm, cm, em, Xm, Ym

  Lambda_Factor   = PS_WENO_Option%Lambda_Factor
  Sound_Speed_Min = PS_WENO_Option%Sound_Speed_Min

  N_x = ZERO ; N_y = ZERO ; N_z = ONE

  S = ZERO ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

#if defined (PARALLEL_MPI)
  if (.NOT. PRESENT(Lambda_Max)) then
    if (Distributed) call Eigenvalues (Q)
  endif
#endif

!--------------------- begin of outer loop in y-dir and z-dir ------------

!$OMP PARALLEL DO PRIVATE (U, V, W, P, C, r, h, F, F_Half, Lambda, MR_Flag_z)  &
!$OMP        FIRSTPRIVATE (E_L, E_R)
  do j = j_0, j_1
!$OMP PARALLEL DO 
  do i = i_0, i_1

    U  = Q(i,j,:,2)/Q(i,j,:,1)
    V  = Q(i,j,:,3)/Q(i,j,:,1)
    W  = Q(i,j,:,4)/Q(i,j,:,1)
    f1 = Q(i,j,:,6)/Q(i,j,:,1)

    f1 = MIN(MAX(ZERO,f1),ONE)
    f2 = ONE - f1

    Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
    Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

    Gamma(i,j,:) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

    P = (Gamma(i,j,:)-ONE)*(Q(i,j,:,5) - HALF*Q(i,j,:,1)*(U**2+V**2+W**2))

#if defined (PARALLEL_MPI)
    if (.NOT. Distributed) C = SQRT(ABS(Gamma(i,j,:)*P/Q(i,j,:,1)))
#else
                           C = SQRT(ABS(Gamma(i,j,:)*P/Q(i,j,:,1)))
#endif

    r = SQRT(ABS(Q(:,j,k,1)))
    r = SQRT(ABS(Q(i,j,:,1)))
    h = (Q(i,j,:,5)+P)/Q(i,j,:,1)

    F(:,1) =   (Q(i,j,:,1)      * W    )
    F(:,2) =   (Q(i,j,:,2)      * W    )
    F(:,3) =   (Q(i,j,:,3)      * W    )
    F(:,4) =   (Q(i,j,:,4)      * W + P)
    F(:,5) =   (Q(i,j,:,5) + P) * W
    F(:,6) =   (Q(i,j,:,6)      * W    )

    if (S(3) /= ZERO) then
      F(:,1) = F(:,1) - S(3)*Q(i,j,:,1)
      F(:,2) = F(:,2) - S(3)*Q(i,j,:,2)
      F(:,3) = F(:,3) - S(3)*Q(i,j,:,3)
      F(:,4) = F(:,4) - S(3)*Q(i,j,:,4)
      F(:,5) = F(:,5) - S(3)*Q(i,j,:,5)
      F(:,6) = F(:,6) - S(3)*Q(i,j,:,6)
    endif

    U_dot_N = U*N_x + V*N_y + W*N_z

  if (.NOT. PRESENT(Lambda_Max)) then
#if defined (PARALLEL_MPI)
    if (Distributed) then
      Lambda = Lambda_i(i,j,:)
    else
#endif
      Lambda(1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(3)+C)))
      Lambda(2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(3)-C)))
      Lambda(3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(3)  )))
      Lambda(4) = Lambda(3)
      Lambda(5) = Lambda(3)
      Lambda(6) = Lambda(3)
#if defined (PARALLEL_MPI)
    endif
#endif
  else
     Lambda = Lambda_Max
  endif
      Lambda = Lambda*Lambda_Factor

! compute left and right eigenvectors of Roe's mean matrix

!$OMP PARALLEL DO PRIVATE (s0, s1, Um, Vm, Wm, hm)           &
!$OMP             PRIVATE (f1m, f2m, gm, gm1, gmp)           &
!$OMP             PRIVATE (qm, cm, em, Xm, Ym, UdotN)        &
!$OMP             PRIVATE (b1, b2, b3, b4, c1, u1, v1, w1)
    do k = k_0-1, k_1

      s0 = r(k)/(r(k) + r(k+1)) ; s1 = ONE-s0

       Um = s0* U(k) + s1* U(k+1)
       Vm = s0* V(k) + s1* V(k+1)
       Wm = s0* W(k) + s1* W(k+1)
       hm = s0* h(k) + s1* h(k+1)
      f1m = s0*f1(k) + s1*f1(k+1)

      f1m = MIN(MAX(ZERO,f1m),ONE)
      f2m = ONE - f1m

!     gm  = f1m*Gamma_Gas(1,1) + f2m*Gamma_Gas(2,1)
!     gmp = Gamma_Gas(1,1) - Gamma_Gas(2,1)

      gm  = (f1m*Cp_1 + f2m*Cp_2)/(f1m*Cv_1 + f2m*Cv_2)
      gmp = (Cp_1*Cv_2-Cp_2*Cv_1)/(f1m*Cv_1 + f2m*Cv_2)**2

      gm1 = gm-ONE

      qm = HALF*(Um*Um + Vm*Vm + Wm*Wm)
      cm = SQRT(ABS(gm1*(hm - qm)))
      em = cm**2/(gm*gm1)

      Xm  = gmp*em ; Ym = Xm/gm1

      UdotN = Um*N_x + Vm*N_y + Wm*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN

      E_R(k,1,1) = ONE
      E_R(k,1,2) = ONE
      E_R(k,1,3) = M_x
      E_R(k,1,4) = M_y
      E_R(k,1,5) = M_z
      E_R(k,1,6) = ZERO

      E_R(k,2,1) =  Um     + cm*N_x
      E_R(k,2,2) =  Um     - cm*N_x
      E_R(k,2,3) = ONE
      E_R(k,2,4) = ZERO
      E_R(k,2,5) = ZERO
      E_R(k,2,6) = ZERO

      E_R(k,3,1) =  Vm     + cm*N_y
      E_R(k,3,2) =  Vm     - cm*N_y
      E_R(k,3,3) = ZERO
      E_R(k,3,4) = ONE
      E_R(k,3,5) = ZERO
      E_R(k,3,6) = ZERO

      E_R(k,4,1) =  Wm     + cm*N_z
      E_R(k,4,2) =  Wm     - cm*N_z
      E_R(k,4,3) = ZERO
      E_R(k,4,4) = ZERO
      E_R(k,4,5) = ONE
      E_R(k,4,6) = ZERO

      E_R(k,5,1) =  hm     + cm*UdotN
      E_R(k,5,2) =  hm     - cm*UdotN
      E_R(k,5,3) = -qm*M_x + Um
      E_R(k,5,4) = -qm*M_y + Vm
      E_R(k,5,5) = -qm*M_z + Wm
      E_R(k,5,6) =         - Ym

      E_R(k,6,1) = f1m
      E_R(k,6,2) = f1m
      E_R(k,6,3) = f1m*M_x
      E_R(k,6,4) = f1m*M_y
      E_R(k,6,5) = f1m*M_z
      E_R(k,6,6) = ONE

      b1 = HALF*gm1/cm**2 ; b2 = qm*b1 ; b3 = Ym*b1 ; b4 = b2 - f1m*b3
      c1 = HALF/cm        ; u1 = Um*b1 ; v1 = Vm*b1 ; w1 = Wm*b1

      E_L(k,1,1) =     - c1*UdotN  + b4
      E_L(k,1,2) =       c1*N_x    - u1
      E_L(k,1,3) =       c1*N_y    - v1
      E_L(k,1,4) =       c1*N_z    - w1
      E_L(k,1,5) =                   b1
      E_L(k,1,6) =                   b3

      E_L(k,2,1) =       c1*UdotN  + b4
      E_L(k,2,2) =     - c1*N_x    - u1
      E_L(k,2,3) =     - c1*N_y    - v1
      E_L(k,2,4) =     - c1*N_z    - w1
      E_L(k,2,5) =                   b1
      E_L(k,2,6) =                   b3

      E_L(k,3,1) =  UdotN*N_x-TWO*Um*b4
      E_L(k,3,2) =            TWO*Um*u1 - N_x*N_x + ONE
      E_L(k,3,3) =            TWO*Um*v1 - N_x*N_y
      E_L(k,3,4) =            TWO*Um*w1 - N_x*N_z
      E_L(k,3,5) =           -TWO*Um*b1
      E_L(k,3,6) =           -TWO*Um*b3

      E_L(k,4,1) =  UdotN*N_y-TWO*Vm*b4
      E_L(k,4,2) =            TWO*Vm*u1 - N_y*N_x
      E_L(k,4,3) =            TWO*Vm*v1 - N_y*N_y + ONE
      E_L(k,4,4) =            TWO*Vm*w1 - N_y*N_z
      E_L(k,4,5) =           -TWO*Vm*b1
      E_L(k,4,6) =           -TWO*Vm*b3

      E_L(k,5,1) =  UdotN*N_z-TWO*Wm*b4
      E_L(k,5,2) =            TWO*Wm*u1 - N_z*N_x
      E_L(k,5,3) =            TWO*Wm*v1 - N_z*N_y 
      E_L(k,5,4) =            TWO*Wm*w1 - N_z*N_z + ONE
      E_L(k,5,5) =           -TWO*Wm*b1
      E_L(k,5,6) =           -TWO*Wm*b3

      E_L(k,6,1) = -f1m
      E_L(k,6,2) =  ZERO
      E_L(k,6,3) =  ZERO
      E_L(k,6,4) =  ZERO
      E_L(k,6,5) =  ZERO
      E_L(k,6,6) =  ONE

    enddo
!$OMP END PARALLEL DO

    call PS_WENO_Euler   (Order, K0, K5, k_0-1, k_1, NV,            &
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

    Subroutine Eigenvalues (Q)

    REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
    REALTYPE, dimension(i_0:i_1,j_0:j_1,NV)   :: Lambda_Local

    integer                             :: MPI_Communicator
    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid


!$OMP PARALLEL DO PRIVATE (U, V, W, P, C, f1, f2, Cp_1, Cp_2, Cv_1, Cv_2)
    do j = j_0, j_1
!$OMP PARALLEL DO
    do i = i_0, i_1

       U = Q(i,j,:,2)/Q(i,j,:,1)
       V = Q(i,j,:,3)/Q(i,j,:,1)
       W = Q(i,j,:,4)/Q(i,j,:,1)
      f1 = Q(i,j,:,6)/Q(i,j,:,1)

      f1 = MIN(MAX(ZERO,f1),ONE)
      f2 = ONE - f1

      Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
      Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

      Gamma(i,j,:) = (f1*Cp_1 + f2*Cp_2)/(f1*Cv_1 + f2*Cv_2)

      P = (Gamma(i,j,:)-ONE)*(Q(i,j,:,5) - HALF*Q(i,j,:,1)*(U**2+V**2+W**2))

      C = SQRT(ABS(Gamma(i,j,:)*P/Q(i,j,:,1)))

      U_dot_N = U*N_x + V*N_y + W*N_z

      Lambda_Local(i,j,1) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(3)+C)))
      Lambda_Local(i,j,2) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(3)-C)))
      Lambda_Local(i,j,3) = MAX(Sound_Speed_Min, MAXVAL(ABS(U_dot_N-S(3)  )))
      Lambda_Local(i,j,4) = Lambda_Local(i,j,3)
      Lambda_Local(i,j,5) = Lambda_Local(i,j,3)
      Lambda_Local(i,j,6) = Lambda_Local(i,j,3)

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

  END Subroutine PS_WENO_Euler_Flux_3D2S_3

END MODULE WENO_Euler_Flux_TWO_Species
