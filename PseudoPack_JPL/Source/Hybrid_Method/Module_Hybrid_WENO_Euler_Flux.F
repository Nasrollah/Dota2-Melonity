#include "pseudopack.h"


MODULE Hybrid_WENO_Euler_Flux

  USE         Euler_Function
  USE    WENO_Euler_Flux_General
  USE Central_Euler_Flux

IMPLICIT NONE

integer , parameter :: Idx = 0

INTERFACE PS_Hybrid_WENO_Euler_Flux
  MODULE PROCEDURE Hybrid_1D_a
  MODULE PROCEDURE Hybrid_2D_a
  MODULE PROCEDURE Hybrid_3D_a

  MODULE PROCEDURE Hybrid_1D_b
  MODULE PROCEDURE Hybrid_2D_b
  MODULE PROCEDURE Hybrid_3D_b
END INTERFACE

PRIVATE
PUBLIC :: PS_Hybrid_WENO_Euler_Flux

CONTAINS
!
!==========================================================================
!
  Subroutine Hybrid_1D_a    (Index, Distributed, WENO_Order,            &
                                    Central_Order, MR_Order, dx,        &
                                    N0, N5, N2, N3, NV,                 &
                             Gamma_Gas, Gamma,                          &
                             Q, U, P, H, C, F, D_F,                     &
                             Lambda_Max,                                &
                             MR_N0, MR_N5,                              &
                             Flag, Flag_WENO,                           &
                             D_F_Epsilon, V_Moving_Frame,               &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5

  REALTYPE, dimension(:,:)                :: Gamma_Gas
  REALTYPE, dimension(:)                  :: Gamma                   

  REALTYPE, dimension(N0:N5,NV)          :: Q, F, D_F
  REALTYPE, dimension(N0:N5)             :: U, P, H, C
  REALTYPE, dimension(NV)                :: Lambda_Max

  REALTYPE, dimension(1)      , OPTIONAL :: V_Moving_Frame

  integer , dimension(MR_N0:MR_N5,1)           :: Flag
  integer , dimension(   N0:   N5,1), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

      call Hybrid_1D_x_a    (Index, Distributed, WENO_Order,         &
                                    Central_Order, MR_Order, dx,     &
                                    N0, N5, N2, N3, NV,              &
                             Gamma_Gas, Gamma,                       &
                             Q, U, P, H, C, F, D_F,                  &
                             Lambda_Max,                             &
                             MR_N0, MR_N5,                           &
                             Flag, Flag_WENO,                        &
                             D_F_Epsilon, V_Moving_Frame,            &
                             MPI_Comm_Type)

  END Subroutine Hybrid_1D_a 
!
!==========================================================================
!
  Subroutine Hybrid_1D_x_a  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, P, H, C, F, D_F,                      &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5,                               &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5

  REALTYPE, dimension(:,:)              :: Gamma_Gas
  REALTYPE, dimension(N0:N5)   , TARGET :: Gamma                   

  REALTYPE, dimension(N0:N5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5)   , TARGET :: U, P, H, C
  REALTYPE, dimension(NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(1)      , OPTIONAL :: V_Moving_Frame

  integer , dimension(   N0:   N5  )           :: Flag_1D
  integer , dimension(MR_N0:MR_N5,1)           :: Flag
  integer , dimension(   N0:   N5,1), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:)  , POINTER :: U_p, P_p, H_p, C_p
  REALTYPE, dimension(:)  , POINTER :: Lambda_p, Gamma_p

  integer  :: i, n
  integer  :: i_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,Index) = 0

  if (MAXVAL(Flag(N2:N3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

    Flag_1D = Flag(N0:N5,Index) ; if (MAXVAL(Flag_1D) == 0) RETURN

    i = N2-1 

    Do WHILE (i <= N3) 
      i = i+1

      if (Flag_1D(i) == 0) CYCLE
    
      do i_End = i+1,N3 ; if (Flag_1D(i_End) == 0) EXIT ; enddo

      call Hybrid_Range (N2, N3, i, i_End, M_s, M_Buffer, &
                              N0_p, N5_p, N2_p, N3_p)

      D_F_0(-M_s:-1 ,:) = D_F(N0_p  :N2_p-1,:)
      D_F_N(   1:M_s,:) = D_F(N3_p+1:N5_p  ,:)

           Q_p =>     Q(N0_p:N5_p,:)
           F_p =>     F(N0_p:N5_p,:)
         D_F_p =>   D_F(N0_p:N5_p,:)

           U_p =>     U(N0_p:N5_p)
           P_p =>     P(N0_p:N5_p)
           H_p =>     H(N0_p:N5_p)
           C_p =>     C(N0_p:N5_p)
       Gamma_p => Gamma(N0_p:N5_p)

      Lambda_p => Lambda_Max(:)

      call PS_WENO_Euler_Flux        (Index, Distributed, WENO_Order, dx,   &
                                             N0_p, N5_p, N2_p, N3_p, NV,    &
                                      Gamma_Gas, Gamma_p,                   &
                                      Q_p, U_p, P_p, H_p, C_p,              &
                                      D_F_p,                                &
                                              Flux_x=F_p,                   &
                                            Lambda_x=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(N0_p  :N2_p-1,:) = D_F_0(-M_s:-1 ,:)
      D_F(N3_p+1:N5_p  ,:) = D_F_N(   1:M_s,:)

      Flag_1D(N2_p:N3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(N2_p:N3_p,Index) = 1

      i = i_End-1
    enddo

  Nullify (Q_p, F_p, D_F_p, U_p, P_p, H_p, C_p)

  END Subroutine Hybrid_1D_x_a
!
!==========================================================================
!
  Subroutine Hybrid_2D_a  (Index, Distributed, WENO_Order,             &
                                  Central_Order, MR_Order, dx,         &
                                  N0, N5, N2, N3,                      &
                                  M0, M5, M2, M3, NV,                  &
                           Gamma_Gas, Gamma,                           &
                           Q, U, V, P, H, C, F, D_F,                   &
                           Lambda_Max,                                 &
                           MR_N0, MR_N5, MR_M0, MR_M5,                 &
                           Flag, Flag_WENO,                            &
                           D_F_Epsilon, V_Moving_Frame,                &
                           MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dx

  integer  :: NV
 
  integer  :: MR_N0, MR_N5, MR_M0, MR_M5

  REALTYPE, dimension(:,:)                     :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5)             :: Gamma                   

  REALTYPE, dimension(N0:N5,M0:M5,NV)         :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5)            :: U, V, P, H, C
  REALTYPE, dimension(:,:)                    :: Lambda_Max

  REALTYPE, dimension(2)            , OPTIONAL :: V_Moving_Frame

  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,2)           :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,2), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  SELECT CASE (Index)
    CASE (1)
      call Hybrid_2D_x_a      (Index, Distributed, WENO_Order,         &
                                      Central_Order, MR_Order, dx,     &
                                      N0, N5, N2, N3,                  &
                                      M0, M5, M2, M3, NV,              &
                               Gamma_Gas, Gamma,                       &
                               Q, U, V, P, H, C, F, D_F,               &
                               Lambda_Max,                             &
                               MR_N0, MR_N5, MR_M0, MR_M5,             &
                               Flag, Flag_WENO,                        &
                               D_F_Epsilon, V_Moving_Frame,            &
                               MPI_Comm_Type)

    CASE (2)
      call Hybrid_2D_y_a      (Index, Distributed, WENO_Order,         &
                                      Central_Order, MR_Order, dx,     &
                                      N0, N5, N2, N3,                  &
                                      M0, M5, M2, M3, NV,              &
                               Gamma_Gas, Gamma,                       &
                               Q, U, V, P, H, C, F, D_F,               &
                               Lambda_Max,                             &
                               MR_N0, MR_N5, MR_M0, MR_M5,             &
                               Flag, Flag_WENO,                        &
                               D_F_Epsilon, V_Moving_Frame,            &
                               MPI_Comm_Type)

  END SELECT

  END Subroutine Hybrid_2D_a 
!
!==========================================================================
!
  Subroutine Hybrid_2D_x_a  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, V, P, H, C, F, D_F,                   &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5,                 &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5

  REALTYPE, dimension(:,:)                    :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5)   , TARGET :: Gamma                   

  REALTYPE, dimension(N0:N5,M0:M5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5)   , TARGET :: U, V, P, H, C
  REALTYPE, dimension(M0:M5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(2)            , OPTIONAL :: V_Moving_Frame

  integer , dimension(   N0:   N5              )           :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,2)           :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,2), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:)  , POINTER :: U_p, V_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:)  , POINTER :: Lambda_p, Gamma_p

  integer  :: i, j, n
  integer  :: i_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do j = M2,M3
    M0_p = j ; M2_p = j ; M3_p = j ; M5_p = j 

    Flag_1D = Flag(N0:N5,j,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    i = N2-1 

    Do WHILE (i <= N3) 
      i = i+1

      if (Flag_1D(i) == 0) CYCLE

      do i_End = i+1,N3 ; if (Flag_1D(i_End) == 0) EXIT ; enddo

      call Hybrid_Range (N2, N3, i, i_End, M_s, M_Buffer, &
                              N0_p, N5_p, N2_p, N3_p)


      D_F_0(-M_s:-1 ,:) = D_F(N0_p  :N2_p-1,j,:)
      D_F_N(   1:M_s,:) = D_F(N3_p+1:N5_p  ,j,:)

           Q_p =>     Q(N0_p:N5_p,j:j,:)
           F_p =>     F(N0_p:N5_p,j:j,:)
         D_F_p =>   D_F(N0_p:N5_p,j:j,:)
  
           U_p =>     U(N0_p:N5_p,j:j)
           V_p =>     V(N0_p:N5_p,j:j)
           P_p =>     P(N0_p:N5_p,j:j)
           H_p =>     H(N0_p:N5_p,j:j)
           C_p =>     C(N0_p:N5_p,j:j)
       Gamma_p => Gamma(N0_p:N5_p,j:j)

      Lambda_p => Lambda_Max(j:j,:)
      
      call PS_WENO_Euler_Flux        (Index, Distributed, WENO_Order, dx,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Idx  , Distributed, WENO_Order, dx,   &
                                             M0_p, M5_p, M2_p, M3_p, NV,    &
                                      Gamma_Gas, Gamma_p,                   &
                                      Q_p, U_p, V_p, P_p, H_p, C_p,         &
                                      D_F_p,                                &
                                              Flux_x=F_p,                   &
                                            Lambda_x=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(N0_p  :N2_p-1,j,:) = D_F_0(-M_s:-1 ,:)
      D_F(N3_p+1:N5_p  ,j,:) = D_F_N(   1:M_s,:)

      Flag_1D(N2_p:N3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(N2_p:N3_p,j,Index) = 1

      i = i_End-1
    enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, P_p, H_p, C_p)

  END Subroutine Hybrid_2D_x_a
!
!==========================================================================
!
  Subroutine Hybrid_2D_y_a  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dy,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, V, P, H, C, F, D_F,                   &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5,                 &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5

  REALTYPE, dimension(:,:)                    :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5)   , TARGET :: Gamma                   

  REALTYPE, dimension(N0:N5,M0:M5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5)   , TARGET :: U, V, P, H, C
  REALTYPE, dimension(N0:N5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(2)            , OPTIONAL :: V_Moving_Frame

  integer , dimension(               M0:   M5  )           :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,2)           :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,2), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:)  , POINTER :: U_p, V_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:)  , POINTER :: Lambda_p, Gamma_p

  integer  :: i, j, n
  integer  :: j_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do i = N2,N3
    N0_p = i ; N2_p = i ; N3_p = i ; N5_p = i 

    Flag_1D = Flag(i,M0:M5,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    j = M2-1  

    Do WHILE (j <= M3) 
      j = j+1

      if (Flag_1D(j) == 0) CYCLE
    
      do j_End = j+1,M3 ; if (Flag_1D(j_End) == 0) EXIT ; enddo

      call Hybrid_Range (M2, M3, j, j_End, M_s, M_Buffer, &
                              M0_p, M5_p, M2_p, M3_p)

      D_F_0(-M_s:-1 ,:) = D_F(i,M0_p  :M2_p-1,:)
      D_F_N(   1:M_s,:) = D_F(i,M3_p+1:M5_p  ,:)

           Q_p =>     Q(i:i,M0_p:M5_p,:)
           F_p =>     F(i:i,M0_p:M5_p,:)
         D_F_p =>   D_F(i:i,M0_p:M5_p,:)
  
           U_p =>     U(i:i,M0_p:M5_p)
           V_p =>     V(i:i,M0_p:M5_p)
           P_p =>     P(i:i,M0_p:M5_p)
           H_p =>     H(i:i,M0_p:M5_p)
           C_p =>     C(i:i,M0_p:M5_p)
       Gamma_p => Gamma(i:i,M0_p:M5_p)

      Lambda_p => Lambda_Max(i:i,:)

      call PS_WENO_Euler_Flux        (Idx  , Distributed, WENO_Order, dy,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Index, Distributed, WENO_Order, dy,   &
                                             M0_p, M5_p, M2_p, M3_p, NV,    &
                                      Gamma_Gas, Gamma_p,                   &
                                      Q_p, U_p, V_p, P_p, H_p, C_p, D_F_p,  &
                                              Flux_y=F_p,                   &
                                            Lambda_y=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(i,M0_p  :M2_p-1,:) = D_F_0(-M_s:-1 ,:)
      D_F(i,M3_p+1:M5_p  ,:) = D_F_N(   1:M_s,:)

      Flag_1D(M2_p:M3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(i,M2_p:M3_p,Index) = 1

      j = j_End - 1
    enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, P_p, H_p, C_p)

  END Subroutine Hybrid_2D_y_a 
!
!==========================================================================
!
  Subroutine Hybrid_3D_a    (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, V, W, P, H, C, F, D_F,                &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,   &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5

  REALTYPE, dimension(:,:)                          :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)            :: Gamma                   

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)         :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)            :: U, V, W, P, H, C
  REALTYPE, dimension(:,:,:)                        :: Lambda_Max

  REALTYPE, dimension(3)                  , OPTIONAL :: V_Moving_Frame

  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,   K0:   K5,3),            &
                                                    OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  SELECT CASE (Index)
    CASE (1)
      call Hybrid_3D_x_a    (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, V, W, P, H, C, F, D_F,                &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,   &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

    CASE (2)
      call Hybrid_3D_y_a    (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, V, W, P, H, C, F, D_F,                &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,   &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

    CASE (3)
      call Hybrid_3D_z_a    (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, V, W, P, H, C, F, D_F,                &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,   &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)
  END SELECT

  END Subroutine Hybrid_3D_a
!
!==========================================================================
!
  Subroutine Hybrid_3D_x_a  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, V, W, P, H, C, F, D_F,                &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,   &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5

  REALTYPE, dimension(:,:)                          :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , TARGET :: Gamma                   

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , TARGET :: U, V, W, P, H, C
  REALTYPE, dimension(M0:M5,K0:K5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(3)                  , OPTIONAL :: V_Moving_Frame

  integer , dimension(   N0:   N5                          ) :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,   K0:   K5,3),            &
                                                    OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:,:)  , POINTER :: U_p, V_p, W_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:,:)  , POINTER :: Lambda_p, Gamma_p

  integer  :: i, j, k, n
  integer  :: i_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p
  integer  :: K0_p, K2_p, K3_p, K5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,K2:K3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do k = K2,K3
    K0_p = k ; K2_p = k ; K3_p = k ; K5_p = k 
  
  do j = M2,M3
    M0_p = j ; M2_p = j ; M3_p = j ; M5_p = j 

    Flag_1D = Flag(N0:N5,j,k,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    i = N2-1 

    Do WHILE (i <= N3) 
      i = i+1

      if (Flag_1D(i) == 0) CYCLE
    
      do i_End = i+1,N3 ; if (Flag_1D(i_End) == 0) EXIT ; enddo

      call Hybrid_Range (N2, N3, i, i_End, M_s, M_Buffer, &
                              N0_p, N5_p, N2_p, N3_p)

      D_F_0(-M_s:-1 ,:) = D_F(N0_p  :N2_p-1,j,k,:)
      D_F_N(   1:M_s,:) = D_F(N3_p+1:N5_p  ,j,k,:)

           Q_p =>     Q(N0_p:N5_p,j:j,k:k,:)
           F_p =>     F(N0_p:N5_p,j:j,k:k,:)
         D_F_p =>   D_F(N0_p:N5_p,j:j,k:k,:)

           U_p =>     U(N0_p:N5_p,j:j,k:k)
           V_p =>     V(N0_p:N5_p,j:j,k:k)
           W_p =>     W(N0_p:N5_p,j:j,k:k)
           P_p =>     P(N0_p:N5_p,j:j,k:k)
           H_p =>     H(N0_p:N5_p,j:j,k:k)
           C_p =>     C(N0_p:N5_p,j:j,k:k)
       Gamma_p => Gamma(N0_p:N5_p,j:j,k:k)

      Lambda_p => Lambda_Max(j:j,k:k,:)

      call PS_WENO_Euler_Flux        (Index, Distributed, WENO_Order, dx,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Idx  , Distributed, WENO_Order, dx,   &
                                             M0_p, M5_p, M2_p, M3_p,        &
                                      Idx  , Distributed, WENO_Order, dx,   &
                                             K0_p, K5_p, K2_p, K3_p, NV,    &
                                      Gamma_Gas, Gamma_p,                   &
                                      Q_p, U_p, V_p, W_p, P_p, H_p, C_p,    &
                                      D_F_p,                                &
                                              Flux_x=F_p,                   &
                                            Lambda_x=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(N0_p  :N2_p-1,j,k,:) = D_F_0(-M_s:-1 ,:)
      D_F(N3_p+1:N5_p  ,j,k,:) = D_F_N(   1:M_s,:)

      Flag_1D(N2_p:N3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(N2_p:N3_p,j,k,Index) = 1

      i = i_End-1
    enddo
  enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  END Subroutine Hybrid_3D_x_a
!
!==========================================================================
!
  Subroutine Hybrid_3D_y_a  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dy,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, V, W, P, H, C, F, D_F,                &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,   &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dy

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5

  REALTYPE, dimension(:,:)                          :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , TARGET :: Gamma                   

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , TARGET :: U, V, W, P, H, C
  REALTYPE, dimension(N0:N5,K0:K5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(3)                  , OPTIONAL :: V_Moving_Frame

  integer , dimension(               M0:   M5              ) :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,   K0:   K5,3),            &
                                                    OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:,:)  , POINTER :: U_p, V_p, W_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:,:)  , POINTER :: Lambda_p, Gamma_p

  integer  :: i, j, k, n
  integer  :: j_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p
  integer  :: K0_p, K2_p, K3_p, K5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,K2:K3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do k = K2,K3
    K0_p = k ; K2_p = k ; K3_p = k ; K5_p = k 
  
  do i = N2,N3
    N0_p = i ; N2_p = i ; N3_p = i ; N5_p = i 

    Flag_1D = Flag(i,M0:M5,k,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    j = M2-1 

    Do WHILE (j <= M3) 
      j = j+1

      if (Flag_1D(j) == 0) CYCLE
    
      do j_End = j+1,M3 ; if (Flag_1D(j_End) == 0) EXIT ; enddo

      call Hybrid_Range (M2, M3, j, j_End, M_s, M_Buffer, &
                              M0_p, M5_p, M2_p, M3_p)

      D_F_0(-M_s:-1 ,:) = D_F(i,M0_p  :M2_p-1,k,:)
      D_F_N(   1:M_s,:) = D_F(i,M3_p+1:M5_p  ,k,:)

           Q_p =>     Q(i:i,M0_p:M5_p,k:k,:)
           F_p =>     F(i:i,M0_p:M5_p,k:k,:)
         D_F_p =>   D_F(i:i,M0_p:M5_p,k:k,:)
  
           U_p =>     U(i:i,M0_p:M5_p,k:k)
           V_p =>     V(i:i,M0_p:M5_p,k:k)
           W_p =>     W(i:i,M0_p:M5_p,k:k)
           P_p =>     P(i:i,M0_p:M5_p,k:k)
           H_p =>     H(i:i,M0_p:M5_p,k:k)
           C_p =>     C(i:i,M0_p:M5_p,k:k)
       Gamma_p => Gamma(i:i,M0_p:M5_p,k:k)

      Lambda_p => Lambda_Max(i:i,k:k,:)

      call PS_WENO_Euler_Flux        (Idx  , Distributed, WENO_Order, dy,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Index, Distributed, WENO_Order, dy,   &
                                             M0_p, M5_p, M2_p, M3_p,        &
                                      Idx  , Distributed, WENO_Order, dy,   &
                                             K0_p, K5_p, K2_p, K3_p, NV,    &
                                      Gamma_Gas, Gamma_p,                   &
                                      Q_p, U_p, V_p, W_p, P_p, H_p, C_p,    &
                                      D_F_p,                                &
                                              Flux_y=F_p,                   &
                                            Lambda_y=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(i,M0_p  :M2_p-1,k,:) = D_F_0(-M_s:-1 ,:)
      D_F(i,M3_p+1:M5_p  ,k,:) = D_F_N(   1:M_s,:)

      Flag_1D(M2_p:M3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(i,M2_p:M3_p,k,Index) = 1

      j = j_End - 1
    enddo
  enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  END Subroutine Hybrid_3D_y_a
!
!==========================================================================
!
  Subroutine Hybrid_3D_z_a  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dz,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma_Gas, Gamma,                           &
                             Q, U, V, W, P, H, C, F, D_F,                &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,   &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dz

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5

  REALTYPE, dimension(:,:)                          :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , TARGET :: Gamma                   

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , TARGET :: U, V, W, P, H, C
  REALTYPE, dimension(N0:N5,M0:M5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(3)                  , OPTIONAL :: V_Moving_Frame

  integer , dimension(                           K0:   K5  ) :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,   K0:   K5,3),            &
                                                    OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:,:)  , POINTER :: U_p, V_p, W_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:,:)  , POINTER :: Lambda_p, Gamma_p

  integer  :: i, j, k, n
  integer  :: k_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p
  integer  :: K0_p, K2_p, K3_p, K5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,K2:K3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do j = M2,M3
    M0_p = j ; M2_p = j ; M3_p = j ; M5_p = j 
  
  do i = N2,N3
    N0_p = i ; N2_p = i ; N3_p = i ; N5_p = i 

    Flag_1D = Flag(i,j,K0:K5,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    k = K2-1  

    Do WHILE (k <= K3) 
      k = k+1

      if (Flag_1D(k) == 0) CYCLE
    
      do k_End = k+1,K3 ; if (Flag_1D(k_End) == 0) EXIT ; enddo

      call Hybrid_Range (K2, K3, k, k_End, M_s, M_Buffer, &
                              K0_p, K5_p, K2_p, K3_p)

      D_F_0(-M_s:-1 ,:) = D_F(i,j,K0_p  :K2_p-1,:)
      D_F_N(   1:M_s,:) = D_F(i,j,K3_p+1:K5_p  ,:)

           Q_p =>     Q(i:i,j:j,K0_p:K5_p,:)
           F_p =>     F(i:i,j:j,K0_p:K5_p,:)
         D_F_p =>   D_F(i:i,j:j,K0_p:K5_p,:)

           U_p =>     U(i:i,j:j,K0_p:K5_p)
           V_p =>     V(i:i,j:j,K0_p:K5_p)
           W_p =>     W(i:i,j:j,K0_p:K5_p)
           P_p =>     P(i:i,j:j,K0_p:K5_p)
           H_p =>     H(i:i,j:j,K0_p:K5_p)
           C_p =>     C(i:i,j:j,K0_p:K5_p)
       Gamma_p => Gamma(i:i,j:j,K0_p:K5_p)

      Lambda_p => Lambda_Max(i:i,j:j,:)

      call PS_WENO_Euler_Flux        (Idx  , Distributed, WENO_Order, dz,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Idx  , Distributed, WENO_Order, dz,   &
                                             M0_p, M5_p, M2_p, M3_p,        &
                                      Index, Distributed, WENO_Order, dz,   &
                                             K0_p, K5_p, K2_p, K3_p, NV,    &
                                      Gamma_Gas, Gamma_p,                   &
                                      Q_p, U_p, V_p, W_p, P_p, H_p, C_p,    &
                                      D_F_p,                                &
                                              Flux_z=F_p,                   &
                                            Lambda_z=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(i,j,K0_p  :K2_p-1,:) = D_F_0(-M_s:-1 ,:)
      D_F(i,j,K3_p+1:K5_p  ,:) = D_F_N(   1:M_s,:)

      Flag_1D(K2_p:K3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(i,j,K2_p:K3_p,Index) = 1

      k = k_End-1
    enddo
  enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  END Subroutine Hybrid_3D_z_a
!
!==========================================================================
!
  Subroutine Hybrid_1D_b  (Index, Distributed, WENO_Order,             &
                                  Central_Order, MR_Order, dx,         &
                                  N0, N5, N2, N3, NV,                  &
                           Gamma, Q, U, P, H, C, F, D_F,               &
                           Lambda_Max,                                 &
                           MR_N0, MR_N5,                               &
                           Flag, Flag_WENO,                            &
                           D_F_Epsilon, V_Moving_Frame,                &
                           MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5

  REALTYPE                               :: Gamma

  REALTYPE, dimension(N0:N5,NV)          :: Q, F, D_F
  REALTYPE, dimension(N0:N5)             :: U, P, H, C
  REALTYPE, dimension(:)                 :: Lambda_Max

  REALTYPE, dimension(1)      , OPTIONAL :: V_Moving_Frame

  integer , dimension(MR_N0:MR_N5,1)           :: Flag
  integer , dimension(   N0:   N5,1), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

      call Hybrid_1D_x_b    (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3, NV,                  &
                             Gamma, Q, U, P, H, C, F, D_F,               &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5,                               &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  END Subroutine Hybrid_1D_b 
!
!==========================================================================
!
  Subroutine Hybrid_1D_x_b  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3, NV,                  &
                             Gamma, Q, U, P, H, C, F, D_F,               &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5,                               &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5

  REALTYPE                              :: Gamma

  REALTYPE, dimension(N0:N5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5)   , TARGET :: U, P, H, C
  REALTYPE, dimension(NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(1)      , OPTIONAL :: V_Moving_Frame

  integer , dimension(   N0:   N5  )           :: Flag_1D
  integer , dimension(MR_N0:MR_N5,1)           :: Flag
  integer , dimension(   N0:   N5,1), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:)  , POINTER :: U_p, P_p, H_p, C_p
  REALTYPE, dimension(:)  , POINTER :: Lambda_p

  integer  :: i, n
  integer  :: i_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,Index) = 0

  if (MAXVAL(Flag(N2:N3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

    Flag_1D = Flag(N0:N5,Index) ; if (MAXVAL(Flag_1D) == 0) RETURN

    i = N2-1 

    Do WHILE (i <= N3) 
      i = i+1

      if (Flag_1D(i) == 0) CYCLE
    
      do i_End = i+1,N3 ; if (Flag_1D(i_End) == 0) EXIT ; enddo

      call Hybrid_Range (N2, N3, i, i_End, M_s, M_Buffer, &
                              N0_p, N5_p, N2_p, N3_p)

      D_F_0(-M_s:-1 ,:) = D_F(N0_p  :N2_p-1,:)
      D_F_N(   1:M_s,:) = D_F(N3_p+1:N5_p  ,:)

           Q_p =>     Q(N0_p:N5_p,:)
           F_p =>     F(N0_p:N5_p,:)
         D_F_p =>   D_F(N0_p:N5_p,:)

           U_p =>     U(N0_p:N5_p)
           P_p =>     P(N0_p:N5_p)
           H_p =>     H(N0_p:N5_p)
           C_p =>     C(N0_p:N5_p)

      Lambda_p => Lambda_Max(:)

      call PS_WENO_Euler_Flux        (Index, Distributed, WENO_Order, dx,   &
                                             N0_p, N5_p, N2_p, N3_p, NV,    &
                                      Gamma, Q_p, U_p,                      &
                                             P_p, H_p, C_p, D_F_p,          &
                                              Flux_x=F_p,                   &
                                            Lambda_x=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(N0_p  :N2_p-1,:) = D_F_0(-M_s:-1 ,:)
      D_F(N3_p+1:N5_p  ,:) = D_F_N(   1:M_s,:)

      Flag_1D(N2_p:N3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(N2_p:N3_p,Index) = 1

      i = i_End-1
    enddo

  Nullify (Q_p, F_p, D_F_p, U_p, P_p, H_p, C_p)

  END Subroutine Hybrid_1D_x_b
!
!==========================================================================
!
  Subroutine Hybrid_2D_b  (Index, Distributed, WENO_Order,             &
                                  Central_Order, MR_Order, dx,         &
                                  N0, N5, N2, N3,                      &
                                  M0, M5, M2, M3, NV,                  &
                           Gamma, Q, U, V, P, H, C, F, D_F,            &
                           Lambda_Max,                                 &
                           MR_N0, MR_N5, MR_M0, MR_M5,                 &
                           Flag, Flag_WENO,                            &
                           D_F_Epsilon, V_Moving_Frame,                &
                           MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5

  REALTYPE                                    :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,NV)         :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5)            :: U, V, P, H, C
  REALTYPE, dimension(:,:)                    :: Lambda_Max

  REALTYPE, dimension(2)            , OPTIONAL :: V_Moving_Frame

  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,2)           :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,2), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  SELECT CASE (Index)
    CASE (1)
      call Hybrid_2D_x_b    (Index, Distributed, WENO_Order,         &
                                    Central_Order, MR_Order, dx,     &
                                    N0, N5, N2, N3,                  &
                                    M0, M5, M2, M3, NV,              &
                             Gamma, Q, U, V, P, H, C, F, D_F,        &
                             Lambda_Max,                             &
                             MR_N0, MR_N5, MR_M0, MR_M5,             &
                             Flag, Flag_WENO,                        &
                             D_F_Epsilon, V_Moving_Frame,            &
                             MPI_Comm_Type)

    CASE (2)
      call Hybrid_2D_y_b    (Index, Distributed, WENO_Order,         &
                                    Central_Order, MR_Order, dx,     &
                                    N0, N5, N2, N3,                  &
                                    M0, M5, M2, M3, NV,              &
                             Gamma, Q, U, V, P, H, C, F, D_F,        &
                             Lambda_Max,                             &
                             MR_N0, MR_N5, MR_M0, MR_M5,             &
                             Flag, Flag_WENO,                        &
                             D_F_Epsilon, V_Moving_Frame,            &
                             MPI_Comm_Type)

  END SELECT

  END Subroutine Hybrid_2D_b 
!
!==========================================================================
!
  Subroutine Hybrid_2D_x_b  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3, NV,                  &
                             Gamma, Q, U, V, P, H, C, F, D_F,            &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5,                 &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5

  REALTYPE                                    :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5)   , TARGET :: U, V, P, H, C
  REALTYPE, dimension(M0:M5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(2)            , OPTIONAL :: V_Moving_Frame

  integer , dimension(   N0:   N5              )           :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,2)           :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,2), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:)  , POINTER :: U_p, V_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:)  , POINTER :: Lambda_p

  integer  :: i, j, n
  integer  :: i_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do j = M2,M3
    M0_p = j ; M2_p = j ; M3_p = j ; M5_p = j 

    Flag_1D = Flag(N0:N5,j,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    i = N2-1 

    Do WHILE (i <= N3) 
      i = i+1

      if (Flag_1D(i) == 0) CYCLE
    
      do i_End = i+1,N3 ; if (Flag_1D(i_End) == 0) EXIT ; enddo

      call Hybrid_Range (N2, N3, i, i_End, M_s, M_Buffer, &
                              N0_p, N5_p, N2_p, N3_p)

      D_F_0(-M_s:-1 ,:) = D_F(N0_p  :N2_p-1,j,:)
      D_F_N(   1:M_s,:) = D_F(N3_p+1:N5_p  ,j,:)

           Q_p =>     Q(N0_p:N5_p,j:j,:)
           F_p =>     F(N0_p:N5_p,j:j,:)
         D_F_p =>   D_F(N0_p:N5_p,j:j,:)

           U_p =>     U(N0_p:N5_p,j:j)
           V_p =>     V(N0_p:N5_p,j:j)
           P_p =>     P(N0_p:N5_p,j:j)
           H_p =>     H(N0_p:N5_p,j:j)
           C_p =>     C(N0_p:N5_p,j:j)

      Lambda_p => Lambda_Max(j:j,:)

      call PS_WENO_Euler_Flux        (Index, Distributed, WENO_Order, dx,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Idx  , Distributed, WENO_Order, dx,   &
                                             M0_p, M5_p, M2_p, M3_p, NV,    &
                                      Gamma, Q_p, U_p, V_p,                 &
                                             P_p, H_p, C_p, D_F_p,          &
                                              Flux_x=F_p,                   &
                                            Lambda_x=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(N0_p  :N2_p-1,j,:) = D_F_0(-M_s:-1 ,:)
      D_F(N3_p+1:N5_p  ,j,:) = D_F_N(   1:M_s,:)

      Flag_1D(N2_p:N3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(N2_p:N3_p,j,Index) = 1

      i = i_End-1
    enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, P_p, H_p, C_p)

  END Subroutine Hybrid_2D_x_b
!
!==========================================================================
!
  Subroutine Hybrid_2D_y_b  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dy,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3, NV,                  &
                             Gamma, Q, U, V, P, H, C, F, D_F,            &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5,                 &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  REALTYPE :: dy

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5

  REALTYPE                                    :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5)   , TARGET :: U, V, P, H, C
  REALTYPE, dimension(N0:N5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(2)            , OPTIONAL :: V_Moving_Frame

  integer , dimension(               M0:   M5  )           :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,2)           :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,2), OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:)  , POINTER :: U_p, V_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:)  , POINTER :: Lambda_p

  integer  :: i, j, n
  integer  :: j_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do i = N2,N3
    N0_p = i ; N2_p = i ; N3_p = i ; N5_p = i 

    Flag_1D = Flag(i,M0:M5,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    j = M2-1 

    Do WHILE (j <= M3) 
      j = j+1

      if (Flag_1D(j) == 0) CYCLE
    
      do j_End = j+1,M3 ; if (Flag_1D(j_End) == 0) EXIT ; enddo

      call Hybrid_Range (M2, M3, j, j_End, M_s, M_Buffer, &
                              M0_p, M5_p, M2_p, M3_p)

      D_F_0(-M_s:-1 ,:) = D_F(i,M0_p  :M2_p-1,:)
      D_F_N(   1:M_s,:) = D_F(i,M3_p+1:M5_p  ,:)

           Q_p =>     Q(i:i,M0_p:M5_p,:)
           F_p =>     F(i:i,M0_p:M5_p,:)
         D_F_p =>   D_F(i:i,M0_p:M5_p,:)

           U_p =>     U(i:i,M0_p:M5_p)
           V_p =>     V(i:i,M0_p:M5_p)
           P_p =>     P(i:i,M0_p:M5_p)
           H_p =>     H(i:i,M0_p:M5_p)
           C_p =>     C(i:i,M0_p:M5_p)

      Lambda_p => Lambda_Max(i:i,:)

      call PS_WENO_Euler_Flux        (Idx  , Distributed, WENO_Order, dy,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Index, Distributed, WENO_Order, dy,   &
                                             M0_p, M5_p, M2_p, M3_p, NV,    &
                                      Gamma, Q_p, U_p, V_p,                 &
                                             P_p, H_p, C_p, D_F_p,          &
                                              Flux_y=F_p,                   &
                                            Lambda_y=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(i,M0_p  :M2_p-1,:) = D_F_0(-M_s:-1 ,:)
      D_F(i,M3_p+1:M5_p  ,:) = D_F_N(   1:M_s,:)

      Flag_1D(M2_p:M3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(i,M2_p:M3_p,Index) = 1

      j = j_End - 1
    enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, P_p, H_p, C_p)

  END Subroutine Hybrid_2D_y_b 
!
!==========================================================================
!
  Subroutine Hybrid_3D_b    (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dx,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma, Q, U, V, W, P, H, C, F, D_F,         &
                             Lambda_Max,                                &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,  &
                             Flag, Flag_WENO,                           &
                             D_F_Epsilon, V_Moving_Frame,               &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5

  REALTYPE                                          :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)         :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)            :: U, V, W, P, H, C
  REALTYPE, dimension(:,:,:)                        :: Lambda_Max

  REALTYPE, dimension(3)                  , OPTIONAL :: V_Moving_Frame

  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,   K0:   K5,3),            &
                                                    OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  SELECT CASE (Index)
    CASE (1)
      call Hybrid_3D_x_b      (Index, Distributed, WENO_Order,            &
                                      Central_Order, MR_Order, dx,        &
                                      N0, N5, N2, N3,                     &
                                      M0, M5, M2, M3,                     &
                                      K0, K5, K2, K3, NV,                 &
                               Gamma, Q, U, V, W, P, H, C, F, D_F,        &
                               Lambda_Max,                                &
                               MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,  &
                               Flag, Flag_WENO,                           &
                               D_F_Epsilon, V_Moving_Frame,               &
                               MPI_Comm_Type)

    CASE (2)
      call Hybrid_3D_y_b      (Index, Distributed, WENO_Order,            &
                                      Central_Order, MR_Order, dx,        &
                                      N0, N5, N2, N3,                     &
                                      M0, M5, M2, M3,                     &
                                      K0, K5, K2, K3, NV,                 &
                               Gamma, Q, U, V, W, P, H, C, F, D_F,        &
                               Lambda_Max,                                &
                               MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,  &
                               Flag, Flag_WENO,                           &
                               D_F_Epsilon, V_Moving_Frame,               &
                               MPI_Comm_Type)

    CASE (3)
      call Hybrid_3D_z_b      (Index, Distributed, WENO_Order,            &
                                      Central_Order, MR_Order, dx,        &
                                      N0, N5, N2, N3,                     &
                                      M0, M5, M2, M3,                     &
                                      K0, K5, K2, K3, NV,                 &
                               Gamma, Q, U, V, W, P, H, C, F, D_F,        &
                               Lambda_Max,                                &
                               MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,  &
                               Flag, Flag_WENO,                           &
                               D_F_Epsilon, V_Moving_Frame,               &
                               MPI_Comm_Type)

  END SELECT

  END Subroutine Hybrid_3D_b
!
!==========================================================================
!
  Subroutine Hybrid_3D_x_b  (Index, Distributed, WENO_Order,            &
                                    Central_Order, MR_Order, dx,        &
                                    N0, N5, N2, N3,                     &
                                    M0, M5, M2, M3,                     &
                                    K0, K5, K2, K3, NV,                 &
                             Gamma, Q, U, V, W, P, H, C, F, D_F,        &
                             Lambda_Max,                                &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,  &
                             Flag, Flag_WENO,                           &
                             D_F_Epsilon, V_Moving_Frame,               &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dx

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5

  REALTYPE                                          :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , TARGET :: U, V, W, P, H, C
  REALTYPE, dimension(M0:M5,K0:K5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(3)                  , OPTIONAL :: V_Moving_Frame

  integer , dimension(   N0:   N5                          ) :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,   K0:   K5,3),            &
                                                    OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:,:)  , POINTER :: U_p, V_p, W_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:,:)  , POINTER :: Lambda_p

  integer  :: i, j, k, n
  integer  :: i_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p
  integer  :: K0_p, K2_p, K3_p, K5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,K2:K3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do k = K2,K3
    K0_p = k ; K2_p = k ; K3_p = k ; K5_p = k 
  
  do j = M2,M3
    M0_p = j ; M2_p = j ; M3_p = j ; M5_p = j 

    Flag_1D = Flag(N0:N5,j,k,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    i = N2-1 

    Do WHILE (i <= N3) 
      i = i+1

      if (Flag_1D(i) == 0) CYCLE
    
      do i_End = i+1,N3 ; if (Flag_1D(i_End) == 0) EXIT ; enddo

      call Hybrid_Range (N2, N3, i, i_End, M_s, M_Buffer, &
                              N0_p, N5_p, N2_p, N3_p)

      D_F_0(-M_s:-1 ,:) = D_F(N0_p  :N2_p-1,j,k,:)
      D_F_N(   1:M_s,:) = D_F(N3_p+1:N5_p  ,j,k,:)

           Q_p =>     Q(N0_p:N5_p,j:j,k:k,:)
           F_p =>     F(N0_p:N5_p,j:j,k:k,:)
         D_F_p =>   D_F(N0_p:N5_p,j:j,k:k,:)

           U_p =>     U(N0_p:N5_p,j:j,k:k)
           V_p =>     V(N0_p:N5_p,j:j,k:k)
           W_p =>     W(N0_p:N5_p,j:j,k:k)
           P_p =>     P(N0_p:N5_p,j:j,k:k)
           H_p =>     H(N0_p:N5_p,j:j,k:k)
           C_p =>     C(N0_p:N5_p,j:j,k:k)

      Lambda_p => Lambda_Max(j:j,k:k,:)

      call PS_WENO_Euler_Flux        (Index, Distributed, WENO_Order, dx,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Idx  , Distributed, WENO_Order, dx,   &
                                             M0_p, M5_p, M2_p, M3_p,        &
                                      Idx  , Distributed, WENO_Order, dx,   &
                                             K0_p, K5_p, K2_p, K3_p, NV,    &
                                      Gamma, Q_p, U_p, V_p, W_p,            &
                                             P_p, H_p, C_p, D_F_p,          &
                                              Flux_x=F_p,                   &
                                            Lambda_x=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(N0_p  :N2_p-1,j,k,:) = D_F_0(-M_s:-1 ,:)
      D_F(N3_p+1:N5_p  ,j,k,:) = D_F_N(   1:M_s,:)

      Flag_1D(N2_p:N3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(N2_p:N3_p,j,k,Index) = 1

      i = i_End-1
    enddo
  enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  END Subroutine Hybrid_3D_x_b
!
!==========================================================================
!
  Subroutine Hybrid_3D_y_b  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dy,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma, Q, U, V, W, P, H, C, F, D_F,         &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,   &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dy

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5

  REALTYPE                                          :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , TARGET :: U, V, W, P, H, C
  REALTYPE, dimension(N0:N5,K0:K5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(3)                  , OPTIONAL :: V_Moving_Frame

  integer , dimension(               M0:   M5              ) :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,   K0:   K5,3),            &
                                                    OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:,:)  , POINTER :: U_p, V_p, W_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:,:)  , POINTER :: Lambda_p

  integer  :: i, j, k, n
  integer  :: j_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p
  integer  :: K0_p, K2_p, K3_p, K5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,K2:K3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do k = K2,K3
    K0_p = k ; K2_p = k ; K3_p = k ; K5_p = k 
  
  do i = N2,N3
    N0_p = i ; N2_p = i ; N3_p = i ; N5_p = i 

    Flag_1D = Flag(i,M0:M5,k,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    j = M2-1 

    Do WHILE (j <= M3) 
      j = j+1

      if (Flag_1D(j) == 0) CYCLE
    
      do j_End = j+1,M3 ; if (Flag_1D(j_End) == 0) EXIT ; enddo

      call Hybrid_Range (M2, M3, j, j_End, M_s, M_Buffer, &
                              M0_p, M5_p, M2_p, M3_p)

      D_F_0(-M_s:-1 ,:) = D_F(i,M0_p  :M2_p-1,k,:)
      D_F_N(   1:M_s,:) = D_F(i,M3_p+1:M5_p  ,k,:)

           Q_p =>     Q(i:i,M0_p:M5_p,k:k,:)
           F_p =>     F(i:i,M0_p:M5_p,k:k,:)
         D_F_p =>   D_F(i:i,M0_p:M5_p,k:k,:)
  
           U_p =>     U(i:i,M0_p:M5_p,k:k)
           V_p =>     V(i:i,M0_p:M5_p,k:k)
           W_p =>     W(i:i,M0_p:M5_p,k:k)
           P_p =>     P(i:i,M0_p:M5_p,k:k)
           H_p =>     H(i:i,M0_p:M5_p,k:k)
           C_p =>     C(i:i,M0_p:M5_p,k:k)

      Lambda_p => Lambda_Max(i:i,k:k,:)

      call PS_WENO_Euler_Flux        (Idx  , Distributed, WENO_Order, dy,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Index, Distributed, WENO_Order, dy,   &
                                             M0_p, M5_p, M2_p, M3_p,        &
                                      Idx  , Distributed, WENO_Order, dy,   &
                                             K0_p, K5_p, K2_p, K3_p, NV,    &
                                      Gamma, Q_p, U_p, V_p, W_p,            &
                                             P_p, H_p, C_p, D_F_p,          &
                                              Flux_y=F_p,                   &
                                            Lambda_y=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(i,M0_p  :M2_p-1,k,:) = D_F_0(-M_s:-1 ,:)
      D_F(i,M3_p+1:M5_p  ,k,:) = D_F_N(   1:M_s,:)

      Flag_1D(M2_p:M3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(i,M2_p:M3_p,k,Index) = 1

      j = j_End - 1
    enddo
  enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  END Subroutine Hybrid_3D_y_b
!
!==========================================================================
!
  Subroutine Hybrid_3D_z_b  (Index, Distributed, WENO_Order,             &
                                    Central_Order, MR_Order, dz,         &
                                    N0, N5, N2, N3,                      &
                                    M0, M5, M2, M3,                      &
                                    K0, K5, K2, K3, NV,                  &
                             Gamma, Q, U, V, W, P, H, C, F, D_F,         &
                             Lambda_Max,                                 &
                             MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5,   &
                             Flag, Flag_WENO,                            &
                             D_F_Epsilon, V_Moving_Frame,                &
                             MPI_Comm_Type)

  integer  :: Index, WENO_Order, Central_Order, MR_Order
  logical  :: Distributed
  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  REALTYPE :: dz

  integer  :: NV

  integer  :: MR_N0, MR_N5, MR_M0, MR_M5, MR_K0, MR_K5

  REALTYPE                                          :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET :: Q, F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , TARGET :: U, V, W, P, H, C
  REALTYPE, dimension(N0:N5,M0:M5,NV)      , TARGET :: Lambda_Max

  REALTYPE, dimension(3)                  , OPTIONAL :: V_Moving_Frame

  integer , dimension(                           K0:   K5  ) :: Flag_1D
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) :: Flag
  integer , dimension(   N0:   N5,   M0:   M5,   K0:   K5,3),            &
                                                    OPTIONAL :: Flag_WENO

  REALTYPE, OPTIONAL :: D_F_Epsilon
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension( -(WENO_Order+1)/2:-1,NV) :: D_F_0
  REALTYPE, dimension(1:(WENO_Order+1)/2   ,NV) :: D_F_N

  REALTYPE, dimension(:,:,:,:), POINTER :: Q_p, F_p, D_F_p
  REALTYPE, dimension(:,:,:)  , POINTER :: U_p, V_p, W_p, P_p, H_p, C_p
  REALTYPE, dimension(:,:,:)  , POINTER :: Lambda_p

  integer  :: i, j, k, n
  integer  :: k_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p
  integer  :: M0_p, M2_p, M3_p, M5_p
  integer  :: K0_p, K2_p, K3_p, K5_p

  if (PRESENT(Flag_WENO)) Flag_WENO(:,:,:,Index) = 0

  if (MAXVAL(Flag(N2:N3,M2:M3,K2:K3,Index)) == 0) RETURN

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  call Hybrid_Buffer (WENO_Order, Central_Order, 0, M_s, M_Buffer)

  do j = M2,M3
    M0_p = j ; M2_p = j ; M3_p = j ; M5_p = j 
  
  do i = N2,N3
    N0_p = i ; N2_p = i ; N3_p = i ; N5_p = i 

    Flag_1D = Flag(i,j,K0:K5,Index) ; if (MAXVAL(Flag_1D) == 0) CYCLE

    k = K2-1  

    Do WHILE (k <= K3) 
      k = k+1

      if (Flag_1D(k) == 0) CYCLE
    
      do k_End = k+1,K3 ; if (Flag_1D(k_End) == 0) EXIT ; enddo

      call Hybrid_Range (K2, K3, k, k_End, M_s, M_Buffer, &
                              K0_p, K5_p, K2_p, K3_p)

      D_F_0(-M_s:-1 ,:) = D_F(i,j,K0_p  :K2_p-1,:)
      D_F_N(   1:M_s,:) = D_F(i,j,K3_p+1:K5_p  ,:)

           Q_p =>     Q(i:i,j:j,K0_p:K5_p,:)
           F_p =>     F(i:i,j:j,K0_p:K5_p,:)
         D_F_p =>   D_F(i:i,j:j,K0_p:K5_p,:)

           U_p =>     U(i:i,j:j,K0_p:K5_p)
           V_p =>     V(i:i,j:j,K0_p:K5_p)
           W_p =>     W(i:i,j:j,K0_p:K5_p)
           P_p =>     P(i:i,j:j,K0_p:K5_p)
           H_p =>     H(i:i,j:j,K0_p:K5_p)
           C_p =>     C(i:i,j:j,K0_p:K5_p)

      Lambda_p => Lambda_Max(i:i,j:j,:)

      call PS_WENO_Euler_Flux        (Idx  , Distributed, WENO_Order, dz,   &
                                             N0_p, N5_p, N2_p, N3_p,        &
                                      Idx  , Distributed, WENO_Order, dz,   &
                                             M0_p, M5_p, M2_p, M3_p,        &
                                      Index, Distributed, WENO_Order, dz,   &
                                             K0_p, K5_p, K2_p, K3_p, NV,    &
                                      Gamma, Q_p, U_p, V_p, W_p,            &
                                             P_p, H_p, C_p, D_F_p,          &
                                              Flux_z=F_p,                   &
                                            Lambda_z=Lambda_p,              &
                                         D_F_Epsilon=D_F_Epsilon,           &
                                      V_Moving_Frame=V_Moving_Frame,        &
                                       MPI_Comm_Type=MPI_Comm_Type)

      D_F(i,j,K0_p  :K2_p-1,:) = D_F_0(-M_s:-1 ,:)
      D_F(i,j,K3_p+1:K5_p  ,:) = D_F_N(   1:M_s,:)

      Flag_1D(K2_p:K3_p) = 1

      if (PRESENT(Flag_WENO)) Flag_WENO(i,j,K2_p:K3_p,Index) = 1

      k = k_End-1
    enddo
  enddo
  enddo

  Nullify (Q_p, F_p, D_F_p, U_p, V_p, W_p, P_p, H_p, C_p)

  END Subroutine Hybrid_3D_z_b
!
!==========================================================================
!
  Subroutine Hybrid_Buffer (WENO_Order, Central_Order, MR_Order, &
                                 M_s, M_Buffer)

  integer  :: WENO_Order, Central_Order, MR_Order, M_Buffer
  integer  :: M_s, M_c, M_m

  M_s      = (   WENO_Order+1)/2
  M_c      = (Central_Order+1)/2
  M_m      =       MR_Order

#if 0
  M_Buffer = MAX(MAX(M_s,M_c),M_m)
#else
  M_Buffer =     MAX(M_s,M_c)
#endif

  END Subroutine Hybrid_Buffer 
!
!==========================================================================
!
  Subroutine Hybrid_Range (N2, N3, i, i_End, M_s, M_Buffer, &
                               N0_p, N5_p, N2_p, N3_p)

  integer  :: N2, N3, i, i_End, M_s, M_Buffer
  integer  :: N0_p, N2_p, N3_p, N5_p

  N2_p = i-M_Buffer    ; N3_p = i_End-1+M_Buffer
  N2_p = MAX(N2,N2_p)  ; N3_p = MIN(N3,N3_p)
  N0_p = N2_p-M_s      ; N5_p = N3_p+M_s

  i_End = N3_p+1

  END Subroutine Hybrid_Range 

END MODULE Hybrid_WENO_Euler_Flux

#if 0
write (6,*) j, i, i_End, Flag_1D(i-1), Flag_1D(i), Flag_1D(i+1)    
write (6,*) j, i, i_End, Flag_1D(i_End-1), Flag_1D(i_End)    
write (6,100) 'x', j, i, i_End, M_Buffer, N0_p, N2_p, N3_p, N5_p
write (6,100) 'y', i, j, j_End, M_Buffer, M0_p, M2_p, M3_p, M5_p

100  format (1x,a,' : ',10(i5,:))
#endif
