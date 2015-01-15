#include "pseudopack.h"


MODULE WENO_Euler_Flux_General

  USE Euler_Function
  USE WENO_Euler_Flux_123D

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
PUBLIC ::    WENO_Euler_Flux_2D_P_a

CONTAINS
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_C_a (Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,  NV,        &
                                     Gamma_Gas, Gamma,                    &
                                     Q, D_Flux,                           &
                                       Flux_x,                            &
                                     Lambda_x,                            &
                                     D_F_Epsilon,                         &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  REALTYPE, dimension(:,:)   :: Gamma_Gas
  REALTYPE, dimension(N0:N5) :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(1)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Flux_x

  REALTYPE, dimension(N0:N5,NV)           :: Q, D_Flux
  REALTYPE, dimension(N0:N5)              :: U, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, NV,                   &
                                    Gamma, Q, U, P, H, C)

  call WENO_Euler_Flux_1D_P_a      (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, P, H, C, D_Flux,               &
                                      Flux_x,                            &
                                    Lambda_x,                            &
                                    D_F_Epsilon,                         &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_1D_C_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_P_a (Index_x, Distributed_x, Order_x, dx, &
                                            N0, N5, N2, N3,  NV,          &
                                     Gamma_Gas, Gamma,                    &
                                     Q, U, P, H, C, D_Flux,               &
                                       Flux_x,                            &
                                     Lambda_x,                            &
                                     D_F_Epsilon,                         &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  REALTYPE, dimension(:,:)   :: Gamma_Gas
  REALTYPE, dimension(N0:N5) :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(1)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Flux_x

  REALTYPE, dimension(N0:N5,NV)           :: Q, Flux, D_Flux
  REALTYPE, dimension(N0:N5)              :: U, P, H, C
!
! ---------
      if (PRESENT(Flux_x)) then
        Flux = Flux_x
      else
        call PS_Euler_Flux         (Index_x, N0, N5, NV,                 &
                                             Q, U, P, Flux,              &
                                    V_Moving_Frame)
      endif

      call PS_WENO_Euler_Flux_1D_x (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, P, H, C, Flux, D_Flux,         &
                                    Lambda_x,                            &
                                    V_Moving_Frame, MPI_Comm_Type)

      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO 
        endif
      endif

  END Subroutine WENO_Euler_Flux_1D_P_a 

!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_C_a (Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,             &
                                     Index_y, Distributed_y, Order_y, dy, &
                                              M0, M5, M2, M3,  NV,        &
                                     Gamma_Gas, Gamma,                    &
                                     Q, D_Flux,                           &
                                       Flux_x,   Flux_y,                  &
                                     Lambda_x, Lambda_y,                  &
                                     D_F_Epsilon,                         &
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

  REALTYPE, dimension(:,:)         :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(2)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(M0:M5,NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,NV)      , OPTIONAL :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Flux_x, Flux_y

  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, NV,           &
                                    Gamma, Q, U, V, P, H, C)

  call WENO_Euler_Flux_2D_P_a      (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                    Index_y, Distributed_y, Order_y, dy, &
                                           M0, M5, M2, M3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, P, H, C, D_Flux,            &
                                      Flux_x,   Flux_y,                  &
                                    Lambda_x, Lambda_y,                  &
                                    D_F_Epsilon,                         &
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
                                     Q, U, V, P, H, C, D_Flux,            &
                                       Flux_x,   Flux_y,                  &
                                     Lambda_x, Lambda_y,                  &
                                     D_F_Epsilon,                         &
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

  REALTYPE, dimension(:,:)         :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5) :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(2)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(M0:M5,NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,NV)      , OPTIONAL :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Flux_x, Flux_y

  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Q, Flux, D_Flux, D_G
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V, P, H, C

  Flux = ZERO ; D_G = ZERO

!
! ---------
  SELECT CASE (Index_x)
    CASE (1,2)
      if (PRESENT(Flux_x)) then
        Flux = Flux_x
      else
        call PS_Euler_Flux         (Index_x, N0, N5, M0, M5, NV,         &
                                             Q, U, V, P, Flux,           &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_x)
    CASE (1)
      call PS_WENO_Euler_Flux_2D_x (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, P, H, C, Flux, D_Flux,      &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_2D_y (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, P, H, C, Flux, D_Flux,      &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)
  END SELECT 
  
#if defined (DEBUG_HOPE)
#include "WENO_Euler_Flux_Debug_x.i"
!RETURN
#endif

  SELECT CASE (Index_x)
    CASE (1,2)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO 
        endif
      endif

  END SELECT

!
! ---------
  SELECT CASE (Index_y)
    CASE (1,2)
      if (PRESENT(Flux_y)) then
        Flux = Flux_y
      else
        call PS_Euler_Flux         (Index_y, N0, N5, M0, M5, NV,         &
                                             Q, U, V, P, Flux,           &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_y)
    CASE (1)
      call PS_WENO_Euler_Flux_2D_x (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, P, H, C, Flux, D_G,         &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_2D_y (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, P, H, C, Flux, D_G,         &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)

  END SELECT

  SELECT CASE (Index_y)
    CASE (1,2)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO 
        endif
      endif

  END SELECT
!
! ---------
  SELECT CASE (Index_y)
    CASE (1,2)
      if (Index_x /= 0) then
        D_Flux = D_Flux + D_G
      else
        D_Flux = D_G
      endif

  END SELECT

#if defined (DEBUG_HOPE)
#include "WENO_Euler_Flux_Debug_y.i"
#endif

  END Subroutine WENO_Euler_Flux_2D_P_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_C_a (Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,             &
                                     Index_y, Distributed_y, Order_y, dy, &
                                              M0, M5, M2, M3,             &
                                     Index_z, Distributed_z, Order_z, dz, &
                                              K0, K5, K2, K3,  NV,        &
                                     Gamma_Gas, Gamma,                    &
                                     Q, D_Flux,                           &
                                       Flux_x,   Flux_y,   Flux_z,        &
                                     Lambda_x, Lambda_y, Lambda_z,        &
                                     D_F_Epsilon,                         &
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

  REALTYPE, dimension(:,:)               :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(3)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(M0:M5,K0:K5,NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,K0:K5,NV)      , OPTIONAL :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV)      , OPTIONAL :: Lambda_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), OPTIONAL :: Flux_x, Flux_y, Flux_z

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, K0, K5, NV,   &
                                    Gamma, Q, U, V, W, P, H, C)

  call WENO_Euler_Flux_3D_P_a      (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                    Index_y, Distributed_y, Order_y, dy, &
                                           M0, M5, M2, M3,               &
                                    Index_z, Distributed_z, Order_z, dz, &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, D_Flux,         &
                                      Flux_x,   Flux_y,   Flux_z,        &
                                    Lambda_x, Lambda_y, Lambda_z,        &
                                    D_F_Epsilon,                         &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_3D_C_a
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_P_a (Index_x, Distributed_x, Order_x, dx, &
                                            N0, N5, N2, N3,               &
                                     Index_y, Distributed_y, Order_y, dy, &
                                            M0, M5, M2, M3,               &
                                     Index_z, Distributed_z, Order_z, dz, &
                                            K0, K5, K2, K3,  NV,          &
                                     Gamma_Gas, Gamma,                    &
                                     Q, U, V, W, P, H, C, D_Flux,         &
                                       Flux_x,   Flux_y,   Flux_z,        &
                                     Lambda_x, Lambda_y, Lambda_z,        &
                                     D_F_Epsilon,                         &
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

  REALTYPE, dimension(:,:)               :: Gamma_Gas
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(3)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(M0:M5,K0:K5,NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,K0:K5,NV)      , OPTIONAL :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV)      , OPTIONAL :: Lambda_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), OPTIONAL :: Flux_x, Flux_y, Flux_z

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, Flux, D_Flux, D_G
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, P, H, C
!
! ---------
  SELECT CASE (Index_x)
    CASE (1,2,3)
      if (PRESENT(Flux_x)) then
        Flux = Flux_x
      else
        call PS_Euler_Flux         (Index_x, N0, N5, M0, M5, K0, K5, NV, &
                                             Q, U, V, W, P, Flux,        &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_x)
    CASE (0)
      D_Flux = ZERO

    CASE (1)
      call PS_WENO_Euler_Flux_3D_x (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_Flux,   &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D_y (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_Flux,   &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D_z (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_Flux,   &
                                    Lambda_z, V_Moving_Frame,            &
                                    MPI_Comm_Type)

  END SELECT 
  
  SELECT CASE (Index_x)
    CASE (1,2,3)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO 
        endif
      endif

  END SELECT

!
! ---------
  SELECT CASE (Index_y)
    CASE (1,2,3)
      if (PRESENT(Flux_y)) then
        Flux = Flux_y
      else
        call PS_Euler_Flux         (Index_y, N0, N5, M0, M5, K0, K5, NV, &
                                             Q, U, V, W, P, Flux,        &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_y)
    CASE (0)
      D_G = ZERO

    CASE (1)
      call PS_WENO_Euler_Flux_3D_x (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D_y (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D_z (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_z, V_Moving_Frame,            &
                                    MPI_Comm_Type)

  END SELECT

  SELECT CASE (Index_y)
    CASE (1,2,3)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO 
        endif
      endif

  END SELECT
!
! ---------
  D_Flux = D_Flux + D_G

!
! ---------
  SELECT CASE (Index_z)
    CASE (1,2,3)
      if (PRESENT(Flux_z)) then
        Flux = Flux_z
      else
        call PS_Euler_Flux         (Index_z, N0, N5, M0, M5, K0, K5, NV, &
                                             Q, U, V, W, P, Flux,        &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_z)
    CASE (0)
      D_G = ZERO
      
    CASE (1)
      call PS_WENO_Euler_Flux_3D_x (Index_z, Distributed_z, Order_z, dz, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D_y (Index_z, Distributed_z, Order_z, dz, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D_z (Index_z, Distributed_z, Order_z, dz, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                    Gamma_Gas, Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_z, V_Moving_Frame,            &
                                    MPI_Comm_Type)

  END SELECT

  SELECT CASE (Index_z)
    CASE (1,2,3)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO 
        endif
      endif

  END SELECT
!
! ---------
  D_Flux = D_Flux + D_G

  END Subroutine WENO_Euler_Flux_3D_P_a

!===================
! Case : Constant Gamma
!===================
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_C_b (Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,  NV,        &
                                                Gamma,                    &
                                     Q, D_Flux,                           &
                                       Flux_x,                            &
                                     Lambda_x,                            &
                                     D_F_Epsilon,                         &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  REALTYPE :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(1)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Flux_x

  REALTYPE, dimension(N0:N5,NV)           :: Q, D_Flux
  REALTYPE, dimension(N0:N5)              :: U, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, NV,                   &
                                    Gamma, Q, U, P, H, C)

  call WENO_Euler_Flux_1D_P_b      (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,  NV,          &
                                               Gamma,                    &
                                    Q, U, P, H, C, D_Flux,               &
                                      Flux_x,                            &
                                    Lambda_x,                            &
                                    D_F_Epsilon,                         &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_1D_C_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_1D_P_b (Index_x, Distributed_x, Order_x, dx, &
                                            N0, N5, N2, N3,  NV,          &
                                                Gamma,                    &
                                     Q, U, P, H, C, D_Flux,               &
                                       Flux_x,                            &
                                     Lambda_x,                            &
                                     D_F_Epsilon,                         &
                                     V_Moving_Frame, MPI_Comm_Type)

  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  REALTYPE :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(1)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,NV), OPTIONAL :: Flux_x

  REALTYPE, dimension(N0:N5,NV)           :: Q, Flux, D_Flux
  REALTYPE, dimension(N0:N5)              :: U, P, H, C
!
! ---------
      if (PRESENT(Flux_x)) then
        Flux = Flux_x
      else
        call PS_Euler_Flux         (Index_x, N0, N5, NV,                 &
                                             Q, U, P, Flux,              &
                                    V_Moving_Frame)
      endif

      call PS_WENO_Euler_Flux_1D_x (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,  NV,          &
                                               Gamma,                    &
                                    Q, U, P, H, C, Flux, D_Flux,         &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO 
        endif
      endif

  END Subroutine WENO_Euler_Flux_1D_P_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_C_b (Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,             &
                                     Index_y, Distributed_y, Order_y, dy, &
                                              M0, M5, M2, M3,  NV,        &
                                                Gamma,                    &
                                     Q, D_Flux,                           &
                                       Flux_x,   Flux_y,                  &
                                     Lambda_x, Lambda_y,                  &
                                     D_F_Epsilon,                         &
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

  REALTYPE :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(2)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(M0:M5,NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,NV)      , OPTIONAL :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Flux_x
  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Flux_y

  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, NV,           &
                                    Gamma, Q, U, V, P, H, C)

  call WENO_Euler_Flux_2D_P_b      (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                    Index_y, Distributed_y, Order_y, dy, &
                                           M0, M5, M2, M3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, P, H, C, D_Flux,            &
                                      Flux_x,   Flux_y,                  &
                                    Lambda_x, Lambda_y,                  &
                                    D_F_Epsilon,                         &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_2D_C_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_2D_P_b (Index_x, Distributed_x, Order_x, dx, &
                                            N0, N5, N2, N3,               &
                                     Index_y, Distributed_y, Order_y, dy, &
                                            M0, M5, M2, M3,  NV,          &
                                                Gamma,                    &
                                     Q, U, V, P, H, C, D_Flux,            &
                                       Flux_x,   Flux_y,                  &
                                     Lambda_x, Lambda_y,                  &
                                     D_F_Epsilon,                         &
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

  REALTYPE :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(2)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(M0:M5,NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,NV)      , OPTIONAL :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV), OPTIONAL :: Flux_x, Flux_y

  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Q, Flux, D_Flux, D_G
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V, P, H, C

!
! ---------
  SELECT CASE (Index_x)
    CASE (1,2)
      if (PRESENT(Flux_x)) then
        Flux = Flux_x
      else
        call PS_Euler_Flux         (Index_x, N0, N5, M0, M5, NV,         &
                                             Q, U, V, P, Flux,           &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_x)
    CASE (1)
      call PS_WENO_Euler_Flux_2D_x (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, P, H, C, Flux, D_Flux,      &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_2D_y (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, P, H, C, Flux, D_Flux,      &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)
  END SELECT 
  
  SELECT CASE (Index_x)
    CASE (1,2)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO 
        endif
      endif

  END SELECT

!
! ---------
  SELECT CASE (Index_y)
    CASE (1,2)
      if (PRESENT(Flux_y)) then
        Flux = Flux_y
      else
        call PS_Euler_Flux         (Index_y, N0, N5, M0, M5, NV,         &
                                             Q, U, V, P, Flux,           &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_y)
    CASE (1)
      call PS_WENO_Euler_Flux_2D_x (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, P, H, C, Flux, D_G,         &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_2D_y (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, P, H, C, Flux, D_G,         &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)

  END SELECT

  SELECT CASE (Index_y)
    CASE (1,2)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO 
        endif
      endif

  END SELECT
!
! ---------
  SELECT CASE (Index_y)
    CASE (1,2)
      if (Index_x /= 0) then
        D_Flux = D_Flux + D_G
      else
        D_Flux = D_G
      endif

  END SELECT

  END Subroutine WENO_Euler_Flux_2D_P_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_C_b (Index_x, Distributed_x, Order_x, dx, &
                                              N0, N5, N2, N3,             &
                                     Index_y, Distributed_y, Order_y, dy, &
                                              M0, M5, M2, M3,             &
                                     Index_z, Distributed_z, Order_z, dz, &
                                              K0, K5, K2, K3,  NV,        &
                                                Gamma,                    &
                                     Q, D_Flux,                           &
                                       Flux_x,   Flux_y,   Flux_z,        &
                                     Lambda_x, Lambda_y, Lambda_z,        &
                                     D_F_Epsilon,                         &
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

  REALTYPE :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(3)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(M0:M5,K0:K5,NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,K0:K5,NV)      , OPTIONAL :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV)      , OPTIONAL :: Lambda_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), OPTIONAL :: Flux_x, Flux_y, Flux_z

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, P, H, C

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, K0, K5, NV,   &
                                    Gamma, Q, U, V, W, P, H, C)

  call WENO_Euler_Flux_3D_P_b      (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                    Index_y, Distributed_y, Order_y, dy, &
                                           M0, M5, M2, M3,               &
                                    Index_z, Distributed_z, Order_z, dz, &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, D_Flux,         &
                                      Flux_x,   Flux_y,   Flux_z,        &
                                    Lambda_x, Lambda_y, Lambda_z,        &
                                    D_F_Epsilon,                         &
                                    V_Moving_Frame, MPI_Comm_Type)

  END Subroutine WENO_Euler_Flux_3D_C_b
!
!===================================================================
!
  Subroutine WENO_Euler_Flux_3D_P_b (Index_x, Distributed_x, Order_x, dx, &
                                            N0, N5, N2, N3,               &
                                     Index_y, Distributed_y, Order_y, dy, &
                                            M0, M5, M2, M3,               &
                                     Index_z, Distributed_z, Order_z, dz, &
                                            K0, K5, K2, K3,  NV,          &
                                                Gamma,                    &
                                     Q, U, V, W, P, H, C, D_Flux,         &
                                       Flux_x,   Flux_y,   Flux_z,        &
                                     Lambda_x, Lambda_y, Lambda_z,        &
                                     D_F_Epsilon,                         &
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

  REALTYPE :: Gamma

  REALTYPE,                 OPTIONAL :: D_F_Epsilon
  REALTYPE, dimension(3)  , OPTIONAL :: V_Moving_Frame
  integer ,                 OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(M0:M5,K0:K5,NV)      , OPTIONAL :: Lambda_x
  REALTYPE, dimension(N0:N5,K0:K5,NV)      , OPTIONAL :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV)      , OPTIONAL :: Lambda_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), OPTIONAL :: Flux_x, Flux_y, Flux_z

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, Flux, D_Flux, D_G
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, P, H, C
!
! ---------
  SELECT CASE (Index_x)
    CASE (1,2,3)
      if (PRESENT(Flux_x)) then
        Flux = Flux_x
      else
        call PS_Euler_Flux         (Index_x, N0, N5, M0, M5, K0, K5, NV, &
                                             Q, U, V, W, P, Flux,        &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_x)
    CASE (0)
      D_Flux = ZERO

    CASE (1)
      call PS_WENO_Euler_Flux_3D_x (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_Flux,   &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D_y (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_Flux,   &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D_z (Index_x, Distributed_x, Order_x, dx, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_Flux,   &
                                    Lambda_z, V_Moving_Frame,            &
                                    MPI_Comm_Type)

  END SELECT 
  
  SELECT CASE (Index_x)
    CASE (1,2,3)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO 
        endif
      endif

  END SELECT
!
! ---------
  SELECT CASE (Index_y)
    CASE (1,2,3)
      if (PRESENT(Flux_y)) then
        Flux = Flux_y
      else
        call PS_Euler_Flux         (Index_y, N0, N5, M0, M5, K0, K5, NV, &
                                             Q, U, V, W, P, Flux,        &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_y)
    CASE (0)
      D_G = ZERO

    CASE (1)
      call PS_WENO_Euler_Flux_3D_x (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D_y (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D_z (Index_y, Distributed_y, Order_y, dy, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_z, V_Moving_Frame,            &
                                    MPI_Comm_Type)

  END SELECT

  SELECT CASE (Index_y)
    CASE (1,2,3)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO 
        endif
      endif

  END SELECT
!
! ---------
  D_Flux = D_Flux + D_G

!
! ---------
  SELECT CASE (Index_z)
    CASE (1,2,3)
      if (PRESENT(Flux_z)) then
        Flux = Flux_z
      else
        call PS_Euler_Flux         (Index_z, N0, N5, M0, M5, K0, K5, NV, &
                                             Q, U, V, W, P, Flux,        &
                                    V_Moving_Frame)
      endif

  END SELECT

  SELECT CASE (Index_z)
    CASE (0)
      D_G = ZERO
      
    CASE (1)
      call PS_WENO_Euler_Flux_3D_x (Index_z, Distributed_z, Order_z, dz, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_x, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Euler_Flux_3D_y (Index_z, Distributed_z, Order_z, dz, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_y, V_Moving_Frame,            &
                                    MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Euler_Flux_3D_z (Index_z, Distributed_z, Order_z, dz, &
                                           N0, N5, N2, N3,               &
                                           M0, M5, M2, M3,               &
                                           K0, K5, K2, K3,  NV,          &
                                               Gamma,                    &
                                    Q, U, V, W, P, H, C, Flux, D_G,      &
                                    Lambda_z, V_Moving_Frame,            &
                                    MPI_Comm_Type)

  END SELECT

  SELECT CASE (Index_z)
    CASE (1,2,3)
      if (PRESENT(D_F_Epsilon)) then
        if (D_F_Epsilon /= ZERO) then
          WHERE (ABS(D_G   ) <  D_F_Epsilon) D_G    = ZERO 
        endif
      endif

  END SELECT
!
! ---------
  D_Flux = D_Flux + D_G

  END Subroutine WENO_Euler_Flux_3D_P_b

END MODULE WENO_Euler_Flux_General
