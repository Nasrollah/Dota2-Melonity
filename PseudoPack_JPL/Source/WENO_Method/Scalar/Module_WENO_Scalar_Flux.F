#include "pseudopack.h"


MODULE WENO_Scalar_Flux

#if defined (PARALLEL_MPI)
  USE Processor
#endif
  USE WENO_Differentiation

implicit NONE

INTERFACE PS_WENO_Scalar_Flux
  MODULE PROCEDURE PS_WENO_Scalar_Flux_1D
  MODULE PROCEDURE PS_WENO_Scalar_Flux_2D
  MODULE PROCEDURE PS_WENO_Scalar_Flux_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_WENO_Scalar_Flux

PUBLIC :: PS_WENO_Scalar_Flux_1D
PUBLIC :: PS_WENO_Scalar_Flux_1D_1

PUBLIC :: PS_WENO_Scalar_Flux_2D
PUBLIC :: PS_WENO_Scalar_Flux_2D_1
PUBLIC :: PS_WENO_Scalar_Flux_2D_2

PUBLIC :: PS_WENO_Scalar_Flux_3D
PUBLIC :: PS_WENO_Scalar_Flux_3D_1
PUBLIC :: PS_WENO_Scalar_Flux_3D_2
PUBLIC :: PS_WENO_Scalar_Flux_3D_3

CONTAINS

!
!===================================================================
!
  Subroutine PS_WENO_Scalar_Flux_1D (Dimensions,                            &
                                     Index_x, Distributed_x, Order_x, dx,   &
                                              N0, N5, N2, N3, NV,           &
                                     Q, D_Flux, MPI_Comm_Type)

  integer  :: Dimensions 
  integer  :: Index_x, Order_x
  logical  :: Distributed_x
  integer  :: N0, N5, N2, N3
  REALTYPE :: dx

  integer  :: NV

  integer ,                        OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,NV) :: Q, D_Flux

  call PS_WENO_Scalar_Flux_1D_1 (Dimensions,                            &
                                 Index_x, Distributed_x, Order_x, dx,   &
                                          N0, N5, N2, N3, NV,           &
                                 Q, D_Flux, MPI_Comm_Type)

  END Subroutine PS_WENO_Scalar_Flux_1D

!
!===================================================================
!
  Subroutine PS_WENO_Scalar_Flux_2D (Dimensions,                            &
                                     Index_x, Distributed_x, Order_x, dx,   &
                                              N0, N5, N2, N3,               &
                                     Index_y, Distributed_y, Order_y, dy,   &
                                              M0, M5, M2, M3, NV,           &
                                     Q, D_Flux, MPI_Comm_Type)

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

  integer ,                        OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, D_Flux

  SELECT CASE (Index_x)
    CASE DEFAULT
      call PS_WENO_Scalar_Flux_2D_1 (Dimensions,                            &
                                     Index_x, Distributed_x, Order_x, dx,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3, NV,           &
                                      Q, D_Flux, MPI_Comm_Type)
    CASE (2)
      call PS_WENO_Scalar_Flux_2D_2 (Dimensions,                            &
                                     Index_x, Distributed_x, Order_x, dx,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3, NV,           &
                                      Q, D_Flux, MPI_Comm_Type)
  END SELECT

  SELECT CASE (Index_y)
    CASE (1)
      call PS_WENO_Scalar_Flux_2D_1 (Dimensions,                            &
                                     Index_y, Distributed_y, Order_y, dy,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3, NV,           &
                                      Q, D_Flux, MPI_Comm_Type)

    CASE DEFAULT
      call PS_WENO_Scalar_Flux_2D_2 (Dimensions,                            &
                                     Index_y, Distributed_y, Order_y, dy,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3, NV,           &
                                      Q, D_Flux, MPI_Comm_Type)
  END SELECT

  END Subroutine PS_WENO_Scalar_Flux_2D
!
!===================================================================
!
  Subroutine PS_WENO_Scalar_Flux_3D (Dimensions,                            &
                                     Index_x, Distributed_x, Order_x, dx,   &
                                              N0, N5, N2, N3,               &
                                     Index_y, Distributed_y, Order_y, dy,   &
                                              M0, M5, M2, M3,               &
                                     Index_z, Distributed_z, Order_z, dz,   &
                                              K0, K5, K2, K3, NV,           &
                                     Q, D_Flux, MPI_Comm_Type)

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

  integer ,                        OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Flux

  SELECT CASE (Index_x)
    CASE DEFAULT
      call PS_WENO_Scalar_Flux_3D_1 (Dimensions,                            &
                                     Index_x, Distributed_x, Order_x, dx,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3,               &
                                              K0, K5, K2, K3, NV,           &
                                      Q, D_Flux, MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Scalar_Flux_3D_2 (Dimensions,                            &
                                     Index_x, Distributed_x, Order_x, dx,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3,               &
                                              K0, K5, K2, K3, NV,           &
                                      Q, D_Flux, MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Scalar_Flux_3D_3 (Dimensions,                            &
                                     Index_x, Distributed_x, Order_x, dx,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3,               &
                                              K0, K5, K2, K3, NV,           &
                                      Q, D_Flux, MPI_Comm_Type)
  END SELECT

  SELECT CASE (Index_y)
    CASE (1)
      call PS_WENO_Scalar_Flux_3D_1 (Dimensions,                            &
                                     Index_y, Distributed_y, Order_y, dy,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3,               &
                                              K0, K5, K2, K3, NV,           &
                                     Q, D_Flux, MPI_Comm_Type)

    CASE DEFAULT
      call PS_WENO_Scalar_Flux_3D_2 (Dimensions,                            &
                                     Index_y, Distributed_y, Order_y, dy,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3,               &
                                              K0, K5, K2, K3, NV,           &
                                     Q, D_Flux, MPI_Comm_Type)

    CASE (3)
      call PS_WENO_Scalar_Flux_3D_3 (Dimensions,                            &
                                     Index_y, Distributed_y, Order_y, dy,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3,               &
                                              K0, K5, K2, K3, NV,           &
                                     Q, D_Flux, MPI_Comm_Type)
  END SELECT

  SELECT CASE (Index_z)
    CASE (1)
      call PS_WENO_Scalar_Flux_3D_1 (Dimensions,                            &
                                     Index_z, Distributed_z, Order_z, dz,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3,               &
                                              K0, K5, K2, K3, NV,           &
                                     Q, D_Flux, MPI_Comm_Type)

    CASE (2)
      call PS_WENO_Scalar_Flux_3D_2 (Dimensions,                            &
                                     Index_z, Distributed_z, Order_z, dz,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3,               &
                                              K0, K5, K2, K3, NV,           &
                                     Q, D_Flux, MPI_Comm_Type)

    CASE DEFAULT
      call PS_WENO_Scalar_Flux_3D_3 (Dimensions,                            &
                                     Index_z, Distributed_z, Order_z, dz,   &
                                              N0, N5, N2, N3,               &
                                              M0, M5, M2, M3,               &
                                              K0, K5, K2, K3, NV,           &
                                     Q, D_Flux, MPI_Comm_Type)
  END SELECT

  END Subroutine PS_WENO_Scalar_Flux_3D 
!
!===================================================================
!
  Subroutine PS_WENO_Scalar_Flux_1D_1 (Dimensions,                    &
                                       Index, Distributed, Order, dx, &
                                       N0, N5, i_0, i_1, NV,          &
                                       Q, D_Flux, MPI_Comm_type)               

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: NV
  REALTYPE :: dx

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,NV) :: Q, D_Flux

  integer  :: n

  do n = 1,NV
    call PS_Diff_WENO   (Index, Order, N0, N5, i_0, i_1, dx,             &
                                Q(:,n), D_Flux(:,n))
  enddo

  END Subroutine PS_WENO_Scalar_Flux_1D_1
!
!===================================================================
!
  Subroutine PS_WENO_Scalar_Flux_2D_1 (Dimensions,                    &
                                       Index, Distributed, Order, dx, &
                                       N0, N5, i_0, i_1,              &
                                       M0, M5, j_0, j_1, NV,          &
                                       Q, D_Flux, MPI_Comm_type)               

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: NV
  REALTYPE :: dx

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, D_Flux

  integer  :: n

  do n = 1,NV
    call PS_Diff_WENO   (Index, Order, N0, N5, i_0, i_1, dx,             &
                                Order, M0, M5, j_0, j_1, dx,             &
                                Q(:,:,n), D_Flux(:,:,n))
  enddo

  END Subroutine PS_WENO_Scalar_Flux_2D_1
!
!===================================================================
!
  Subroutine PS_WENO_Scalar_Flux_2D_2 (Dimensions,                    &
                                       Index, Distributed, Order, dy, &
                                       N0, N5, i_0, i_1,              &
                                       M0, M5, j_0, j_1, NV,          &
                                       Q, D_Flux, MPI_Comm_type)               

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: NV
  REALTYPE :: dy

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, D_Flux

  integer  :: n

  do n = 1,NV
    call PS_Diff_WENO   (Index, Order, N0, N5, i_0, i_1, dy,             &
                                Order, M0, M5, j_0, j_1, dy,             &
                                Q(:,:,n), D_Flux(:,:,n))
  enddo

  END Subroutine PS_WENO_Scalar_Flux_2D_2
!
!===================================================================
!
  Subroutine PS_WENO_Scalar_Flux_3D_1 (Dimensions,                    &
                                       Index, Distributed, Order, dx, &
                                       N0, N5, i_0, i_1,              &
                                       M0, M5, j_0, j_1,              &
                                       K0, K5, k_0, k_1, NV,          &
                                       Q, D_Flux, MPI_Comm_type)               

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1
  integer  :: NV
  REALTYPE :: dx

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Flux

  integer  :: n

  do n = 1,NV
    call PS_Diff_WENO   (Index, Order, N0, N5, i_0, i_1, dx,             &
                                Order, M0, M5, j_0, j_1, dx,             &
                                Order, K0, K5, k_0, k_1, dx,             &
                                Q(:,:,:,n), D_Flux(:,:,:,n))
  enddo

  END Subroutine PS_WENO_Scalar_Flux_3D_1
!
!===================================================================
!
  Subroutine PS_WENO_Scalar_Flux_3D_2 (Dimensions,                    &
                                       Index, Distributed, Order, dy, &
                                       N0, N5, i_0, i_1,              &
                                       M0, M5, j_0, j_1,              &
                                       K0, K5, k_0, k_1, NV,          &
                                       Q, D_Flux, MPI_Comm_type)               

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1
  integer  :: NV
  REALTYPE :: dy

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Flux

  integer  :: n

  do n = 1,NV
    call PS_Diff_WENO   (Index, Order, N0, N5, i_0, i_1, dy,             &
                                Order, M0, M5, j_0, j_1, dy,             &
                                Order, K0, K5, k_0, k_1, dy,             &
                                Q(:,:,:,n), D_Flux(:,:,:,n))
  enddo

  END Subroutine PS_WENO_Scalar_Flux_3D_2
!
!===================================================================
!
  Subroutine PS_WENO_Scalar_Flux_3D_3 (Dimensions,                    &
                                       Index, Distributed, Order, dz, &
                                       N0, N5, i_0, i_1,              &
                                       M0, M5, j_0, j_1,              &
                                       K0, K5, k_0, k_1, NV,          &
                                       Q, D_Flux, MPI_Comm_type)               

  integer  :: Dimensions, Index, Order
  logical  :: Distributed
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1
  integer  :: NV
  REALTYPE :: dz

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Flux

  integer  :: n

  do n = 1,NV
    call PS_Diff_WENO   (Index, Order, N0, N5, i_0, i_1, dz,             &
                                Order, M0, M5, j_0, j_1, dz,             &
                                Order, K0, K5, k_0, k_1, dz,             &
                                Q(:,:,:,n), D_Flux(:,:,:,n))
  enddo

  END Subroutine PS_WENO_Scalar_Flux_3D_3

END MODULE WENO_Scalar_Flux
