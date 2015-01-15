#include "pseudopack.h"

! ----------------------------------------------------------------------
! FUNCTION NAME: PS_Scalar_Differentiation
! AUTHOR       : Wai-Sun Don
! DESCRIPTION  : Hybrid Differentiation of a scalar data.
!               
! Index      The index of the multi-dimensional array F to be worked on 
!            (Input)
! Order      Order of the WENO/Hybrid scheme (3, 5)    (Input)
! N0         The Lower bound of the conservative variables Q (Input)
! N5         The Upper bound of the conservative variables Q (Input)
! i_0         The starting index of the computed Hybrid Flux e.g. (0)
!               (Input)
! i_1         The ending   index of the computed Hybrid Flux e.g. (31)
!               (Input)
! dx         Grid spacing in x
!               (Input)
!              the Reconstructed function at the cell boundaries. (Output)
! F          1D array of size (N0:N5) of the function at the cell center.
!              (Input)
! D_F        1D array of size (N0:N5) containing
!              the Reconstructed function at the cell boundaries. (Output)
! MR_Flag    1D array of size (N0:N5) containing
!              Multi-resolution Flag
!            (Input, OPTIONAL)
!*********************************************************************


MODULE Scalar_Differentiation

  USE Hybrid_Scalar_Reconstruction
  USE   WENO_Scalar_Reconstruction

implicit NONE

INTERFACE PS_Diff_Scalar
  MODULE PROCEDURE Scalar_Diff_1D
  MODULE PROCEDURE Scalar_Diff_2D
  MODULE PROCEDURE Scalar_Diff_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_Diff_Scalar

CONTAINS

  Subroutine Scalar_Diff_1D (Index, Order, N0, N5, i_0, i_1, dx,   &
                                    F, D_F, MR_Flag)

  integer  :: Index, Order, N0, N5, i_0, i_1
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5)           :: F, D_F, F_Half
  integer , dimension(N0:N5), OPTIONAL :: MR_Flag

  integer  :: i

  if (PRESENT(MR_Flag)) then
        call PS_Hybrid_Scalar_Reconstruction (Order  , N0, N5, i_0-1, i_1, &
                                              F, F_Half, MR_Flag           )
  else
        call PS_WENO_Scalar_Reconstruction   (Order  , N0, N5, i_0-1, i_1, &
                                              F, F_Half                    )
  endif

!$OMP PARALLEL DO
        do i = i_0,i_1
          D_F(i) = (F_Half(i)-F_Half(i-1))/dx
        enddo
!$OMP END PARALLEL DO

  END Subroutine Scalar_Diff_1D
!
!===================================================================
!
  Subroutine Scalar_Diff_2D (Index, Order_x, N0, N5, i_0, i_1, dx,   &
                                    Order_y, M0, M5, j_0, j_1, dy,   &
                                    F, D_F, MR_Flag)

  integer  :: Index, Order_x, N0, N5, i_0, i_1
  integer  ::        Order_y, M0, M5, j_0, j_1
  REALTYPE :: dx, dy

  REALTYPE, dimension(N0:N5,M0:M5)           :: F, D_F, F_Half
  integer , dimension(N0:N5,M0:M5), OPTIONAL :: MR_Flag

  integer  :: i, j

  SELECT CASE (Index)
    CASE (1)
      if (PRESENT(MR_Flag)) then
!$OMP PARALLEL DO
        do j = j_0,j_1
          call PS_Hybrid_Scalar_Reconstruction (Order_x, N0, N5, i_0-1, i_1,   &
                                                F(:,j), F_Half(:,j),           &
                                                       MR_Flag(:,j)            )
        enddo
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO
        do j = j_0,j_1
          call PS_WENO_Scalar_Reconstruction   (Order_x, N0, N5, i_0-1, i_1,   &
                                                F(:,j), F_Half(:,j)            )
        enddo
!$OMP END PARALLEL DO
      endif

!$OMP PARALLEL DO
      do j = j_0,j_1
        do i = i_0,i_1
          D_F(i,j) = (F_Half(i,j)-F_Half(i-1,j))/dx
        enddo
      enddo
!$OMP END PARALLEL DO
  
    CASE (2)
      if (PRESENT(MR_Flag)) then
!$OMP PARALLEL DO
        do i = i_0,i_1
          call PS_Hybrid_Scalar_Reconstruction (Order_y, M0, M5, j_0-1, j_1,   &
                                                F(i,:), F_Half(i,:),           &
                                                       MR_Flag(i,:)            )
        enddo
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO
        do i = i_0,i_1
          call PS_WENO_Scalar_Reconstruction   (Order_y, M0, M5, j_0-1, j_1,   &
                                                F(i,:), F_Half(i,:)            )
        enddo
!$OMP END PARALLEL DO
      endif

!$OMP PARALLEL DO
      do j = j_0,j_1
        do i = i_0,i_1
          D_F(i,j) = (F_Half(i,j)-F_Half(i,j-1))/dy
        enddo
      enddo
!$OMP END PARALLEL DO
 
  END SELECT

  END Subroutine Scalar_Diff_2D
!
!===================================================================
!
  Subroutine Scalar_Diff_3D (Index, Order_x, N0, N5, i_0, i_1, dx,   &
                                    Order_y, M0, M5, j_0, j_1, dy,   &
                                    Order_z, K0, K5, k_0, k_1, dz,   &
                                    F, D_F, MR_Flag)

  integer  :: Index, Order_x, N0, N5, i_0, i_1
  integer  ::        Order_y, M0, M5, j_0, j_1
  integer  ::        Order_z, K0, K5, k_0, k_1
  REALTYPE :: dx, dy, dz

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)           :: F, D_F, F_Half
  integer , dimension(N0:N5,M0:M5,K0:K5), OPTIONAL :: MR_Flag

  integer  :: i, j, k

  SELECT CASE (Index)
    CASE (1)
      if (PRESENT(MR_Flag)) then
!$OMP PARALLEL DO
        do k = k_0,k_1
          do j = j_0,j_1
            call PS_Hybrid_Scalar_Reconstruction (Order_x, N0, N5, i_0-1, i_1, &
                                                  F(:,j,k), F_Half(:,j,k),     &
                                                           MR_Flag(:,j,k)      )
          enddo
        enddo
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO
        do k = k_0,k_1
          do j = j_0,j_1
            call PS_Hybrid_Scalar_Reconstruction (Order_x, N0, N5, i_0-1, i_1, &
                                                  F(:,j,k), F_Half(:,j,k)      )
          enddo
        enddo
!$OMP END PARALLEL DO
      endif

!$OMP PARALLEL DO
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            D_F(i,j,k) = (F_Half(i,j,k)-F_Half(i-1,j,k))/dx
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    CASE (2)
      if (PRESENT(MR_Flag)) then
!$OMP PARALLEL DO
        do k = k_0,k_1
          do i = i_0,i_1
            call PS_Hybrid_Scalar_Reconstruction (Order_y, M0, M5, j_0-1, j_1, &
                                                  F(i,:,k), F_Half(i,:,k),     &
                                                           MR_Flag(i,:,k)      )
          enddo
        enddo
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO
        do k = k_0,k_1
          do i = i_0,i_1
            call PS_Hybrid_Scalar_Reconstruction (Order_y, M0, M5, j_0-1, j_1, &
                                                  F(i,:,k), F_Half(i,:,k)      )
          enddo
        enddo
!$OMP END PARALLEL DO
      endif

!$OMP PARALLEL DO
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            D_F(i,j,k) = (F_Half(i,j,k)-F_Half(i,j-1,k))/dy
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    CASE (3)
      if (PRESENT(MR_Flag)) then
!$OMP PARALLEL DO
        do j = j_0,j_1
          do i = i_0,i_1
            call PS_Hybrid_Scalar_Reconstruction (Order_z, K0, K5, k_0-1, k_1, &
                                                  F(i,j,:), F_Half(i,j,:),     &
                                                           MR_Flag(i,j,:)      )
          enddo
        enddo
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO
        do j = j_0,j_1
          do i = i_0,i_1
            call PS_Hybrid_Scalar_Reconstruction (Order_z, K0, K5, k_0-1, k_1, &
                                                  F(i,j,:), F_Half(i,j,:)      )
          enddo
        enddo
!$OMP END PARALLEL DO
      endif

!$OMP PARALLEL DO
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            D_F(i,j,k) = (F_Half(i,j,k)-F_Half(i,j,k-1))/dz
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

  END SELECT

  END Subroutine Scalar_Diff_3D

END MODULE Scalar_Differentiation
