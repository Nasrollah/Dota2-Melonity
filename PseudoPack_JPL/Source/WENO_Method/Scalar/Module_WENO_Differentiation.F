#include "pseudopack.h"

! ----------------------------------------------------------------------
! FUNCTION NAME: PS_WENO_Scalar_Reconstruction
! AUTHOR       : Wai-Sun Don
! DESCRIPTION  : WENO Reconstruction
!               
! Dimension  The dimension of the PDE (2, 3, 4) (Input)
! Order      Order of the WENO scheme (3, 5)    (Input)
! N0         The Lower bound of the conservative variables Q (Input)
! N5         The Upper bound of the conservative variables Q (Input)
! i_0         The starting index of the computed WENO Flux e.g. (-1) 
!               (Input)
! i_1         The ending   index of the computed WENO Flux e.g. (N3-N2)
!               (Input)
! F          1D array of size (N0:N5) of the function
!              (Input)
! F_Half     1D array of size (N0:N5) containing
!              the Reconstructed function at the cell boundaries. (Output)
!*********************************************************************


MODULE WENO_Differentiation

  USE WENO_Scalar_Reconstruction

implicit NONE

INTERFACE PS_Diff_WENO        
  MODULE PROCEDURE WENO_Diff_1D
  MODULE PROCEDURE WENO_Diff_2D
  MODULE PROCEDURE WENO_Diff_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_Diff_WENO

CONTAINS

  Subroutine WENO_Diff_1D (Index, Order, N0, N5, i_0, i_1, dx, F, D_F)

  integer  :: Index, Order, N0, N5, i_0, i_1
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5)  :: F, D_F, F_Half

  integer  :: i

  call PS_WENO_Scalar_Reconstruction (Order, N0, N5, i_0-1, i_1, F, F_Half)

!$OMP PARALLEL DO
  do i = i_0,i_1
    D_F(i) = (F_Half(i)-F_Half(i-1))/dx
  enddo
!$OMP END PARALLEL DO

  END Subroutine WENO_Diff_1D
!
!===================================================================
!
  Subroutine WENO_Diff_2D (Index, Order_x, N0, N5, i_0, i_1, dx,   &
                                  Order_y, M0, M5, j_0, j_1, dy,   &
                                  F, D_F)

  integer  :: Index, Order_x, N0, N5, i_0, i_1
  integer  ::        Order_y, M0, M5, j_0, j_1
  REALTYPE :: dx, dy

  REALTYPE, dimension(N0:N5,M0:M5)  :: F, D_F, F_Half

  integer  :: i, j

  SELECT CASE (Index)
    CASE (1)
!$OMP PARALLEL DO
      do j = j_0,j_1
        call PS_WENO_Scalar_Reconstruction (Order_x, N0, N5, i_0-1, i_1, &
                                            F(:,j), F_Half(:,j)          )
      enddo
!$OMP END PARALLEL DO
       
!$OMP PARALLEL DO
      do j = j_0,j_1
        do i = i_0,i_1
          D_F(i,j) = (F_Half(i,j)-F_Half(i-1,j))/dx
        enddo
      enddo
!$OMP END PARALLEL DO
  
    CASE (2)
!$OMP PARALLEL DO
      do i = i_0,i_1
        call PS_WENO_Scalar_Reconstruction (Order_y, M0, M5, j_0-1, j_1, &
                                            F(i,:), F_Half(i,:)          )
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
      do j = j_0,j_1
        do i = i_0,i_1
          D_F(i,j) = (F_Half(i,j)-F_Half(i,j-1))/dy
        enddo
      enddo
!$OMP END PARALLEL DO
 
  END SELECT

  END Subroutine WENO_Diff_2D
!
!===================================================================
!
  Subroutine WENO_Diff_3D (Index, Order_x, N0, N5, i_0, i_1, dx,   &
                                  Order_y, M0, M5, j_0, j_1, dy,   &
                                  Order_z, K0, K5, k_0, k_1, dz,   &
                                  F, D_F)

  integer  :: Index, Order_x, N0, N5, i_0, i_1
  integer  ::        Order_y, M0, M5, j_0, j_1
  integer  ::        Order_z, K0, K5, k_0, k_1
  REALTYPE :: dx, dy, dz

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)  :: F, D_F, F_Half

  integer  :: i, j, k

  SELECT CASE (Index)
    CASE (1)
!$OMP PARALLEL DO
      do k = k_0,k_1
        do j = j_0,j_1
          call PS_WENO_Scalar_Reconstruction (Order_x, N0, N5, i_0-1, i_1, &
                                              F(:,j,k), F_Half(:,j,k)      )
        enddo
      enddo
!$OMP END PARALLEL DO

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
!$OMP PARALLEL DO
      do k = k_0,k_1
        do i = i_0,i_1
          call PS_WENO_Scalar_Reconstruction (Order_y, M0, M5, j_0-1, j_1, &
                                              F(i,:,k), F_Half(i,:,k)      )
        enddo
      enddo
!$OMP END PARALLEL DO

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
!$OMP PARALLEL DO
      do j = j_0,j_1
        do i = i_0,i_1
          call PS_WENO_Scalar_Reconstruction (Order_z, K0, K5, k_0-1, k_1, &
                                              F(i,j,:), F_Half(i,j,:)      )
        enddo
      enddo
!$OMP END PARALLEL DO

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

  END Subroutine WENO_Diff_3D

END MODULE WENO_Differentiation
