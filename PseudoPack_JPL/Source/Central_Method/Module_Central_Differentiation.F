#include "pseudopack.h"

! ----------------------------------------------------------------------
! FUNCTION NAME: PS_Central_Diff
! AUTHOR       : Wai-Sun Don
! DESCRIPTION  : Central Differencing
!               
! Dimension  The dimension of the PDE (2, 3, 4) (Input)
! Order      Order of the Central scheme (3, 5)    (Input)
! Degree     Degree of polynomial interpolation    (Input)
! N          Total number of collocation points used (Odd)  (Input)
!              (Degree <= N-1)
! N0         The Lower bound of the conservative variables Q (Input)
! N5         The Upper bound of the conservative variables Q (Input)
! i0         The starting index of the computed Central Flux e.g. (-1) 
!               (Input)
! i1         The ending   index of the computed Central Flux e.g. (N3-N2)
!               (Input)
! F          1D array of size (N0:N5) of the function
!              (Input)
! D_F        1D array of size (N0:N5) of the derivative of F
!              (Input)
!*********************************************************************


MODULE Central_Differentiation

  USE  Central_Difference_Coefficient

IMPLICIT NONE

INTERFACE PS_Central_Diff        
  MODULE PROCEDURE Central_Diff_1D
  MODULE PROCEDURE Central_Diff_2D
  MODULE PROCEDURE Central_Diff_3D
END INTERFACE

INTERFACE PS_Central_Derivative
  MODULE PROCEDURE Central_Diff_1D
  MODULE PROCEDURE Central_Diff_2D
  MODULE PROCEDURE Central_Diff_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_Central_Diff
PUBLIC  :: PS_Central_Derivative

CONTAINS

  Subroutine Central_Diff_1D (Index,                                         &
                              Order_x, Degree_x, Nx, N0, N5, i0, i1, dx,   &
                              F, D_F)

  integer  :: Index, Order_x, Degree_x, Nx, N0, N5, i0, i1
  REALTYPE :: dx

  REALTYPE, dimension(N0:N5)  :: F, D_F

  integer                               :: i, M_s, M, N
  REALTYPE, dimension(:,:), ALLOCATABLE :: w

  M_s = (Nx+1)/2 ; M = MIN(M_s,(Degree_x+1)/2) ; N = 2*M+1
  
  ALLOCATE (w(-M:M,0:Order_x))

  call PS_CD_Coefficient (Order_x, N, Degree_x, dx, w)

!$OMP PARALLEL DO
  do i = i0,i1
    D_F(i) = SUM(w(:,Order_x)*F(-M+i:i+M))
  enddo
!$OMP END PARALLEL DO

  DEALLOCATE (w)

  END Subroutine Central_Diff_1D
!
!===================================================================
!
  Subroutine Central_Diff_2D (Index,                                         &
                              Order_x, Degree_x, Nx, N0, N5, i0, i1, dx,   &
                              Order_y, Degree_y, Ny, M0, M5, j0, j1, dy,   &
                              F, D_F)

  integer  :: Index, Order_x, Degree_x, Nx, N0, N5, i0, i1
  integer  ::        Order_y, Degree_y, Ny, M0, M5, j0, j1
  REALTYPE :: dx, dy

  REALTYPE, dimension(N0:N5,M0:M5)  :: F, D_F

  integer                               :: i, j, M_s, M, N
  REALTYPE, dimension(:,:), ALLOCATABLE :: w

  SELECT CASE (Index)
    CASE (1)
      M_s = (Nx+1)/2 ; M = MIN(M_s,(Degree_x+1)/2) ; N = 2*M+1
  
      ALLOCATE (w(-M:M,0:Order_x))

      call PS_CD_Coefficient (Order_x, N, Degree_x, dx, w)

!$OMP PARALLEL DO
      do j = j0,j1
        do i = i0,i1
          D_F(i,j) = SUM(w(:,Order_x)*F(-M+i:i+M,j))
        enddo
      enddo
!$OMP END PARALLEL DO
  
    CASE (2)
      M_s = (Ny+1)/2 ; M = MIN(M_s,(Degree_y+1)/2) ; N = 2*M+1
  
      ALLOCATE (w(-M:M,0:Order_y))

      call PS_CD_Coefficient (Order_y, N, Degree_y, dy, w)

!$OMP PARALLEL DO
      do j = j0,j1
        do i = i0,i1
          D_F(i,j) = SUM(w(:,Order_y)*F(i,-M+j:j+M))
        enddo
      enddo
!$OMP END PARALLEL DO
 
  END SELECT

  DEALLOCATE (w)

  END Subroutine Central_Diff_2D
!
!===================================================================
!
  Subroutine Central_Diff_3D (Index,                                         &
                              Order_x, Degree_x, Nx, N0, N5, i0, i1, dx,   &
                              Order_y, Degree_y, Ny, M0, M5, j0, j1, dy,   &
                              Order_z, Degree_z, Nz, K0, K5, l0, l1, dz,   &
                              F, D_F)

  integer  :: Index, Order_x, Degree_x, Nx, N0, N5, i0, i1
  integer  ::        Order_y, Degree_y, Ny, M0, M5, j0, j1
  integer  ::        Order_z, Degree_z, Nz, K0, K5, l0, l1
  REALTYPE :: dx, dy, dz

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)  :: F, D_F

  integer                               :: i, j, k, M_s, M, N
  REALTYPE, dimension(:,:), ALLOCATABLE :: w

  SELECT CASE (Index)
    CASE (1)
      M_s = (Nx+1)/2 ; M = MIN(M_s,(Degree_x+1)/2) ; N = 2*M+1
  
      write (6,*) M_s, M, Order_x, N, Degree_x, dx

      ALLOCATE (w(-M:M,0:Order_x))

      call PS_CD_Coefficient (Order_x, N, Degree_x, dx, w)

!$OMP PARALLEL DO
      do k = l0,l1
        do j = j0,j1
          do i = i0,i1
            D_F(i,j,k) = SUM(w(:,Order_x)*F(-M+i:i+M,j,k))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
  
    CASE (2)
      M_s = (Ny+1)/2 ; M = MIN(M_s,(Degree_y+1)/2) ; N = 2*M+1
  
      ALLOCATE (w(-M:M,0:Order_y))

      call PS_CD_Coefficient (Order_y, N, Degree_y, dy, w)

!$OMP PARALLEL DO
      do k = l0,l1
        do j = j0,j1
          do i = i0,i1
            D_F(i,j,k) = SUM(w(:,Order_y)*F(i,-M+j:j+M,k))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
  
    CASE (3)
      M_s = (Nz+1)/2 ; M = MIN(M_s,(Degree_z+1)/2) ; N = 2*M+1
  
      ALLOCATE (w(-M:M,0:Order_z))

      call PS_CD_Coefficient (Order_z, N, Degree_z, dz, w)

!$OMP PARALLEL DO
      do k = l0,l1
        do j = j0,j1
          do i = i0,i1
            D_F(i,j,k) = SUM(w(:,Order_z)*F(i,j,-M+k:k+M))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
  
  END SELECT

  DEALLOCATE (w)

  END Subroutine Central_Diff_3D

END MODULE Central_Differentiation
