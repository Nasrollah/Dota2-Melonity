#include "pseudopack.h"

! ----------------------------------------------------------------------
! FUNCTION NAME: PS_CD_Sum
! AUTHOR       : Wai-Sun Don
! DESCRIPTION  : Central Differencing
!               
! Order      Order of the Central scheme                     (Input)
! N0         The Lower bound of the conservative variables Q (Input)
! N5         The Upper bound of the conservative variables Q (Input)
! i_0         The starting index of the computed Central Flux e.g. (-1) 
!               (Input)
! i_1         The ending   index of the computed Central Flux e.g. (N3-N2)
!               (Input)
! W          1D array of size (-20:20,-1:4) of the weight function 
!              (Input)
! F          1D array of size (N0:N5) of the function
!              (Input)
! D_F        1D array of size (N0:N5) of the derivative of F
!              (Input)
! alpha      post-multipication of D_F 
!*********************************************************************


MODULE Central_Difference_Function

! USE               Full_Dot_Product, PS_CD_Dot_Product=>PS_Dot_Product
  USE Central_Difference_Dot_Product

IMPLICIT NONE

INTERFACE PS_CD_Diff
  MODULE PROCEDURE CD_Sum_10
  MODULE PROCEDURE CD_Sum_20
  MODULE PROCEDURE CD_Sum_30

  MODULE PROCEDURE CD_Sum_11
  MODULE PROCEDURE CD_Sum_21
  MODULE PROCEDURE CD_Sum_31
END INTERFACE

PRIVATE
PUBLIC  :: PS_CD_Diff

CONTAINS

  Subroutine CD_Sum_10 (Index,                                       &
                        N0,N5, i_0,i_1, Order_x, W_x,                &
                        F, D_F, alpha)

  integer  :: Index, N0,N5, i_0,i_1, Order_x

  REALTYPE, dimension(-20:20,-1:4)                  :: W_x
  REALTYPE, dimension(N0:N5)                        :: F, D_F

  REALTYPE, OPTIONAL :: alpha

  integer  :: M, Order

      M = W_x(0,-1) ; Order = Order_x

  call PS_CD_Dot_Product (N0,N5, i_0,i_1, M, W_x(-M:M,Order), F, D_F, alpha)

  END Subroutine CD_Sum_10
!
!===================================================================
!
  Subroutine CD_Sum_20 (Index,                                       &
                        N0,N5, i_0,i_1, Order_x, W_x,                &
                        M0,M5, j_0,j_1, Order_y, W_y,                &
                        F, D_F, alpha)

  integer  :: Index, N0,N5, i_0,i_1, Order_x
  integer  ::        M0,M5, j_0,j_1, Order_y

  REALTYPE, dimension(-20:20,-1:4)                  :: W_x, W_y
  REALTYPE, dimension(N0:N5,M0:M5)                  :: F, D_F

  REALTYPE, OPTIONAL :: alpha

  integer  :: M, Order, i, j

  SELECT CASE (Index)
    CASE DEFAULT
      M = W_x(0,-1) ; Order = Order_x

!$OMP PARALLEL DO
      do j = j_0,j_1
        call PS_CD_Dot_Product (N0,N5, i_0,i_1, M, W_x(-M:M,Order),    &
                                F(:,j), D_F(:,j), alpha)
      enddo
!$OMP END PARALLEL DO
  
    CASE (2)
      M = W_y(0,-1) ; Order = Order_y

!$OMP PARALLEL DO
      do i = i_0,i_1
        call PS_CD_Dot_Product (M0,M5, j_0,j_1, M, W_y(-M:M,Order),    &
                                F(i,:), D_F(i,:), alpha)
      enddo
!$OMP END PARALLEL DO
 
  END SELECT

  END Subroutine CD_Sum_20
!
!===================================================================
!
  Subroutine CD_Sum_30 (Index,                                       &
                        N0,N5, i_0,i_1, Order_x, W_x,                &
                        M0,M5, j_0,j_1, Order_y, W_y,                &
                        K0,K5, k_0,k_1, Order_z, W_z,                &
                        F, D_F, alpha)

  integer  :: Index, N0,N5, i_0,i_1, Order_x
  integer  ::        M0,M5, j_0,j_1, Order_y
  integer  ::        K0,K5, k_0,k_1, Order_z

  REALTYPE, dimension(-20:20,-1:4)                  :: W_x, W_y, W_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)            :: F, D_F

  REALTYPE, OPTIONAL :: alpha

  integer  :: M, Order, i, j, k

  SELECT CASE (Index)
    CASE (1)
      M = W_x(0,-1) ; Order = Order_x

!$OMP PARALLEL DO
      do k = k_0,k_1
        do j = j_0,j_1
          call PS_CD_Dot_Product (N0,N5, i_0,i_1, M, W_x(-M:M,Order),   &
                                  F(:,j,k), D_F(:,j,k), alpha)
        enddo
      enddo
!$OMP END PARALLEL DO
  
    CASE (2)
      M = W_y(0,-1) ; Order = Order_y

!$OMP PARALLEL DO
      do k = k_0,k_1
        do i = i_0,i_1
          call PS_CD_Dot_Product (M0,M5, j_0,j_1, M, W_y(-M:M,Order),   &
                                  F(i,:,k), D_F(i,:,k), alpha)
        enddo
      enddo
!$OMP END PARALLEL DO
  
    CASE (3)
      M = W_z(0,-1) ; Order = Order_z

!$OMP PARALLEL DO
      do j = j_0,j_1
        do i = i_0,i_1
          call PS_CD_Dot_Product (K0,K5, k_0,k_1, M, W_z(-M:M,Order),   &
                                  F(i,j,:), D_F(i,j,:), alpha)
        enddo
      enddo
!$OMP END PARALLEL DO
  
  END SELECT

  END Subroutine CD_Sum_30
!
!===================================================================
!
  Subroutine CD_Sum_11 (Index,                                       &
                        N0,N5, i_0,i_1, Order_x, W_x,                &
                        NV, F, D_F, alpha)

  integer  :: Index, N0,N5, i_0,i_1, Order_x, NV

  REALTYPE, dimension(-20:20,-1:4)                  :: W_x
  REALTYPE, dimension(N0:N5,NV)                     :: F, D_F

  REALTYPE, OPTIONAL :: alpha

  integer  :: M, Order

      M = W_x(0,-1) ; Order = Order_x

  call PS_CD_Dot_Product (N0,N5, i_0,i_1, NV, M, W_x(-M:M,Order), F, D_F, alpha)

  END Subroutine CD_Sum_11
!
!===================================================================
!
  Subroutine CD_Sum_21 (Index,                                       &
                        N0,N5, i_0,i_1, Order_x, W_x,                &
                        M0,M5, j_0,j_1, Order_y, W_y,                &
                        NV, F, D_F, alpha)

  integer  :: Index, N0,N5, i_0,i_1, Order_x
  integer  ::        M0,M5, j_0,j_1, Order_y, NV

  REALTYPE, dimension(-20:20,-1:4)                  :: W_x, W_y
  REALTYPE, dimension(N0:N5,M0:M5,NV)               :: F, D_F

  REALTYPE, OPTIONAL :: alpha

  integer  :: M, Order, i, j

  SELECT CASE (Index)
    CASE DEFAULT
      M = W_x(0,-1) ; Order = Order_x

!$OMP PARALLEL DO
      do j = j_0,j_1
        call PS_CD_Dot_Product (N0,N5, i_0,i_1, NV, M, W_x(-M:M,Order),     &
                                F(:,j,:), D_F(:,j,:), alpha)
      enddo
!$OMP END PARALLEL DO
  
    CASE (2)
      M = W_y(0,-1) ; Order = Order_y

!$OMP PARALLEL DO
      do i = i_0,i_1
        call PS_CD_Dot_Product (M0,M5, j_0,j_1, NV, M, W_y(-M:M,Order),     &
                                F(i,:,:), D_F(i,:,:), alpha)
      enddo
!$OMP END PARALLEL DO
 
  END SELECT

  END Subroutine CD_Sum_21
!
!===================================================================
!
  Subroutine CD_Sum_31 (Index,                                       &
                        N0,N5, i_0,i_1, Order_x, W_x,                &
                        M0,M5, j_0,j_1, Order_y, W_y,                &
                        K0,K5, k_0,k_1, Order_z, W_z,                &
                        NV, F, D_F, alpha)

  integer  :: Index, N0,N5, i_0,i_1, Order_x
  integer  ::        M0,M5, j_0,j_1, Order_y
  integer  ::        K0,K5, k_0,k_1, Order_z, NV

  REALTYPE, dimension(-20:20,-1:4)                  :: W_x, W_y, W_z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)         :: F, D_F

  REALTYPE, OPTIONAL :: alpha

  integer  :: M, Order, i, j, k

  SELECT CASE (Index)
    CASE DEFAULT
      M = W_x(0,-1) ; Order = Order_x

!$OMP PARALLEL DO
      do k = k_0,k_1
        do j = j_0,j_1
          call PS_CD_Dot_Product (N0,N5, i_0,i_1, NV, M, W_x(-M:M,Order),     &
                                  F(:,j,k,:), D_F(:,j,k,:), alpha)
        enddo
      enddo
!$OMP END PARALLEL DO
  
    CASE (2)
      M = W_y(0,-1) ; Order = Order_y

!$OMP PARALLEL DO
      do k = k_0,k_1
        do i = i_0,i_1
          call PS_CD_Dot_Product (M0,M5, j_0,j_1, NV, M, W_y(-M:M,Order),     &
                                  F(i,:,k,:), D_F(i,:,k,:), alpha)
        enddo
      enddo
!$OMP END PARALLEL DO
  
    CASE (3)
      M = W_z(0,-1) ; Order = Order_z

!$OMP PARALLEL DO
      do j = j_0,j_1
        do i = i_0,i_1
          call PS_CD_Dot_Product (K0,K5, k_0,k_1, NV, M, W_z(-M:M,Order),     &
                                  F(i,j,:,:), D_F(i,j,:,:), alpha)
        enddo
      enddo
!$OMP END PARALLEL DO
  
  END SELECT

  END Subroutine CD_Sum_31

END MODULE Central_Difference_Function
