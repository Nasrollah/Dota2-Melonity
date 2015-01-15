#include "pseudopack.h"


MODULE Backward_Differencing

implicit NONE

INTERFACE PS_Backward_Difference
  MODULE PROCEDURE Backward_Difference_1D
END INTERFACE

PRIVATE
PUBLIC  :: PS_Backward_Difference

CONTAINS

  Subroutine Backward_Difference_1D (N0, N5, i0, i1, NV, dx, F, D_F, alpha, &
                                     Sum_D_F)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: dx
  
  REALTYPE, OPTIONAL  :: alpha
  logical , OPTIONAL  :: Sum_D_F
  logical             :: Do_Sum_D_F

  REALTYPE, dimension(N0:N5,1:NV) :: F, D_F
  REALTYPE, dimension(N0:N5,1:NV) :: D_F_0

  integer  :: i, n
  REALTYPE :: c

  if (PRESENT(alpha)) then
    if (alpha == ZERO) then ; D_F = ZERO ; RETURN ; endif
  endif

  Do_Sum_D_F = .FALSE. ; if (PRESENT(Sum_D_F)) Do_Sum_D_F = Sum_D_F

  if (Do_Sum_D_F) D_F_0 = D_F

  SELECT CASE (NV)
    CASE (1)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
      enddo
 
    CASE (2)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
      enddo
 
    CASE (3)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
      enddo
 
    CASE (4)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
        D_F(i, 4) = F(i, 4) - F(i-1, 4)
      enddo
 
    CASE (5)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
        D_F(i, 4) = F(i, 4) - F(i-1, 4)
        D_F(i, 5) = F(i, 5) - F(i-1, 5)
      enddo

    CASE (6)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
        D_F(i, 4) = F(i, 4) - F(i-1, 4)
        D_F(i, 5) = F(i, 5) - F(i-1, 5)
        D_F(i, 6) = F(i, 6) - F(i-1, 6)
     enddo

    CASE (7)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
        D_F(i, 4) = F(i, 4) - F(i-1, 4)
        D_F(i, 5) = F(i, 5) - F(i-1, 5)
        D_F(i, 6) = F(i, 6) - F(i-1, 6)
        D_F(i, 7) = F(i, 7) - F(i-1, 7)
     enddo

    CASE (8)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
        D_F(i, 4) = F(i, 4) - F(i-1, 4)
        D_F(i, 5) = F(i, 5) - F(i-1, 5)
        D_F(i, 6) = F(i, 6) - F(i-1, 6)
        D_F(i, 7) = F(i, 7) - F(i-1, 7)
        D_F(i, 8) = F(i, 8) - F(i-1, 8)
     enddo

    CASE (9)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
        D_F(i, 4) = F(i, 4) - F(i-1, 4)
        D_F(i, 5) = F(i, 5) - F(i-1, 5)
        D_F(i, 6) = F(i, 6) - F(i-1, 6)
        D_F(i, 7) = F(i, 7) - F(i-1, 7)
        D_F(i, 8) = F(i, 8) - F(i-1, 8)
        D_F(i, 9) = F(i, 9) - F(i-1, 9)
      enddo

    CASE (10)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
        D_F(i, 4) = F(i, 4) - F(i-1, 4)
        D_F(i, 5) = F(i, 5) - F(i-1, 5)
        D_F(i, 6) = F(i, 6) - F(i-1, 6)
        D_F(i, 7) = F(i, 7) - F(i-1, 7)
        D_F(i, 8) = F(i, 8) - F(i-1, 8)
        D_F(i, 9) = F(i, 9) - F(i-1, 9)
        D_F(i,10) = F(i,10) - F(i-1,10)
      enddo

    CASE (11)
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
        D_F(i, 4) = F(i, 4) - F(i-1, 4)
        D_F(i, 5) = F(i, 5) - F(i-1, 5)
        D_F(i, 6) = F(i, 6) - F(i-1, 6)
        D_F(i, 7) = F(i, 7) - F(i-1, 7)
        D_F(i, 8) = F(i, 8) - F(i-1, 8)
        D_F(i, 9) = F(i, 9) - F(i-1, 9)
        D_F(i,10) = F(i,10) - F(i-1,10)
        D_F(i,11) = F(i,11) - F(i-1,11)
      enddo

    CASE DEFAULT
      do i = i0,i1
        D_F(i, 1) = F(i, 1) - F(i-1, 1)
        D_F(i, 2) = F(i, 2) - F(i-1, 2)
        D_F(i, 3) = F(i, 3) - F(i-1, 3)
        D_F(i, 4) = F(i, 4) - F(i-1, 4)
        D_F(i, 5) = F(i, 5) - F(i-1, 5)
        D_F(i, 6) = F(i, 6) - F(i-1, 6)
        D_F(i, 7) = F(i, 7) - F(i-1, 7)
        D_F(i, 8) = F(i, 8) - F(i-1, 8)
        D_F(i, 9) = F(i, 9) - F(i-1, 9)
        D_F(i,10) = F(i,10) - F(i-1,10)
        D_F(i,11) = F(i,11) - F(i-1,11)
        D_F(i,12) = F(i,12) - F(i-1,12)
      enddo

  END SELECT

  if (NV > 12) then
    do n = 13,NV
      do i = i0,i1
        D_F(i,n) = F(i,n) - F(i-1,n)
      enddo
    enddo
  endif
 
  D_F(N0:i0-1,:) = ZERO ; D_F(i1+1:N5,:) = ZERO

  if (.NOT. PRESENT(alpha)) then
    c = ONE/dx 
  else
    c = alpha/dx
  endif

    D_F = c*D_F

  if (Do_Sum_D_F) D_F = D_F_0 + D_F

  END Subroutine Backward_Difference_1D

END MODULE Backward_Differencing
