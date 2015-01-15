#include "pseudopack.h"


MODULE Full_Dot_Product

implicit NONE

INTERFACE PS_Dot_Product
  MODULE PROCEDURE Dot_Product_0
  MODULE PROCEDURE Dot_Product_NV
END INTERFACE

PRIVATE
PUBLIC  :: PS_Dot_Product

CONTAINS

  Subroutine Dot_Product_NV (N0,N5, i0,i1, NV, M, W, F, D_F, alpha)

  integer  :: N0,N5, i0,i1, NV, M

  REALTYPE, OPTIONAL :: alpha

  REALTYPE, dimension(N0:N5,NV) :: F, D_F
  REALTYPE, dimension(-M:M)     :: W

  integer  :: i, n

  if (PRESENT(alpha)) then
    if (alpha == ZERO) then ; D_F = ZERO ; RETURN ; endif
  endif

  SELECT CASE (M)
    CASE (1)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) =                                   W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n)
        enddo
      enddo

    CASE (2)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) =                  W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n)
        enddo
      enddo

    CASE (3)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n)
        enddo
      enddo

    CASE (4)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) =                                   W(-4)*F(i-4,n) + &
                     W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n) + &
                     W( 4)*F(i+4,n) 
        enddo
      enddo

    CASE (5)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) =                  W(-5)*F(i-5,n) + W(-4)*F(i-4,n) + &
                     W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n) + &
                     W( 4)*F(i+4,n) + W( 5)*F(i+5,n)
        enddo
      enddo

    CASE (6)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W(-6)*F(i-6,n) + W(-5)*F(i-5,n) + W(-4)*F(i-4,n) + &
                     W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n) + &
                     W( 4)*F(i+4,n) + W( 5)*F(i+5,n) + W( 6)*F(i+6,n)
        enddo
      enddo

    CASE (7)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) =                                   W(-7)*F(i-7,n) + &
                     W(-6)*F(i-6,n) + W(-5)*F(i-5,n) + W(-4)*F(i-4,n) + &
                     W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n) + &
                     W( 4)*F(i+4,n) + W( 5)*F(i+5,n) + W( 6)*F(i+6,n) + &
                     W( 7)*F(i+7,n)                                         
        enddo
      enddo

    CASE (8)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) =                  W(-8)*F(i-8,n) + W(-7)*F(i-7,n) + &
                     W(-6)*F(i-6,n) + W(-5)*F(i-5,n) + W(-4)*F(i-4,n) + &
                     W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n) + &
                     W( 4)*F(i+4,n) + W( 5)*F(i+5,n) + W( 6)*F(i+6,n) + &
                     W( 7)*F(i+7,n) + W( 8)*F(i+8,n)                      
        enddo
      enddo

    CASE (9)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W(-9)*F(i-9,n) + W(-8)*F(i-8,n) + W(-7)*F(i-7,n) + &
                     W(-6)*F(i-6,n) + W(-5)*F(i-5,n) + W(-4)*F(i-4,n) + &
                     W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n) + &
                     W( 4)*F(i+4,n) + W( 5)*F(i+5,n) + W( 6)*F(i+6,n) + &
                     W( 7)*F(i+7,n) + W( 8)*F(i+8,n) + W( 9)*F(i+9,n) 
        enddo
      enddo

    CASE (10)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W(-10)*F(i-10,n) +                                     &
                     W(-9)*F(i-9,n) + W(-8)*F(i-8,n) + W(-7)*F(i-7,n) + &
                     W(-6)*F(i-6,n) + W(-5)*F(i-5,n) + W(-4)*F(i-4,n) + &
                     W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n) + &
                     W( 4)*F(i+4,n) + W( 5)*F(i+5,n) + W( 6)*F(i+6,n) + &
                     W( 7)*F(i+7,n) + W( 8)*F(i+8,n) + W( 9)*F(i+9,n) + &
                     W( 10)*F(i+10,n) 
        enddo
      enddo

    CASE (11)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W(-11)*F(i-11,n) +                                     &
                     W(-10)*F(i-10,n) +                                     &
                     W(-9)*F(i-9,n) + W(-8)*F(i-8,n) + W(-7)*F(i-7,n) + &
                     W(-6)*F(i-6,n) + W(-5)*F(i-5,n) + W(-4)*F(i-4,n) + &
                     W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n) + &
                     W( 4)*F(i+4,n) + W( 5)*F(i+5,n) + W( 6)*F(i+6,n) + &
                     W( 7)*F(i+7,n) + W( 8)*F(i+8,n) + W( 9)*F(i+9,n) + &
                     W( 10)*F(i+10,n) +                                     &
                     W( 11)*F(i+11,n) 
        enddo
      enddo

    CASE DEFAULT
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W(-12)*F(i-12,n) +                                     &
                     W(-11)*F(i-11,n) +                                     &
                     W(-10)*F(i-10,n) +                                     &
                     W(-9)*F(i-9,n) + W(-8)*F(i-8,n) + W(-7)*F(i-7,n) + &
                     W(-6)*F(i-6,n) + W(-5)*F(i-5,n) + W(-4)*F(i-4,n) + &
                     W(-3)*F(i-3,n) + W(-2)*F(i-2,n) + W(-1)*F(i-1,n) + &
                     W( 0)*F(i  ,n) +                                       &
                     W( 1)*F(i+1,n) + W( 2)*F(i+2,n) + W( 3)*F(i+3,n) + &
                     W( 4)*F(i+4,n) + W( 5)*F(i+5,n) + W( 6)*F(i+6,n) + &
                     W( 7)*F(i+7,n) + W( 8)*F(i+8,n) + W( 9)*F(i+9,n) + &
                     W( 10)*F(i+10,n) +                                     &
                     W( 11)*F(i+11,n) +                                     &
                     W( 12)*F(i+12,n) 
        enddo
      enddo

  END SELECT

  if (M >  12) then
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = D_F(i,n) + SUM(W(-M:-13)*F(i-M :i-13,n)) +   &
                                SUM(W(13:M  )*F(i+13:i+M ,n))
        enddo
      enddo
  endif

  D_F(N0:i0-1,:) = ZERO ; D_F(i1+1:N5,:) = ZERO

  if (.NOT. PRESENT(alpha)) RETURN
  if (alpha == ONE)         RETURN 

  D_F = alpha*D_F 

  END Subroutine Dot_Product_NV
!
! ==========================================================================
!
  Subroutine Dot_Product_0 (N0,N5, i0,i1, M, W, F, D_F, alpha)

  integer  :: N0,N5, i0,i1, M

  REALTYPE, OPTIONAL :: alpha

  REALTYPE, dimension(N0:N5) :: F, D_F
  REALTYPE, dimension(-M:M)  :: W

  integer  :: i

  if (PRESENT(alpha)) then
    if (alpha == ZERO) then ; D_F = ZERO ; RETURN ; endif
  endif

  SELECT CASE (M)
    CASE (1)
        do i = i0,i1
          D_F(i) =                                 W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1)
        enddo

    CASE (2)
        do i = i0,i1
          D_F(i) =                W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2)
        enddo

    CASE (3)
        do i = i0,i1
          D_F(i) = W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3)
        enddo

    CASE (4)
        do i = i0,i1
          D_F(i) =                               W(-4)*F(i-4) + &
                   W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3) + &
                   W( 4)*F(i+4) 
        enddo

    CASE (5)
        do i = i0,i1
          D_F(i) =                W(-5)*F(i-5) + W(-4)*F(i-4) + &
                   W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3) + &
                   W( 4)*F(i+4) + W( 5)*F(i+5)
        enddo

    CASE (6)
        do i = i0,i1
          D_F(i) = W(-6)*F(i-6) + W(-5)*F(i-5) + W(-4)*F(i-4) + &
                   W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3) + &
                   W( 4)*F(i+4) + W( 5)*F(i+5) + W( 6)*F(i+6)
        enddo

    CASE (7)
        do i = i0,i1
          D_F(i) =                               W(-7)*F(i-7) + &
                   W(-6)*F(i-6) + W(-5)*F(i-5) + W(-4)*F(i-4) + &
                   W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3) + &
                   W( 4)*F(i+4) + W( 5)*F(i+5) + W( 6)*F(i+6) + &
                   W( 7)*F(i+7)                                         
        enddo

    CASE (8)
        do i = i0,i1
          D_F(i) =                W(-8)*F(i-8) + W(-7)*F(i-7) + &
                   W(-6)*F(i-6) + W(-5)*F(i-5) + W(-4)*F(i-4) + &
                   W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3) + &
                   W( 4)*F(i+4) + W( 5)*F(i+5) + W( 6)*F(i+6) + &
                   W( 7)*F(i+7) + W( 8)*F(i+8)                      
        enddo

    CASE (9)
        do i = i0,i1
          D_F(i) = W(-9)*F(i-9) + W(-8)*F(i-8) + W(-7)*F(i-7) + &
                   W(-6)*F(i-6) + W(-5)*F(i-5) + W(-4)*F(i-4) + &
                   W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3) + &
                   W( 4)*F(i+4) + W( 5)*F(i+5) + W( 6)*F(i+6) + &
                   W( 7)*F(i+7) + W( 8)*F(i+8) + W( 9)*F(i+9) 
        enddo

    CASE (10)
        do i = i0,i1
          D_F(i) = W(-10)*F(i-10) +                                     &
                   W(-9)*F(i-9) + W(-8)*F(i-8) + W(-7)*F(i-7) + &
                   W(-6)*F(i-6) + W(-5)*F(i-5) + W(-4)*F(i-4) + &
                   W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3) + &
                   W( 4)*F(i+4) + W( 5)*F(i+5) + W( 6)*F(i+6) + &
                   W( 7)*F(i+7) + W( 8)*F(i+8) + W( 9)*F(i+9) + &
                   W( 10)*F(i+10) 
        enddo

    CASE (11)
        do i = i0,i1
          D_F(i) = W(-11)*F(i-11) +                                     &
                   W(-10)*F(i-10) +                                     &
                   W(-9)*F(i-9) + W(-8)*F(i-8) + W(-7)*F(i-7) + &
                   W(-6)*F(i-6) + W(-5)*F(i-5) + W(-4)*F(i-4) + &
                   W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3) + &
                   W( 4)*F(i+4) + W( 5)*F(i+5) + W( 6)*F(i+6) + &
                   W( 7)*F(i+7) + W( 8)*F(i+8) + W( 9)*F(i+9) + &
                   W( 10)*F(i+10) +                                     &
                   W( 11)*F(i+11) 
        enddo

    CASE DEFAULT
        do i = i0,i1
          D_F(i) = W(-12)*F(i-12) +                                     &
                   W(-11)*F(i-11) +                                     &
                   W(-10)*F(i-10) +                                     &
                   W(-9)*F(i-9) + W(-8)*F(i-8) + W(-7)*F(i-7) + &
                   W(-6)*F(i-6) + W(-5)*F(i-5) + W(-4)*F(i-4) + &
                   W(-3)*F(i-3) + W(-2)*F(i-2) + W(-1)*F(i-1) + &
                   W( 0)*F(i  ) +                                       &
                   W( 1)*F(i+1) + W( 2)*F(i+2) + W( 3)*F(i+3) + &
                   W( 4)*F(i+4) + W( 5)*F(i+5) + W( 6)*F(i+6) + &
                   W( 7)*F(i+7) + W( 8)*F(i+8) + W( 9)*F(i+9) + &
                   W( 10)*F(i+10) +                                     &
                   W( 11)*F(i+11) +                                     &
                   W( 12)*F(i+12) 
        enddo

  END SELECT

  if (M >  12) then
        do i = i0,i1
          D_F(i) = D_F(i) + SUM(W(-M:-13)*F(i-M :i-13)) +   &
                            SUM(W(13:M  )*F(i+13:i+M ))
        enddo
  endif

  D_F(N0:i0-1) = ZERO ; D_F(i1+1:N5) = ZERO

  if (.NOT. PRESENT(alpha)) RETURN
  if (alpha == ONE)         RETURN 

  D_F = alpha*D_F 

  END Subroutine Dot_Product_0

END MODULE Full_Dot_Product
!
! =======================================================================
!
MODULE Central_Difference_Dot_Product

implicit NONE

INTERFACE PS_CD_Dot_Product
  MODULE PROCEDURE Dot_Product_0
  MODULE PROCEDURE Dot_Product_NV
END INTERFACE

PRIVATE
PUBLIC  :: PS_CD_Dot_Product

CONTAINS

  Subroutine Dot_Product_NV (N0,N5, i0,i1, NV, M, W, F, D_F, alpha)

  integer  :: N0,N5, i0,i1, NV, M

  REALTYPE, OPTIONAL :: alpha

  REALTYPE, dimension(N0:N5,NV) :: F, D_F
  REALTYPE, dimension(-M:M)     :: W

  integer  :: i, n

  if (PRESENT(alpha)) then
    if (alpha == ZERO) then ; D_F = ZERO ; RETURN ; endif
  endif

  SELECT CASE (M)
    CASE (1)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) 
        enddo
      enddo

    CASE (2)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) 
        enddo
      enddo

    CASE (3)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) 
        enddo
      enddo

    CASE (4)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) +                       &
                     W( 4)*(F(i+ 4,n) - F(i- 4,n)) 
        enddo
      enddo

    CASE (5)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) +                       &
                     W( 4)*(F(i+ 4,n) - F(i- 4,n)) +                       &
                     W( 5)*(F(i+ 5,n) - F(i- 5,n)) 
        enddo
      enddo

    CASE (6)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) +                       &
                     W( 4)*(F(i+ 4,n) - F(i- 4,n)) +                       &
                     W( 5)*(F(i+ 5,n) - F(i- 5,n)) +                       &
                     W( 6)*(F(i+ 6,n) - F(i- 6,n)) 
        enddo
      enddo

    CASE (7)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) +                       &
                     W( 4)*(F(i+ 4,n) - F(i- 4,n)) +                       &
                     W( 5)*(F(i+ 5,n) - F(i- 5,n)) +                       &
                     W( 6)*(F(i+ 6,n) - F(i- 6,n)) +                       &
                     W( 7)*(F(i+ 7,n) - F(i- 7,n)) 
        enddo
      enddo

    CASE (8)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) +                       &
                     W( 4)*(F(i+ 4,n) - F(i- 4,n)) +                       &
                     W( 5)*(F(i+ 5,n) - F(i- 5,n)) +                       &
                     W( 6)*(F(i+ 6,n) - F(i- 6,n)) +                       &
                     W( 7)*(F(i+ 7,n) - F(i- 7,n)) +                       &
                     W( 8)*(F(i+ 8,n) - F(i- 8,n)) 
        enddo
      enddo

    CASE (9)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) +                       &
                     W( 4)*(F(i+ 4,n) - F(i- 4,n)) +                       &
                     W( 5)*(F(i+ 5,n) - F(i- 5,n)) +                       &
                     W( 6)*(F(i+ 6,n) - F(i- 6,n)) +                       &
                     W( 7)*(F(i+ 7,n) - F(i- 7,n)) +                       &
                     W( 8)*(F(i+ 8,n) - F(i- 8,n)) +                       &
                     W( 9)*(F(i+ 9,n) - F(i- 9,n)) 
        enddo
      enddo

    CASE (10)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) +                       &
                     W( 4)*(F(i+ 4,n) - F(i- 4,n)) +                       &
                     W( 5)*(F(i+ 5,n) - F(i- 5,n)) +                       &
                     W( 6)*(F(i+ 6,n) - F(i- 6,n)) +                       &
                     W( 7)*(F(i+ 7,n) - F(i- 7,n)) +                       &
                     W( 8)*(F(i+ 8,n) - F(i- 8,n)) +                       &
                     W( 9)*(F(i+ 9,n) - F(i- 9,n)) +                       &
                     W(10)*(F(i+10,n) - F(i-10,n)) 
        enddo
      enddo

    CASE (11)
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) +                       &
                     W( 4)*(F(i+ 4,n) - F(i- 4,n)) +                       &
                     W( 5)*(F(i+ 5,n) - F(i- 5,n)) +                       &
                     W( 6)*(F(i+ 6,n) - F(i- 6,n)) +                       &
                     W( 7)*(F(i+ 7,n) - F(i- 7,n)) +                       &
                     W( 8)*(F(i+ 8,n) - F(i- 8,n)) +                       &
                     W( 9)*(F(i+ 9,n) - F(i- 9,n)) +                       &
                     W(10)*(F(i+10,n) - F(i-10,n)) +                       &
                     W(11)*(F(i+11,n) - F(i-11,n)) 
        enddo
      enddo

    CASE DEFAULT
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = W( 0)* F(i   ,n)              +                       &
                     W( 1)*(F(i+ 1,n) - F(i- 1,n)) +                       &
                     W( 2)*(F(i+ 2,n) - F(i- 2,n)) +                       &
                     W( 3)*(F(i+ 3,n) - F(i- 3,n)) +                       &
                     W( 4)*(F(i+ 4,n) - F(i- 4,n)) +                       &
                     W( 5)*(F(i+ 5,n) - F(i- 5,n)) +                       &
                     W( 6)*(F(i+ 6,n) - F(i- 6,n)) +                       &
                     W( 7)*(F(i+ 7,n) - F(i- 7,n)) +                       &
                     W( 8)*(F(i+ 8,n) - F(i- 8,n)) +                       &
                     W( 9)*(F(i+ 9,n) - F(i- 9,n)) +                       &
                     W(10)*(F(i+10,n) - F(i-10,n)) +                       &
                     W(11)*(F(i+11,n) - F(i-11,n)) +                       &
                     W(12)*(F(i+12,n) - F(i-12,n))
        enddo
      enddo

  END SELECT

  if (M >  12) then
      do n = 1,NV
        do i = i0,i1
          D_F(i,n) = D_F(i,n) + SUM(W(13:M)*(F(i+13:i+M:1,n)-F(i-13:i-M:-1,n)))
        enddo
      enddo
  endif

  D_F(N0:i0-1,:) = ZERO ; D_F(i1+1:N5,:) = ZERO

  if (.NOT. PRESENT(alpha)) RETURN
  if (alpha == ONE)         RETURN 

  D_F = alpha*D_F 

  END Subroutine Dot_Product_NV
!
! ==========================================================================
!
  Subroutine Dot_Product_0 (N0,N5, i0,i1, M, W, F, D_F, alpha)

  integer  :: N0,N5, i0,i1, M

  REALTYPE, OPTIONAL :: alpha

  REALTYPE, dimension(N0:N5)    :: F, D_F
  REALTYPE, dimension(-M:M)     :: W

  integer  :: i

  if (PRESENT(alpha)) then
    if (alpha == ZERO) then ; D_F = ZERO ; RETURN ; endif
  endif

  SELECT CASE (M)
    CASE (1)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) 
        enddo

    CASE (2)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) 
        enddo

    CASE (3)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) 
        enddo

    CASE (4)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) +                       &
                   W( 4)*(F(i+ 4) - F(i- 4)) 
        enddo

    CASE (5)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) +                       &
                   W( 4)*(F(i+ 4) - F(i- 4)) +                       &
                   W( 5)*(F(i+ 5) - F(i- 5)) 
        enddo

    CASE (6)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) +                       &
                   W( 4)*(F(i+ 4) - F(i- 4)) +                       &
                   W( 5)*(F(i+ 5) - F(i- 5)) +                       &
                   W( 6)*(F(i+ 6) - F(i- 6)) 
        enddo

    CASE (7)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) +                       &
                   W( 4)*(F(i+ 4) - F(i- 4)) +                       &
                   W( 5)*(F(i+ 5) - F(i- 5)) +                       &
                   W( 6)*(F(i+ 6) - F(i- 6)) +                       &
                   W( 7)*(F(i+ 7) - F(i- 7)) 
        enddo

    CASE (8)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) +                       &
                   W( 4)*(F(i+ 4) - F(i- 4)) +                       &
                   W( 5)*(F(i+ 5) - F(i- 5)) +                       &
                   W( 6)*(F(i+ 6) - F(i- 6)) +                       &
                   W( 7)*(F(i+ 7) - F(i- 7)) +                       &
                   W( 8)*(F(i+ 8) - F(i- 8)) 
        enddo

    CASE (9)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) +                       &
                   W( 4)*(F(i+ 4) - F(i- 4)) +                       &
                   W( 5)*(F(i+ 5) - F(i- 5)) +                       &
                   W( 6)*(F(i+ 6) - F(i- 6)) +                       &
                   W( 7)*(F(i+ 7) - F(i- 7)) +                       &
                   W( 8)*(F(i+ 8) - F(i- 8)) +                       &
                   W( 9)*(F(i+ 9) - F(i- 9)) 
        enddo

    CASE (10)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) +                       &
                   W( 4)*(F(i+ 4) - F(i- 4)) +                       &
                   W( 5)*(F(i+ 5) - F(i- 5)) +                       &
                   W( 6)*(F(i+ 6) - F(i- 6)) +                       &
                   W( 7)*(F(i+ 7) - F(i- 7)) +                       &
                   W( 8)*(F(i+ 8) - F(i- 8)) +                       &
                   W( 9)*(F(i+ 9) - F(i- 9)) +                       &
                   W(10)*(F(i+10) - F(i-10)) 
        enddo

    CASE (11)
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) +                       &
                   W( 4)*(F(i+ 4) - F(i- 4)) +                       &
                   W( 5)*(F(i+ 5) - F(i- 5)) +                       &
                   W( 6)*(F(i+ 6) - F(i- 6)) +                       &
                   W( 7)*(F(i+ 7) - F(i- 7)) +                       &
                   W( 8)*(F(i+ 8) - F(i- 8)) +                       &
                   W( 9)*(F(i+ 9) - F(i- 9)) +                       &
                   W(10)*(F(i+10) - F(i-10)) +                       &
                   W(11)*(F(i+11) - F(i-11)) 
        enddo

    CASE DEFAULT
        do i = i0,i1
          D_F(i) = W( 0)* F(i   )            +                       &
                   W( 1)*(F(i+ 1) - F(i- 1)) +                       &
                   W( 2)*(F(i+ 2) - F(i- 2)) +                       &
                   W( 3)*(F(i+ 3) - F(i- 3)) +                       &
                   W( 4)*(F(i+ 4) - F(i- 4)) +                       &
                   W( 5)*(F(i+ 5) - F(i- 5)) +                       &
                   W( 6)*(F(i+ 6) - F(i- 6)) +                       &
                   W( 7)*(F(i+ 7) - F(i- 7)) +                       &
                   W( 8)*(F(i+ 8) - F(i- 8)) +                       &
                   W( 9)*(F(i+ 9) - F(i- 9)) +                       &
                   W(10)*(F(i+10) - F(i-10)) +                       &
                   W(11)*(F(i+11) - F(i-11)) +                       &
                   W(12)*(F(i+12) - F(i-12))
        enddo

  END SELECT

  if (M >  12) then
        do i = i0,i1
          D_F(i) = D_F(i) + SUM(W(13:M)*(F(i+13:i+M:1)-F(i-13:i-M:-1)))
        enddo
  endif

  D_F(N0:i0-1) = ZERO ; D_F(i1+1:N5) = ZERO

  if (.NOT. PRESENT(alpha)) RETURN
  if (alpha == ONE)         RETURN 

  D_F = alpha*D_F 

  END Subroutine Dot_Product_0

END MODULE Central_Difference_Dot_Product
