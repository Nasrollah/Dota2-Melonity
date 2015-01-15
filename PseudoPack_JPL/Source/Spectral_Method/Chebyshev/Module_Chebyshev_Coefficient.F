#include "pseudopack.h"

! ==================================================================
!  Perform the interpolation at the point z_l
!  Along the first  dimension,   x_0 <= z_l <= x_N for all  j=0,..,M
!     g(z_l, y_j) = Sum(a_kj * T_kl)  
!
!  Along the second dimension,   y_0 <= z_l <= y_N for all i=0,..,N
!     g(x_i, z_l) = Sum(a_ki * T_kl)
!
!  where T_kl = COS(k*ACOS(z_l))
!        a_kj = the kth Chebyshev Lobatto coefficients for each j
!        a_ki = the kth Chebyshev Lobatto coefficients for each i
!
! ==================================================================
 
MODULE Chebyshev_Coefficient

implicit NONE

PRIVATE

PUBLIC  :: PS_Interpolation_Sum_a_T_1D
PUBLIC  :: PS_Interpolation_Sum_a_T_2D
PUBLIC  :: PS_Chebyshev_Lobatto_T_k_at_z 
PUBLIC  :: PS_Chebyshev_Lobatto_a_k 

CONTAINS

  Subroutine PS_Interpolation_Sum_a_T_1D (a, T, g)

  integer                    :: i
  REALTYPE, dimension(0:)    :: a, g
  REALTYPE, dimension(0:,0:) :: T

    if (SIZE(a,DIM=1) /= SIZE(T,DIM=1)) STOP 'PS_Interpolation_Sum_a_T_1D -- 1'

    if (SIZE(g,DIM=1) <= SIZE(T,DIM=2)) STOP 'PS_Interpolation_Sum_a_T_1D -- 2'

  do i = 0,SIZE(g)-1
    g(i) = SUM(a*T(:,i))
  enddo

  END Subroutine PS_Interpolation_Sum_a_T_1D
!
! ==================================================================
!
  Subroutine PS_Interpolation_Sum_a_T_2D (a, T, g, Transpose)

  integer                    :: i, j, l
  logical                    :: Transpose
  REALTYPE, dimension(0:,0:) :: a, T, g

    if (SIZE(a,DIM=1) /= SIZE(T,DIM=1)) STOP 'PS_Interpolation_Sum_a_T_2D -- 1'

  if (Transpose) then
    if (SIZE(g,DIM=1) <= SIZE(T,DIM=2)) STOP 'PS_Interpolation_Sum_a_T_2D -- 2'
    if (SIZE(g,DIM=2) <= SIZE(a,DIM=2)) STOP 'PS_Interpolation_Sum_a_T_2D -- 3'
  else
    if (SIZE(g,DIM=1) <= SIZE(a,DIM=2)) STOP 'PS_Interpolation_Sum_a_T_2D -- 2'
    if (SIZE(g,DIM=2) <= SIZE(T,DIM=2)) STOP 'PS_Interpolation_Sum_a_T_2D -- 3'
  endif

  if (Transpose) then

    do l = 0,SIZE(T,DIM=2)-1          ! Interpolate along the second dimension
      do i = 0,SIZE(a,DIM=2)-1
        g(i,l) = SUM(a(:,i)*T(:,l))
      enddo
    enddo

  else

    do j = 0,SIZE(a,DIM=2)-1          ! Interpolate along the first dimension
      do l = 0,SIZE(T,DIM=2)-1       
        g(l,j) = SUM(a(:,j)*T(:,l))
      enddo
    enddo

  endif
 
  END Subroutine PS_Interpolation_Sum_a_T_2D
!
! ==================================================================
!
  Subroutine PS_Chebyshev_Lobatto_T_k_at_z (z, T)

  REALTYPE, dimension(0:)    :: z
  REALTYPE, dimension(0:,0:) :: T

  integer  :: k, M

    if (SIZE(z,DIM=1) >= SIZE(T,DIM=2)) STOP 'PS_Chebyshev_Lobatto_T_k_at_z'

  M = SIZE(T,DIM=1)-1

  do k = 0,M
    T(k,:) = COS(k*ACOS(z))
  enddo

  END Subroutine PS_Chebyshev_Lobatto_T_k_at_z 
!
! ==================================================================
! Should pre-compute the cos(i*k*Pi/N) for multiple j
!
  Subroutine PS_Chebyshev_Lobatto_a_k (f, a, Even_Or_Odd, N_Total)

  REALTYPE, dimension(0:) :: f, a
  integer ,      OPTIONAL :: N_Total

  integer ,      OPTIONAL :: Even_Or_Odd
  integer                 :: EOD
               
  REALTYPE, dimension(0:2*SIZE(f)-1) :: b
  REALTYPE, dimension(0:  SIZE(f)-1) :: d

  integer  :: N, M, i, k
  REALTYPE :: PiN, TwoN, delta

  EOD = 0 ; if (PRESENT(Even_Or_Odd)) EOD = Even_Or_Odd

                               M = SIZE(f)-1
  if (PRESENT(N_Total)) then ; N = N_Total-1 ; else ; N = SIZE(f)-1 ; endif
  
  d(0) = TWO ; d(1:M-1) = ONE ; d(M) = TWO

  if ((MOD(N,2) /= 0) .AND. (EOD /= 0)) d(M) = ONE

  PiN = PI/N ; TwoN = TWO/N

#if 0
  SELECT CASE (EOD)
    CASE (-1)
      do k = 1,N,2

        b(k) = f(0)/d(0)
        do i = 1,M
          b(k) = b(k) + f(i)*COS(i*k*PiN)/d(i)
        enddo

      enddo

      b(1:N:2) = TWO*b(1:N:2)

    CASE ( 0)
      do k = 0,N

        b(k) = f(0)/d(0)
        do i = 1,M
          b(k) = b(k) + f(i)*COS(i*k*PiN)/d(i)
        enddo

      enddo

    CASE ( 1)
      do k = 0,N,2

        b(k) = f(0)/d(0)
        do i = 1,M
          b(k) = b(k) + f(i)*COS(i*k*PiN)/d(i)
        enddo

      enddo

      b(0:N:2) = TWO*b(0:N:2)

    END SELECT
#else
  b = ZERO

  do k = 0,N

    delta = (ONE+EOD*(-1)**k)

    b(k) = delta*f(0)/d(0)
    do i = 1,M
      b(k) = b(k) + delta*f(i)*COS(i*k*PiN)/d(i)
    enddo

  enddo
#endif

  b(0) = HALF*b(0)
  b(N) = HALF*b(N)

  if (EOD == -1) a = b(1:N:2)*TwoN
  if (EOD ==  0) a = b(0:N:1)*TwoN
  if (EOD ==  1) a = b(0:N:2)*TwoN

  END Subroutine PS_Chebyshev_Lobatto_a_k

END MODULE Chebyshev_Coefficient
