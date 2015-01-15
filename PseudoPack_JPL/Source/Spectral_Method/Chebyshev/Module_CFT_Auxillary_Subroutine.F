#include "pseudopack.h"


MODULE CFT_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC  :: CFT_1_Copy
PUBLIC  :: CFT_1_Smooth
PUBLIC  :: CFT_1_Recursive
PUBLIC  :: CFT_2_Copy
PUBLIC  :: CFT_2_Smooth
PUBLIC  :: CFT_2_Recursive

PUBLIC  :: CFT_1_Copy_Even
PUBLIC  :: CFT_1_Copy_Odd
PUBLIC  :: CFT_1_Smooth_Even
PUBLIC  :: CFT_1_Smooth_Odd
PUBLIC  :: CFT_1_Recursive_Even
PUBLIC  :: CFT_1_Recursive_Odd

PUBLIC  :: CFT_2_Copy_Even
PUBLIC  :: CFT_2_Copy_Odd
PUBLIC  :: CFT_2_Smooth_Even
PUBLIC  :: CFT_2_Smooth_Odd
PUBLIC  :: CFT_2_Recursive_Even
PUBLIC  :: CFT_2_Recursive_Odd

CONTAINS

  Subroutine CFT_1_Copy   (N, M, a, b)

  integer  :: N, M

  REALTYPE, dimension(0:,:) :: a, b

  b(1:N,1:M) = a(1:N,1:M)

  END Subroutine CFT_1_Copy   

! -----------------------------------------------------------------------

  Subroutine CFT_1_Smooth (N, M, Sigma, a)

  integer  :: N, M, j

  REALTYPE, dimension(0:,:) :: a
  REALTYPE, dimension(0:)   :: Sigma

  do j = 1,M
    a(0:N,j) = a(0:N,j)*Sigma(0:N)
  enddo

  END Subroutine CFT_1_Smooth 

! -----------------------------------------------------------------------

  Subroutine CFT_1_Recursive (N, M, a, b)

  integer  :: N, M, i

  REALTYPE, dimension(0:,:) :: a, b

    b(N  ,1:M) = ZERO
    b(N-1,1:M) = N*a(N,1:M)

  do i = N-2,1,-2
    b(i  ,1:M) = b(i+2,1:M) + 2*(i+1)*a(i+1,1:M)
    b(i-1,1:M) = b(i+1,1:M) + 2*(i  )*a(i  ,1:M)
  enddo

  if (MOD(N,2) == 0) then
    b(0,1:M) = b(2,1:M) + TWO*a(1,1:M)
  endif

  END Subroutine CFT_1_Recursive 

! ----------------------------------------------------------------------

  Subroutine CFT_2_Copy   (N, M, a, b)

  integer  :: N, M

  REALTYPE, dimension(:,0:) :: a, b

  b(1:M,1:N) = a(1:M,1:N)

  END Subroutine CFT_2_Copy   

! -----------------------------------------------------------------------

  Subroutine CFT_2_Smooth (N, M, Sigma, a)

  integer  :: N, M, i

  REALTYPE, dimension(:,0:) :: a
  REALTYPE, dimension(0:)   :: Sigma

  do i = 0,N
    a(1:M,i) = a(1:M,i)*Sigma(i)
  enddo

  END Subroutine CFT_2_Smooth 

! -----------------------------------------------------------------------

  Subroutine CFT_2_Recursive (N, M, a, b)

  integer  :: N, M, i

  REALTYPE, dimension(:,0:) :: a, b

    b(1:M,N  ) = ZERO
    b(1:M,N-1) = N*a(1:M,N)

#if defined (ESSL) || defined (CRAY)
  do i = N-2,1,-2
    b(1:M,i  ) = b(1:M,i+2) + 2*(i+1)*a(1:M,i+1)
    b(1:M,i-1) = b(1:M,i+1) + 2*(i  )*a(1:M,i  )
  enddo

  if (MOD(N,2) == 0) then
    b(1:M,0) = b(1:M,2) + TWO*a(1:M,1)
  endif
#endif

#if defined (VFFT)
  do i = N-2,0,-1
    b(1:M,i  ) = b(1:M,i+2) + 2*(i+1)*a(1:M,i+1)
  enddo
#endif

  END Subroutine CFT_2_Recursive 

! -----------------------------------------------------------------------
  Subroutine CFT_1_Copy_Even (N, M, a, b)

  integer  :: N, M

  REALTYPE, dimension(0:,:) :: a, b

  b(0:N,1:M) = a(0:N,1:M)

  END Subroutine CFT_1_Copy_Even

! -----------------------------------------------------------------------

  Subroutine CFT_1_Copy_Odd  (N, M, a, b)

  integer  :: N, M

  REALTYPE, dimension(0:,:) :: a, b

  b(N,1:M)     = ZERO
  b(0:N-1,1:M) = a(0:N-1,1:M)

  END Subroutine CFT_1_Copy_Odd

! -----------------------------------------------------------------------

  Subroutine CFT_1_Smooth_Even (N, M, Sigma, a)

  integer  :: N, M, j

  REALTYPE, dimension(0:,:) :: a
  REALTYPE, dimension(0:)   :: Sigma

  do j = 1,M
    a(0:N,j) = a(0:N,j)*Sigma(0:2*N:2)
  enddo

  END Subroutine CFT_1_Smooth_Even 

! -----------------------------------------------------------------------

  Subroutine CFT_1_Smooth_Odd  (N, M, Sigma, a)

  integer  :: N, M, j

  REALTYPE, dimension(0:,:) :: a
  REALTYPE, dimension(0:)   :: Sigma

  do j = 1,M
    a(0:N-1,j) = a(0:N-1,j)*Sigma(1:2*(N-1)+1:2)
  enddo

  END Subroutine CFT_1_Smooth_Odd

! -----------------------------------------------------------------------

  Subroutine CFT_1_Recursive_Even (N, M, a, b)

  integer  :: N, M

  REALTYPE, dimension(0:,:) :: a, b

#if defined (ESSL) || defined (CRAY)
  integer :: i

    b(N-1,1:M) = 2*N*a(N,1:M)

  do i = N-2,0,-1
    b(i,1:M) = b(i+1,1:M) + 4*(i+1)*a(i+1,1:M)
  enddo
#endif

  END Subroutine CFT_1_Recursive_Even 

! -----------------------------------------------------------------------

  Subroutine CFT_1_Recursive_Odd  (N, M, a, b)

  integer  :: N, M, i

  REALTYPE, dimension(0:,:) :: a, b

    b(N,1:M) = ZERO

  do i = N-1,0,-1
    b(i,1:M) = b(i+1,1:M) + 2*(2*i+1)*a(i,1:M)
  enddo

  END Subroutine CFT_1_Recursive_Odd

! -----------------------------------------------------------------------
  Subroutine CFT_2_Copy_Even (N, M, a, b)

  integer  :: N, M

  REALTYPE, dimension(:,0:) :: a, b

  b(1:M,0:N) = a(1:M,0:N)

  END Subroutine CFT_2_Copy_Even

! -----------------------------------------------------------------------

  Subroutine CFT_2_Copy_Odd  (N, M, a, b)

  integer  :: N, M

  REALTYPE, dimension(:,0:) :: a, b

  b(1:M,0:N-1) = a(1:M,0:N-1)

  END Subroutine CFT_2_Copy_Odd

! -----------------------------------------------------------------------

  Subroutine CFT_2_Smooth_Even (N, M, Sigma, a)

  integer  :: N, M, i

  REALTYPE, dimension(:,0:) :: a
  REALTYPE, dimension(0:)   :: Sigma

  do i = 0,N
    a(1:M,i) = a(1:M,i)*Sigma(2*i)
  enddo

  END Subroutine CFT_2_Smooth_Even 

! -----------------------------------------------------------------------

  Subroutine CFT_2_Smooth_Odd  (N, M, Sigma, a)

  integer  :: N, M, i

  REALTYPE, dimension(:,0:) :: a
  REALTYPE, dimension(0:)   :: Sigma

  do i = 0,N-1
    a(1:M,i) = a(1:M,i)*Sigma(2*i+1)
  enddo

  END Subroutine CFT_2_Smooth_Odd

! -----------------------------------------------------------------------

  Subroutine CFT_2_Recursive_Even (N, M, a, b)

  integer  :: N, M, i

  REALTYPE, dimension(:,0:) :: a, b

    b(1:M,N-1) = 2*N*a(1:M,N)

  do i = N-2,0,-1
    b(1:M,i) = b(1:M,i+1) + 4*(i+1)*a(1:M,i+1)
  enddo

  END Subroutine CFT_2_Recursive_Even 

! -----------------------------------------------------------------------

  Subroutine CFT_2_Recursive_Odd  (N, M, a, b)

  integer  :: N, M, i

  REALTYPE, dimension(:,0:) :: a, b

    b(1:M,N) = ZERO

  do i = N-1,0,-1
    b(1:M,i) = b(1:M,i+1) + 2*(2*i+1)*a(1:M,i)
  enddo

  END Subroutine CFT_2_Recursive_odd

END MODULE CFT_Auxillary_Subroutine
