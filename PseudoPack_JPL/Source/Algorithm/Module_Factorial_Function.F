#include "pseudopack.h"

MODULE Factorial_Function

implicit NONE

PRIVATE
PUBLIC  :: Factorial
PUBLIC  :: Factorial_HALF
PUBLIC  :: Factorial_Omega

CONTAINS

  RECURSIVE Function Factorial (n) RESULT (Factorial_n)

  REALTYPE              :: Factorial_n
  integer , INTENT (IN) :: n

  SELECT CASE (n)
    CASE (0)
      Factorial_n = ONE

    CASE (1:)
      Factorial_n = n*Factorial(n-1)

    CASE DEFAULT
      Factorial_n = ZERO

  END SELECT

  END Function Factorial
!
!=============================================================
!
  REALTYPE Function Factorial_HALF (n) 

  REALTYPE              :: F
  integer               :: i
  integer , INTENT (IN) :: n

  F = SQRT(PI)
  do i = 1,n
    F = (i-HALF)*F
  enddo

  Factorial_HALF = F

  END Function Factorial_HALF
!
!=============================================================
!
  REALTYPE Function Factorial_Omega (n, m) 

  integer , INTENT (IN) :: n, m

  if (m == 0) then ; Factorial_Omega = ONE ; RETURN ; endif

  SELECT CASE (n)
    CASE (0)
      Factorial_Omega = Factorial(m-1)

    CASE (1:)
      Factorial_Omega = Factorial(n+m-1)*ONE/Factorial(n-1)

    CASE DEFAULT
      Factorial_Omega = ZERO

  END SELECT

  END Function Factorial_Omega

END MODULE Factorial_Function
