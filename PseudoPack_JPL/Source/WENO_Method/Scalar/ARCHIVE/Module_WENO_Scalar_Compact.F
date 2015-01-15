#include "pseudopack.h"


Module WENO_Scalar_Compact

  USE WENO_Coefficients
  USE WENO_Option

implicit NONE

PRIVATE
PUBLIC  :: PS_WENO_Scalar_Compact

CONTAINS

  Subroutine PS_WENO_Scalar_Compact (Order, N0, N5, i0, i1, F, F_Half, &
                                            WENO_Epsilon)

  integer                    :: Order, N0, N5, i0, i1
  REALTYPE, dimension(N0:N5) :: F, F_Half

  REALTYPE, OPTIONAL         :: WENO_Epsilon
  REALTYPE                   :: Epsilon_WENO

  integer                    :: i, j, k, m, r

  REALTYPE, dimension(0:(Order-1)/2                            ) :: V, a, b
  REALTYPE, dimension(0:(Order-1)/2                            ) :: d
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2              ) :: C
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2,0:(Order-1)/2) :: S

                             Epsilon_WENO = PS_WENO_Option%Epsilon
  if (PRESENT(WENO_Epsilon)) Epsilon_WENO = WENO_Epsilon

  call PS_WENO_Coefficient (Order, d, C, S)

  k = (Order-1)/2

!$OMP PARALLEL DO
  do i = i0,i1

    do r = 0,k
      j = i+r 

        V(r) = SUM(C(r,:)*F(j-k:j))

        b(r) = F(j-0)*SUM(S(k-0,k-0:k,r)*F(j-0:j))
      do m = 1,k
        b(r) = F(j-m)*SUM(S(k-m,k-m:k,r)*F(j-m:j)) + b(r)
      enddo
    enddo

    a = d/(Epsilon_WENO + b)**2

#if defined (MAPPED_WENO)
      if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif      

    F_Half(i) = SUM(a*V)/SUM(a)

  enddo
!$OMP END PARALLEL DO

  END Subroutine PS_WENO_Scalar_Compact 
!
! ========================================================
!
  FUNCTION W(d, W_k)
 
  REALTYPE                     :: S_W_k
  REALTYPE, dimension(:)       :: d, W_k
  REALTYPE, dimension(SIZE(d)) :: W, Omega

  S_W_k = SUM(W_k)
  Omega = W_k/S_W_k
  W     = Omega*(d*(d+ONE)+Omega*(Omega-THREE*d))/                 &
                (d*d      +Omega*(ONE  -  TWO*d))

  END FUNCTION 

END Module WENO_Scalar_Compact
