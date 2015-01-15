  call PS_WENO_Coefficient (Order, d, C, S)

  k = (Order-1)/2

  F_Half = ZERO

  do l = 1,2
!$OMP PARALLEL DO
    do i = i0,i1

#include "WENO_Euler_Compact_Beta.i"

#include "WENO_Zico_Tau.i"
#include "WENO_Euler_Zico.i"

      F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
    enddo
!$OMP END PARALLEL DO
  enddo
