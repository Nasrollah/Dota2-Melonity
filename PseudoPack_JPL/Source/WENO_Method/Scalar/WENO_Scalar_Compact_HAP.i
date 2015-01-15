  call PS_WENO_Coefficient (Order, d, C, S)

  k = (Order-1)/2

!$OMP PARALLEL DO
  do i = i0,i1
#include "WENO_Scalar_Compact_Beta.i"

#include "WENO_Euler_JS.i"

#include "WENO_Euler_HAP.i"

    F_Half(i) = SUM(a*V)/SUM(a)

  enddo
!$OMP END PARALLEL DO
