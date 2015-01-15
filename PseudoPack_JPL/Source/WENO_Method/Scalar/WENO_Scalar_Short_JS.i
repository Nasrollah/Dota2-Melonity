  call PS_WENO_Coefficient (Order, d, C, S)

  k = (Order-1)/2

  SELECT CASE (Order)
    CASE (1)
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Short_Beta_1.i"

#include "WENO_Euler_JS.i"

        F_Half(i) = SUM(a*V)/Sum(a)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Short_Beta_3.i"

#include "WENO_Euler_JS.i"

        F_Half(i) = SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO

    CASE (5)
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Short_Beta_5.i"

#include "WENO_Euler_JS.i"

        F_Half(i) = SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO

    CASE (7) 
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Short_Beta_7.i"

#include "WENO_Euler_JS.i"

        F_Half(i) = SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO

    CASE (9) 
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Short_Beta_9.i"

#include "WENO_Euler_JS.i"

        F_Half(i) = SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO

    CASE (11)
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Short_Beta_11.i"

#include "WENO_Euler_JS.i"

        F_Half(i) = SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO
  END SELECT
