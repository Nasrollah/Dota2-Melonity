  call PS_WENO_Coefficient (Order, d, C, S)

  SELECT CASE (Order)
    CASE (1)
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Long_Beta_1.i"

#include "WENO_Euler_Long_JS_1.i"

        Sum_a = a(0) 

        F_Half(i) = (a(0)*V(0))/Sum_a 
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Long_Beta_3.i"

#include "WENO_Euler_Long_JS_3.i"

        Sum_a = a(0) + a(1) 

        F_Half(i) = (a(0)*V(0) + a(1)*V(1))/Sum_a 
      enddo
!$OMP END PARALLEL DO

    CASE (5)
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Long_Beta_5.i"

#include "WENO_Euler_Long_JS_5.i"

        Sum_a = a(0) + a(1) + a(2)

        F_Half(i) = (a(0)*V(0) + a(1)*V(1) + a(2)*V(2))/Sum_a
      enddo
!$OMP END PARALLEL DO

    CASE (7) 
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Long_Beta_7.i"

#include "WENO_Euler_Long_JS_7.i"

        Sum_a = a(0) + a(1) + a(2) + a(3) 

        F_Half(i) = (a(0)*V(0) + a(1)*V(1) + a(2)*V(2) + a(3)*V(3))/Sum_a
      enddo
!$OMP END PARALLEL DO

    CASE (9) 
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Long_Beta_9.i"

#include "WENO_Euler_Long_JS_9.i"

        Sum_a = a(0) + a(1) + a(2) + a(3) + a(4) 

        F_Half(i) = (a(0)*V(0) + a(1)*V(1) + a(2)*V(2) + &
                     a(3)*V(3) + a(4)*V(4)            )/Sum_a

    enddo
!$OMP END PARALLEL DO

    CASE (11)
!$OMP PARALLEL DO
      do i = i0,i1  
#include "WENO_Scalar_Long_Beta_11.i"

#include "WENO_Euler_Long_JS_11.i"

        Sum_a = a(0) + a(1) + a(2) + a(3) + a(4) + a(5)

        F_Half(i) = (a(0)*V(0) + a(1)*V(1) + a(2)*V(2) + &
                     a(3)*V(3) + a(4)*V(4) + a(5)*V(5))/Sum_a
      enddo
!$OMP END PARALLEL DO

  END SELECT
