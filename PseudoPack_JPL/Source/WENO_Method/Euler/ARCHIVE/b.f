#include "pseudopack.h"


Module WENO_Euler_Long_JS

  USE WENO_Coefficients
  USE WENO_Option

implicit NONE

PRIVATE
PUBLIC  :: PS_WENO_Euler_Long_JS

CONTAINS

  Subroutine PS_WENO_Euler_Long_JS  (Order, N0, N5, i0, i1, F, F_Half,  &
                                            WENO_Epsilon, WENO_Power)

#include "WENO_Euler_Declaration.i"
 
#if 0
#ifdef WENO_NON_INTEGER_POWER
                             Power        = PS_WENO_Option%Power
  if (PRESENT(WENO_Power  )) Power        = WENO_Power
#else
                             Power        = INT(PS_WENO_Option%Power)
  if (PRESENT(WENO_Power  )) Power        = INT(WENO_Power)
#endif
#else
                             Power        = PS_WENO_Option%Power
  if (PRESENT(WENO_Power  )) Power        = WENO_Power

  p = INT(Power)
#endif

                             Epsilon_WENO = PS_WENO_Option%Epsilon
  if (PRESENT(WENO_Epsilon)) Epsilon_WENO = WENO_Epsilon

  call PS_WENO_Coefficient (Order, d, C, S)

  F_Half = ZERO

  do l = 1,2

  SELECT CASE (Order)
    CASE (1)
      do i = i0,i1  
#include "WENO_Euler_Long_Order_1_JS.i"

#ifdef NO_ABS_BETA_JS
        a(0) = d(0)/(    b(0)  + Epsilon_WENO)**p
#else
        a(0) = d(0)/(ABS(b(0)) + Epsilon_WENO)**p
#endif

        Sum_a = a(0) 

        F_Half(i) = F_Half(i) + (a(0)*V(0))/Sum_a 
      enddo

    CASE (3)
      do i = i0,i1  
#include "WENO_Euler_Long_Order_3_JS.i"

#ifdef NO_ABS_BETA_JS
        a(0) = d(0)/(    b(0)  + Epsilon_WENO)**p
        a(1) = d(1)/(    b(1)  + Epsilon_WENO)**p
#else
        a(0) = d(0)/(ABS(b(0)) + Epsilon_WENO)**p
        a(1) = d(1)/(ABS(b(1)) + Epsilon_WENO)**p
#endif

        Sum_a = a(0) + a(1) 

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1))/Sum_a 
      enddo

    CASE (5)
      do i = i0,i1  
#include "WENO_Euler_Long_Order_5_JS.i"

#ifdef NO_ABS_BETA_JS
        a(0) = d(0)/(    b(0)  + Epsilon_WENO)**p
        a(1) = d(1)/(    b(1)  + Epsilon_WENO)**p
        a(2) = d(2)/(    b(2)  + Epsilon_WENO)**p
#else
        a(0) = d(0)/(ABS(b(0)) + Epsilon_WENO)**p
        a(1) = d(1)/(ABS(b(1)) + Epsilon_WENO)**p
        a(2) = d(2)/(ABS(b(2)) + Epsilon_WENO)**p
#endif

        Sum_a = a(0) + a(1) + a(2)

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1) + a(2)*V(2))/Sum_a
      enddo

    CASE (7) 
      do i = i0,i1  
#include "WENO_Euler_Long_Order_7_JS.i"

#ifdef NO_ABS_BETA_JS
        a(0) = d(0)/(    b(0)  + Epsilon_WENO)**p
        a(1) = d(1)/(    b(1)  + Epsilon_WENO)**p
        a(2) = d(2)/(    b(2)  + Epsilon_WENO)**p
        a(3) = d(3)/(    b(3)  + Epsilon_WENO)**p
#else
        a(0) = d(0)/(ABS(b(0)) + Epsilon_WENO)**p
        a(1) = d(1)/(ABS(b(1)) + Epsilon_WENO)**p
        a(2) = d(2)/(ABS(b(2)) + Epsilon_WENO)**p
        a(3) = d(3)/(ABS(b(3)) + Epsilon_WENO)**p
#endif

        Sum_a = a(0) + a(1) + a(2) + a(3) 

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1) + a(2)*V(2) + &
                                 a(3)*V(3)                        )/Sum_a
      enddo

    CASE (9) 
      do i = i0,i1  
#include "WENO_Euler_Long_Order_9_JS.i"

#ifdef NO_ABS_BETA_JS
        a(0) = d(0)/(    b(0)  + Epsilon_WENO)**p
        a(1) = d(1)/(    b(1)  + Epsilon_WENO)**p
        a(2) = d(2)/(    b(2)  + Epsilon_WENO)**p
        a(3) = d(3)/(    b(3)  + Epsilon_WENO)**p
        a(4) = d(4)/(    b(4)  + Epsilon_WENO)**p
#else
        a(0) = d(0)/(ABS(b(0)) + Epsilon_WENO)**p
        a(1) = d(1)/(ABS(b(1)) + Epsilon_WENO)**p
        a(2) = d(2)/(ABS(b(2)) + Epsilon_WENO)**p
        a(3) = d(3)/(ABS(b(3)) + Epsilon_WENO)**p
        a(4) = d(4)/(ABS(b(4)) + Epsilon_WENO)**p
#endif

        Sum_a = a(0) + a(1) + a(2) + a(3) + a(4) 

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1) + a(2)*V(2) + &
                                 a(3)*V(3) + a(4)*V(4)            )/Sum_a
            
    enddo

    CASE (11)
      do i = i0,i1  
#include "WENO_Euler_Long_Order_11_JS.i"

#ifdef NO_ABS_BETA_JS
        a(0) = d(0)/(    b(0)  + Epsilon_WENO)**p
        a(1) = d(1)/(    b(1)  + Epsilon_WENO)**p
        a(2) = d(2)/(    b(2)  + Epsilon_WENO)**p
        a(3) = d(3)/(    b(3)  + Epsilon_WENO)**p
        a(4) = d(4)/(    b(4)  + Epsilon_WENO)**p
        a(5) = d(5)/(    b(5)  + Epsilon_WENO)**p
#else
        a(0) = d(0)/(ABS(b(0)) + Epsilon_WENO)**p
        a(1) = d(1)/(ABS(b(1)) + Epsilon_WENO)**p
        a(2) = d(2)/(ABS(b(2)) + Epsilon_WENO)**p
        a(3) = d(3)/(ABS(b(3)) + Epsilon_WENO)**p
        a(4) = d(4)/(ABS(b(4)) + Epsilon_WENO)**p
        a(5) = d(5)/(ABS(b(5)) + Epsilon_WENO)**p
#endif

        Sum_a = a(0) + a(1) + a(2) + a(3) + a(4) + a(5)

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1) + a(2)*V(2) + &
                                 a(3)*V(3) + a(4)*V(4) + a(5)*V(5))/Sum_a
      enddo

  END SELECT

  enddo

  END Subroutine PS_WENO_Euler_Long_JS

END Module WENO_Euler_Long_JS
