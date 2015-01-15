  call PS_WENO_Coefficient (Order, d, C, S)

  k = (Order-1)/2

  F_Half = ZERO

  do l = 1,2

  SELECT CASE (Order)
    CASE (1)
      do i = i0,i1  
#include "WENO_Euler_Short_Beta_1.i"

#include "WENO_Euler_JS.i"

#include "WENO_Euler_HAP.i"

        F_Half(i) = F_Half(i) + SUM(a*V)/Sum(a)
      enddo

    CASE (3)
      do i = i0,i1  
#include "WENO_Euler_Short_Beta_3.i"

#include "WENO_Euler_JS.i"

#include "WENO_Euler_HAP.i"

        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo

    CASE (5)
      do i = i0,i1  
#include "WENO_Euler_Short_Beta_5.i"

#include "WENO_Euler_JS.i"

#include "WENO_Euler_HAP.i"

        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo

    CASE (7) 
      do i = i0,i1  
#include "WENO_Euler_Short_Beta_7.i"

#include "WENO_Euler_JS.i"

#include "WENO_Euler_HAP.i"

        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo

    CASE (9) 
      do i = i0,i1  
#include "WENO_Euler_Short_Beta_9.i"

#include "WENO_Euler_JS.i"

#include "WENO_Euler_HAP.i"

        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo

    CASE (11)
      do i = i0,i1  
#include "WENO_Euler_Short_Beta_11.i"

#include "WENO_Euler_JS.i"

#include "WENO_Euler_HAP.i"

        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo

  END SELECT

  enddo
