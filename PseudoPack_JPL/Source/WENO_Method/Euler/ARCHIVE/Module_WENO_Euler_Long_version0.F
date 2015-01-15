#include "pseudopack.h"


Module WENO_Euler_Long

  USE WENO_Coefficients
  USE WENO_Option

implicit NONE

PRIVATE
PUBLIC  :: PS_WENO_Euler_Long

CONTAINS

  Subroutine PS_WENO_Euler_Long (Order, N0, N5, i0, i1, F, F_Half,  &
                                        WENO_Epsilon)

  integer                    :: Order, N0, N5, i0, i1

  REALTYPE, dimension(N0:N5,-(Order-1)/2:(Order-1)/2,2) :: F
  REALTYPE, dimension(N0:N5                           ) :: F_Half

  REALTYPE, OPTIONAL         :: WENO_Epsilon
  REALTYPE                   :: Epsilon_WENO

  integer                    :: i, l
  REALTYPE                   :: Sum_a

  REALTYPE, dimension(0:(Order-1)/2                            ) :: V, a, b
  REALTYPE, dimension(0:(Order-1)/2                            ) :: d
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2              ) :: C
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2,0:(Order-1)/2) :: S

                             Epsilon_WENO = PS_WENO_Option%Epsilon
  if (PRESENT(WENO_Epsilon)) Epsilon_WENO = WENO_Epsilon

  call PS_WENO_Coefficient (Order, d, C, S)

  F_Half = ZERO

  do l = 1,2

  SELECT CASE (Order)
    CASE (1)
!$OMP PARALLEL DO
      do i = i0,i1  
          V(0) = C(0,0)*F(i, 0,l)
          b(0) = F(i, 0,l)*(S(0,0,0)*F(i, 0,l))

        a(0) = d(0)/(Epsilon_WENO + b(0))**2

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        Sum_a = a(0) 

        F_Half(i) = F_Half(i) + (a(0)*V(0))/Sum_a 
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = i0,i1  
          V(0) = C(0,0)*F(i,-1,l) + C(0,1)*F(i, 0,l) 
          V(1) = C(1,0)*F(i, 0,l) + C(1,1)*F(i, 1,l) 

          b(0) = F(i,-1,l)*(S(0,0,0)*F(i,-1,l) + S(0,1,0)*F(i, 0,l))+ &
                 F(i, 0,l)*(S(1,1,0)*F(i, 0,l)                     )

          b(1) = F(i, 0,l)*(S(0,0,1)*F(i, 0,l) + S(0,1,1)*F(i, 1,l))+ &
                 F(i, 1,l)*(S(1,1,1)*F(i, 1,l)                     )

        a(0) = d(0)/(Epsilon_WENO + b(0))**2
        a(1) = d(1)/(Epsilon_WENO + b(1))**2

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        Sum_a = a(0) + a(1) 

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1))/Sum_a 
      enddo
!$OMP END PARALLEL DO

    CASE (5)
!$OMP PARALLEL DO
      do i = i0,i1  
          V(0) = C(0,0)*F(i,-2,l) + C(0,1)*F(i,-1,l) + C(0,2)*F(i, 0,l)
          V(1) = C(1,0)*F(i,-1,l) + C(1,1)*F(i, 0,l) + C(1,2)*F(i, 1,l)
          V(2) = C(2,0)*F(i, 0,l) + C(2,1)*F(i, 1,l) + C(2,2)*F(i, 2,l)

          b(0) = F(i,-2,l)*(S(0,0,0)*F(i,-2,l) + S(0,1,0)*F(i,-1,l) + &
                            S(0,2,0)*F(i, 0,l)                        ) + &
                 F(i,-1,l)*(S(1,1,0)*F(i,-1,l) + S(1,2,0)*F(i, 0,l)   ) + &
                 F(i, 0,l)*(S(2,2,0)*F(i, 0,l)                        )

          b(1) = F(i,-1,l)*(S(0,0,1)*F(i,-1,l) + S(0,1,1)*F(i, 0,l) + &
                            S(0,2,1)*F(i, 1,l)                        ) + &
                 F(i, 0,l)*(S(1,1,1)*F(i, 0,l) + S(1,2,1)*F(i, 1,l)   ) + &
                 F(i, 1,l)*(S(2,2,1)*F(i, 1,l)                        )

          b(2) = F(i, 0,l)*(S(0,0,2)*F(i, 0,l) + S(0,1,2)*F(i, 1,l) + &
                            S(0,2,2)*F(i, 2,l)                        ) + &
                 F(i, 1,l)*(S(1,1,2)*F(i, 1,l) + S(1,2,2)*F(i, 2,l)   ) + &
                 F(i, 2,l)*(S(2,2,2)*F(i, 2,l)                        )

        a(0) = d(0)/(Epsilon_WENO + b(0))**2
        a(1) = d(1)/(Epsilon_WENO + b(1))**2
        a(2) = d(2)/(Epsilon_WENO + b(2))**2

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        Sum_a = a(0) + a(1) + a(2)

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1) + a(2)*V(2))/Sum_a
      enddo
!$OMP END PARALLEL DO

    CASE (7) 
!$OMP PARALLEL DO
      do i = i0,i1  
          V(0) = C(0,0)*F(i,-3,l) + C(0,1)*F(i,-2,l) + C(0,2)*F(i,-1,l) + &
                 C(0,3)*F(i, 0,l) 
          V(1) = C(1,0)*F(i,-2,l) + C(1,1)*F(i,-1,l) + C(1,2)*F(i, 0,l) + &
                 C(1,3)*F(i, 1,l) 
          V(2) = C(2,0)*F(i,-1,l) + C(2,1)*F(i, 0,l) + C(2,2)*F(i, 1,l) + &
                 C(2,3)*F(i, 2,l) 
          V(3) = C(3,0)*F(i, 0,l) + C(3,1)*F(i, 1,l) + C(3,2)*F(i, 2,l) + &
                 C(3,3)*F(i, 3,l) 

          b(0) = F(i,-3,l)*(S(0,0,0)*F(i,-3,l) + S(0,1,0)*F(i,-2,l) + &
                            S(0,2,0)*F(i,-1,l) + S(0,3,0)*F(i, 0,l)   )+ &
                 F(i,-2,l)*(S(1,1,0)*F(i,-2,l) + S(1,2,0)*F(i,-1,l) + &
                            S(1,3,0)*F(i, 0,l)                        ) + &
                 F(i,-1,l)*(S(2,2,0)*F(i,-1,l) + S(2,3,0)*F(i, 0,l)   )+ &
                 F(i, 0,l)*(S(3,3,0)*F(i, 0,l)                        ) 

          b(1) = F(i,-2,l)*(S(0,0,1)*F(i,-2,l) + S(0,1,1)*F(i,-1,l) + &
                            S(0,2,1)*F(i, 0,l) + S(0,3,1)*F(i, 1,l)   )+ &
                 F(i,-1,l)*(S(1,1,1)*F(i,-1,l) + S(1,2,1)*F(i, 0,l) + &
                            S(1,3,1)*F(i, 1,l)                        )+ &
                 F(i, 0,l)*(S(2,2,1)*F(i, 0,l) + S(2,3,1)*F(i, 1,l)   )+ &
                 F(i, 1,l)*(S(3,3,1)*F(i, 1,l)                        )

          b(2) = F(i,-1,l)*(S(0,0,2)*F(i,-1,l) + S(0,1,2)*F(i, 0,l) + &
                            S(0,2,2)*F(i, 1,l) + S(0,3,2)*F(i, 2,l)   )+ &
                 F(i, 0,l)*(S(1,1,2)*F(i, 0,l) + S(1,2,2)*F(i, 1,l) + &
                            S(1,3,2)*F(i, 2,l)                        )+ &
                 F(i, 1,l)*(S(2,2,2)*F(i, 1,l) + S(2,3,2)*F(i, 2,l)   )+ &
                 F(i, 2,l)*(S(3,3,2)*F(i, 2,l)                        )

          b(3) = F(i, 0,l)*(S(0,0,3)*F(i, 0,l) + S(0,1,3)*F(i, 1,l) + &
                            S(0,2,3)*F(i, 2,l) + S(0,3,3)*F(i, 3,l)   )+ &
                 F(i, 1,l)*(S(1,1,3)*F(i, 1,l) + S(1,2,3)*F(i, 2,l) + &
                            S(1,3,3)*F(i, 3,l)                        )+ &
                 F(i, 2,l)*(S(2,2,3)*F(i, 2,l) + S(2,3,3)*F(i, 3,l)   )+ &
                 F(i, 3,l)*(S(3,3,3)*F(i, 3,l)                        )

        a(0) = d(0)/(Epsilon_WENO + b(0))**2
        a(1) = d(1)/(Epsilon_WENO + b(1))**2
        a(2) = d(2)/(Epsilon_WENO + b(2))**2
        a(3) = d(3)/(Epsilon_WENO + b(3))**2

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        Sum_a = a(0) + a(1) + a(2) + a(3) 

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1) + a(2)*V(2) + &
                                 a(3)*V(3)                        )/Sum_a
      enddo
!$OMP END PARALLEL DO

    CASE (9) 
!$OMP PARALLEL DO
      do i = i0,i1  
          V(0) = C(0,0)*F(i,-4,l) + C(0,1)*F(i,-3,l) + C(0,2)*F(i,-2,l) + &
                 C(0,3)*F(i,-1,l) + C(0,4)*F(i, 0,l) 
          V(1) = C(1,0)*F(i,-3,l) + C(1,1)*F(i,-2,l) + C(1,2)*F(i,-1,l) + &
                 C(1,3)*F(i, 0,l) + C(1,4)*F(i, 1,l) 
          V(2) = C(2,0)*F(i,-2,l) + C(2,1)*F(i,-1,l) + C(2,2)*F(i, 0,l) + &
                 C(2,3)*F(i, 1,l) + C(2,4)*F(i, 2,l) 
          V(3) = C(3,0)*F(i,-1,l) + C(3,1)*F(i, 0,l) + C(3,2)*F(i, 1,l) + &
                 C(3,3)*F(i, 2,l) + C(3,4)*F(i, 3,l) 
          V(4) = C(4,0)*F(i, 0,l) + C(4,1)*F(i, 1,l) + C(4,2)*F(i, 2,l) + &
                 C(4,3)*F(i, 3,l) + C(4,4)*F(i, 4,l) 

          b(0) = F(i,-4,l)*(S(0,0,0)*F(i,-4,l) + S(0,1,0)*F(i,-3,l) + &
                            S(0,2,0)*F(i,-2,l) + S(0,3,0)*F(i,-1,l) + &
                            S(0,4,0)*F(i, 0,l)                        ) + &
                 F(i,-3,l)*(S(1,1,0)*F(i,-3,l) + S(1,2,0)*F(i,-2,l) + &
                            S(1,3,0)*F(i,-1,l) + S(1,4,0)*F(i, 0,l)   ) + &
                 F(i,-2,l)*(S(2,2,0)*F(i,-2,l) + S(2,3,0)*F(i,-1,l) + &
                            S(2,4,0)*F(i, 0,l)                        ) + &
                 F(i,-1,l)*(S(3,3,0)*F(i,-1,l) + S(3,4,0)*F(i, 0,l)   ) + &
                 F(i, 0,l)*(S(4,4,0)*F(i, 0,l)                        ) 

          b(1) = F(i,-3,l)*(S(0,0,1)*F(i,-3,l) + S(0,1,1)*F(i,-2,l) + &
                            S(0,2,1)*F(i,-1,l) + S(0,3,1)*F(i, 0,l) + &
                            S(0,4,1)*F(i, 1,l)                        ) + &
                 F(i,-2,l)*(S(1,1,1)*F(i,-2,l) + S(1,2,1)*F(i,-1,l) + &
                            S(1,3,1)*F(i, 0,l) + S(1,4,1)*F(i, 1,l)   ) + &
                 F(i,-1,l)*(S(2,2,1)*F(i,-1,l) + S(2,3,1)*F(i, 0,l) + &
                            S(2,4,1)*F(i, 1,l)                        ) + &
                 F(i, 0,l)*(S(3,3,1)*F(i, 0,l) + S(3,4,1)*F(i, 1,l)   ) + &
                 F(i, 1,l)*(S(4,4,1)*F(i, 1,l)                        ) 

          b(2) = F(i,-2,l)*(S(0,0,2)*F(i,-2,l) + S(0,1,2)*F(i,-1,l) + &
                            S(0,2,2)*F(i, 0,l) + S(0,3,2)*F(i, 1,l) + &
                            S(0,4,2)*F(i, 2,l)                        ) + &
                 F(i,-1,l)*(S(1,1,2)*F(i,-1,l) + S(1,2,2)*F(i, 0,l) + &
                            S(1,3,2)*F(i, 1,l) + S(1,4,2)*F(i, 2,l)   ) + &
                 F(i, 0,l)*(S(2,2,2)*F(i, 0,l) + S(2,3,2)*F(i, 1,l) + &
                            S(2,4,2)*F(i, 2,l)                        ) + &
                 F(i, 1,l)*(S(3,3,2)*F(i, 1,l) + S(3,4,2)*F(i, 2,l)   ) + &
                 F(i, 2,l)*(S(4,4,2)*F(i, 2,l)                        ) 

          b(3) = F(i,-1,l)*(S(0,0,3)*F(i,-1,l) + S(0,1,3)*F(i, 0,l) + &
                            S(0,2,3)*F(i, 1,l) + S(0,3,3)*F(i, 2,l) + &
                            S(0,4,3)*F(i, 3,l)                        ) + &
                 F(i, 0,l)*(S(1,1,3)*F(i, 0,l) + S(1,2,3)*F(i, 1,l) + &
                            S(1,3,3)*F(i, 2,l) + S(1,4,3)*F(i, 3,l)   ) + &
                 F(i, 1,l)*(S(2,2,3)*F(i, 1,l) + S(2,3,3)*F(i, 2,l) + &
                            S(2,4,3)*F(i, 3,l)                        ) + &
                 F(i, 2,l)*(S(3,3,3)*F(i, 2,l) + S(3,4,3)*F(i, 3,l)   ) + &
                 F(i, 3,l)*(S(4,4,3)*F(i, 3,l)                        ) 

          b(4) = F(i, 0,l)*(S(0,0,4)*F(i, 0,l) + S(0,1,4)*F(i, 1,l) + &
                            S(0,2,4)*F(i, 2,l) + S(0,3,4)*F(i, 3,l) + &
                            S(0,4,4)*F(i, 4,l)                        ) + &
                 F(i, 1,l)*(S(1,1,4)*F(i, 1,l) + S(1,2,4)*F(i, 2,l) + &
                            S(1,3,4)*F(i, 3,l) + S(1,4,4)*F(i, 4,l)   ) + &
                 F(i, 2,l)*(S(2,2,4)*F(i, 2,l) + S(2,3,4)*F(i, 3,l) + &
                            S(2,4,4)*F(i, 4,l)                        ) + &
                 F(i, 3,l)*(S(3,3,4)*F(i, 3,l) + S(3,4,4)*F(i, 4,l)   ) + &
                 F(i, 4,l)*(S(4,4,4)*F(i, 4,l)                        ) 

        a(0) = d(0)/(Epsilon_WENO + b(0))**2
        a(1) = d(1)/(Epsilon_WENO + b(1))**2
        a(2) = d(2)/(Epsilon_WENO + b(2))**2
        a(3) = d(3)/(Epsilon_WENO + b(3))**2
        a(4) = d(4)/(Epsilon_WENO + b(4))**2

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        Sum_a = a(0) + a(1) + a(2) + a(3) + a(4) 

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1) + a(2)*V(2) + &
                                 a(3)*V(3) + a(4)*V(4)            )/Sum_a
            
    enddo
!$OMP END PARALLEL DO

    CASE (11)
!$OMP PARALLEL DO
      do i = i0,i1  
          V(0) = C(0,0)*F(i,-5,l) + C(0,1)*F(i,-4,l) + C(0,2)*F(i,-3,l) + &
                 C(0,3)*F(i,-2,l) + C(0,4)*F(i,-1,l) + C(0,5)*F(i, 0,l) 
          V(1) = C(1,0)*F(i,-4,l) + C(1,1)*F(i,-3,l) + C(1,2)*F(i,-2,l) + &
                 C(1,3)*F(i,-1,l) + C(1,4)*F(i, 0,l) + C(1,5)*F(i, 1,l) 
          V(2) = C(2,0)*F(i,-3,l) + C(2,1)*F(i,-2,l) + C(2,2)*F(i,-1,l) + &
                 C(2,3)*F(i, 0,l) + C(2,4)*F(i, 1,l) + C(2,5)*F(i, 2,l) 
          V(3) = C(3,0)*F(i,-2,l) + C(3,1)*F(i,-1,l) + C(3,2)*F(i, 0,l) + &
                 C(3,3)*F(i, 1,l) + C(3,4)*F(i, 2,l) + C(3,5)*F(i, 3,l) 
          V(4) = C(4,0)*F(i,-1,l) + C(4,1)*F(i, 0,l) + C(4,2)*F(i, 1,l) + &
                 C(4,3)*F(i, 2,l) + C(4,4)*F(i, 3,l) + C(4,5)*F(i, 4,l)
          V(5) = C(5,0)*F(i, 0,l) + C(5,1)*F(i, 1,l) + C(5,2)*F(i, 2,l) + &
                 C(5,3)*F(i, 3,l) + C(5,4)*F(i, 4,l) + C(5,5)*F(i, 5,l)

          b(0) = F(i,-5,l)*(S(0,0,0)*F(i,-5,l) + S(0,1,0)*F(i,-4,l) + &
                            S(0,2,0)*F(i,-3,l) + S(0,3,0)*F(i,-2,l) + &
                            S(0,4,0)*F(i,-1,l) + S(0,5,0)*F(i, 0,l)   ) + &
                 F(i,-4,l)*(S(1,1,0)*F(i,-4,l) + S(1,2,0)*F(i,-3,l) + &
                            S(1,3,0)*F(i,-2,l) + S(1,4,0)*F(i,-1,l) + &
                            S(1,5,0)*F(i, 0,l)                        ) + &
                 F(i,-3,l)*(S(2,2,0)*F(i,-3,l) + S(2,3,0)*F(i,-2,l) + &
                            S(2,4,0)*F(i,-1,l) + S(2,5,0)*F(i, 0,l)   ) + &
                 F(i,-2,l)*(S(3,3,0)*F(i,-2,l) + S(3,4,0)*F(i,-1,l) + &
                            S(3,5,0)*F(i, 0,l)                        ) + &
                 F(i,-1,l)*(S(4,4,0)*F(i,-1,l) + S(4,5,0)*F(i, 0,l)   ) + &
                 F(i, 0,l)*(S(5,5,0)*F(i, 0,l)                        ) 

          b(1) = F(i,-4,l)*(S(0,0,1)*F(i,-4,l) + S(0,1,1)*F(i,-3,l) + &
                            S(0,2,1)*F(i,-2,l) + S(0,3,1)*F(i,-1,l) + &
                            S(0,4,1)*F(i, 0,l) + S(0,5,1)*F(i, 1,l)   ) + &
                 F(i,-3,l)*(S(1,1,1)*F(i,-3,l) + S(1,2,1)*F(i,-2,l) + &
                            S(1,3,1)*F(i,-1,l) + S(1,4,1)*F(i, 0,l) + &
                            S(1,5,1)*F(i, 1,l)                        ) + &
                 F(i,-2,l)*(S(2,2,1)*F(i,-2,l) + S(2,3,1)*F(i,-1,l) + &
                            S(2,4,1)*F(i, 0,l) + S(2,5,1)*F(i, 1,l)   ) + &
                 F(i,-1,l)*(S(3,3,1)*F(i,-1,l) + S(3,4,1)*F(i, 0,l) + &
                            S(3,5,1)*F(i, 1,l)                        ) + &
                 F(i, 0,l)*(S(4,4,1)*F(i, 0,l) + S(4,5,1)*F(i, 1,l)   ) + &
                 F(i, 1,l)*(S(5,5,1)*F(i, 1,l)                        ) 

          b(2) = F(i,-3,l)*(S(0,0,2)*F(i,-3,l) + S(0,1,2)*F(i,-2,l) + &
                            S(0,2,2)*F(i,-1,l) + S(0,3,2)*F(i, 0,l) + &
                            S(0,4,2)*F(i, 1,l) + S(0,5,2)*F(i, 2,l)   ) + &
                 F(i,-2,l)*(S(1,1,2)*F(i,-2,l) + S(1,2,2)*F(i,-1,l) + &
                            S(1,3,2)*F(i, 0,l) + S(1,4,2)*F(i, 1,l) + &
                            S(1,5,2)*F(i, 2,l)                        ) + &
                 F(i,-1,l)*(S(2,2,2)*F(i,-1,l) + S(2,3,2)*F(i, 0,l) + &
                            S(2,4,2)*F(i, 1,l) + S(2,5,2)*F(i, 2,l)   ) + &
                 F(i, 0,l)*(S(3,3,2)*F(i, 0,l) + S(3,4,2)*F(i, 1,l) + &
                            S(3,5,2)*F(i, 2,l)                        ) + &
                 F(i, 1,l)*(S(4,4,2)*F(i, 1,l) + S(4,5,2)*F(i, 2,l)   ) + &
                 F(i, 2,l)*(S(5,5,2)*F(i, 2,l)                        ) 

          b(3) = F(i,-2,l)*(S(0,0,3)*F(i,-2,l) + S(0,1,3)*F(i,-1,l) + &
                            S(0,2,3)*F(i, 0,l) + S(0,3,3)*F(i, 1,l) + &
                            S(0,4,3)*F(i, 2,l) + S(0,5,3)*F(i, 3,l)   ) + &
                 F(i,-1,l)*(S(1,1,3)*F(i,-1,l) + S(1,2,3)*F(i, 0,l) + &
                            S(1,3,3)*F(i, 1,l) + S(1,4,3)*F(i, 2,l) + &
                            S(1,5,3)*F(i, 3,l)                        ) + &
                 F(i, 0,l)*(S(2,2,3)*F(i, 0,l) + S(2,3,3)*F(i, 1,l) + &
                            S(2,4,3)*F(i, 2,l) + S(2,5,3)*F(i, 3,l)   ) + &
                 F(i, 1,l)*(S(3,3,3)*F(i, 1,l) + S(3,4,3)*F(i, 2,l) + &
                            S(3,5,3)*F(i, 3,l)                        ) + &
                 F(i, 2,l)*(S(4,4,3)*F(i, 2,l) + S(4,5,3)*F(i, 3,l)   ) + &
                 F(i, 3,l)*(S(5,5,3)*F(i, 3,l)                        ) 

          b(4) = F(i,-1,l)*(S(0,0,4)*F(i,-1,l) + S(0,1,4)*F(i, 0,l) + &
                            S(0,2,4)*F(i, 1,l) + S(0,3,4)*F(i, 2,l) + &
                            S(0,4,4)*F(i, 3,l) + S(0,5,4)*F(i, 4,l)   ) + &
                 F(i, 0,l)*(S(1,1,4)*F(i, 0,l) + S(1,2,4)*F(i, 1,l) + &
                            S(1,3,4)*F(i, 2,l) + S(1,4,4)*F(i, 3,l) + &
                            S(1,5,4)*F(i, 4,l)                        ) + &
                 F(i, 1,l)*(S(2,2,4)*F(i, 1,l) + S(2,3,4)*F(i, 2,l) + &
                            S(2,4,4)*F(i, 3,l) + S(2,5,4)*F(i, 4,l)   ) + &
                 F(i, 2,l)*(S(3,3,4)*F(i, 2,l) + S(3,4,4)*F(i, 3,l) + &
                            S(3,5,4)*F(i, 4,l)                        ) + &
                 F(i, 3,l)*(S(4,4,4)*F(i, 3,l) + S(4,5,4)*F(i, 4,l)   ) + &
                 F(i, 4,l)*(S(5,5,4)*F(i, 4,l)                        ) 

          b(5) = F(i, 0,l)*(S(0,0,5)*F(i, 0,l) + S(0,1,5)*F(i, 1,l) + &
                            S(0,2,5)*F(i, 2,l) + S(0,3,5)*F(i, 3,l) + &
                            S(0,4,5)*F(i, 4,l) + S(0,5,5)*F(i, 5,l)   ) + &
                 F(i, 1,l)*(S(1,1,5)*F(i, 1,l) + S(1,2,5)*F(i, 2,l) + &
                            S(1,3,5)*F(i, 3,l) + S(1,4,5)*F(i, 4,l) + &
                            S(1,5,5)*F(i, 5,l)                        ) + &
                 F(i, 2,l)*(S(2,2,5)*F(i, 2,l) + S(2,3,5)*F(i, 3,l) + &
                            S(2,4,5)*F(i, 4,l) + S(2,5,5)*F(i, 5,l)   ) + &
                 F(i, 3,l)*(S(3,3,5)*F(i, 3,l) + S(3,4,5)*F(i, 4,l) + &
                            S(3,5,5)*F(i, 5,l)                        ) + &
                 F(i, 4,l)*(S(4,4,5)*F(i, 4,l) + S(4,5,5)*F(i, 5,l)   ) + &
                 F(i, 5,l)*(S(5,5,5)*F(i, 5,l)                        ) 

        a(0) = d(0)/(Epsilon_WENO + b(0))**2
        a(1) = d(1)/(Epsilon_WENO + b(1))**2
        a(2) = d(2)/(Epsilon_WENO + b(2))**2
        a(3) = d(3)/(Epsilon_WENO + b(3))**2
        a(4) = d(4)/(Epsilon_WENO + b(4))**2
        a(5) = d(5)/(Epsilon_WENO + b(5))**2

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        Sum_a = a(0) + a(1) + a(2) + a(3) + a(4) + a(5)

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1) + a(2)*V(2) + &
                                 a(3)*V(3) + a(4)*V(4) + a(5)*V(5))/Sum_a
      enddo
!$OMP END PARALLEL DO

  END SELECT

  enddo

  END Subroutine PS_WENO_Euler_Long
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

END Module WENO_Euler_Long
