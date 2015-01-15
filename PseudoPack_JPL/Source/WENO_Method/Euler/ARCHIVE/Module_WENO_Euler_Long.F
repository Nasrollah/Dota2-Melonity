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
  REALTYPE                   :: Sum_a, Sum_of_a

  REALTYPE, dimension(0:(Order-1)/2                            ) :: V, a, b, bb
  REALTYPE, dimension(0:(Order-1)/2                            ) :: d
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2              ) :: C
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2,0:(Order-1)/2) :: S

  REALTYPE, dimension(0:(Order-1)/2                            ) :: g

  REALTYPE                   :: h, ch
  REALTYPE, dimension(3, 3)  :: A1, A2

  !                           Epsilon_WENO = PS_WENO_Option%Epsilon
  !if (PRESENT(WENO_Epsilon)) Epsilon_WENO = WENO_Epsilon


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

#ifdef MAPPED_WENO
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

#ifdef MAPPED_WENO
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        Sum_a = a(0) + a(1) 

        F_Half(i) = F_Half(i) + (a(0)*V(0) + a(1)*V(1))/Sum_a 
      enddo
!$OMP END PARALLEL DO

    CASE (5)
!$OMP PARALLEL DO

      !A1(1,1) = -NINE/TEN;   A1(1,2) =  THREE/FOUR; A1(1,3) =  THREE/(TEN+TEN)
      !A1(2,1) =  THREE/FIVE; A1(2,2) = -THREE/FOUR; A1(2,3) =  THREE/(TEN+TEN)
      !A1(3,1) =  THREE/TEN;  A1(3,2) =  ZERO;       A1(3,3) = -THREE/TEN
      !A2(1,1) = -ONE/TEN;    A2(1,2) =  ZERO;       A2(1,3) =  ONE/TEN
      !A2(2,1) = -ONE/TEN;    A2(2,2) = -HALF;       A2(2,3) =  THREE/FIVE
      !A2(3,1) =  ONE/FIVE;   A2(3,2) =  HALF;       A2(3,3) = -SEVEN/TEN


      do i = i0,i1  
        V(0) = C(0,0)*F(i,-2,l) + C(0,1)*F(i,-1,l) + C(0,2)*F(i, 0,l)
        V(1) = C(1,0)*F(i,-1,l) + C(1,1)*F(i, 0,l) + C(1,2)*F(i, 1,l)
        V(2) = C(2,0)*F(i, 0,l) + C(2,1)*F(i, 1,l) + C(2,2)*F(i, 2,l)

! zico = new_weno + weno_zico 
! min  = new_weno + weno_min  
! js = sem new_weno 
! hap = sem new_weno + mapped_weno
#ifdef NEW_WENO
! new weno nao presta 
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


          if (b(0) < 1d-16) b(0) = 1d-16
          if (b(1) < 1d-16) b(1) = 1d-16
          if (b(2) < 1d-16) b(2) = 1d-16
          h = max(abs(b(0)-b(2)), 1d-16)  ! h = beta5

#ifdef WENO_ZICO
          ! Zico
          a(0) = d(0)*((b(0)+h)/b(0))/2
          a(1) = d(1)*((b(1)+h)/b(1))/2
          a(2) = d(2)*((b(2)+h)/b(2))/2
#endif 
#ifdef WENO_MIN
          ! Minimum
          a(0) = d(0)*((b(0)+h)/min(b(0),h))/2
          a(1) = d(1)*((b(1)+h)/min(b(1),h))/2
          a(2) = d(2)*((b(2)+h)/min(b(2),h))/2
#endif 
#else
!JS 
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
!js   epsilon  is 10^-6
!hap  epsilon  is 10^-16
!      Epsilon_WENO = 1d-6     
      Epsilon_WENO = 1.0d-10     

#ifdef MAPPED_WENO 
!     Epsilon_WENO = 1d-16     
      Epsilon_WENO = 1d-10
#endif 

          a(0) = d(0)/(Epsilon_WENO + b(0))**2
          a(1) = d(1)/(Epsilon_WENO + b(1))**2
          a(2) = d(2)/(Epsilon_WENO + b(2))**2
!hap 
! erase this. Rafael said so
#ifdef MAPPED_WENO
!            if (PS_WENO_Option%Weight_Type == 1) 
#if 0
            a = W(d,a)
#else
  Sum_of_a = a(0) + a(1) + a(2)

  g(0) = a(0)/Sum_of_a
  g(1) = a(1)/Sum_of_a
  g(2) = a(2)/Sum_of_a

  a    = g*(d*(d+ONE)+g*(g   - THREE*d))/                 &
           (d* d     +g*(ONE -   TWO*d))
#endif

#endif

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

#ifdef NEW_WENO
          b = abs(b)/240
          where (b < 1d-16) b = 1d-16
          h = min(b(0), b(1), b(2), b(3))
          ch = log(h)/log(1.0/10)         

          if (ch >= 7) then
            h = 25*h
          elseif (ch >= 0) then
            h = (3*ch + 4)*h
          else
            h = 4*h
          end if

          where (b > h)
            b = 1
          elsewhere
            b = 0
          end where

          h = 8*b(0) + 4*b(1) + 2*b(2) + b(3)
          select case (int(h))
            case (0) ! b = (0,0,0,0)
              a(0) = 1.0/35; a(1) = 12.0/35; a(2) = 18.0/35; a(3) = 4.0/35;
            case (1) ! b = (0,0,0,1)
              a(0) = 1.0/15; a(1) = 8.0/15; a(2) = 2.0/5; a(3) = 0.0;
            case (2,3) ! b = (0,0,1,0), (0,0,1,1)
              a(0) = 1.0/5; a(1) = 4.0/5; a(2) = 0.0; a(3) = 0.0;
            case (4,12) ! b = (0,1,0,0), (1,1,0,0)
              a(0) = 0.0; a(1) = 0.0; a(2) = 3.0/5; a(3) = 2.0/5;
            case (7) ! b = (0,1,1,1)
              a(0) = 1.0; a(1) = 0.0; a(2) = 0.0; a(3) = 0.0;
            case (8) ! b = (1,0,0,0)
              a(0) = 0.0; a(1) = 1.0/5; a(2) = 3.0/5; a(3) = 1.0/5;
            case (9) ! b = (1,0,0,1)
              a(0) = 0.0; a(1) = 2.0/5; a(2) = 3.0/5; a(3) = 0.0;
            case (11) ! b = (1,0,1,1)
              a(0) = 0.0; a(1) = 1.0; a(2) = 0.0; a(3) = 0.0;
            case (5,13) ! b = (0,1,0,1), (1,1,0,1)
              a(0) = 0.0; a(1) = 0.0; a(2) = 1.0; a(3) = 0.0;
            case (6,10,14) ! b = (0,1,1,0), (1,0,1,0), (1,1,1,0)
              a(0) = 0.0; a(1) = 0.0; a(2) = 0.0; a(3) = 1.0;
          end select
#else
          a(0) = d(0)/(Epsilon_WENO + b(0))**2
          a(1) = d(1)/(Epsilon_WENO + b(1))**2
          a(2) = d(2)/(Epsilon_WENO + b(2))**2
          a(3) = d(3)/(Epsilon_WENO + b(3))**2

#ifdef MAPPED_WENO
            if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
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

#ifdef NEW_WENO
          b = abs(b)/5040
          where (b < 1d-16) b = 1d-16
          h = min(b(0), b(1), b(2), b(3), b(4))
          ch = log(h)/log(1.0/10)         

          if (ch >= 9) then
            h = 140*h
          elseif (ch >= 0) then
            h = (15*ch + 5)*h
          else
            h = 5*h
          end if

          where (b > h)
            b = 1
          elsewhere
            b = 0
          end where

          h = 16*b(0) + 8*b(1) + 4*b(2) + 2*b(3) + b(4)
          select case (int(h))
            case (0) ! (00000)
              a(0) = 1.0/126; a(1) = 10.0/63; a(2) = 10.0/21; a(3) = 20.0/63; a(4) = 5.0/126;
            case (1) ! (00001)
              a(0) = 1.0/56; a(1) = 15.0/56; a(2) = 15.0/28; a(3) = 5.0/28; a(4) = 0.0;
            case (16) ! (10000)
              a(0) = 0.0; a(1) = 1.0/14; a(2) = 3.0/7; a(3) = 3.0/7; a(4) = 1.0/14;
            case (3, 2) ! (00011)
              a(0) = 1.0/21; a(1) = 10.0/21; a(2) = 10.0/21; a(3) = 0.0; a(4) = 0.0;
            case (17) ! (10001)
              a(0) = 0.0; a(1) = 1.0/7; a(2) = 4.0/7; a(3) = 2.0/7; a(4) = 0.0;
            case (24, 8) ! (11000)
              a(0) = 0.0; a(1) = 0.0; a(2) = 2.0/7; a(3) = 4.0/7; a(4) = 1.0/7;
            case (7, 5, 6) ! (00111)
              a(0) = 1.0/6; a(1) = 5.0/6; a(2) = 0.0; a(3) = 0.0; a(4) = 0.0;
            case (19, 18) ! (10011)
              a(0) = 0.0; a(1) = 1.0/3; a(2) = 2.0/3; a(3) = 0.0; a(4) = 0.0;
            case (25, 9) ! (11001)
              a(0) = 0.0; a(1) = 0.0; a(2) = 1.0/2; a(3) = 1.0/2; a(4) = 0.0;
            case (28, 4, 12, 20) ! (11100)
              a(0) = 0.0; a(1) = 0.0; a(2) = 0.0; a(3) = 2.0/3; a(4) = 1.0/3;
            case (15) ! (01111)
              a(0) = 1.0; a(1) = 0.0; a(2) = 0.0; a(3) = 0.0; a(4) = 0.0;
            case (23) ! (10111)
              a(0) = 0.0; a(1) = 1.0; a(2) = 0.0; a(3) = 0.0; a(4) = 0.0;
            case (27, 10, 11, 26) ! (11011)
              a(0) = 0.0; a(1) = 0.0; a(2) = 1.0; a(3) = 0.0; a(4) = 0.0;
            case (29, 13, 21) ! (11101)
              a(0) = 0.0; a(1) = 0.0; a(2) = 0.0; a(3) = 1.0; a(4) = 0.0;
            case (30, 14, 22) ! (11110)
              a(0) = 0.0; a(1) = 0.0; a(2) = 0.0; a(3) = 0.0; a(4) = 1.0;
          end select
#else
          a(0) = d(0)/(Epsilon_WENO + b(0))**2
          a(1) = d(1)/(Epsilon_WENO + b(1))**2
          a(2) = d(2)/(Epsilon_WENO + b(2))**2
          a(3) = d(3)/(Epsilon_WENO + b(3))**2
          a(4) = d(4)/(Epsilon_WENO + b(4))**2

#ifdef MAPPED_WENO
            if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
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

#ifdef NEW_WENO
          b = abs(b)/120960
          where (b < 1d-16) b = 1d-16
          h = min(b(0), b(1), b(2), b(3), b(4), b(5))
          ch = log(h)/log(1.0/10)         

          if (ch >= 11) then
            h = 72*h
          elseif (ch >= 0) then
            h = (6*ch + 6)*h
          else
            h = 6*h
          end if

          where (b > h)
            b = 1
          elsewhere
            b = 0
          end where

          h = 32*b(0) + 16*b(1) + 8*b(2) + 4*b(3) + 2*b(4) + b(5)
          select case (int(h))
            case (0) ! (000000)
              a(0)=1.0/462; a(1)=5.0/77; a(2)=25.0/77; a(3)=100.0/231; a(4)=25.0/154; a(5)=1.0/77;
            case (1) ! (000001)
              a(0)=1.0/210; a(1)=4.0/35; a(2)=3.0/7; a(3)=8.0/21; a(4)=1.0/14; a(5)=0.0;
            case (32) ! (100000)
              a(0)=0.0; a(1)=1.0/42; a(2)=5.0/21; a(3)=10.0/21; a(4)=5.0/21; a(5)=1.0/42;
            case (3,2) ! (000011)
              a(0)=1.0/84; a(1)=3.0/14; a(2)=15.0/28; a(3)=5.0/21; a(4)=0.0; a(5)=0.0;
            case (33) ! (100001)
              a(0)=0.0; a(1)=1.0/21; a(2)=5.0/14; a(3)=10.0/21; a(4)=5.0/42; a(5)=0.0;
            case (48,16) ! (110000)
              a(0)=0.0; a(1)=0.0; a(2)=5.0/42; a(3)=10.0/21; a(4)=5.0/14; a(5)=1.0/21;
            case (7,4,5,6) ! (000111)
              a(0)=1.0/28; a(1)=3.0/7; a(2)=15.0/28; a(3)=0.0; a(4)=0.0; a(5)=0.0;
            case (35,34) ! (100011)
              a(0)=0.0; a(1)=3.0/28; a(2)=15.0/28; a(3)=5.0/14; a(4)=0.0; a(5)=0.0;
            case (49,17) ! (110001)
              a(0)=0.0; a(1)=0.0; a(2)=3.0/14; a(3)=4.0/7; a(4)=3.0/14; a(5)=0.0;
            case (56,8,24,40) ! (111000)
              a(0)=0.0; a(1)=0.0; a(2)=0.0; a(3)=5.0/14; a(4)=15.0/28; a(5)=3.0/28;
            case (15,10,11,13,14) ! (001111)
              a(0)=1.0/7; a(1)=6.0/7; a(2)=0.0; a(3)=0.0; a(4)=0.0; a(5)=0.0;
            case (39,37,38) ! (100111)
              a(0)=0.0; a(1)=2.0/7; a(2)=5.0/7; a(3)=0.0; a(4)=0.0; a(5)=0.0;
            case (51,18,19,50) ! (110011)
              a(0)=0.0; a(1)=0.0; a(2)=3.0/7; a(3)=4.0/7; a(4)=0.0; a(5)=0.0;
            case (57,9,25,41) ! (111001)
              a(0)=0.0; a(1)=0.0; a(2)=0.0; a(3)=4.0/7; a(4)=3.0/7; a(5)=0.0;
            case (60,12,20,28,36,44,52) ! (111100)
              a(0)=0.0; a(1)=0.0; a(2)=0.0; a(3)=0.0; a(4)=5.0/7; a(5)=2.0/7;
            case (31) ! (011111)
              a(0)=1.0; a(1)=0.0; a(2)=0.0; a(3)=0.0; a(4)=0.0; a(5)=0.0;
            case (47) ! (101111)
              a(0)=0.0; a(1)=1.0; a(2)=0.0; a(3)=0.0; a(4)=0.0; a(5)=0.0;
            case (55,21,22,23,54) ! (110111)
              a(0)=0.0; a(1)=0.0; a(2)=1.0; a(3)=0.0; a(4)=0.0; a(5)=0.0;
            case (59,26,27,42,43,58) ! (111011)
              a(0)=0.0; a(1)=0.0; a(2)=0.0; a(3)=1.0; a(4)=0.0; a(5)=0.0;
            case (61,29,45,53) ! (111101)
              a(0)=0.0; a(1)=0.0; a(2)=0.0; a(3)=0.0; a(4)=1.0; a(5)=0.0;
            case (62,30,46) ! (111110)
              a(0)=0.0; a(1)=0.0; a(2)=0.0; a(3)=0.0; a(4)=0.0; a(5)=1.0;
          end select
#else
          a(0) = d(0)/(Epsilon_WENO + b(0))**2
          a(1) = d(1)/(Epsilon_WENO + b(1))**2
          a(2) = d(2)/(Epsilon_WENO + b(2))**2
          a(3) = d(3)/(Epsilon_WENO + b(3))**2
          a(4) = d(4)/(Epsilon_WENO + b(4))**2
          a(5) = d(5)/(Epsilon_WENO + b(5))**2

#ifdef MAPPED_WENO
            if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
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
