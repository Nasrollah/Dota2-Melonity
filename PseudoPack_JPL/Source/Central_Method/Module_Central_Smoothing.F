#include "pseudopack.h"

! ----------------------------------------------------------------------------
! JCP 194 (2004) 194-214 
! A family of low dispersive and low dissipative explicit schemes for flow
! and noise computations
! Christophe Bogey, Christophe Bailly
!
! JCP 146 (1998) 82-104
! A General Class of Commutative Filters for LES in Complex Geometries
! Oleg V. Vasilyev, Thomas S. Lund, Parviz Moin
! ----------------------------------------------------------------------------

MODULE Central_Difference_Smoothing

IMPLICIT NONE

INTERFACE PS_Central_Smooth
  MODULE PROCEDURE Central_Smooth_1
  MODULE PROCEDURE Central_Smooth_2
  MODULE PROCEDURE Central_Smooth_3

  MODULE PROCEDURE Central_Smooth_1N
  MODULE PROCEDURE Central_Smooth_2N
  MODULE PROCEDURE Central_Smooth_3N
END INTERFACE

PRIVATE
PUBLIC  :: PS_Central_Smooth
PUBLIC  :: PS_Central_Smooth_Coefficient 

CONTAINS

  Subroutine Central_Smooth_1  (Direction,                         &
                                N0, N5, i_0, i_1,                  & 
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: Direction
  integer  :: N0, N5, i_0, i_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  REALTYPE, dimension(N0:N5)                       :: Q
  integer , dimension(N0:N5)            , OPTIONAL :: Flag

      if (Choice == 0) RETURN

      call Central_Smooth_11   (N0, N5, i_0, i_1,                  & 
                                Q, Choice, Order, Viscosity, Flag)

  END Subroutine Central_Smooth_1 
!
! ================================================================
!
  Subroutine Central_Smooth_2  (Direction,                         &
                                N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: Direction
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  REALTYPE, dimension(N0:N5,M0:M5)                 :: Q
  integer , dimension(N0:N5,M0:M5)      , OPTIONAL :: Flag

  if (Choice == 0) RETURN

  SELECT CASE (Direction)
    CASE (1)
      call Central_Smooth_21   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

    CASE (2)
      call Central_Smooth_22   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

  END SELECT

  END Subroutine Central_Smooth_2 
!
! ================================================================
!
  Subroutine Central_Smooth_3  (Direction,                         &
                                N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: Direction
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)           :: Q
  integer , dimension(N0:N5,M0:M5,K0:K5), OPTIONAL :: Flag

  if (Choice == 0) RETURN

  SELECT CASE (Direction)
    CASE (1)
      call Central_Smooth_31   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

    CASE (2)
      call Central_Smooth_32   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

    CASE (3)
      call Central_Smooth_33   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

  END SELECT

  END Subroutine Central_Smooth_3 
!
! ================================================================
!
  Subroutine Central_Smooth_1N (Direction,                         &
                                N0, N5, i_0, i_1, NV,              & 
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: Direction, NV, n
  integer  :: N0, N5, i_0, i_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  REALTYPE, dimension(N0:N5,NV)                       :: Q
  integer , dimension(N0:N5)               , OPTIONAL :: Flag

  if (Choice == 0) RETURN

  do n = 1,NV
      call Central_Smooth_11   (N0, N5, i_0, i_1,                  & 
                                Q(:,n)    , Choice, Order, Viscosity, Flag)
  enddo

  END Subroutine Central_Smooth_1N
!
! ================================================================
!
  Subroutine Central_Smooth_2N (Direction,                         &
                                N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1, NV,              &
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: Direction, NV, n
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  REALTYPE, dimension(N0:N5,M0:M5,NV)                 :: Q
  integer , dimension(N0:N5,M0:M5)         , OPTIONAL :: Flag

  if (Choice == 0) RETURN

  do n = 1,NV

  SELECT CASE (Direction)
    CASE (1)
      call Central_Smooth_21   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                Q(:,:,n)  , Choice, Order, Viscosity, Flag)

    CASE (2)
      call Central_Smooth_22   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                Q(:,:,n)  , Choice, Order, Viscosity, Flag)

  END SELECT

  enddo

  END Subroutine Central_Smooth_2N
!
! ================================================================
!
  Subroutine Central_Smooth_3N (Direction,                         &
                                N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1, NV,              &
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: Direction, NV, n
  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q
  integer , dimension(N0:N5,M0:M5,K0:K5)   , OPTIONAL :: Flag

  if (Choice == 0) RETURN

  do n = 1,NV

  SELECT CASE (Direction)
    CASE (1)
      call Central_Smooth_31   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q(:,:,:,n), Choice, Order, Viscosity, Flag)

    CASE (2)
      call Central_Smooth_32   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q(:,:,:,n), Choice, Order, Viscosity, Flag)

    CASE (3)
      call Central_Smooth_33   (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q(:,:,:,n), Choice, Order, Viscosity, Flag)

  END SELECT

  enddo

  END Subroutine Central_Smooth_3N
!
! ================================================================
!
  Subroutine Central_Smooth_11 (N0, N5, i_0, i_1,                  & 
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: N0, N5, i_0, i_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  integer            :: Order_S
  REALTYPE           :: Epsilon
  integer            :: Choice_S

  integer                       :: i
  REALTYPE, dimension(0:20)     :: a
  REALTYPE, dimension(N0:N5)    :: Q, QQ

  integer , dimension(N0:N5), OPTIONAL :: Flag

  Choice_S = 0      ; if (PRESENT(Choice   )) Choice_S = Choice
   Order_S = 4      ; if (PRESENT(Order    ))  Order_S = Order
  Epsilon  = 2.0d-1 ; if (PRESENT(Viscosity)) Epsilon  = Viscosity
  
  if (Choice_S     == 0   ) RETURN
  if (ABS(Epsilon) == ZERO) RETURN

#if defined (DEBUG)
  if ((Epsilon <  ZERO) .OR. (Epsilon >  ONE)) then
    write (6 ,*) 'STOP!  Viscosity should be between 0 and 1'
    write (99,*) 'STOP!  Viscosity should be between 0 and 1'
    STOP
  endif
#endif

  Order_S = Order_S/2

  call PS_Central_Smooth_Coefficient (Order_S, a, Choice_S)

  QQ = ZERO

  SELECT CASE (Order_S)
    CASE (1)
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) + a(1)*(Q(i-1) + Q(i+1)) 
      enddo

    CASE (2)
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) + a(1)*(Q(i-1) + Q(i+1)) + &
                            a(2)*(Q(i-2) + Q(i+2)) 
      enddo

    CASE (3)
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) + a(1)*(Q(i-1) + Q(i+1)) + &
                            a(2)*(Q(i-2) + Q(i+2)) + &
                            a(3)*(Q(i-3) + Q(i+3)) 
      enddo

    CASE (4)
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) + a(1)*(Q(i-1) + Q(i+1)) + &
                            a(2)*(Q(i-2) + Q(i+2)) + &
                            a(3)*(Q(i-3) + Q(i+3)) + &
                            a(4)*(Q(i-4) + Q(i+4)) 
      enddo

    CASE (5)
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) + a(1)*(Q(i-1) + Q(i+1)) + &
                            a(2)*(Q(i-2) + Q(i+2)) + &
                            a(3)*(Q(i-3) + Q(i+3)) + &
                            a(4)*(Q(i-4) + Q(i+4)) + &
                            a(5)*(Q(i-5) + Q(i+5)) 
      enddo

    CASE (6)
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) + a(1)*(Q(i-1) + Q(i+1)) + &
                            a(2)*(Q(i-2) + Q(i+2)) + &
                            a(3)*(Q(i-3) + Q(i+3)) + &
                            a(4)*(Q(i-4) + Q(i+4)) + &
                            a(5)*(Q(i-5) + Q(i+5)) + &
                            a(6)*(Q(i-6) + Q(i+6)) 
      enddo

    CASE (7)
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) + a(1)*(Q(i-1) + Q(i+1)) + &
                            a(2)*(Q(i-2) + Q(i+2)) + &
                            a(3)*(Q(i-3) + Q(i+3)) + &
                            a(4)*(Q(i-4) + Q(i+4)) + &
                            a(5)*(Q(i-5) + Q(i+5)) + &
                            a(6)*(Q(i-6) + Q(i+6)) + &
                            a(7)*(Q(i-7) + Q(i+7)) 
      enddo

    CASE (8)
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) + a(1)*(Q(i-1) + Q(i+1)) + &
                            a(2)*(Q(i-2) + Q(i+2)) + &
                            a(3)*(Q(i-3) + Q(i+3)) + &
                            a(4)*(Q(i-4) + Q(i+4)) + &
                            a(5)*(Q(i-5) + Q(i+5)) + &
                            a(6)*(Q(i-6) + Q(i+6)) + &
                            a(7)*(Q(i-7) + Q(i+7)) + &
                            a(8)*(Q(i-8) + Q(i+8)) 
      enddo

    CASE (9)
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) + a(1)*(Q(i-1) + Q(i+1)) + &
                            a(2)*(Q(i-2) + Q(i+2)) + &
                            a(3)*(Q(i-3) + Q(i+3)) + &
                            a(4)*(Q(i-4) + Q(i+4)) + &
                            a(5)*(Q(i-5) + Q(i+5)) + &
                            a(6)*(Q(i-6) + Q(i+6)) + &
                            a(7)*(Q(i-7) + Q(i+7)) + &
                            a(8)*(Q(i-8) + Q(i+8)) + & 
                            a(9)*(Q(i-9) + Q(i+9)) 
      enddo

    CASE DEFAULT
      do i = i_0,i_1
        QQ(i) = a(0)*Q(i) +                           &
                SUM( a(1:Order_S)*(Q(i-1:i-Order_S) + & 
                                   Q(i+1:i+Order_S)) )
      enddo

  END SELECT

  if (PRESENT(Flag)) then
    WHERE (Flag == 0) Q = Q - Epsilon*QQ
  else
                      Q = Q - Epsilon*QQ
  endif

  END Subroutine Central_Smooth_11
!
! ================================================================
!
  Subroutine Central_Smooth_21 (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  integer            :: Order_S
  REALTYPE           :: Epsilon
  integer            :: Choice_S

  integer                             :: i, j
  REALTYPE, dimension(0:20)           :: a
  REALTYPE, dimension(N0:N5,M0:M5)    :: Q, QQ

  integer , dimension(N0:N5,M0:M5), OPTIONAL :: Flag

  Choice_S = 0      ; if (PRESENT(Choice   )) Choice_S = Choice
   Order_S = 4      ; if (PRESENT(Order    ))  Order_S = Order
  Epsilon  = 2.0d-1 ; if (PRESENT(Viscosity)) Epsilon  = Viscosity
  
  if (Choice_S     == 0   ) RETURN
  if (ABS(Epsilon) == ZERO) RETURN

#if defined (DEBUG)
  if ((Epsilon <  ZERO) .OR. (Epsilon >  ONE)) then
    write (6 ,*) 'STOP!  Viscosity should be between 0 and 1'
    write (99,*) 'STOP!  Viscosity should be between 0 and 1'
    STOP
  endif
#endif

  Order_S = Order_S/2

  call PS_Central_Smooth_Coefficient (Order_S, a, Choice_S)

  QQ = ZERO

  SELECT CASE (Order_S)
    CASE (1)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i-1,j) + Q(i+1,j)) 
        enddo
      enddo

    CASE (2)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i-1,j) + Q(i+1,j)) + &
                                  a(2)*(Q(i-2,j) + Q(i+2,j)) 
        enddo
      enddo

    CASE (3)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i-1,j) + Q(i+1,j)) + &
                                  a(2)*(Q(i-2,j) + Q(i+2,j)) + &
                                  a(3)*(Q(i-3,j) + Q(i+3,j)) 
        enddo
      enddo

    CASE (4)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i-1,j) + Q(i+1,j)) + &
                                  a(2)*(Q(i-2,j) + Q(i+2,j)) + &
                                  a(3)*(Q(i-3,j) + Q(i+3,j)) + &
                                  a(4)*(Q(i-4,j) + Q(i+4,j)) 
        enddo
      enddo

    CASE (5)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i-1,j) + Q(i+1,j)) + &
                                  a(2)*(Q(i-2,j) + Q(i+2,j)) + &
                                  a(3)*(Q(i-3,j) + Q(i+3,j)) + &
                                  a(4)*(Q(i-4,j) + Q(i+4,j)) + &
                                  a(5)*(Q(i-5,j) + Q(i+5,j)) 
        enddo
      enddo

    CASE (6)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i-1,j) + Q(i+1,j)) + &
                                  a(2)*(Q(i-2,j) + Q(i+2,j)) + &
                                  a(3)*(Q(i-3,j) + Q(i+3,j)) + &
                                  a(4)*(Q(i-4,j) + Q(i+4,j)) + &
                                  a(5)*(Q(i-5,j) + Q(i+5,j)) + &
                                  a(6)*(Q(i-6,j) + Q(i+6,j)) 
        enddo
      enddo

    CASE (7)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i-1,j) + Q(i+1,j)) + &
                                  a(2)*(Q(i-2,j) + Q(i+2,j)) + &
                                  a(3)*(Q(i-3,j) + Q(i+3,j)) + &
                                  a(4)*(Q(i-4,j) + Q(i+4,j)) + &
                                  a(5)*(Q(i-5,j) + Q(i+5,j)) + &
                                  a(6)*(Q(i-6,j) + Q(i+6,j)) + &
                                  a(7)*(Q(i-7,j) + Q(i+7,j)) 
        enddo
      enddo

    CASE (8)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i-1,j) + Q(i+1,j)) + &
                                  a(2)*(Q(i-2,j) + Q(i+2,j)) + &
                                  a(3)*(Q(i-3,j) + Q(i+3,j)) + &
                                  a(4)*(Q(i-4,j) + Q(i+4,j)) + &
                                  a(5)*(Q(i-5,j) + Q(i+5,j)) + &
                                  a(6)*(Q(i-6,j) + Q(i+6,j)) + &
                                  a(7)*(Q(i-7,j) + Q(i+7,j)) + &
                                  a(8)*(Q(i-8,j) + Q(i+8,j)) 
        enddo
      enddo

    CASE (9)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i-1,j) + Q(i+1,j)) + &
                                  a(2)*(Q(i-2,j) + Q(i+2,j)) + &
                                  a(3)*(Q(i-3,j) + Q(i+3,j)) + &
                                  a(4)*(Q(i-4,j) + Q(i+4,j)) + &
                                  a(5)*(Q(i-5,j) + Q(i+5,j)) + &
                                  a(6)*(Q(i-6,j) + Q(i+6,j)) + &
                                  a(7)*(Q(i-7,j) + Q(i+7,j)) + &
                                  a(8)*(Q(i-8,j) + Q(i+8,j)) + &
                                  a(9)*(Q(i-9,j) + Q(i+9,j)) 
        enddo
      enddo

    CASE DEFAULT
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) +                           &
                    SUM( a(1:Order_S)*(Q(i-1:i-Order_S,j) + &
                                       Q(i+1:i+Order_S,j)) )
        enddo
      enddo

  END SELECT

  if (PRESENT(Flag)) then
    WHERE (Flag == 0) Q = Q - Epsilon*QQ
  else
                      Q = Q - Epsilon*QQ
  endif

  END Subroutine Central_Smooth_21
!
! ================================================================
!
  Subroutine Central_Smooth_22 (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  integer            :: Order_S
  REALTYPE           :: Epsilon
  integer            :: Choice_S

  integer                             :: i, j
  REALTYPE, dimension(0:20)           :: a
  REALTYPE, dimension(N0:N5,M0:M5)    :: Q, QQ

  integer , dimension(N0:N5,M0:M5), OPTIONAL :: Flag

  Choice_S = 0      ; if (PRESENT(Choice   )) Choice_S = Choice
   Order_S = 4      ; if (PRESENT(Order    ))  Order_S = Order
  Epsilon  = 2.0d-1 ; if (PRESENT(Viscosity)) Epsilon  = Viscosity
  
  if (Choice_S     == 0   ) RETURN
  if (ABS(Epsilon) == ZERO) RETURN

#if defined (DEBUG)
  if ((Epsilon <  ZERO) .OR. (Epsilon >  ONE)) then
    write (6 ,*) 'STOP!  Viscosity should be between 0 and 1'
    write (99,*) 'STOP!  Viscosity should be between 0 and 1'
    STOP
  endif
#endif

  Order_S = Order_S/2

  call PS_Central_Smooth_Coefficient (Order_S, a, Choice_S)

  QQ = ZERO

  SELECT CASE (Order_S)
    CASE (1)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i,j-1) + Q(i,j+1)) 
        enddo
      enddo

    CASE (2)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i,j-1) + Q(i,j+1)) + &
                                  a(2)*(Q(i,j-2) + Q(i,j+2)) 
        enddo
      enddo

    CASE (3)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i,j-1) + Q(i,j+1)) + &
                                  a(2)*(Q(i,j-2) + Q(i,j+2)) + &
                                  a(3)*(Q(i,j-3) + Q(i,j+3)) 
        enddo
      enddo

    CASE (4)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i,j-1) + Q(i,j+1)) + &
                                  a(2)*(Q(i,j-2) + Q(i,j+2)) + &
                                  a(3)*(Q(i,j-3) + Q(i,j+3)) + &
                                  a(4)*(Q(i,j-4) + Q(i,j+4)) 
        enddo
      enddo

    CASE (5)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i,j-1) + Q(i,j+1)) + &
                                  a(2)*(Q(i,j-2) + Q(i,j+2)) + &
                                  a(3)*(Q(i,j-3) + Q(i,j+3)) + &
                                  a(4)*(Q(i,j-4) + Q(i,j+4)) + &
                                  a(5)*(Q(i,j-5) + Q(i,j+5)) 
        enddo
      enddo

    CASE (6)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i,j-1) + Q(i,j+1)) + &
                                  a(2)*(Q(i,j-2) + Q(i,j+2)) + &
                                  a(3)*(Q(i,j-3) + Q(i,j+3)) + &
                                  a(4)*(Q(i,j-4) + Q(i,j+4)) + &
                                  a(5)*(Q(i,j-5) + Q(i,j+5)) + &
                                  a(6)*(Q(i,j-6) + Q(i,j+6)) 
        enddo
      enddo

    CASE (7)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i,j-1) + Q(i,j+1)) + &
                                  a(2)*(Q(i,j-2) + Q(i,j+2)) + &
                                  a(3)*(Q(i,j-3) + Q(i,j+3)) + &
                                  a(4)*(Q(i,j-4) + Q(i,j+4)) + &
                                  a(5)*(Q(i,j-5) + Q(i,j+5)) + &
                                  a(6)*(Q(i,j-6) + Q(i,j+6)) + &
                                  a(7)*(Q(i,j-7) + Q(i,j+7)) 
        enddo
      enddo

    CASE (8)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i,j-1) + Q(i,j+1)) + &
                                  a(2)*(Q(i,j-2) + Q(i,j+2)) + &
                                  a(3)*(Q(i,j-3) + Q(i,j+3)) + &
                                  a(4)*(Q(i,j-4) + Q(i,j+4)) + &
                                  a(5)*(Q(i,j-5) + Q(i,j+5)) + &
                                  a(6)*(Q(i,j-6) + Q(i,j+6)) + &
                                  a(7)*(Q(i,j-7) + Q(i,j+7)) + &
                                  a(8)*(Q(i,j-8) + Q(i,j+8)) 
        enddo
      enddo

    CASE (9)
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) + a(1)*(Q(i,j-1) + Q(i,j+1)) + &
                                  a(2)*(Q(i,j-2) + Q(i,j+2)) + &
                                  a(3)*(Q(i,j-3) + Q(i,j+3)) + &
                                  a(4)*(Q(i,j-4) + Q(i,j+4)) + &
                                  a(5)*(Q(i,j-5) + Q(i,j+5)) + &
                                  a(6)*(Q(i,j-6) + Q(i,j+6)) + &
                                  a(7)*(Q(i,j-7) + Q(i,j+7)) + &
                                  a(8)*(Q(i,j-8) + Q(i,j+8)) + &
                                  a(9)*(Q(i,j-9) + Q(i,j+9)) 
        enddo
      enddo

    CASE DEFAULT
      do j = j_0,j_1
        do i = i_0,i_1
          QQ(i,j) = a(0)*Q(i,j) +                            &
                    SUM( a(1:Order_S)*(Q(i,j-1:j-Order_S) +  &
                                       Q(i,j+1:j+Order_S)) )
        enddo
      enddo

  END SELECT

  if (PRESENT(Flag)) then
    WHERE (Flag == 0) Q = Q - Epsilon*QQ
  else
                      Q = Q - Epsilon*QQ
  endif

  END Subroutine Central_Smooth_22
!
! ================================================================
!
  Subroutine Central_Smooth_31 (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  integer            :: Order_S
  REALTYPE           :: Epsilon
  integer            :: Choice_S

  integer                                :: i, j, k
  REALTYPE, dimension(0:20)              :: a
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Q, QQ

  integer , dimension(N0:N5,M0:M5,K0:K5), OPTIONAL :: Flag

  Choice_S = 0      ; if (PRESENT(Choice   )) Choice_S = Choice
   Order_S = 4      ; if (PRESENT(Order    ))  Order_S = Order
  Epsilon  = 2.0d-1 ; if (PRESENT(Viscosity)) Epsilon  = Viscosity
  
  if (Choice_S     == 0   ) RETURN
  if (ABS(Epsilon) == ZERO) RETURN

#if defined (DEBUG)
  if ((Epsilon <  ZERO) .OR. (Epsilon >  ONE)) then
    write (6 ,*) 'STOP!  Viscosity should be between 0 and 1'
    write (99,*) 'STOP!  Viscosity should be between 0 and 1'
    STOP
  endif
#endif

  Order_S = Order_S/2

  call PS_Central_Smooth_Coefficient (Order_S, a, Choice_S)

  QQ = ZERO

  SELECT CASE (Order_S)
    CASE (1)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i-1,j,k) + Q(i+1,j,k)) 
          enddo
        enddo
      enddo

    CASE (2)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i-1,j,k) + Q(i+1,j,k)) + &
                                        a(2)*(Q(i-2,j,k) + Q(i+2,j,k)) 
          enddo
        enddo
      enddo

    CASE (3)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i-1,j,k) + Q(i+1,j,k)) + &
                                        a(2)*(Q(i-2,j,k) + Q(i+2,j,k)) + &
                                        a(3)*(Q(i-3,j,k) + Q(i+3,j,k)) 
          enddo
        enddo
      enddo

    CASE (4)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i-1,j,k) + Q(i+1,j,k)) + &
                                        a(2)*(Q(i-2,j,k) + Q(i+2,j,k)) + &
                                        a(3)*(Q(i-3,j,k) + Q(i+3,j,k)) + &
                                        a(4)*(Q(i-4,j,k) + Q(i+4,j,k)) 
          enddo
        enddo
      enddo

    CASE (5)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i-1,j,k) + Q(i+1,j,k)) + &
                                        a(2)*(Q(i-2,j,k) + Q(i+2,j,k)) + &
                                        a(3)*(Q(i-3,j,k) + Q(i+3,j,k)) + &
                                        a(4)*(Q(i-4,j,k) + Q(i+4,j,k)) + &
                                        a(5)*(Q(i-5,j,k) + Q(i+5,j,k)) 
          enddo
        enddo
      enddo

    CASE (6)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i-1,j,k) + Q(i+1,j,k)) + &
                                        a(2)*(Q(i-2,j,k) + Q(i+2,j,k)) + &
                                        a(3)*(Q(i-3,j,k) + Q(i+3,j,k)) + &
                                        a(4)*(Q(i-4,j,k) + Q(i+4,j,k)) + &
                                        a(5)*(Q(i-5,j,k) + Q(i+5,j,k)) + &
                                        a(6)*(Q(i-6,j,k) + Q(i+6,j,k)) 
          enddo
        enddo
      enddo

    CASE (7)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i-1,j,k) + Q(i+1,j,k)) + &
                                        a(2)*(Q(i-2,j,k) + Q(i+2,j,k)) + &
                                        a(3)*(Q(i-3,j,k) + Q(i+3,j,k)) + &
                                        a(4)*(Q(i-4,j,k) + Q(i+4,j,k)) + &
                                        a(5)*(Q(i-5,j,k) + Q(i+5,j,k)) + &
                                        a(6)*(Q(i-6,j,k) + Q(i+6,j,k)) + &
                                        a(7)*(Q(i-7,j,k) + Q(i+7,j,k)) 
          enddo
        enddo
      enddo

    CASE (8)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i-1,j,k) + Q(i+1,j,k)) + &
                                        a(2)*(Q(i-2,j,k) + Q(i+2,j,k)) + &
                                        a(3)*(Q(i-3,j,k) + Q(i+3,j,k)) + &
                                        a(4)*(Q(i-4,j,k) + Q(i+4,j,k)) + &
                                        a(5)*(Q(i-5,j,k) + Q(i+5,j,k)) + &
                                        a(6)*(Q(i-6,j,k) + Q(i+6,j,k)) + &
                                        a(7)*(Q(i-7,j,k) + Q(i+7,j,k)) + &
                                        a(8)*(Q(i-8,j,k) + Q(i+8,j,k)) 
          enddo
        enddo
      enddo

   CASE (9)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i-1,j,k) + Q(i+1,j,k)) + &
                                        a(2)*(Q(i-2,j,k) + Q(i+2,j,k)) + &
                                        a(3)*(Q(i-3,j,k) + Q(i+3,j,k)) + &
                                        a(4)*(Q(i-4,j,k) + Q(i+4,j,k)) + &
                                        a(5)*(Q(i-5,j,k) + Q(i+5,j,k)) + &
                                        a(6)*(Q(i-6,j,k) + Q(i+6,j,k)) + &
                                        a(7)*(Q(i-7,j,k) + Q(i+7,j,k)) + &
                                        a(8)*(Q(i-8,j,k) + Q(i+8,j,k)) + &
                                        a(9)*(Q(i-9,j,k) + Q(i+9,j,k)) 
          enddo
        enddo
      enddo

    CASE DEFAULT
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) +                            &
                        SUM( a(1:Order_S)*(Q(i-1:i-Order_S,j,k) +  &
                                           Q(i+1:i+Order_S,j,k)) )
          enddo
        enddo
      enddo

  END SELECT

  if (PRESENT(Flag)) then
    WHERE (Flag == 0) Q = Q - Epsilon*QQ
  else
                      Q = Q - Epsilon*QQ
  endif

  END Subroutine Central_Smooth_31
!
! ================================================================
!
  Subroutine Central_Smooth_32 (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  integer            :: Order_S
  REALTYPE           :: Epsilon
  integer            :: Choice_S

  integer                                :: i, j, k
  REALTYPE, dimension(0:20)              :: a
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Q, QQ

  integer , dimension(N0:N5,M0:M5,K0:K5), OPTIONAL :: Flag

  Choice_S = 0      ; if (PRESENT(Choice   )) Choice_S = Choice
   Order_S = 4      ; if (PRESENT(Order    ))  Order_S = Order
  Epsilon  = 2.0d-1 ; if (PRESENT(Viscosity)) Epsilon  = Viscosity
  
  if (Choice_S     == 0   ) RETURN
  if (ABS(Epsilon) == ZERO) RETURN

#if defined (DEBUG)
  if ((Epsilon <  ZERO) .OR. (Epsilon >  ONE)) then
    write (6 ,*) 'STOP!  Viscosity should be between 0 and 1'
    write (99,*) 'STOP!  Viscosity should be between 0 and 1'
    STOP
  endif
#endif

  Order_S = Order_S/2

  call PS_Central_Smooth_Coefficient (Order_S, a, Choice_S)

  QQ = ZERO

  SELECT CASE (Order_S)
    CASE (1)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j-1,k) + Q(i,j+1,k)) 
          enddo
        enddo
      enddo

    CASE (2)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j-1,k) + Q(i,j+1,k)) + &
                                        a(2)*(Q(i,j-2,k) + Q(i,j+2,k)) 
          enddo
        enddo
      enddo

    CASE (3)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j-1,k) + Q(i,j+1,k)) + &
                                        a(2)*(Q(i,j-2,k) + Q(i,j+2,k)) + &
                                        a(3)*(Q(i,j-3,k) + Q(i,j+3,k)) 
          enddo
        enddo
      enddo

    CASE (4)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j-1,k) + Q(i,j+1,k)) + &
                                        a(2)*(Q(i,j-2,k) + Q(i,j+2,k)) + &
                                        a(3)*(Q(i,j-3,k) + Q(i,j+3,k)) + &
                                        a(4)*(Q(i,j-4,k) + Q(i,j+4,k)) 
          enddo
        enddo
      enddo

    CASE (5)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j-1,k) + Q(i,j+1,k)) + &
                                        a(2)*(Q(i,j-2,k) + Q(i,j+2,k)) + &
                                        a(3)*(Q(i,j-3,k) + Q(i,j+3,k)) + &
                                        a(4)*(Q(i,j-4,k) + Q(i,j+4,k)) + &
                                        a(5)*(Q(i,j-5,k) + Q(i,j+5,k)) 
          enddo
        enddo
      enddo

    CASE (6)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j-1,k) + Q(i,j+1,k)) + &
                                        a(2)*(Q(i,j-2,k) + Q(i,j+2,k)) + &
                                        a(3)*(Q(i,j-3,k) + Q(i,j+3,k)) + &
                                        a(4)*(Q(i,j-4,k) + Q(i,j+4,k)) + &
                                        a(5)*(Q(i,j-5,k) + Q(i,j+5,k)) + &
                                        a(6)*(Q(i,j-6,k) + Q(i,j+6,k)) 
          enddo
        enddo
      enddo

    CASE (7)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j-1,k) + Q(i,j+1,k)) + &
                                        a(2)*(Q(i,j-2,k) + Q(i,j+2,k)) + &
                                        a(3)*(Q(i,j-3,k) + Q(i,j+3,k)) + &
                                        a(4)*(Q(i,j-4,k) + Q(i,j+4,k)) + &
                                        a(5)*(Q(i,j-5,k) + Q(i,j+5,k)) + &
                                        a(6)*(Q(i,j-6,k) + Q(i,j+6,k)) + &
                                        a(7)*(Q(i,j-7,k) + Q(i,j+7,k)) 
          enddo
        enddo
      enddo

    CASE (8)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j-1,k) + Q(i,j+1,k)) + &
                                        a(2)*(Q(i,j-2,k) + Q(i,j+2,k)) + &
                                        a(3)*(Q(i,j-3,k) + Q(i,j+3,k)) + &
                                        a(4)*(Q(i,j-4,k) + Q(i,j+4,k)) + &
                                        a(5)*(Q(i,j-5,k) + Q(i,j+5,k)) + &
                                        a(6)*(Q(i,j-6,k) + Q(i,j+6,k)) + &
                                        a(7)*(Q(i,j-7,k) + Q(i,j+7,k)) + &
                                        a(8)*(Q(i,j-8,k) + Q(i,j+8,k)) 
          enddo
        enddo
      enddo

    CASE (9)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j-1,k) + Q(i,j+1,k)) + &
                                        a(2)*(Q(i,j-2,k) + Q(i,j+2,k)) + &
                                        a(3)*(Q(i,j-3,k) + Q(i,j+3,k)) + &
                                        a(4)*(Q(i,j-4,k) + Q(i,j+4,k)) + &
                                        a(5)*(Q(i,j-5,k) + Q(i,j+5,k)) + &
                                        a(6)*(Q(i,j-6,k) + Q(i,j+6,k)) + &
                                        a(7)*(Q(i,j-7,k) + Q(i,j+7,k)) + &
                                        a(8)*(Q(i,j-8,k) + Q(i,j+8,k)) + &
                                        a(9)*(Q(i,j-9,k) + Q(i,j+9,k)) 
          enddo
        enddo
      enddo

    CASE DEFAULT
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) +                            &
                        SUM( a(1:Order_S)*(Q(i,j-1:j-Order_S,k) +  &
                                           Q(i,j+1:j+Order_S,k)) )
          enddo
        enddo
      enddo

  END SELECT

  if (PRESENT(Flag)) then
    WHERE (Flag == 0) Q = Q - Epsilon*QQ
  else
                      Q = Q - Epsilon*QQ
  endif

  END Subroutine Central_Smooth_32
!
! ================================================================
!
  Subroutine Central_Smooth_33 (N0, N5, i_0, i_1,                  & 
                                M0, M5, j_0, j_1,                  &
                                K0, K5, k_0, k_1,                  &
                                Q, Choice, Order, Viscosity, Flag)

  integer  :: N0, N5, i_0, i_1
  integer  :: M0, M5, j_0, j_1
  integer  :: K0, K5, k_0, k_1

  integer , OPTIONAL :: Order
  REALTYPE, OPTIONAL :: Viscosity
  integer , OPTIONAL :: Choice

  integer            :: Order_S
  REALTYPE           :: Epsilon
  integer            :: Choice_S

  integer                                :: i, j, k
  REALTYPE, dimension(0:20)              :: a
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Q, QQ

  integer , dimension(N0:N5,M0:M5,K0:K5), OPTIONAL :: Flag

  Choice_S = 0      ; if (PRESENT(Choice   )) Choice_S = Choice
   Order_S = 4      ; if (PRESENT(Order    ))  Order_S = Order
  Epsilon  = 2.0d-1 ; if (PRESENT(Viscosity)) Epsilon  = Viscosity
  
  if (Choice_S     == 0   ) RETURN
  if (ABS(Epsilon) == ZERO) RETURN

#if defined (DEBUG)
  if ((Epsilon <  ZERO) .OR. (Epsilon >  ONE)) then
    write (6 ,*) 'STOP!  Viscosity should be between 0 and 1'
    write (99,*) 'STOP!  Viscosity should be between 0 and 1'
    STOP
  endif
#endif

  Order_S = Order_S/2

  call PS_Central_Smooth_Coefficient (Order_S, a, Choice_S)

  QQ = ZERO

  SELECT CASE (Order_S)
    CASE (1)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j,k-1) + Q(i,j,k+1)) 
          enddo
        enddo
      enddo

    CASE (2)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j,k-1) + Q(i,j,k+1)) + &
                                        a(2)*(Q(i,j,k-2) + Q(i,j,k+2)) 
          enddo
        enddo
      enddo

    CASE (3)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j,k-1) + Q(i,j,k+1)) + &
                                        a(2)*(Q(i,j,k-2) + Q(i,j,k+2)) + &
                                        a(3)*(Q(i,j,k-3) + Q(i,j,k+3)) 
          enddo
        enddo
      enddo

    CASE (4)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j,k-1) + Q(i,j,k+1)) + &
                                        a(2)*(Q(i,j,k-2) + Q(i,j,k+2)) + &
                                        a(3)*(Q(i,j,k-3) + Q(i,j,k+3)) + &
                                        a(4)*(Q(i,j,k-4) + Q(i,j,k+4)) 
          enddo
        enddo
      enddo

    CASE (5)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j,k-1) + Q(i,j,k+1)) + &
                                        a(2)*(Q(i,j,k-2) + Q(i,j,k+2)) + &
                                        a(3)*(Q(i,j,k-3) + Q(i,j,k+3)) + &
                                        a(4)*(Q(i,j,k-4) + Q(i,j,k+4)) + &
                                        a(5)*(Q(i,j,k-5) + Q(i,j,k+5)) 
          enddo
        enddo
      enddo

    CASE (6)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j,k-1) + Q(i,j,k+1)) + &
                                        a(2)*(Q(i,j,k-2) + Q(i,j,k+2)) + &
                                        a(3)*(Q(i,j,k-3) + Q(i,j,k+3)) + &
                                        a(4)*(Q(i,j,k-4) + Q(i,j,k+4)) + &
                                        a(5)*(Q(i,j,k-5) + Q(i,j,k+5)) + &
                                        a(6)*(Q(i,j,k-6) + Q(i,j,k+6)) 
          enddo
        enddo
      enddo

    CASE (7)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j,k-1) + Q(i,j,k+1)) + &
                                        a(2)*(Q(i,j,k-2) + Q(i,j,k+2)) + &
                                        a(3)*(Q(i,j,k-3) + Q(i,j,k+3)) + &
                                        a(4)*(Q(i,j,k-4) + Q(i,j,k+4)) + &
                                        a(5)*(Q(i,j,k-5) + Q(i,j,k+5)) + &
                                        a(6)*(Q(i,j,k-6) + Q(i,j,k+6)) + &
                                        a(7)*(Q(i,j,k-7) + Q(i,j,k+7)) 
          enddo
        enddo
      enddo

    CASE (8)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j,k-1) + Q(i,j,k+1)) + &
                                        a(2)*(Q(i,j,k-2) + Q(i,j,k+2)) + &
                                        a(3)*(Q(i,j,k-3) + Q(i,j,k+3)) + &
                                        a(4)*(Q(i,j,k-4) + Q(i,j,k+4)) + &
                                        a(5)*(Q(i,j,k-5) + Q(i,j,k+5)) + &
                                        a(6)*(Q(i,j,k-6) + Q(i,j,k+6)) + &
                                        a(7)*(Q(i,j,k-7) + Q(i,j,k+7)) + &
                                        a(8)*(Q(i,j,k-8) + Q(i,j,k+8)) 
          enddo
        enddo
      enddo

    CASE (9)
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) + a(1)*(Q(i,j,k-1) + Q(i,j,k+1)) + &
                                        a(2)*(Q(i,j,k-2) + Q(i,j,k+2)) + &
                                        a(3)*(Q(i,j,k-3) + Q(i,j,k+3)) + &
                                        a(4)*(Q(i,j,k-4) + Q(i,j,k+4)) + &
                                        a(5)*(Q(i,j,k-5) + Q(i,j,k+5)) + &
                                        a(6)*(Q(i,j,k-6) + Q(i,j,k+6)) + &
                                        a(7)*(Q(i,j,k-7) + Q(i,j,k+7)) + &
                                        a(8)*(Q(i,j,k-8) + Q(i,j,k+8)) + &
                                        a(9)*(Q(i,j,k-9) + Q(i,j,k+9)) 
          enddo
        enddo
      enddo

    CASE DEFAULT
      do k = k_0,k_1
        do j = j_0,j_1
          do i = i_0,i_1
            QQ(i,j,k) = a(0)*Q(i,j,k) +                            &
                        SUM( a(1:Order_S)*(Q(i,j,k-1:k-Order_S) +  &
                                           Q(i,j,k+1:k+Order_S)) )
          enddo
        enddo
      enddo

  END SELECT

  if (PRESENT(Flag)) then
    WHERE (Flag == 0) Q = Q - Epsilon*QQ
  else
                      Q = Q - Epsilon*QQ
  endif

  END Subroutine Central_Smooth_33
!
! ================================================================
!
  Subroutine PS_Central_Smooth_Coefficient (Order, a, Choice)
  
  integer                 :: Order
  REALTYPE, dimension(0:) :: a
  integer , OPTIONAL      :: Choice

  integer                 :: Choice_S

  Choice_S = 1 ; if (PRESENT(Choice)) Choice_S = Choice

  SELECT CASE (Choice_S)
    CASE (1)
      call CD_Smooth_Coefficient_0 (Order, a)

    CASE (2)                                          ! Optimized
      call CD_Smooth_Coefficient_1 (Order, a)

  END SELECT

  END Subroutine PS_Central_Smooth_Coefficient 
!
! ================================================================
!
  Subroutine CD_Smooth_Coefficient_0 (Order, a)
  
  integer                 :: Order
  REALTYPE, dimension(0:) :: a

  SELECT CASE (Order)
    CASE (1) 
      a(0) =  1.0d0/2.0d0
      a(1) = -1.0d0/4.0d0

    CASE (2)
      a(0) =  3.0d0/ 8.0d0
      a(1) = -1.0d0/ 4.0d0
      a(2) =  1.0d0/16.0d0

    CASE (3)
      a(0) =   5.0d0/16.0d0
      a(1) = -15.0d0/64.0d0
      a(2) =   3.0d0/32.0d0
      a(3) = - 1.0d0/64.0d0

    CASE (4)
      a(0) =  35.0d0/128.0d0
      a(1) = - 7.0d0/ 32.0d0
      a(2) =   7.0d0/ 64.0d0
      a(3) = - 1.0d0/ 32.0d0
      a(4) =   1.0d0/256.0d0

    CASE (5)
      a(0) =   63.0d0/ 256.0d0
      a(1) = -105.0d0/ 512.0d0
      a(2) =   15.0d0/ 128.0d0
      a(3) = - 45.0d0/1024.0d0
      a(4) =    5.0d0/ 512.0d0
      a(5) = -  1.0d0/1024.0d0

    CASE (6)
      a(0) =  231.0d0/1024.0d0
      a(1) = - 99.0d0/ 512.0d0
      a(2) =  495.0d0/4096.0d0
      a(3) = - 55.0d0/1024.0d0
      a(4) =   33.0d0/2048.0d0
      a(5) = -  3.0d0/1024.0d0
      a(6) =    1.0d0/4096.0d0

    CASE (7)
      a(0) =   429.0d0/ 2048.0d0
      a(1) = -3003.0d0/16384.0d0
      a(2) =  1001.0d0/ 8192.0d0
      a(3) = -1001.0d0/16384.0d0
      a(4) =    91.0d0/ 4096.0d0
      a(5) = -  91.0d0/16384.0d0
      a(6) =     7.0d0/ 8192.0d0
      a(7) = -   1.0d0/16384.0d0

    CASE (8)
      a(0) =  6435.0d0/32768.0d0
      a(1) = - 715.0d0/ 4096.0d0
      a(2) =  1001.0d0/ 8192.0d0
      a(3) = - 273.0d0/ 4096.0d0
      a(4) =   455.0d0/16384.0d0
      a(5) = -  35.0d0/ 4096.0d0
      a(6) =    15.0d0/ 8192.0d0
      a(7) = -   1.0d0/ 4096.0d0
      a(8) =     1.0d0/65536.0d0

    CASE (9)
      a(0) =  12155.0d0/ 65536.0d0 
      a(1) = -21879.0d0/131072.0d0
      a(2) =   1989.0d0/ 16384.0d0
      a(3) = - 4641.0d0/ 65536.0d0
      a(4) =   1071.0d0/ 32768.0d0
      a(5) = -  765.0d0/ 65536.0d0
      a(6) =     51.0d0/ 16384.0d0
      a(7) = -  153.0d0/262144.0d0
      a(8) =      9.0d0/131072.0d0
      a(9) = -    1.0d0/262144.0d0

    CASE (10)
      a(0)  =   46189.0d0/ 262144.0d0
      a(1)  = - 20995.0d0/ 131072.0d0
      a(2)  =   62985.0d0/ 524288.0d0
      a(3)  = -  4845.0d0/  65536.0d0
      a(4)  =    4845.0d0/ 131072.0d0
      a(5)  = -   969.0d0/  65536.0d0
      a(6)  =    4845.0d0/1048576.0d0
      a(7)  = -   285.0d0/ 262144.0d0
      a(8)  =      95.0d0/ 524288.0d0
      a(9)  = -     5.0d0/ 262144.0d0
      a(10) =       1.0d0/1048576.0d0

  END SELECT

  END Subroutine CD_Smooth_Coefficient_0
!
! ================================================================
!
  Subroutine CD_Smooth_Coefficient_1 (Order, a)
  
  integer                 :: Order
  REALTYPE, dimension(0:) :: a

  SELECT CASE (Order)
    CASE (1) 
      a(0) =  1.0d0/2.0d0
      a(1) = -1.0d0/4.0d0

    CASE (2)
      a(0) =  3.0d0/ 8.0d0
      a(1) = -1.0d0/ 4.0d0
      a(2) =  1.0d0/16.0d0

    CASE (3)
      a(0) =   5.0d0/16.0d0
      a(1) = -15.0d0/64.0d0
      a(2) =   3.0d0/32.0d0
      a(3) = - 1.0d0/64.0d0

    CASE (4)
      a(0) =  0.243527493120
      a(1) = -0.204788880640 
      a(2) =  0.120007591680  
      a(3) = -0.045211119360 
      a(4) =  0.008228661760 

    CASE (5)
      a(0) =  0.215044884112  
      a(1) = -0.187772883589 
      a(2) =  0.123755948787 
      a(3) = -0.059227575576 
      a(4) =  0.018721609157 
      a(5) = -0.002999540835 

    CASE (6)
      a(0) =  0.190899511506
      a(1) = -0.171503832236
      a(2) =  0.123632891797
      a(3) = -0.069975429105
      a(4) =  0.029662754736
      a(5) = -0.008520738659
      a(6) =  0.001254597714

  END SELECT

  END Subroutine CD_Smooth_Coefficient_1

END MODULE Central_Difference_Smoothing
