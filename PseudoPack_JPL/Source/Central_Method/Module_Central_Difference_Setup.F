#include "pseudopack.h"


MODULE Central_Difference_Setup

  USE PS_IO_Unit, ONLY: lid6, lid99
  USE Interpolation

IMPLICIT NONE

INTERFACE PS_CD_Operator_Setup
  MODULE PROCEDURE Setup_CD_Operator
END INTERFACE

PRIVATE

PUBLIC  :: PS_CD_Operator_Setup
PUBLIC  :: PS_Central_Domain_Index
PUBLIC  :: PS_Grid_Index_Update

CONTAINS

  Subroutine Setup_CD_Operator (Degree, W, dx, Order_Of_D, N_Ghost_Cell)

  integer  , parameter :: Max_Order_Of_D = 4

  integer              :: N, Degree
  REALTYPE , OPTIONAL  :: dx
  integer  , OPTIONAL  :: Order_Of_D
  integer  , OPTIONAL  :: N_Ghost_Cell

  REALTYPE, dimension(-20:20,-1:4)              :: W
  REALTYPE, dimension(:)          , ALLOCATABLE :: x

  integer  :: Order, M, i
  REALTYPE :: ddx

  w = HUGE(ONE)

  M = Degree/2 ; N = 2*M+1

  Order = Max_Order_Of_D ; if (PRESENT(Order_Of_D)) Order = Order_Of_D
  ddx   = ONE            ; if (PRESENT(dx))         ddx   = dx

  W(0,-1) = M               ! Number of Ghost points or (W(1,-1)-1)/2
  W(1,-1) = Degree          ! Desired Degree of polynomial interoplation
  W(2,-1) = N               ! Desired number of collocation points used
  W(3,-1) = Order           ! Max order of differentiations
  W(4,-1) = ddx             ! dx

#if defined (DEBUG)
  call Data_Check
#endif

  ALLOCATE (x(-M:M))

  do i = -M,M ; x(i) = i*ddx ; enddo

  call PS_Lagrangian_Weight (x(0), Order, Degree, x, N, w(-M:M,0:Order))

  DEALLOCATE (x)

  if (PRESENT(N_Ghost_Cell)) N_Ghost_Cell = M

  if (M >  0) then
      w(0,0) = ONE-( SUM(w(-M:-1,0)) + SUM(w(1:M,0)) )
    do i = 1,Order
      w(0,i) =    -( SUM(w(-M:-1,i)) + SUM(w(1:M,i)) )
    enddo
  endif

#if defined (DEBUG_OUTPUT)
  do i = 0,Order
    write (6,101) i, w(-M:M,i)*dx**i, SUM(w(-M:M,i))
  enddo
  101 format (i5, 10(1x,g12.5,:))
#endif

#if defined (DEBUG)
  CONTAINS
!
! =================================================================
!
    Subroutine Data_Check

    if (Degree >  N-1) then
      write (lid6 ,*) 'STOP ! PS_CD_Coefficient'
      write (lid6 ,*) '       Degree=', Degree, ' is greater than N = ',N
      write (lid99,*) 'STOP ! PS_CD_Coefficient'
      write (lid99,*) '       Degree=', Degree, ' is greater than N = ',N

      STOP
    endif

    END Subroutine Data_Check
#endif

  END Subroutine Setup_CD_Operator 
!
! ==================================================================
!
  Subroutine PS_Central_Domain_Index (Method, Symmetry, Order, N_Pnts,   &
                                      N, N0, N1, N2, N3, N4, N5          )

  integer  :: Method, Symmetry, Order
  integer  :: N_Pnts, N, N0, N1, N2, N3, N4, N5, M_s

  M_s = Order/2

  N = N_Pnts-2 ; if (Symmetry == 1) N = (N-1)/2

  N0 = -M_s ; N1 = -1 ; N2 = 0 ; N3 = N ; N4 = N3+1 ; N5 = N3+M_s

  END Subroutine PS_Central_Domain_Index
!
! ==================================================================
!
  Subroutine PS_Grid_Index_Update (WENO_Order, Central_Order, MR_Order,    &
                                   N0, N2, N3, N5,                         &
                                   Symmetry, N_Ghost)

  integer                            :: WENO_Order, Central_Order, MR_Order
  integer                            :: N0, N2, N3, N5
  integer ,                 OPTIONAL :: Symmetry
  integer , dimension(0:1), OPTIONAL :: N_Ghost
  integer                            :: N_Ghost_Point

  if (MOD(MR_Order,2) == 0) then
    N_Ghost_Point = MAX(MR_Order  , (WENO_Order+1)/2, Central_Order/2)
  else
    N_Ghost_Point = MAX(MR_Order+1, (WENO_Order+1)/2, Central_Order/2)
  endif

  N0 = N2 - N_Ghost_Point
  N5 = N3 + N_Ghost_Point

  if (MOD(N3-N2+1,2) == 1) then
    N0 = N0-1 ; N5 = N5+1

    if (PRESENT(Symmetry)) then
      if (Symmetry /= 1) N0 = N0+1
    endif
  endif

  if (PRESENT(N_Ghost)) then
    N_Ghost(0) = N2-N0 ; N_Ghost(1) = N5-N3
  endif

  END Subroutine PS_Grid_Index_Update 

END MODULE Central_Difference_Setup
