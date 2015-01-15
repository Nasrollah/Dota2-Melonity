#include "pseudopack.h"

  
Module Central_Difference_Coefficient
  
  USE PS_IO_Unit, ONLY: lid6, lid99
  USE Interpolation
  
IMPLICIT NONE
  
PRIVATE
PUBLIC  :: PS_CD_Coefficient
  
CONTAINS
  
  Subroutine PS_CD_Coefficient (Order, N, Degree, dx, W)
  
  integer  :: Order, N, Degree
  REALTYPE :: dx
  
  REALTYPE, dimension(-(Degree+1)/2:(Degree+1)/2,0:Order) :: W
  REALTYPE, dimension(:), ALLOCATABLE                     :: x
  
  integer  :: M_s, N_Degree, i
  REALTYPE :: Scale
  
#if defined (DEBUG)
  call Data_Check
#endif

  M_s = MIN((N-1)/2, (Degree+1)/2) ; N_Degree = 2*M_s
  
  ALLOCATE (x(-M_s:M_s))
  
  do i = -M_s,M_s ; x(i) = i*dx ; enddo
  
  w = HUGE(ONE)
  
  call PS_Lagrangian_Weight (x(0), Order, N_Degree, x, SIZE(w,DIM=1), w)
  
  DEALLOCATE (x)
  
#if 1
  if (M_s >  0) then
      w(0,0) = ONE-( SUM(w(-M_s:-1,0)) + SUM(w(1:M_s,0)) )
    do i = 1,Order
      w(0,i) =    -( SUM(w(-M_s:-1,i)) + SUM(w(1:M_s,i)) )
    enddo
  endif
#endif

#if defined (DEBUG_OUTPUT)
  do i = 0,Order
    write (6,101) i, w(-M_s:M_s,i)*dx**i, SUM(w(-M_s:M_s,i))
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

  END Subroutine PS_CD_Coefficient 

END Module Central_Difference_Coefficient
