#include "pseudopack.h"


MODULE Fourier_Grid_Point

implicit NONE

PRIVATE

PUBLIC  :: F_GL_Point_Full
PUBLIC  :: F_GL_Point_Half
PUBLIC  :: F_GR_Point_Full

CONTAINS

! ----------------------------------------------------------------------
! FUNCTION NAME: F_GL_Point_Full
! DESCRIPTION  : Fourier collocation point for Full Domain
! ---------------------------------------------------------------------
  Subroutine F_GL_Point_Full (N, x)
 
  integer  :: N, i
  REALTYPE :: Angle

  REALTYPE, dimension(0:N-1) :: x

  Angle = TWO*PI/N

  do i = 0,N-1
    x(i) = Angle*i
  enddo
 
  END Subroutine F_GL_Point_Full 

! ----------------------------------------------------------------------
! FUNCTION NAME: F_GL_Point_Half
! DESCRIPTION  : Fourier collocation point for Half Domain
! ---------------------------------------------------------------------
  Subroutine F_GL_Point_Half (N, x)
 
  integer  :: N, i
  REALTYPE :: Angle

  REALTYPE, dimension(0:N/2) :: x
 
  Angle = TWO*PI/N

  do i = 0,N/2
    x(i) = Angle*i
  enddo
 
  END Subroutine F_GL_Point_Half 

! ----------------------------------------------------------------------
! FUNCTION NAME: F_GR_Point_Full
! DESCRIPTION  : Fourier collocation point for Full Domain
! ---------------------------------------------------------------------
  Subroutine F_GR_Point_Full (N, x)

  integer  :: N, i
  REALTYPE :: Angle

  REALTYPE, dimension(0:N-1) :: x

  Angle = TWO*PI/N

  do i = 0,N-1
    x(i) = Angle*i
  enddo

  END Subroutine F_GR_Point_Full 

END MODULE Fourier_Grid_Point
