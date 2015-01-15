#include "pseudopack.h"


MODULE Chebyshev_Grid_Point

implicit NONE

PRIVATE

PUBLIC  :: C_GL_Point_Full
PUBLIC  :: C_GL_Point_Half
PUBLIC  :: C_GR_Point_Full
PUBLIC  :: C_Ga_Point_Full
PUBLIC  :: C_Ga_Point_Half

CONTAINS

! ----------------------------------------------------------------------
! FUNCTION NAME: C_GL_Point_Full
! DESCRIPTION  : computes the Un-Mapped Chebyshev-Guass_Lobatto points
!                without assumption of symmetry.
! ----------------------------------------------------------------------
  Subroutine C_GL_Point_Full (N, x)

  integer  :: N, N2, i
  REALTYPE :: PiN

  REALTYPE, dimension(0:N) :: x

#ifdef DEBUG
  ASSERT(N>0, "C_GL_Point_Full: N<=0")
#endif

   N2 = N/2
  PiN = PI/N

    x(0) = ONE
  do i = 1,N2
    x(i) = COS(PiN*i)
  enddo

  if (MOD(N,2) == 0) then
    do i = 0,N2-1
      x(N-i) = -x(i)
    enddo
    x(N2) = ZERO
  else
    do i = 0,N2
      x(N-i) = -x(i)
    enddo
  endif

  do i = 0,N
    if (i*3 == N) then
      x(  i) =  HALF
      x(N-i) = -HALF
      RETURN
    endif
  enddo

  END Subroutine C_GL_Point_Full 

! ----------------------------------------------------------------------
! FUNCTION NAME: C_GL_Point_Half
! DESCRIPTION  : computes the Un-Mapped Chebyshev-Guass_Lobatto points
!                with assumption of symmetry.
! ----------------------------------------------------------------------
  Subroutine C_GL_Point_Half (N, x)

  integer  :: N, N2, i
  REALTYPE :: PiN

  REALTYPE, dimension(0:N/2) :: x

#ifdef DEBUG
  ASSERT(N>0, "C_GL_Point_Half: N<=0")
#endif

   N2 = N/2
  PiN = PI/N

    x(0) = ONE
  do i = 1,N2
    x(i) = COS(PiN*i)
  enddo

  if (MOD(N,2) == 0) x(N2) = ZERO

  do i = 0,N2
    if (i*3 == N) then
      x(i) = HALF
      RETURN
    endif
  enddo

  END Subroutine C_GL_Point_Half 

! ----------------------------------------------------------------------
! FUNCTION NAME: C_GR_Point_Full
! DESCRIPTION  : computes the Un-Mapped Chebyshev-Guass_Radau points
! ----------------------------------------------------------------------
  Subroutine C_GR_Point_Full (N, x)

  integer  :: N, i
  REALTYPE :: PiN

  REALTYPE, dimension(0:N) :: x

#ifdef DEBUG
  ASSERT(N>0, "C_GR_Point_Full: N<=0")
#endif

  PiN = TWO*PI/(2*N+1)

  x(0)   = ONE
  do i = 1,N
    x(i) = COS(PiN*i)
  enddo

  END Subroutine C_GR_Point_Full 

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Ga_Point_Full
! DESCRIPTION  : computes the Un-Mapped Chebyshev-Guass points
!                without assumption of symmetry.
! ----------------------------------------------------------------------
  Subroutine C_Ga_Point_Full (N, x)

  integer  :: N, N2, i
  REALTYPE :: PiN

  REALTYPE, dimension(0:N) :: x

#ifdef DEBUG
  ASSERT(N>0, "C_Ga_Point_Full: N<=0")
#endif

  N2  = N/2
  PiN = PI/(2*N+2)

  do i = 0,N2
    x(i) = COS(PiN*(2*i+1))
  enddo

  if (MOD(N,2) == 0) then
    do i = 0,N2-1
      x(N-i) = -x(i)
    enddo

    x(N2) = ZERO
  else
    do i = 0,N2
      x(N-i) = -x(i)
    enddo
  endif

  END Subroutine C_Ga_Point_Full 

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Ga_Point_Half
! DESCRIPTION  : computes the Un-Mapped Chebyshev-Guass points
!                with assumption of symmetry.
! ----------------------------------------------------------------------
  Subroutine C_Ga_Point_Half (N, x)

  integer  :: N, N2, i
  REALTYPE :: PiN

  REALTYPE, dimension(0:N/2) :: x

#ifdef DEBUG
  ASSERT(N>0, "C_Ga_Point_Half: N<=0")
#endif

  N2  = N/2
  PiN = PI/(2*N+2)

  do i = 0,N2
    x(i) = COS(PiN*(2*i+1))
  enddo

  if (MOD(N,2) == 0) x(N2) = ZERO

  END Subroutine C_Ga_Point_Half 

END MODULE Chebyshev_Grid_Point
