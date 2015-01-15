#include "pseudopack.h"


MODULE Legendre_Grid_Point
  USE Gauss_Quadrature

implicit NONE

PRIVATE

PUBLIC  :: L_GL_Point_Full
PUBLIC  :: L_GL_Point_Half
PUBLIC  :: L_GR_Point_Full
PUBLIC  :: L_Ga_Point_Full
PUBLIC  :: L_Ga_Point_Half

CONTAINS

! ----------------------------------------------------------------------
! FUNCTION NAME: L_GL_Point_Full
! DESCRIPTION  : computes the Un-Mapped Legendre-Guass_Lobatto points
!                without assumption of symmetry.
! ----------------------------------------------------------------------
  Subroutine L_GL_Point_Full (N, x)

  integer                    :: N, Nfix, IERR
  REALTYPE, dimension(2)     :: qxfix

  REALTYPE, dimension(0:N)   :: x
  REALTYPE, dimension(0:N)   :: w, Work 

#ifdef DEBUG
  ASSERT(N>0, "L_GL_Point_Full: N<=0")
#endif

  Nfix = 2 ; qxfix(1) = -ONE ; qxfix(2) = ONE

  call gaussq (1, N+1, ZERO, ZERO, Nfix, qxfix, work, x, w, IERR)

  if (IERR /= 0) then
    write (6 ,*) "L_GL_Point_Full: Gauss Labatto Quadrature failed. Quitting."
    write (99,*) "L_GL_Point_Full: Gauss Labatto Quadrature failed. Quitting."
    STOP
  endif

  x = -x

                     x(0  ) =  ONE
  if (MOD(N,2) == 0) x(N/2) = ZERO
                     x(N  ) = -ONE

  END Subroutine L_GL_Point_Full 

! ----------------------------------------------------------------------
! FUNCTION NAME: L_GL_Point_Half
! DESCRIPTION  : computes the Un-Mapped Legendre-Guass_Lobatto points
!                with assumption of symmetry.
! ----------------------------------------------------------------------
  Subroutine L_GL_Point_Half (N, x)

  integer                    :: N, Nfix, IERR
  REALTYPE, dimension(2)     :: qxfix

  REALTYPE, dimension(0:N/2) :: x
  REALTYPE, dimension(0:N)   :: w, Work, y

#ifdef DEBUG
  ASSERT(N>0, "L_GL_Point_Half: N<=0")
#endif

  Nfix = 2 ; qxfix(1) = -ONE ; qxfix(2) = ONE

  call gaussq (1, N+1, ZERO, ZERO, Nfix, qxfix, work, y, w, IERR)

  if (IERR /= 0) then
    write (6 ,*) "L_GL_Point_Half: Gauss Labatto Quadrature failed. Quitting."
    write (99,*) "L_GL_Point_Half: Gauss Labatto Quadrature failed. Quitting."
    Stop
  endif

  x = -y(0:N/2)

                     x(0  ) =  ONE
  if (MOD(N,2) == 0) x(N/2) = ZERO

  END Subroutine L_GL_Point_Half 

! ----------------------------------------------------------------------
! FUNCTION NAME: L_GR_Point_Full
! DESCRIPTION  : computes the Un-Mapped Legendre-Guass_Radau points
! ----------------------------------------------------------------------
  Subroutine L_GR_Point_Full (N, x)

  integer                    :: N, Nfix, IERR
  REALTYPE, dimension(2)     :: qxfix

  REALTYPE, dimension(0:N)   :: x
  REALTYPE, dimension(0:N)   :: w, Work

#ifdef DEBUG
  ASSERT(N>0, "L_GR_Point_Full: N<=0")
#endif

  Nfix = 1 ; qxfix(1) = -ONE 

  call gaussq (1, N+1, ZERO, ZERO, Nfix, qxfix, work, x, w, IERR)

  if (IERR /= 0) then
    write (6 ,*) "L_GR_Point_Full: Gauss Radau Quadrature failed. Quitting."
    write (99,*) "L_GR_Point_Full: Gauss Radau Quadrature failed. Quitting."
    STOP
  endif

  x = -x 

  x(0) = ONE

  END Subroutine L_GR_Point_Full 

! ----------------------------------------------------------------------
! FUNCTION NAME: L_Ga_Point_Full
! DESCRIPTION  : computes the Un-Mapped Legendre-Guass points
!                without assumption of symmetry.
! ----------------------------------------------------------------------
  Subroutine L_Ga_Point_Full (N, x)

  integer                    :: N, Nfix, IERR
  REALTYPE, dimension(2)     :: qxfix

  REALTYPE, dimension(0:N)   :: x
  REALTYPE, dimension(0:N)   :: w, Work

#ifdef DEBUG
  ASSERT(N>0, "L_Ga_Point_Full: N<=0")
#endif

  Nfix =  0

  call gaussq (1, N+1, ZERO, ZERO, Nfix, qxfix, work, x, w, IERR)

  if (IERR /= 0) then
    write (6 ,*) "L_Ga_Point_Full: Gauss Quadrature failed. Quitting."
    write (99,*) "L_Ga_Point_Full: Gauss Quadrature failed. Quitting."
    STOP
  endif

  x = -x 

  if (MOD(N,2) == 0) x(N/2) = ZERO

  END Subroutine L_Ga_Point_Full 

! ----------------------------------------------------------------------
! FUNCTION NAME: L_Ga_Point_Half
! DESCRIPTION  : computes the Un-Mapped Legendre-Guass points
!                with assumption of symmetry.
! ----------------------------------------------------------------------
  Subroutine L_Ga_Point_Half (N, x)

  integer                    :: N, Nfix, IERR
  REALTYPE, dimension(2)     :: qxfix

  REALTYPE, dimension(0:N/2) :: x
  REALTYPE, dimension(0:N)   :: w, Work, y

#ifdef DEBUG
  ASSERT(N>0, "L_Ga_Point_Half: N<=0")
#endif

  Nfix =  0

  call gaussq (1, N+1, ZERO, ZERO, Nfix, qxfix, work, y, w, IERR)

  if (IERR /= 0) then
    write (6 ,*) "L_Ga_Point_Half: Gauss Quadrature failed. Quitting."
    write (99,*) "L_Ga_Point_Half: Gauss Quadrature failed. Quitting."
    STOP
  endif

  x = -y(0:N/2)

  if (MOD(N,2) == 0) x(N/2) = ZERO

  END Subroutine L_Ga_Point_Half 

END MODULE Legendre_Grid_Point
