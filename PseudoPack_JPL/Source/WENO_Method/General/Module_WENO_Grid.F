#include "pseudopack.h"

! ----------------------------------------------------------------------
! FUNCTION NAME: WENO_Uniform_Grid_1D
! AUTHOR       : Wai-Sun Don
! DESCRIPTION  : Setup 1D Uniform Grid for the WENO Methods.
!               
! N          Number of collocation Grid for x(-M:N+M) (Input)
! Order      Order of the WENO scheme (3, 5, 7, 9, 11)    (Input)
! x0         The left domain  (Input)
! x1         The right domain (Input)
! dx         Grid spacing  (Input)
! x          Array of size (-M:N+M) containing the collocation grid (Output)
!            where M =(Order+1)/2
! Direction   1 if Grid index 0 and N corresponding with x0 and x1 respectively.
!            -1 if Grid index 0 and N corresponding with x1 and x0 respectively
!                  and dx = -dx.
! ISign       Shifting Factor for the grid spacing  (Input)
!                x_i = x0 + (i-ISign/2)*dx  if (Direction =  1)
!             or x_i = x1 - (i+ISign/2)*dx  if (Direction = -1)
!*********************************************************************


MODULE WENO_Grid

implicit NONE

INTERFACE PS_WENO_Uniform_Grid
  MODULE PROCEDURE WENO_Uniform_Grid_1D
END INTERFACE

INTERFACE PS_WENO_Domain_Index  
  MODULE PROCEDURE WENO_Domain_Index_1D
END INTERFACE

PRIVATE

PUBLIC  :: PS_WENO_Uniform_Grid
PUBLIC  :: PS_WENO_Domain_Index  
PUBLIC  :: WENO_Uniform_Grid_1D
PUBLIC  :: WENO_Domain_Index_1D

CONTAINS

  Subroutine WENO_Uniform_Grid_1D (N0, N2, N3, N5, &
                                   x0, x1, dx, x, Symmetry, Direction, ISign)

  integer  :: N0, N2, N3, N5, i
  REALTYPE :: x0, x1, x2, dx

  REALTYPE, dimension(N0:N5) :: x

  integer , OPTIONAL :: Symmetry
  integer , OPTIONAL :: Direction
  REALTYPE, OPTIONAL :: ISign

  integer  :: Symmetric
  integer  :: Wind
  REALTYPE :: Sign

  Symmetric = 0    ; if (PRESENT(Symmetry) ) Symmetric = Symmetry
  Wind      = 1    ; if (PRESENT(Direction)) Wind      = Direction
  Sign      = -ONE ; if (PRESENT(ISign)    ) Sign      = ISign

  if (Symmetric == 1) x1 = HALF*(x0+x1)

  dx = Wind*ABS(x1-x0)/(N3-N2+1)

  SELECT CASE (Wind)
    CASE DEFAULT
      do i = N0,N5
        x(i) = x0 + (i-Sign*HALF)*dx
      enddo

    CASE (-1)
      do i = N0,N5
        x(i) = x1 + (i+Sign*HALF)*dx
      enddo

  END SELECT

  END Subroutine WENO_Uniform_Grid_1D
!
!======================================================================
!
  Subroutine WENO_Domain_Index_1D (Method, Symmetry, Order, N_Pnts, &
                                   N, N0, N1, N2, N3, N4, N5, Tight_Bound)

  integer  :: Method, Symmetry, Order, N_Pnts, N, N0, N1, N2, N3, N4, N5, M_s

  integer , OPTIONAL :: Tight_Bound

  M_s = (Order+1)/2
  N   = N_Pnts-2 ; if (Symmetry == 1) N = (N-1)/2

  N0 = -M_s ; N1 = -1 ; N2 = 0 ; N3 = N ; N4 = N3+1 ; N5 = N3+M_s

  if (.NOT. PRESENT(Tight_Bound)) RETURN

  if (Tight_Bound == 1) N5 = N5-1

  END Subroutine WENO_Domain_Index_1D 

END MODULE WENO_Grid
