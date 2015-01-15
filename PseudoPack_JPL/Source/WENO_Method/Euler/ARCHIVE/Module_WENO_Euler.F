#include "pseudopack.h"

! ----------------------------------------------------------------------
! FUNCTION NAME: PS_WENO_Euler
! AUTHOR       : Wai-Sun Don
! DESCRIPTION  : WENO Flux Reconstruction based on the Global-Lax-Friedrichs
!                Flux Splitting for the Euler equations.
!                The conservative variables Q is an 2D array of size
!                (N0:N5,NV)
!
! Order      Order of the WENO scheme (3, 5)    (Input)
! N0         The Lower bound of the conservative variable Q (Input)
! N5         The Upper bound of the conservative variable Q (Input)
! i0         The starting index of the computed WENO Flux e.g. (-1)
!               (Input)
! i1         The ending   index of the computed WENO Flux e.g. (N3-N2)
!               (Input)
! NV         The Number of Conservative variable Q (Input)
! Q          2D array of size (N0:N5,NV) of the
!              conservative variables. (Input)
! F          2D array of size (N0:N5,NV) of the Euler Flux.
!              (Input)
! lambda     1D array of size (NV) containing the
!              maximum value of 'm' eigenvalues.  The eigenvalues
!              of the Euler equations is u-c, u, u, u+c, where
!              c is the local sound speed. (Input)
! E_L        3D array of size (N0:N5,NV,NV)
!              containing the LEFT  Eigenvectors of the Euler
!              equations. (Input)
! E_R        3D array of size (N0:N5,NV,NV)
!              containing the RIGHT Eigenvectors of the Euler
!              equations. (Input)
! F_Half     2D array of size (N0:N5,NV) containing
!              the Reconstructed Flux at the cell center. (Output)
!*********************************************************************


MODULE WENO_Euler

  USE Euler_Projection
  USE Lax_Friedrichs_Flux_Splitting
  USE WENO_Euler_Reconstruction
  USE WENO_Option

implicit NONE

INTERFACE PS_WENO_Euler
  MODULE PROCEDURE WENO_Euler_Flux
END INTERFACE

PRIVATE
PUBLIC  :: PS_WENO_Euler

CONTAINS
!
!===================================================================
!
  Subroutine WENO_Euler_Flux (Order, N0, N5, i0, i1,                   &
                              NV, Q, F, lambda_G,                      &
                              E_L, E_R, F_Half)

  integer  :: Order, N0, N5, i0, i1, NV

  REALTYPE, dimension(N0:N5,NV)    :: Q
  REALTYPE, dimension(N0:N5,NV)    :: F, F_Half
  REALTYPE, dimension(N0:N5,NV,NV) :: E_L, E_R
  REALTYPE, dimension(      NV)    :: Lambda_G

  REALTYPE, dimension(N0:N5,-(Order-1)/2:(Order-1)/2,2) :: h
  REALTYPE, dimension(N0:N5,NV                      ,2) :: g
  REALTYPE, dimension(N0:N5,NV)                         :: Fm

  integer  :: m 

  do m = 1,NV
    call PS_Global_Lax_Friedrichs     (       N0,N5, NV, Lambda_G(m), F, Q, g)

    call PS_Euler_Left_Projection     (Order, N0,N5, i0,i1, NV, m, E_L, g, h) 


    call PS_WENO_Euler_Reconstruction (Order, N0,N5, i0,i1, h, Fm(:,m))
  enddo

  call PS_Euler_Right_Projection      (       N0,N5, i0,i1, NV, E_R, Fm, F_Half)

  END Subroutine WENO_Euler_Flux

END MODULE WENO_Euler
