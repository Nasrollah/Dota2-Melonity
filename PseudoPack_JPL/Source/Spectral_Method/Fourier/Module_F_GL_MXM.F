#include "pseudopack.h"


MODULE Fourier_GL_MXM
  USE Filtering

implicit NONE

PRIVATE

PUBLIC  :: F_GL_MXM_D
PUBLIC  :: F_GL_MXM_DS

CONTAINS

! --------------------------------------------------------------------
!  (1/N) Sin[N/2 (x-x_j)] Cot[(x-x_j)/2], x_j=2*pi*j/N, N=2n, j=0,..,2n-1
! --------------------------------------------------------------------
  Subroutine F_GL_MXM_D (N, O_Of_D, D)
 
  integer  :: N, O_Of_D, N2, NN, M, k
  REALTYPE :: PiN2

  REALTYPE, dimension(0:N-1,0:N-1,O_Of_D) :: D
  REALTYPE, dimension(1:N/2)              :: CotX, CscX, CosX
  REALTYPE, dimension(-N+1:N-1)           :: DD

  M  = O_Of_D

  N2   = N/2
  NN   = N*N
  PiN2 = TWO*PI/N

    CosX(1) = -HALF
  do k = 2,N2
    CosX(k) = -CosX(k-1)
  enddo

  do k = 1,N2-1
    CotX(k) = ONE/TAN(HALF*k*PiN2)
  enddo
    CotX(N2) = ZERO

  if (M >  1) then
    do k = 1,N2
      CscX(k) = ONE/SIN(HALF*k*PiN2)
    enddo
  endif
    
  if (M >= 1) then

      DD(0) = ZERO
    do k = 1,N2
      DD(k) =  CosX(k)*CotX(k)
    enddo

    call F_GL_MXM_Full (1, N, DD(-N+1), D(0,0,1))
  endif

  if (M >= 2) then

      DD(0) = -(NN+2)/(6*TWO)
    do k = 1,N2
      DD(k) = -CosX(k)*CscX(k)*CscX(k)
    enddo

    call F_GL_MXM_Full (2, N, DD(-N+1), D(0,0,2))
  endif
 
  if (M >= 3) then

      DD(0) = ZERO
    do k = 1,N2
      DD(k) =  CosX(k)*CotX(k)*QUARTER*(SIX*CscX(k)**2-NN)
    enddo

    call F_GL_MXM_Full (3, N, DD(-N+1), D(0,0,3))
  endif
 
  if (M >= 4) then

      DD(0) = (3*NN*NN+20*NN-8)/(24*TEN)
    do k = 1,N2
      DD(k) = -CosX(k)*CscX(k)*CscX(k)*HALF*(FOUR*CotX(k)**2+TWO*CscX(k)**2-NN)
    enddo

    call F_GL_MXM_Full (4, N, DD(-N+1), D(0,0,4))
  endif

  END Subroutine F_GL_MXM_D 
!
! ---------------------------------------------------------------------
!
  Subroutine F_GL_MXM_DS (O_Of_D, Filter_Choice, N, Mode_CutOff, &
                                  Omega, Order, D)
 
  integer  :: O_Of_D, Filter_Choice, N, Mode_CutOff
  integer  :: N2, M, k, l
  REALTYPE :: x, XN2, PiN2, D_k, Omega, Order

  REALTYPE, dimension(0:N-1,0:N-1,O_Of_D) :: D
  REALTYPE, dimension(0:N/2)              :: Sigma
  REALTYPE, dimension(-N+1:N-1)           :: DD

  M = O_Of_D

  N2   = N/2
  XN2  = TWO/N
  PiN2 = XN2*PI

  call Filter (Filter_Choice, N2, Mode_CutOff, Omega, Order, Sigma)

  if (M >= 1) then

    DD(0) = ZERO
    do k = 0,N2
      x   = k*PiN2
      D_k = ZERO
      do l = 1,N2
        D_k = D_k + l   *Sigma(l)*SIN(l*x)
      enddo
      DD(k) = -D_k*XN2
    enddo

    call F_GL_MXM_Full (1, N, DD(-N+1), D(0,0,1))
  endif
 
  if (M >= 2) then

    DD(0) = ZERO
    do k = 0,N2
      x   = k*PiN2
      D_k = ZERO
      do l = 1,N2
        D_k = D_k + l**2*Sigma(l)*COS(l*x)
      enddo
      DD(k) = -D_k*XN2
    enddo

    call F_GL_MXM_Full (2, N, DD(-N+1), D(0,0,2))
  endif

  if (M >= 3) then

    DD(0) = ZERO
    do k = 0,N2
      x   = k*PiN2
      D_k = ZERO
      do l = 1,N2
        D_k = D_k + l**3*Sigma(l)*SIN(l*x)
      enddo
      DD(k) =  D_k*XN2
    enddo

    call F_GL_MXM_Full (3, N, DD(-N+1), D(0,0,3))
  endif

  if (M >= 4) then

    DD(0) = ZERO
    do k = 0,N2
      x   = k*PiN2
      D_k = ZERO
      do l = 1,N2
        D_k = D_k + l**4*Sigma(l)*COS(l*x)
      enddo
      DD(k) =  D_k*XN2
    enddo

    call F_GL_MXM_Full (4, N, DD(-N+1), D(0,0,4))
  endif

  END Subroutine F_GL_MXM_DS 
!
! ---------------------------------------------------------------------
!
  Subroutine F_GL_MXM_Full (O_Of_D, N, DD, D)
 
  integer  :: O_Of_D, N, k, j
  REALTYPE :: Factor

  REALTYPE, dimension(-N+1:N-1)    :: DD
  REALTYPE, dimension(0:N-1,0:N-1) :: D

  Factor = (-1)**O_Of_D

  do k = 1,N/2-1
    DD(N-k) = Factor*DD(k)
  enddo

  do k = 1,N-1
    DD( -k) = Factor*DD(k)
  enddo
     
  do j = 0,N-1
    do k = 0,N-1
      D(k,j) = DD(k-j)
    enddo
  enddo
   
  END Subroutine F_GL_MXM_Full 

END MODULE Fourier_GL_MXM
