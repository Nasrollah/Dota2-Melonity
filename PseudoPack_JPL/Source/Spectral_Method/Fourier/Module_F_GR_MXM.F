#include "pseudopack.h"


MODULE Fourier_GR_MXM
  USE Filtering

implicit NONE

PRIVATE

PUBLIC  :: F_GR_MXM_D
PUBLIC  :: F_GR_MXM_DS

CONTAINS

! --------------------------------------------------------------------
!  (1/N) Sin[N/2 (x-x_j)] / Sin[(x-x_j)/2], x_j=2*pi*j/N, N=2n+1, j=0,..,2n
! --------------------------------------------------------------------
  Subroutine F_GR_MXM_D (N, O_Of_D, D)
 
  integer  :: N, O_Of_D, NN, M, k
  REALTYPE :: PiN2

  REALTYPE, dimension(0:N-1,0:N-1,O_Of_D) :: D
  REALTYPE, dimension(1:N-1)              :: CotX, CscX, CosX
  REALTYPE, dimension(-N+1:N-1)           :: DD

  M  = O_Of_D

  NN   = N*N
  PiN2 = TWO*PI/N

    CosX(1) = -HALF
  do k = 2,N-1
    CosX(k) = -CosX(k-1)
  enddo

  do k = 1,N-1
    CscX(k) = ONE/SIN(HALF*k*PiN2)
  enddo
    
  if (M >  1) then
    do k = 1,N-1
      CotX(k) = ONE/TAN(HALF*k*PiN2)
    enddo
  endif
  
  if (M >= 1) then

      DD(0) = ZERO
    do k = 1,N-1
      DD(k) =  CosX(k)*CscX(k)
    enddo

    call F_GR_MXM_Full (1, N, DD(-N+1), D(0,0,1))
  endif

  if (M >= 2) then

      DD(0) = (1-NN)/(THREE*FOUR)
    do k = 1,N-1
      DD(k) = -CosX(k)*CscX(k)*CotX(k)
    enddo

    call F_GR_MXM_Full (2, N, DD(-N+1), D(0,0,2))
  endif
 
  if (M >= 3) then

      DD(0) = ZERO
    do k = 1,N-1
      DD(k) =  CosX(k)*CscX(k)*QUARTER*(THREE*CotX(k)**2+THREE*CscX(k)**2-NN)
    enddo

    call F_GR_MXM_Full (3, N, DD(-N+1), D(0,0,3))
  endif
 
  if (M >= 4) then

      DD(0) = (3*NN*NN-10*NN+7)/(24*TEN)
    do k = 1,N-1
      DD(k) = -CosX(k)*CscX(k)*CotX(k)*HALF*(CotX(k)**2+ FIVE*CscX(k)**2-NN)
    enddo

    call F_GR_MXM_Full (4, N, DD(-N+1), D(0,0,4))
  endif

  END Subroutine F_GR_MXM_D 
!
! ---------------------------------------------------------------------
!
  Subroutine F_GR_MXM_DS (O_Of_D, Filter_Choice, N, Mode_CutOff, &
                                  Omega, Order, D)
 
  integer  :: O_Of_D, Filter_Choice, N, Mode_CutOff
  integer  :: N2, M, k, l
  REALTYPE :: XN2, PiN2, x_k, D_k, Omega, Order

  REALTYPE, dimension(0:N-1,0:N-1,O_Of_D) :: D
  REALTYPE, dimension(0:(N-1)/2)          :: Sigma
  REALTYPE, dimension(-N+1:N-1)           :: DD

  M = O_Of_D

  N2   = (N-1)/2
  XN2  = TWO/N
  PiN2 = XN2*PI

  call Filter (Filter_Choice, N2, Mode_CutOff, Omega, Order, Sigma)

  if (M >= 1) then 

      DD(0) = ZERO
    do k = 0,N-1
      x_k = k*PiN2
      D_k = ZERO
      do l = 1,(N-1)/2
        D_k = D_k + l   *Sigma(l)*SIN(l*x_k)
      enddo
      DD(k) = -D_k*XN2
    enddo

    call F_GR_MXM_Full (1, N, DD(-N+1), D(0,0,1))
  endif
 
  if (M >= 2) then

    DD(0) = ZERO
    do k = 0,N-1
      x_k = k*PiN2
      D_k = ZERO
      do l = 1,(N-1)/2
        D_k = D_k + l**2*Sigma(l)*COS(l*x_k)
      enddo
      DD(k) = -D_k*XN2
    enddo

    call F_GR_MXM_Full (2, N, DD(-N+1), D(0,0,2))
  endif

  if (M >= 3) then

      DD(0) = ZERO
    do k = 0,N-1
      x_k = k*PiN2
      D_k = ZERO
      do l = 1,(N-1)/2
        D_k = D_k + l**3*Sigma(l)*SIN(l*x_k)
      enddo
      DD(k) =  D_k*XN2
    enddo

    call F_GR_MXM_Full (3, N, DD(-N+1), D(0,0,3))
  endif

  if (M >= 4) then

      DD(0) = ZERO
    do k = 0,N-1
      x_k = k*PiN2
      D_k = ZERO
      do l = 1,(N-1)/2
        D_k = D_k + l**4*Sigma(l)*COS(l*x_k)
      enddo
      DD(k) =  D_k*XN2
    enddo

    call F_GR_MXM_Full (4, N, DD(-N+1), D(0,0,4))
  endif

  END Subroutine F_GR_MXM_DS 
!
! ---------------------------------------------------------------------
!
  Subroutine F_GR_MXM_Full (O_Of_D, N, DD, D)
 
  integer  :: O_Of_D, N, k, j
  REALTYPE :: Factor

  REALTYPE, dimension(-N+1:N-1)    :: DD
  REALTYPE, dimension(0:N-1,0:N-1) :: D

  Factor = (-1)**O_Of_D

  do k = 1,N-1
    DD( -k) = Factor*DD(k)
  enddo

  do j = 0,N-1
    do k = 0,N-1
      D(k,j) = DD(k-j)
    enddo
  enddo

  END Subroutine F_GR_MXM_Full 

END MODULE Fourier_GR_MXM
