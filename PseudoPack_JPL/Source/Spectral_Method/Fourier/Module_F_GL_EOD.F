#include "pseudopack.h"


MODULE Fourier_GL_EOD
  USE Filtering

implicit NONE

PRIVATE

PUBLIC  :: F_GL_EOD_D
PUBLIC  :: F_GL_EOD_DS

CONTAINS

  Subroutine F_GL_EOD_D  (N, O_OF_D, D_Even, D_Odd)

  integer  :: N, O_Of_D
  integer  :: N2, NN, M, j, k
  REALTYPE :: PiN2, Factor

  REALTYPE, dimension(0:N/2,0:N/2,O_Of_D)  :: D_Even, D_Odd
  REALTYPE, dimension(1:N/2)               :: CotX, CscX, CosX
  REALTYPE, dimension(-N:N,MAX_DERIVATIVE) :: D

  M  = O_Of_D

  N2 = N/2
  NN = N*N

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
      D(0,1) = ZERO
    do k = 1,N2
      D(k,1) =  CosX(k)*CotX(k)
    enddo
  endif

  if (M >= 2) then
      D(0,2) = -(NN+2)/(6*TWO)
    do k = 1,N2
      D(k,2) = -CosX(k)*CscX(k)*CscX(k)
    enddo
  endif

  if (M >= 3) then
      D(0,3) = ZERO
    do k = 1,N2
      D(k,3) =  CosX(k)*CotX(k)*QUARTER*(SIX*CscX(k)**2-NN)
    enddo
  endif

  if (M >= 4) then
      D(0,4) = (3*NN*NN+20*NN-8)/(24*TEN)
    do k = 1,N2
      D(k,4) = -CosX(k)*CscX(k)*CscX(k)*HALF*(FOUR*CotX(k)**2+TWO*CscX(k)**2-NN)
    enddo
  endif

  do M = 1,O_Of_D
    Factor = (-1)**M

    do k = 1,N2-1
      D(N-k,M) = Factor*D(k,M)
    enddo

    do k = 1,N-1
      D( -k,M) = Factor*D(k,M)
    enddo
  enddo
     
  do M = 1, O_Of_D
    D_Even(:,0,M) = D(0:N2,M)
    D_Odd (:,0,M) = ZERO

    do j = 1,N2
      do k = 0,N2
        D_Even(k,j,M) = D(k-j,M)+D(k+j-N,M)
        D_Odd (k,j,M) = D(k-j,M)-D(k+j-N,M)
      enddo
    enddo

    if (Mod(N,2) == 0) then
      D_Even(:,N2,M) = HALF*D_Even(:,N2,M)
      D_Odd (:,N2,M) = ZERO
    endif
  enddo
 
  END Subroutine F_GL_EOD_D  
!
! ----------------------------------------------------------------------
!
  Subroutine F_GL_EOD_DS (O_Of_D, Filter_Choice, N, Mode_CutOff, Omega, Order, &
                                  D_Even, D_Odd)

  integer  :: O_Of_D, Filter_Choice, N, Mode_CutOff
  integer  :: N2, M, j, k, l
  REALTYPE :: Omega, Order
  REALTYPE :: XN2, PiN2, D_k, Factor, x

  REALTYPE, dimension(0:N/2,0:N/2,O_Of_D)  :: D_Even, D_Odd
  REALTYPE, dimension(0:N)                 :: Sigma
  REALTYPE, dimension(-N:N,MAX_DERIVATIVE) :: D

  M = O_Of_D

  N2 = N/2

  call Filter (Filter_Choice, N2, Mode_CutOff, Omega, Order, Sigma)

  XN2  = TWO/N
  PiN2 = XN2*PI

  if (M >= 1) then
    D(0,1) = ZERO
    do k = 0,N2
      x = k*PiN2
      D_k = ZERO
      do l = 1,N2
        D_k = D_k + l   *Sigma(l)*SIN(l*x)
      enddo
      D(k,1) = -D_k*XN2
    enddo
  endif

  if (M >= 2) then
    D(0,2) = ZERO
    do k = 0,N2
      x = k*PiN2
      D_k = ZERO
      do l = 1,N2-1
        D_k = D_k + l**2*Sigma(l)*COS(l*x)
      enddo
      D(k,2) = -D_k*XN2
    enddo
  endif

  if (M >= 3) then
    D(0,3) = ZERO
    do k = 0,N2
      x = k*PiN2
      D_k = ZERO
      do l = 1,N2-1
        D_k = D_k + l**3*Sigma(l)*SIN(l*x)
      enddo
      D(k,3) =  D_k*XN2
    enddo
  endif

  if (M >= 4) then
    D(0,4) = ZERO
    do k = 0,N2
      x = k*PiN2
      D_k = ZERO
      do l = 1,N2-1
        D_k = D_k + l**4*Sigma(l)*COS(l*x)
      enddo
      D(k,4) =  D_k*XN2
    enddo
  endif

  do M = 1, O_Of_D
    Factor = (-1)**M

    do k = 1,N2
      D(N-k,M) = Factor*D(k,M)
    enddo

    do k = 1,N-1
      D( -k,M) = Factor*D(k,M)
    enddo
  enddo

  do M = 1, O_Of_D
    D_Even(:,0,M) = D(0:N2,M)
    D_Odd (:,0,M) = ZERO

    do j = 1,N2
      do k = 0,N2
        D_Even(k,j,M) = D(k-j,M)+D(k+j-N,M)
        D_Odd (k,j,M) = D(k-j,M)-D(k+j-N,M)
      enddo
    enddo

    if (Mod(N,2) == 0) then
      D_Even(:,N2,M) = HALF*D_Even(:,N2,M)
      D_Odd (:,N2,M) = ZERO
    endif
  enddo

  END Subroutine F_GL_EOD_DS 

END MODULE Fourier_GL_EOD
