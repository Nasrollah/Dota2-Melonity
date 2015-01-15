#include "pseudopack.h"


MODULE Chebyshev_GL_EOD
  USE Filtering

implicit NONE

PRIVATE

PUBLIC  :: C_GL_EOD_D
PUBLIC  :: C_GL_EOD_S_Even
PUBLIC  :: C_GL_EOD_S_Odd

CONTAINS

  Subroutine C_GL_EOD_D (N, D_Even, D_Odd)

  integer  :: N

  REALTYPE, dimension(0:N/2,0:N/2) :: D_Even, D_Odd
  REALTYPE, dimension(0:N)         :: Sin2

  integer  :: j, k, N2
  REALTYPE :: PiN, PiN2
  REALTYPE :: ckj, S_kj, S_km, S_jpk, S_jmk

  N2 = N/2

  PiN  = PI/N
  PiN2 = HALF*PiN

  Sin2(0) = ZERO
  do k = 1,N-1
    Sin2(k) = ONE/SIN(PiN2*k)
  enddo
  Sin2(N) = ONE

! k = 0
  do j = 1,N
    ckj = ((-1)**j)*Sin2(j)*Sin2(j)

    if (j <= N2) then
      D_Even(0,  j) = ckj
    else
      D_Odd (0,N-j) = ckj
    endif
  enddo
  D_Odd(0,0) = HALF*D_Odd(0,0)

! k = 1,N2 : j = 0, j= N
  do k = 1,N2
    D_Even(k,0) = -((-1)**(  k))*QUARTER*Sin2(  k)*Sin2(  k)
    D_Odd (k,0) =  ((-1)**(N+k))*QUARTER*Sin2(N-k)*Sin2(N-k)
  enddo

! k = 1,N2 : j = 1,N-1
  do k = 1,N2
    do j = 1,N-1
      if (k == j) CYCLE

      if (j+k <= N) then 
        S_jpk = Sin2(j+k) 
      else
        S_jpk = Sin2(N-j+N-k)
      endif

      if (j-k >= 0) then
        S_jmk =  Sin2(j-k)
      else
        S_jmk = -Sin2(k-j)
      endif

      ckj = ((-1)**(k+j))*HALF*S_jpk*S_jmk

      if (j <= N2) then
        D_Even(k,  j) = ckj
      else
        D_Odd (k,N-j) = ckj
      endif

    enddo
  enddo

  D_Even(0,0) = SIXTH*(2*N*N+1)
  do k = 1,N2
    D_Even(k,k) = -HALF*COS(k*PIN)*Sin2(2*k)*Sin2(2*k)
  enddo

  if (MOD(N,2) == 0) then
    D_Odd(:,N2) = D_Even(:,N2)
  endif

  do j = 0,N2
    do k = 0,N2
             S_kj = D_Even(k,j)
             S_km = D_Odd (k,j)

      D_Even(k,j) = S_kj + S_km
      D_Odd (k,j) = S_kj - S_km
    enddo
  enddo

  if (MOD(N,2) == 0) then
    D_Even(:,N2) = HALF*D_Even(:,N2)
    D_Odd (:,N2) = ZERO
  endif

  END Subroutine C_GL_EOD_D 
!
!=======================================================================
!
  Subroutine C_GL_EOD_S_Even (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

  integer  :: Filter_Choice, N, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N/2,0:N/2) :: S
  REALTYPE, dimension(0:N)         :: Sigma

  integer  :: i, j, k, m, N2
  REALTYPE :: PiN, TWON
  REALTYPE :: x_j, x_k, x_m, S_kj, S_km

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  N2 = N/2

  PiN  = PI/N
  TWON = TWO/N

  do k = 0,N2
    x_k  = PiN*k
    do j = 0,N2

      m    = j
      x_j  = PiN*m
      S_kj = HALF*(ONE+Sigma(N)*(-1)**(k+m))
      do i = 1,N-1
        S_kj = S_kj + Sigma(i)*COS(i*x_k)*COS(i*x_j)
      enddo

      m    = N-j
      x_m  = PiN*m
      S_km = HALF*(ONE+Sigma(N)*(-1)**(k+m))
      do i = 1,N-1
        S_km = S_km + Sigma(i)*COS(i*x_k)*COS(i*x_m)
      enddo

      S(k,j) = (S_kj+S_km)*TWON
    enddo
  enddo

    S(:,0) = HALF*S(:,0)

  if (MOD(N,2) == 0) then
    S(:,N2) = HALF*S(:,N2)
  endif

  END Subroutine C_GL_EOD_S_Even 
!
!=======================================================================
!
  Subroutine C_GL_EOD_S_Odd  (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

  integer  :: Filter_Choice, N, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N/2,0:N/2) :: S
  REALTYPE, dimension(0:N)         :: Sigma

  integer  :: i, j, k, m, N2
  REALTYPE :: PiN, TWON
  REALTYPE :: x_j, x_k, x_m, S_kj, S_km

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  N2 = N/2

  PiN  = PI/N
  TWON = TWO/N

  do k = 0,N2
    x_k  = PiN*k
    do j = 0,N2

      m    = j
      x_j  = PiN*m
      S_kj = HALF*(ONE+Sigma(N)*(-1)**(k+m))
      do i = 1,N-1
        S_kj = S_kj + Sigma(i)*COS(i*x_k)*COS(i*x_j)
      enddo

      m    = N-j
      x_m  = PiN*m
      S_km = HALF*(ONE+Sigma(N)*(-1)**(k+m))
      do i = 1,N-1
        S_km = S_km + Sigma(i)*COS(i*x_k)*COS(i*x_m)
      enddo

      S(k,j) = (S_kj-S_km)*TWON
    enddo
  enddo

    S(:,0) = HALF*S(:,0)

  END Subroutine C_GL_EOD_S_Odd  

END MODULE Chebyshev_GL_EOD
