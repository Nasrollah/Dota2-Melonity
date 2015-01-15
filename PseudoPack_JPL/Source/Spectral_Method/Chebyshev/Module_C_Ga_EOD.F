#include "pseudopack.h"


MODULE Chebyshev_Ga_EOD
  USE Filtering

implicit NONE

PRIVATE

PUBLIC  :: C_Ga_EOD_D
PUBLIC  :: C_Ga_EOD_S_Even
PUBLIC  :: C_Ga_EOD_S_Odd

CONTAINS

  Subroutine C_Ga_EOD_D (N, D_Even, D_Odd)

  integer  :: N

  REALTYPE, dimension(0:N/2,0:N/2) :: D_Even, D_Odd
  REALTYPE, dimension(-N:2*N+1)    :: Sin2
  REALTYPE, dimension(0:N)         :: c

  integer  :: j, k, l, N2
  REALTYPE :: PiN
  REALTYPE :: D_kj, D_kl

  N2 = N/2

  PiN  = PI/(2*(N+1))

  do k = -N,-1
    Sin2(k) = ONE/SIN(PiN*k)
  enddo
    Sin2(0) = ZERO
  do k = 1,2*N+1
    Sin2(k) = ONE/SIN(PiN*k)
  enddo

  do k = 0,N
    c(k) = COS(PiN*(2*k+1))
  enddo

  do j = 0,N2
    do k = 0,N2

      if (k /= j) then
        D_kj = HALF*(-1)**(j+k)* Sin2(j+k+1)*Sin2(j-k)*Sin2(2*k+1)/Sin2(2*j+1)
      else
        D_kj = HALF*c(k)       * Sin2(2*k+1)          *Sin2(2*k+1)
      endif

           l = N-j
        D_kl = HALF*(-1)**(l+k)* Sin2(l+k+1)*Sin2(l-k)*Sin2(2*k+1)/Sin2(2*l+1)

      D_Even(k,j) = (D_kj + D_kl)
      D_Odd (k,j) = (D_kj - D_kl)
    enddo
  enddo

  if (MOD(N,2) == 0) then
    D_Even(:,N2) = HALF*D_Even(:,N2)
    D_Odd (:,N2) = ZERO
  endif

  END Subroutine C_Ga_EOD_D 
!
!=======================================================================
!
  Subroutine C_Ga_EOD_S_Even (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

  integer  :: Filter_Choice, N, Mode_CutOff, N2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N/2,0:N/2) :: S
  REALTYPE, dimension(0:N)         :: Sigma

  integer  :: i, j, k, m
  REALTYPE :: PiN, TWON
  REALTYPE :: x_j, x_k, x_m, S_kj, S_km

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  N2 = N/2

  PiN  = PI/(2*(N+1))
  TWON = TWO/(N+1)

  do k = 0,N2
    x_k  = PiN*(2*k+1)
    do j = 0,N2

      m    = j
      x_j  = PiN*(2*m+1)
      S_kj = HALF
      do i = 1,N
        S_kj = S_kj + Sigma(i)*COS(i*x_k)*COS(i*x_j)
      enddo

      m    = N-j
      x_m  = PiN*(2*m+1)
      S_km = HALF
      do i = 1,N
        S_km = S_km + Sigma(i)*COS(i*x_k)*COS(i*x_m)
      enddo

      S(k,j) = (S_kj+S_km)*TWON
    enddo
  enddo

  if (MOD(N,2) == 0) then
    S(:,N2) = HALF*S(:,N2)
  endif

  END Subroutine C_Ga_EOD_S_Even 
!
!=======================================================================
!
  Subroutine C_Ga_EOD_S_Odd  (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

  integer  :: Filter_Choice, N, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N/2,0:N/2) :: S
  REALTYPE, dimension(0:N)         :: Sigma

  integer  :: i, j, k, m, N2
  REALTYPE :: PiN, TWON
  REALTYPE :: x_j, x_k, x_m, S_kj, S_km

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  N2 = N/2

  PiN  = PI/(2*(N+1))
  TWON = TWO/(N+1)

  do k = 0,N2
    x_k  = PiN*(2*k+1)
    do j = 0,N2

      m    = j
      x_j  = PiN*(2*m+1)
      S_kj = HALF
      do i = 1,N
        S_kj = S_kj + Sigma(i)*COS(i*x_k)*COS(i*x_j)
      enddo

      m    = N-j
      x_m  = PiN*(2*m+1)
      S_km = HALF
      do i = 1,N
        S_km = S_km + Sigma(i)*COS(i*x_k)*COS(i*x_m)
      enddo

      S(k,j) = (S_kj-S_km)*TWON
    enddo
  enddo

  END Subroutine C_Ga_EOD_S_Odd  

END MODULE Chebyshev_Ga_EOD
