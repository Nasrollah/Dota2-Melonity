#include "pseudopack.h"


MODULE Chebyshev_GL_MXM
  USE Filtering
  USE MXM_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC  :: C_GL_MXM_D
PUBLIC  :: C_GL_MXM_S
PUBLIC  :: C_GL_MXM_SDS

CONTAINS

  Subroutine C_GL_MXM_D (N, D)

  integer  :: N

  REALTYPE, dimension(0:N,0:N) :: D
  REALTYPE, dimension(0:N)     :: Sin2

  integer  :: j, k, N2
  REALTYPE :: PiN, PiN2
  REALTYPE :: S_jpk

! Setup :

  N2 = N/2

  PiN  = PI/N
  PiN2 = HALF*PiN

  Sin2(0) = ZERO
  do k = 1,N-1
    Sin2(k) = ONE/SIN(PiN2*k)
  enddo
  Sin2(N) = ONE

! K=0 :
  do j = 1,N-1
    D(0,j) = ((-1)**j)*Sin2(j)*Sin2(j)
  enddo

! J=0 :
  do k = 1,N2
    D(k,0) = -QUARTER*D(0,k)
  enddo

! J=N :
  D(0,N) = ((-1)**N)*HALF
  do k = 1,N2
    D(k,N) = ((-1)**(N+k))*QUARTER*Sin2(N-k)*Sin2(N-k)
  enddo

! Upper Triangular :
  do k = 1,N2
    do j = k+1,N-1
      if (j+k >  N) then
        S_jpk = Sin2(N-j+N-k)
      else
        S_jpk = Sin2(j+k)
      endif

      D(k,j) = ((-1)**(j+k))*HALF*Sin2(j-k)*S_jpk
    enddo
  enddo

! Lower Triangular :
  do k = 1,N2-1
    do j = k+1,N2
      D(j,k) = -D(k,j)
    enddo
  enddo

! Diagonal :
#if 1
  call Diagonal_Element_By_Sum_MXM (.TRUE., N, 1, D)
#else
  D(0,0) = SIXTH*(2*N*N+1)
  do k = 1,N2
    D(k,k) = -HALF*COS(PiN*k)*Sin2(2*k)*Sin2(2*k)
  enddo
#endif

  call Centro_Symmetry_Property    (.TRUE., -1, N, 1, D)

  END Subroutine C_GL_MXM_D

!
! ----------------------------------------------------------------------
!
  Subroutine C_GL_MXM_S (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

  integer  :: Filter_Choice, N, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N,0:N) :: S
  REALTYPE, dimension(0:N)     :: Sigma

  integer  :: i, j, k
  REALTYPE :: PiN, TWON
  REALTYPE :: x_j, x_k, S_kj

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  PiN  = PI/N
  TWON = TWO/N

  do k = 0,N
    x_k  = PiN*k
    do j = 0,N
      x_j  = PiN*j
      S_kj = HALF*(ONE+Sigma(N)*(-1)**(k+j))
      do i = 1,N-1
        S_kj = S_kj + Sigma(i)*COS(i*x_k)*COS(i*x_j)
      enddo
      S(k,j) = S_kj*TWON
    enddo
  enddo

  S(:,0) = HALF*S(:,0)
  S(:,N) = HALF*S(:,N)

  END Subroutine C_GL_MXM_S

!
! ----------------------------------------------------------------------
!
  Subroutine C_GL_MXM_SDS (O_Of_D, Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2, D)

  integer  :: O_Of_D
  integer  :: Filter_Choice, N, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N,0:N) :: D, S
  REALTYPE, dimension(0:N)     :: F_k

  call C_GL_MXM_S (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

! Compute DS

  if ((Smooth_1 == 1) .AND. (Smooth_2 == 0)) then
    call C_MXM_DS  ; RETURN
  endif

! Compute SD

  if ((Smooth_1 == 0) .AND. (Smooth_2 == 1)) then
    call C_MXM_SD  ; RETURN
  endif

! Compute SDS

  if ((Smooth_1 == 1) .AND. (Smooth_2 == 1)) then
    call C_MXM_DS  
    call C_MXM_SD  ; RETURN
  endif

  CONTAINS
!
! ----------------------------------------------------------------------
!
    Subroutine C_MXM_DS  

    integer  :: j, k, l

    do k = 0,N

      do j = 0,N
        F_k(j) = ZERO
        do l = 0,N
          F_k(j) = F_k(j) + D(k,l)*S(l,j)
        enddo
      enddo

      D(k,:) = F_k
    enddo

    END Subroutine C_MXM_DS
!
! ----------------------------------------------------------------------
!
    Subroutine C_MXM_SD  

    integer  :: j, k, l

    do k = 0,N

      do j = 0,N
        F_k(j) = ZERO
        do l = 0,N
          F_k(j) = F_k(j)+S(k,l)*D(l,j)
        enddo
      enddo

      S(k,:) = F_k
    enddo

    D = S

    END Subroutine C_MXM_SD

  END Subroutine C_GL_MXM_SDS

END MODULE Chebyshev_GL_MXM
