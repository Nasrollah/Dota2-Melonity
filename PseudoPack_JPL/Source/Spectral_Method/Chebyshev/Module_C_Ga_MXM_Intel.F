#include "pseudopack.h"


MODULE Chebyshev_Ga_MXM
  USE Filtering
  USE MXM_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC  :: C_Ga_MXM_D
PUBLIC  :: C_Ga_MXM_S
PUBLIC  :: C_Ga_MXM_SDS

CONTAINS

  Subroutine C_Ga_MXM_D (N, D)

  integer  :: N

  REALTYPE, dimension(0:N,0:N)  :: D
  REALTYPE, dimension(-N:2*N+1) :: Sin2

  integer  :: j, k, N2
  REALTYPE :: PiN

  N2  = N/2
  PiN = PI/(2*(N+1))

  do k = -N,-1
    Sin2(k) = ONE/SIN(PiN*k)
  enddo
    Sin2(0) = ZERO
  do k = 1,2*N+1
    Sin2(k) = ONE/SIN(PiN*k)
  enddo

  do j = 0,N
    do k = 0,N2
      D(k,j) = -HALF*(-1)**(j+k)*Sin2(k+j+1)*Sin2(k-j)*Sin2(2*k+1)/Sin2(2*j+1)
    enddo
  enddo

! Diagonal :
#if 1
  call Diagonal_Element_By_Sum_MXM (.TRUE., N, 1, D)
#else
  do k = 0,N2
    D(k,k) = HALF*COS(PiN*(2*k+1))*(Sin2(2*k+1)*Sin2(2*k+1))
  enddo
#endif

#if 1
  call Centro_Symmetry_Property    (.TRUE., -1, N, 1, D)
#else
  do j = 0,N
    do k = 0,N2
      D(N-k,N-j) = -D(k,j)
    enddo
  enddo
#endif

  END Subroutine C_Ga_MXM_D
!
! ----------------------------------------------------------------------
!
  Subroutine C_Ga_MXM_S (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

  integer  :: Filter_Choice, N, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N,0:N) :: S
  REALTYPE, dimension(0:N)     :: Sigma

  integer  :: i, j, k
  REALTYPE :: PiN, TWON
  REALTYPE :: x_j, x_k, S_kj

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  PiN  = PI/(2*(N+1))
  TWON = TWO/(N+1)

  do k = 0,N
    x_k  = PiN*(2*k+1)
    do j = 0,N
      x_j  = PiN*(2*j+1)
      S_kj = HALF
      do i = 1,N
        S_kj = S_kj + Sigma(i)*COS(i*x_k)*COS(i*x_j)
      enddo
      S(k,j) = S_kj*TWON
    enddo
  enddo

  END Subroutine C_Ga_MXM_S

!
! ----------------------------------------------------------------------
!
  Subroutine C_Ga_MXM_SDS (O_Of_D, Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2, D)

  integer  :: O_Of_D
  integer  :: Filter_Choice, N, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N,0:N) :: S, D
  
  call C_Ga_MXM_S (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

! Compute DS

  if ((Smooth_1 == 1) .AND. (Smooth_2 == 0)) then
    call C_MXM_DS (N, S, D) ; RETURN
  endif

! Compute SD

  if ((Smooth_1 == 0) .AND. (Smooth_2 == 1)) then
    call C_MXM_SD (N, S, D) ; RETURN
  endif

! Compute SDS

  if ((Smooth_1 == 1) .AND. (Smooth_2 == 1)) then
    call C_MXM_DS (N, S, D) 
    call C_MXM_SD (N, S, D) ; RETURN
  endif

  CONTAINS

    Subroutine C_MXM_DS (N, S, D)

    integer                      :: N
    REALTYPE, dimension(0:N,0:N) :: S, D

    integer                      :: j, k, l
    REALTYPE, dimension(0:N)     :: F_k

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
    Subroutine C_MXM_SD (N, S, D) 

    integer                      :: N
    REALTYPE, dimension(0:N,0:N) :: S, D

    integer                      :: j, k, l
    REALTYPE, dimension(0:N)     :: F_k

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

  END Subroutine C_Ga_MXM_SDS

END MODULE Chebyshev_Ga_MXM
