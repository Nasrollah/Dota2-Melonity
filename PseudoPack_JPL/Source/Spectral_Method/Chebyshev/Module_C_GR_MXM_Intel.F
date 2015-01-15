#include "pseudopack.h"


MODULE Chebyshev_GR_MXM
  USE Filtering
  USE MXM_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC  :: C_GR_MXM_D
PUBLIC  :: C_GR_MXM_S
PUBLIC  :: C_GR_MXM_SDS

CONTAINS

  Subroutine C_GR_MXM_D (N, D)

  integer  :: N

  REALTYPE, dimension(0:N,0:N) :: D
  REALTYPE, dimension(0:2*N)   :: Sin2
  REALTYPE, dimension(0:N)     :: c

  integer  :: j, k
  REALTYPE :: PiN

! Setup :

  PiN = PI/(2*N+1)

  Sin2(0) = ZERO
  do k = 1,2*N
    Sin2(k) = ONE/SIN(PiN*k)
  enddo

  c(0) = HALF
  do j = 1,N
    c(j) = (-1)**j*COS(PiN*j)
  enddo

! Off-Diagonal :
  do k = 0,N
    do j = 0,k-1
      D(k,j) = -HALF*Sin2(k+j)*Sin2(k-j)*c(j)/c(k)
    enddo

    do j = k+1,N
      D(k,j) =  HALF*Sin2(k+j)*Sin2(j-k)*c(j)/c(k)
    enddo
  enddo

! Diagonal :
#if 0
  call Diagonal_Element_By_Sum_MXM (.FALSE., N, 1, D)
#else
  D(0,0) = N*(N+1)/THREE
  do k = 1,N
    D(k,k) = -HALF/SIN(TWO*PiN*k)**2
  enddo
#endif

  END Subroutine C_GR_MXM_D

!
! ----------------------------------------------------------------------
!
  Subroutine C_GR_MXM_S (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

  integer  :: Filter_Choice, N, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N,0:N) :: S
  REALTYPE, dimension(0:N)     :: Sigma

  integer  :: i, j, k
  REALTYPE :: PiN, FOURN
  REALTYPE :: x_j, x_k, S_kj

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  PiN   = TWO*PI/(2*N+1)
  FOURN = FOUR/(2*N+1)

  do k = 0,N
    x_k  = PiN*k
    do j = 0,N
      x_j  = PiN*j
      S_kj = HALF
      do i = 1,N
        S_kj = S_kj + Sigma(i)*COS(i*x_k)*COS(i*x_j)
      enddo
      S(k,j) = S_kj*FOURN
    enddo
  enddo

    S(:,0) = HALF*S(:,0)

  END Subroutine C_GR_MXM_S

!
! ----------------------------------------------------------------------
!
  Subroutine C_GR_MXM_SDS (O_Of_D, Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2, D)

  integer  :: O_Of_D
  integer  :: Filter_Choice, N, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N,0:N) :: D, S

  call C_GR_MXM_S (Filter_Choice, N, Mode_CutOff, Omega, Order, S)

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

  END Subroutine C_GR_MXM_SDS

END MODULE Chebyshev_GR_MXM
