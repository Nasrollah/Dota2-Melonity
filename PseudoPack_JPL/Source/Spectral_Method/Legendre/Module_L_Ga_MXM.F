#include "pseudopack.h"


MODULE Legendre_Ga_MXM
  USE Legendre_Grid_Point

  USE Filtering
  USE MXM_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC  :: L_Ga_MXM_D
PUBLIC  :: L_Ga_MXM_SDS

CONTAINS

  Subroutine L_Ga_MXM_D (N, D)

  integer  :: N

  REALTYPE, dimension(0:N,0:N) :: D
  REALTYPE, dimension(0:N)     :: x, c, LN

  integer  :: j, k, l, N2
  REALTYPE :: Product

! Setup :
  N2 = N/2

  call L_Ga_Point_Full (N, x)

#if 1
  do k = 0,N

    Product = ZERO
    do l = 0,N
      if (l /= k) Product = Product + LOG(ABS(x(k)-x(l)))
    enddo
    c(k) = Product

  enddo

! Off-Diagonal :
  do k = 0,N2
    do j = 0,N
      if (j /= k) D(k,j) = (-1)**(k+j)*EXP(c(k)-c(j))/(x(k)-x(j))
    enddo
  enddo

#else
  c(0) = HALF*(N+1)*(N+2)
  do j = 1,N2
    LN(0) = ONE
    LN(1) = x(j)
    do k = 1,N
      LN(k+1) = ((2*k+1)*x(j)*LN(k)-k*LN(k-1))/(k+1)
    enddo
    c(j) = (N+1)*(LN(N)-x(j)*LN(N+1))/(ONE-x(j)*x(j))
  enddo

  do j = 0,N2
    c(N-j) = (-1)**(N+1)*c(j)
  enddo

! Off-Diagonal :
  do k = 0,N2
    do j = 0,N
      if (k /= j) D(k,j) = (c(k)/c(j))/(x(k)-x(j))
    enddo
  enddo
#endif

! Diagonal :
#if 1
  call Diagonal_Element_By_Sum_MXM (.TRUE., N, 1, D)
#else
  do k = 0,N2
    D(k,k) = x(k)/(ONE-x(k)*x(k))
  enddo
#endif

  call Centro_Symmetry_Property    (.TRUE., -1, N, 1, D)

  END Subroutine L_Ga_MXM_D
!
! ----------------------------------------------------------------------
!
  Subroutine L_Ga_MXM_SDS (O_Of_D, Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2, D)

  integer  :: O_Of_D
  integer  :: Filter_Choice, N, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N,0:N) :: D, S
  REALTYPE, dimension(0:N)     :: Sigma, F_k
  
  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  call Smoothing_Matrix (2, 3, 2, N, N, N, Sigma, S)

! Compute DS

  if ((Smooth_1 == 1) .AND. (Smooth_2 == 0)) then
    call L_MXM_DS  ; RETURN
  endif

! Compute SD

  if ((Smooth_1 == 0) .AND. (Smooth_2 == 1)) then
    call L_MXM_SD  ; RETURN
  endif

! Compute SDS

  if ((Smooth_1 == 1) .AND. (Smooth_2 == 1)) then
    call L_MXM_DS  
    call L_MXM_SD  ; RETURN
  endif

  CONTAINS

    Subroutine L_MXM_DS  

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

    END Subroutine L_MXM_DS
!
! ----------------------------------------------------------------------
!
    Subroutine L_MXM_SD  

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

    END Subroutine L_MXM_SD

  END Subroutine L_Ga_MXM_SDS

END MODULE Legendre_Ga_MXM
