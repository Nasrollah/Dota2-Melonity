#include "pseudopack.h"


MODULE Legendre_GR_MXM
  USE Legendre_Grid_Point

  USE Filtering
  USE MXM_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC  :: L_GR_MXM_D
PUBLIC  :: L_GR_MXM_SDS

CONTAINS

  Subroutine L_GR_MXM_D (N, D)

  integer  :: N

  REALTYPE, dimension(0:N,0:N) :: D
  REALTYPE, dimension(0:N)     :: x, c

  integer  :: j, k, l
  REALTYPE :: Product

! Setup :
  call L_GR_Point_Full (N, x)

! Evaluate the legendre polynomial at the Grid Points.
  do k = 0,N

    Product = ZERO
    do l = 0,N
      if (l /= k) Product = Product + LOG(ABS(x(k)-x(l)))
    enddo
    c(k) = Product

  enddo

! Off-Diagonal :
  do k = 0,N
    do j = 0,N
      if (j /= k) D(k,j) = (-1)**(k+j)*EXP(c(k)-c(j))/(x(k)-x(j))
    enddo
  enddo

! Diagonal :
#if 1
  call Diagonal_Element_By_Sum_MXM (.FALSE., N, 1, D)
#else
  D(0,0) = QUARTER*N*(N+1)
  do k = 1,N
    D(k,k) = ZERO
  enddo
#endif

  END Subroutine L_GR_MXM_D

!
! ----------------------------------------------------------------------
!
  Subroutine L_GR_MXM_SDS (O_Of_D, Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2, D)

  integer  :: O_Of_D
  integer  :: Filter_Choice, N, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N,0:N) :: D, S
  REALTYPE, dimension(0:N)     :: Sigma, F_k

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  call Smoothing_Matrix (2, 2, 2, N, N, N, Sigma, S)

! Compute DS

  if ((Smooth_1 == 1) .AND. (Smooth_2 == 0)) then
    call L_MXM_DS   ; RETURN
  endif

! Compute SD

  if ((Smooth_1 == 0) .AND. (Smooth_2 == 1)) then
    call L_MXM_SD   ; RETURN
  endif

! Compute SDS

  if ((Smooth_1 == 1) .AND. (Smooth_2 == 1)) then
    call L_MXM_DS  
    call L_MXM_SD   ; RETURN
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

  END Subroutine L_GR_MXM_SDS

END MODULE Legendre_GR_MXM
