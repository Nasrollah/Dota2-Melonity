#include "pseudopack.h"


MODULE Legendre_GL_MXM
  USE Legendre_Grid_Point

  USE Filtering
  USE MXM_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC  :: L_GL_MXM_D
PUBLIC  :: L_GL_MXM_SDS

CONTAINS

  Subroutine L_GL_MXM_D (N, D)

  integer  :: N

  REALTYPE, dimension(0:N,0:N) :: D
  REALTYPE, dimension(0:N)     :: x, c, LN

  integer  :: j, k, l, N2
  REALTYPE :: Product

! Setup :
  N2 = N/2

  call L_GL_Point_Full (N, x)

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
! Evaluate the legendre polynomial at the Grid Points.
  c(0) = ONE
  do j = 1,N2
    LN(0) = ONE
    LN(1) = x(j)
    do k = 1,N-1
      LN(k+1) = ((2*k+1)*x(j)*LN(k)-k*LN(k-1))/(k+1)
    enddo
    c(j) = LN(N)
  enddo

  do j = 0,N2
    c(N-j) = (-1)**N*c(j)
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
  D(0,0) = QUARTER*N*(N+1)
  do k = 1,N2
    D(k,k) = ZERO
  enddo
#endif

  call Centro_Symmetry_Property    (.TRUE., -1, N, 1, D)

  END Subroutine L_GL_MXM_D

!
! ----------------------------------------------------------------------
!
  Subroutine L_GL_MXM_SDS (O_Of_D, Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2, D)

  integer  :: O_Of_D
  integer  :: Filter_Choice, N, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N  ,0:N) :: D
  REALTYPE, dimension(0:N/2,0:N) :: S
  REALTYPE, dimension(0:N)       :: Sigma

  integer  :: N2

  N2 = N/2

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

  call Smoothing_Matrix (2, 1, 2, N, N2, N, Sigma, S)

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
!
! ----------------------------------------------------------------------
!
    Subroutine L_MXM_DS  

    integer  :: j, k, l, N2
    REALTYPE :: Factor, DS_kj, SD_kj

    N2 = N/2

    Factor = (-1)**O_Of_D

    if (MOD(N,2) == 0) then
      do k = 0,N2-1
        do j = 0,N

          DS_kj = ZERO
          do l = 0,N2-1
            DS_kj = DS_kj+D(k,l)*S(l,j)+D(k,N-l)*S(l,N-j)
          enddo
          D(N-k,N-j) = Factor*(DS_kj+D(k,N2)*S(N2,j))

        enddo
      enddo

      k = N2
        do j = 0,N2

          DS_kj = ZERO
          do l = 0,N2-1
            DS_kj = DS_kj+D(k,l)*S(l,j)+D(k,N-l)*S(l,N-j)
          enddo
          D(0,N-j) = Factor*(DS_kj+D(k,N2)*S(N2,j))

        enddo

        do j = 0,N2
          D(k,N-j) = D(0,N-j)
        enddo
    else
      do k = 0,N2
        do j = 0,N

          DS_kj = ZERO
          do l = 0,N2
            DS_kj = DS_kj+D(k,l)*S(l,j)+D(k,N-l)*S(l,N-j)
          enddo
          D(N-k,N-j) = Factor*DS_kj

        enddo
      enddo
    endif

    call Flip_Up   (O_Of_D, N, D)

    END Subroutine L_MXM_DS
!
! ----------------------------------------------------------------------
!
    Subroutine L_MXM_SD  

    integer  :: j, k, l, N2
    REALTYPE :: Factor, DS_kj, SD_kj

    N2 = N/2

    Factor = (-1)**O_Of_D

    if (MOD(N,2) == 0) then
      do k = 0,N2-1
        do j = 0,N

          SD_kj = ZERO
          do l = 0,N2-1
            SD_kj = SD_kj+S(k,l)*D(l,j)+Factor*S(k,N-l)*D(l,N-j)
          enddo
          D(N-k,N-j) = Factor*(SD_kj+S(k,N2)*D(N2,j))

        enddo
      enddo

      k = N2
        do j = 0,N2

          SD_kj = ZERO
          do l = 0,N2-1
            SD_kj = SD_kj+S(k,l)*D(l,j)+Factor*S(k,N-l)*D(l,N-j)
          enddo
          D(k,N-j) = Factor*(SD_kj+S(k,N2)*D(N2,j))

        enddo
    else
      do k = 0,N2
        do j = 0,N

          SD_kj = ZERO
          do l = 0,N2
            SD_kj = SD_kj+S(k,l)*D(l,j)+Factor*S(k,N-l)*D(l,N-j)
          enddo
          D(N-k,N-j) = Factor*SD_kj

       enddo
      enddo
    endif

    call Flip_Up   (O_Of_D, N, D)

    END Subroutine L_MXM_SD

  END Subroutine L_GL_MXM_SDS

END MODULE Legendre_GL_MXM
