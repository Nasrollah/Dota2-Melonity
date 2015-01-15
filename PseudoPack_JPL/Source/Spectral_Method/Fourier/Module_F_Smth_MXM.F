#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: F_Setup_Aux_Smth_MXM
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : sets up the aux Smooth structure for the (possibly Mapped)
!                Chebyshev Solution Smoothing using the plain MXM algorithm.
! 
!  Index      = 1 if Smoothing along the first index, 
!               2 if Smoothing along the second index. (NOT USED)
!
!  N_Pnts     = length of the vectors
!  M          = total number of vectors.
!  ldy        = the leading dimension of the array to be processed. (NOT USED)
!  BlockSize  = number of vectors to be processed at one time.
!               if BlockSize=0, BlockSize=M  (NOT USED HERE)
!
!  Point_Type = 1 Lobatto Point
!             = 2 Radau   Point
!             = 3 Gauss   Point  (Same as Lobatto Point)
!
!
! For Smoothing the Smoothative...
!
!     call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)
!
!  Filter_Choice = Choice of Filtering function
!                    1. Exponential Filter (Mode_CutOff, Omega, Order)
!                    2. Lanczos Filter (2nd order)
!                    3. Raised Cosine Filter (2nd order)
!                    4. 8th order Sharpened Rasied Cosine Filter
!  Mode_CutOff   = Cut-Off Frequency for the Exponential Filter
!                  if equal N_Pnts-1, Filter_Choice=0
!  Omega         = constant in which the value of the Nth mode = Exp(-Omega)
!  Order         = Order of exponential filter = Exp(-a(k/N)**Order)
!  Sigma         = Filter function
! ----------------------------------------------------------------------

MODULE Fourier_Smoothing_MXM
  USE Filtering

  USE Smoothing_MXM_1, F_Smooth_MXM_1=>Smooth_MXM_1
  USE Smoothing_MXM_2, F_Smooth_MXM_2=>Smooth_MXM_2

implicit NONE

PRIVATE
PUBLIC  :: F_Setup_Aux_Smth_MXM
PUBLIC  :: F_Smooth_MXM_1
PUBLIC  :: F_Smooth_MXM_2
 
CONTAINS

  Subroutine F_Setup_Aux_Smth_MXM (Index, Method, Point_Type,   &
                                   Algorithm, Symmetry,         &
                                   N_Pnts, M, ldy, BlockSize,   &
                                   Filter_Choice, Mode_CutOff,  &
                                   Omega, Order, Aux_Smth)

  integer  :: Index, Method, Point_Type, Algorithm, Symmetry
  integer  :: N_Pnts, M, ldy, BlockSize
  integer  :: Filter_Choice, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(*) :: Aux_Smth

  integer  :: N, N1
  integer  :: i_d_s

  N = N_Pnts ; N1 = N

  i_d_s = MXM_MATRIX_START(N1)

  SELECT CASE (Point_Type)
    CASE (1,3)
      call F_GL_Smooth_Matrix_MXM (Aux_Smth(i_d_s))

    CASE (2)
      call F_GR_Smooth_Matrix_MXM (Aux_Smth(i_d_s))

  END SELECT

  CONTAINS

! ----------------------------------------------------------------------
! FUNCTION : Gauss-Lobatto point
! ----------------------------------------------------------------------
    Subroutine F_GL_Smooth_Matrix_MXM (S)

    integer  :: N2, k, l
    REALTYPE :: XN2, PiN2, x_k, S_k

    REALTYPE, dimension(0:N-1,0:N-1) :: S
    REALTYPE, dimension(-N+1:N-1)    :: SS
    REALTYPE, dimension(0:N/2)       :: Sigma

    N2 = N/2

    call Filter (Filter_Choice, N2, Mode_CutOff, Omega, Order, Sigma)

    XN2  = TWO/N
    PiN2 = XN2*PI

    do k = 0,N2
      x_k = k*PiN2
      S_k = HALF*(ONE+Sigma(N2)*COS(N2*x_k))
      do l = 1,N2-1
        S_k = S_k + Sigma(l)*COS(l*x_k)
      enddo
        SS(k) = S_k*XN2
    enddo

    do k = 1,N2
      SS(N-k) = SS(k)
    enddo

    do k = 1,N-1
      SS( -k) = SS(k)
    enddo

    do l = 0,N-1
      do k = 0,N-1
        S(k,l) = SS(k-l)
      enddo
    enddo

    END Subroutine F_GL_Smooth_Matrix_MXM 
! ----------------------------------------------------------------------
! FUNCTION : Gauss-Radau point
! ----------------------------------------------------------------------
    Subroutine F_GR_Smooth_Matrix_MXM (S)

    integer  :: N2, k, l
    REALTYPE :: XN2, PiN2, x_k, S_k

    REALTYPE, dimension(0:N-1,0:N-1) :: S
    REALTYPE, dimension(-N+1:N-1)    :: SS
    REALTYPE, dimension(0:N/2)       :: Sigma

    N2 = (N-1)/2

    call Filter (Filter_Choice, N2, Mode_CutOff, Omega, Order, Sigma)

    XN2  = TWO/N
    PiN2 = XN2*PI

    do k = 0,N-1
      x_k = k*PiN2
      S_k = HALF
      do l = 1,N2
        S_k = S_k + Sigma(l)*COS(l*x_k)
      enddo
        SS(k) = S_k*XN2
    enddo

    do k = 1,N-1
      SS(-k) = SS(k)
    enddo

    do l = 0,N-1
      do k = 0,N-1
        S(k,l) = SS(k-l)
      enddo
    enddo

    END Subroutine F_GR_Smooth_Matrix_MXM 

  END Subroutine F_Setup_Aux_Smth_MXM 

END MODULE Fourier_Smoothing_MXM
