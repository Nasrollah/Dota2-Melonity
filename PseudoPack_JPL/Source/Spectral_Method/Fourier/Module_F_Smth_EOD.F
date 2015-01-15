#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: F_Setup_Aux_Smth_EOD
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : sets up the aux Smooth structure for the (possibly Mapped)
!                Chebyshev Solution Smoothing using the Even-Odd algorithm.
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
! For Smoothing the Smoothative...
!
!     call Filter (Filter_Choice, N, MODe_CutOff, Omega, Order, Sigma)
!
!  Filter_Choice = Choice of Filtering function
!                    1. Exponential Filter (MODe_CutOff, Omega, Order)
!                    2. Lanczos Filter (2nd order)
!                    3. Raised Cosine Filter (2nd order)
!                    4. 8th order Sharpened Rasied Cosine Filter
!  MODe_CutOff   = Cut-Off Frequency for the Exponential Filter
!  Omega         = constant in which the value of the Nth mode = Exp(-Omega)
!  Order         = Order of exponential filter = Exp(-a(k/N)**Order)
!  Sigma         = Filter function
! ----------------------------------------------------------------------

MODULE Fourier_Smoothing_EOD
  USE EOD_Auxillary_Subroutine

  USE Filtering

  USE Smoothing_EOD_1, F_Smooth_EOD_1=>Smooth_EOD_1
  USE Smoothing_EOD_2, F_Smooth_EOD_2=>Smooth_EOD_2

implicit NONE

PRIVATE
PUBLIC  :: F_Setup_Aux_Smth_EOD
PUBLIC  :: F_Smooth_EOD_1
PUBLIC  :: F_Smooth_EOD_2
 
CONTAINS

  Subroutine F_Setup_Aux_Smth_EOD (Index, Method, Point_Type,   &
                                   Algorithm, Symmetry,         &
                                   N_Pnts, M, ldy, BlockSize,   &
                                   Filter_Choice, Mode_CutOff,  &
                                   Omega, Order, Aux_Smth)

  integer  :: Index, Method, Point_Type, Algorithm, Symmetry
  integer  :: N_Pnts, M, ldy, BlockSize
  integer  :: Filter_Choice, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(*) :: Aux_Smth

  integer  :: N, N1, N2
  integer  :: i_e_s, i_o_s

  N = N_Pnts   ; N1 = N+1 ; N2 = N/2

  i_e_s = EOD_MATRIX_START_EVEN(N1)
  i_o_s = EOD_MATRIX_START_EVEN(N1) + EOD_HALF_SIZE(N1)

  if (Point_Type == 2) then
    write (6,100) ; write (99,100) ; STOP
  endif

  call F_GL_Smooth_Matrix_EOD (Aux_Smth(i_e_s), Aux_Smth(i_o_s))

  if (Symmetry == 0)  &
    call Adjust_EOD_Matrix (N2, Aux_Smth(i_e_s), Aux_Smth(i_o_s))

 100 format (1x/1x,'FATAL Error : F_Setup_Aux_Smth_EOD'/  &
                1x,'              Gauss Radau point type is not applicable'/ &
                1x,'              for Even Odd Decomposition Algorithm'/)

  CONTAINS
! ----------------------------------------------------------------------
! FUNCTION NAME: F_GL_Smooth_Matrix
! DESCRIPTION  : computes the elements of the Chebyshev 
!                Smoothing Matrix in Even_Odd format.
!
!  N = Total number of grid points before Symmetry
! ----------------------------------------------------------------------
    Subroutine F_GL_Smooth_Matrix_EOD (S_Even, S_Odd)

    integer  :: N2, k, l
    REALTYPE :: XN2, PiN2, x_k, S_k

    REALTYPE, dimension(-N:N)        :: S
    REALTYPE, dimension(0:N/2)       :: Sigma
    REALTYPE, dimension(0:N/2,0:N/2) :: S_Even, S_Odd

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
      S(k) = S_k*XN2
    enddo

    do k = 0,N2
      S(N-k) = S(k)
    enddo

    do k = 1,N
      S( -k) = S(k)
    enddo

    S_Even(:,0) = S(0:N2)
    S_Odd (:,0) = ZERO

    do l = 1,N2
      do k = 0,N2
        S_Even(k,l) = S(k-l)+S(k+l-N)
        S_Odd (k,l) = S(k-l)-S(k+l-N)
      enddo
    enddo

    if (MOD(N,2) == 0) then
      S_Even(:,N2) = HALF*S_Even(:,N2)
      S_Odd (:,N2) = ZERO
    endif

    END Subroutine F_GL_Smooth_Matrix_EOD 

  END Subroutine F_Setup_Aux_Smth_EOD 

END MODULE Fourier_Smoothing_EOD
