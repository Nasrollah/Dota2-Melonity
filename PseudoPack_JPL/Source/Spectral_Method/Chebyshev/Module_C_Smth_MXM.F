#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: C_Setup_Aux_Smth_MXM
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
! For Smoothing the function   ...
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

MODULE Chebyshev_Smoothing_MXM
  USE Chebyshev_GL_MXM, ONLY: C_GL_MXM_S
  USE Chebyshev_GR_MXM, ONLY: C_GR_MXM_S
  USE Chebyshev_Ga_MXM, ONLY: C_Ga_MXM_S

  USE Smoothing_MXM_1, C_Smooth_MXM_1=>Smooth_MXM_1
  USE Smoothing_MXM_2, C_Smooth_MXM_2=>Smooth_MXM_2

implicit NONE

PRIVATE

PUBLIC  :: C_Setup_Aux_Smth_MXM
PUBLIC  :: C_Smooth_MXM_1
PUBLIC  :: C_Smooth_MXM_2

CONTAINS

  Subroutine C_Setup_Aux_Smth_MXM (Index, Method, Point_Type,   &
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

  N = N_Pnts-1 ; N1 = N+1

  i_d_s = MXM_MATRIX_START(N1)

  SELECT CASE (Point_Type)
    CASE (1)
      call C_GL_MXM_S (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                                      Aux_Smth(i_d_s))

    CASE (2)
      call C_GR_MXM_S (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                                      Aux_Smth(i_d_s))

    CASE (3)
      call C_Ga_MXM_S (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                                      Aux_Smth(i_d_s))

  END SELECT

  END Subroutine C_Setup_Aux_Smth_MXM 

END MODULE Chebyshev_Smoothing_MXM
