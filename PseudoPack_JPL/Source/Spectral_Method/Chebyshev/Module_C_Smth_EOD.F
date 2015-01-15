#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Smth_EOD
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : sets up the aux Smooth structure for the (possibly Mapped)
!                Chebyshev Solution Smoothing using the Even-Odd algorithm.
!
!  Index      = 1 if Smoothing along the first index,
!               2 if Smoothing along the second index.
!
!  N_Pnts     = length of the vectors
!  M          = total number of vectors.
!  ldy        = the leading dimension of the array to be processed. (NOT USED)
!  BlockSize  = number of vectors to be processed at one time.
!               if BlockSize=0, BlockSize=M  (NOT USED HERE)
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
!  Omega         = constant in which the value of the Nth mode = Exp(-Omega)
!  Order         = Order of exponential filter = Exp(-a(k/N)**Order)
!  Sigma         = Filter function
! ----------------------------------------------------------------------

MODULE Chebyshev_Smoothing_EOD
  USE EOD_Auxillary_Subroutine

  USE Chebyshev_GL_EOD, ONLY: C_GL_EOD_S_Even, C_GL_EOD_S_Odd
  USE Chebyshev_Ga_EOD, ONLY: C_Ga_EOD_S_Even, C_Ga_EOD_S_Odd

  USE Smoothing_EOD_1, C_Smooth_EOD_1=>Smooth_EOD_1 
  USE Smoothing_EOD_2, C_Smooth_EOD_2=>Smooth_EOD_2 

implicit NONE

PRIVATE
PUBLIC  :: C_Setup_Aux_Smth_EOD
PUBLIC  :: C_Smooth_EOD_1
PUBLIC  :: C_Smooth_EOD_2

CONTAINS

  Subroutine C_Setup_Aux_Smth_EOD (Index, Method, Point_Type,   &
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

  N = N_Pnts-1 ; N1 = N+1 ; N2 = N/2

  i_e_s = EOD_MATRIX_START_EVEN(N1)
  i_o_s = EOD_MATRIX_START_EVEN(N1) + EOD_HALF_SIZE(N1)

  SELECT CASE (Point_Type)
    CASE (1)
      call C_GL_EOD_S_Even (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                                           Aux_Smth(i_e_s))
      call C_GL_EOD_S_Odd  (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                                           Aux_Smth(i_o_s))

    CASE (2)
      write(6 ,*) 'FATAL ERROR! Chebyshev Radau pt is not allowed for EOD'
      write(99,*) 'FATAL ERROR! Chebyshev Radau pt is not allowed for EOD'
      STOP

    CASE (3)
      call C_Ga_EOD_S_Even (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                                           Aux_Smth(i_e_s))
      call C_Ga_EOD_S_Odd  (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                                           Aux_Smth(i_o_s))

  END SELECT

  if (Symmetry == 0)  &
    call Adjust_EOD_Matrix (N2, Aux_Smth(i_e_s), Aux_Smth(i_o_s))

  END Subroutine C_Setup_Aux_Smth_EOD

END MODULE Chebyshev_Smoothing_EOD
