#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Setup_Aux_Diff_CFT
! DESCRIPTION  : sets up the Aux Diff Data structure for using
!                the cosine transform algorithm.
!
!  Index      = 1 if Differentiating along the first index,
!               2 if Differentiating along the second index.
!
!  x_Left     = Left  Physical Domain Point
!  x_Right    = Right Physical Domain Point
!
!  N_Pnts     = length of the vectors
!  M          = total number of vectors.
!  ldy        = the leading dimension of the array to be processed.
!  BlockSize  = number of vectors to be processed at one time.
!               if BlockSize=0, BlockSize=M
!  Symmetry   = 0 if no symmetry assumed
!               1 if data set is symmetric/anti-symmetric about mid-points.
!
! For Smoothing the Derivative...
!
!  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)
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
!  Smooth_1      = Smooth the Flux BEFORE Recursion
!                    0 = NO
!                    1 = Yes
!  Smooth_2      = Smooth the Flux AFTER  Recursion
!                    0 = NO
!                    1 = Yes
! ----------------------------------------------------------------------

MODULE Chebyshev_Differentiation_CFT
  USE Grid_Generation, ONLY: Grid_And_Mapping
  USE Filtering

  USE Chebyshev_Differentiation_CFT_1
  USE Chebyshev_Differentiation_CFT_2

#if defined (ESSL)
  USE ESSL_FFT
#endif
#if defined (VFFT)
  USE VFFT_FFT
#endif
#if defined (CRAY)
  USE CRAY_FFT
#endif

implicit NONE

PRIVATE
PUBLIC  :: C_Setup_Aux_Diff_CFT
PUBLIC  :: C_Diff_CFT_1
PUBLIC  :: C_Diff_CFT_2

CONTAINS

  Subroutine C_Setup_Aux_Diff_CFT (Index, Method, Point_Type,          &
                                   Max_Order_Of_D,                     &
                                   Algorithm, Symmetry,                &
                                   N_Pnts, M, ldy, BlockSize,          &
                                   x, x_Left,x_Right,                  &
                                   Map, Manual, alpha, beta,           &
                                   Filter_Choice, Mode_CutOff,         &
                                   Omega, Order, Smooth_1, Smooth_2,   &
                                   Aux_Diff, EOD_Code)

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry
  integer  :: N_Pnts, M, ldy, BlockSize, Map, Manual
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  integer  :: EOD_Code
  REALTYPE :: x_Left, x_Right
  REALTYPE :: alpha, beta, Omega, Order

  REALTYPE, dimension(:) :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: N, N0, N1, N3
  integer  :: i_m_s, i_s_s, i_w_s, i_q_s

#if defined (DEBUG)
  SELECT CASE (Point_Type)
    CASE (2)
      write (6 ,*) 'CFT for Chebyshev Gauss Radau point is not supported!'
      write (99,*) 'CFT for Chebyshev Gauss Radau point is not supported!'
      STOP

    CASE (3)
      write (6 ,*) 'CFT for Chebyshev Gauss point is not supported!'
      write (99,*) 'CFT for Chebyshev Gauss point is not supported!'
      STOP

  END SELECT
#endif

  N  = N_Pnts-1 ; N0 = N ; N1 = N+1

  if (Symmetry == 1) N0 = N/2 ; N3 = N0+1

  if (Symmetry == 0) then
    i_m_s = CFT_METRIC_START_FULL(N1)
  else
    i_m_s = CFT_METRIC_START_HALF(N1)
  endif

  call Grid_And_Mapping (Method, Point_Type,  &
                         N, N0, Map, Manual, alpha, beta,  &
                         x_Left, x_Right, x, Aux_Diff(i_m_s))

  if (Filter_Choice == 0) then
    Smooth_1 = 0 ; Smooth_2 = 0

    i_s_s = 0
  else
    if (Symmetry == 0) then
      i_s_s = CFT_SIGMA_START_FULL(N1)
    else
      i_s_s = CFT_SIGMA_START_HALF(N1)
    endif

    call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Aux_Diff(i_s_s))
  endif

  Aux_Diff(_I_MATRIX_START)     = 0
  Aux_Diff(_I_ODD_MATRIX_START) = 0
  Aux_Diff(_I_SIGMA_START)      = i_s_s
  Aux_Diff(_I_METRIC_START)     = i_m_s

  EOD_Code = 0

  if ((Symmetry == 1) .AND. (MOD(N,4) /= 0)) then
    write (6,100) N ; write (99,100) N ; STOP
  endif

  if (Symmetry == 0) then
    i_w_s = CFT_WSAVE_START_FULL(N1)
  else
    i_w_s = CFT_WSAVE_START_HALF(N1)
    i_q_s = CFT_QWAVE_START_HALF(N1)
  endif

#if defined (ESSL)

# if defined (SINGLE_PRECISION)
  if (Symmetry == 1) then
    write (6,110) ; write (99,110) ; STOP
  endif
# endif

                     call ESSL_CFT_Init (N0, ldy, BlockSize, Index,  &
                                             Aux_Diff(i_w_s))
  if (Symmetry == 1) call ESSL_QWV_Init (N3, Aux_Diff(i_q_s))
#endif

#if defined (VFFT)

                     call VCosTI        (N3, Aux_Diff(i_w_s))
  if (Symmetry == 1) call VFFT_QWV_Init (N3, Aux_Diff(i_q_s))
#endif

#if defined (CRAY)

                     call Cray_CFT_Init (N3, Aux_Diff(i_w_s))
  if (Symmetry == 1) call Cray_QWV_Init (N3, Aux_Diff(i_q_s))
#endif

 100 format (1x/1x,'Error : Invalid Transform Length, N-1=',i5,' !'/ &
                1x,'        For Symmetry case, MOD(N-1,4) = 0'/)

#if defined (ESSL)
 110 format (1x/1x,'Error : ESSL, Single Precision, Symmetry Case'/ &
                1x,'        is NOT Supported!'/)
#endif

  END Subroutine C_Setup_Aux_Diff_CFT 

END MODULE Chebyshev_Differentiation_CFT
