#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Smth_CFT
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : sets up the Aux Filt Data structure for Filttion
!                Smoothing using the cosine transform algorithm.
!
!  Index      = 1 if Filterentiating along the first index,
!               2 if Filterentiating along the second index.
!
!  N_Pnts     = length of the vectors
!  M          = total number of vectors.
!  ldy        = the leading dimension of the array to be processed.
!  BlockSize  = number of vectors to be processed at one time.
!               if BlockSize=0, BlockSize=M
!
! For Smoothing the Derivative...
!
!  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)
!
!  Filter_Choice = Choice of Filtering function
!                    0. Turn off the Smoothing.
!                    1. Exponential Filter (Mode_CutOff, Omega, Order)
!                    2. Lanczos Filter (2nd order)
!                    3. Raised Cosine Filter (2nd order)
!                    4. 8th order Sharpened Rasied Cosine Filter
!  Mode_CutOff   = Cut-Off Frequency for the Exponential Filter
!  Omega         = constant in which the value of the Nth mode = Exp(-Omega)
!  Order         = Order of exponential filter = Exp(-a(k/N)**Order)
!  Sigma         = Filter function
! ----------------------------------------------------------------------
 
MODULE Chebyshev_Smoothing_CFT
  USE Filtering

  USE Chebyshev_Smoothing_CFT_1
  USE Chebyshev_Smoothing_CFT_2

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
PUBLIC  :: C_Setup_Aux_Smth_CFT
PUBLIC  :: C_Smooth_CFT_1
PUBLIC  :: C_Smooth_CFT_2

CONTAINS

  Subroutine C_Setup_Aux_Smth_CFT (Index, Method, Point_Type,   &
                                   Algorithm, Symmetry,         &
                                   N_Pnts, M, ldy, BlockSize,   &
                                   Filter_Choice, Mode_CutOff,  &
                                   Omega, Order, Aux_Smth)

  integer  :: Index, Method, Point_Type, Algorithm, Symmetry
  integer  :: N_Pnts, M, ldy, BlockSize
  integer  :: Filter_Choice, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(*) :: Aux_Smth

  integer  :: N, N0, N1, N3
  integer  :: i_s_s, i_w_s, i_q_s

  N  = N_Pnts-1 ; N0 = N ; N1 = N+1

  if (Symmetry == 1) N0 = N/2 ; N3 = N0+1

  if (Symmetry == 0) then
    i_s_s = CFT_SIGMA_START_FULL(N1)
  else
    i_s_s = CFT_SIGMA_START_HALF(N1)
  endif

  call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Aux_Smth(i_s_s))

  if ((Symmetry == 1) .AND. (MOD(N,4) /= 0)) then
    write (6,100) N ; write (99,100) N ; STOP
  endif

  if (Symmetry == 0) then
    i_w_s = CFT_WSAVE_START_FULL(N1)
  else
    i_w_s = CFT_WSAVE_START_HALF(N1)
    i_q_s = CFT_QWAVE_START_HALF(N1)
  endif

#if defined(ESSL)

# if defined (SINGLE_PRECISION)
  if (Symmetry == 1) then
    write (6,110) ; write (99,110) ; STOP
  endif
# endif

                     call ESSL_CFT_Init (N0, ldy, BlockSize, Index,  &
                                             Aux_Smth(i_w_s))
  if (Symmetry == 1) call ESSL_QWV_Init (N3, Aux_Smth(i_q_s))
#endif

#if defined(VFFT)

                     call VCosTI        (N3, Aux_Smth(i_w_s))
  if (Symmetry == 1) call VFFT_QWV_Init (N3, Aux_Smth(i_q_s))
#endif

#if defined(CRAY)

                     call Cray_CFT_Init (N3, Aux_Smth(i_w_s))
  if (Symmetry == 1) call Cray_QWV_Init (N3, Aux_Smth(i_q_s))
#endif

 100 format (1x/1x,'Error : Invalid Transform Length, N-1=',i5,' !'/ &
                1x,'        For Symmetry case, Mod(N-1,4) = 0'/)

#if defined(ESSL)
 110 format (1x/1x,'Error : ESSL, Single Precision, Symmetry Case'/ &
                1x,'        is NOT Supported!'/)
#endif

  END Subroutine C_Setup_Aux_Smth_CFT 

END MODULE Chebyshev_Smoothing_CFT
