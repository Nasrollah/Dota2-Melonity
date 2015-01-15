#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: F_Setup_Aux_Smth_FFT
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : sets up the Aux Smth Data structure for the
!                Smoothing using the transform algorithm.
!
!  Index      = 1 if Filtering along the first index,
!               2 if Filtering along the second index.
!
!  N_Pnts     = length of the vectors
!  M          = total number of vectors.
!  ldy        = the leading dimension of the array to be processed.
!  BlockSize  = number of vectors to be processed at one time. 
!               if BlockSize=0, BlockSize=M
!
! For Smoothing the Derivative...
!
!  call Filter (Filter_Choice, N, MODe_CutOff, Omega, Order, Sigma) 
!
!  Filter_Choice = Choice of Filtering function
!                    0. Turn off the Smoothing.
!                    1. Exponential Filter (MODe_CutOff, Omega, Order)
!                    2. Lanczos Filter (2nd order)
!                    3. Raised Cosine Filter (2nd order)
!                    4. 8th order Sharpened Rasied Cosine Filter
!  MODe_CutOff   = Cut-Off Frequency for the Exponential Filter
!  Omega         = constant in which the value of the Nth mode = Exp(-Omega)
!  Order         = Order of exponential filter = Exp(-a(k/N)**Order)
!  Sigma         = Filter function
! ----------------------------------------------------------------------

MODULE Fourier_Smoothing_FFT
  USE Filtering

  USE Fourier_Smoothing_FFT_1
  USE Fourier_Smoothing_FFT_2

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
PUBLIC  :: F_Setup_Aux_Smth_FFT
PUBLIC  :: F_Smooth_FFT_1
PUBLIC  :: F_Smooth_FFT_2

CONTAINS

  Subroutine F_Setup_Aux_Smth_FFT (Index, Method, Point_Type,   &
                                   Algorithm, Symmetry,         &
                                   N_Pnts, M, ldy, BlockSize,   &
                                   Filter_Choice, Mode_CutOff,  &
                                   Omega, Order, Aux_Smth)

  integer  :: Index, Method, Point_Type, Algorithm, Symmetry
  integer  :: N_Pnts, M, ldy, BlockSize
  integer  :: Filter_Choice, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(*) :: Aux_Smth

  integer  :: N, N0, N1, N2
  integer  :: i_fft, i_bft, i_cos, i_sin, i_s_s 

#if defined(DEBUG)
  SELECT CASE (Point_Type)
    CASE (1)
      ASSERT(MOD(N_Pnts,2)==0,"F_set_a_s_f:N Odd")

    CASE (2)
      write (6 ,100) ; write (99,100) ; STOP

  END SELECT
#endif

    N  = N_Pnts

  if (Symmetry == 0) then
    N0 = N   ; N1 = N      ; N2 = N/2
  else
    N0 = N/2 ; N1 = 2*N0+1 ; N2 = N0
  endif

  if (Symmetry == 0) then
    i_s_s = FFT_SIGMA_START_FULL(N1)
  else
    i_s_s = FFT_SIGMA_START_HALF(N1)
  endif

  call Filter (Filter_Choice, N2, Mode_CutOff, Omega, Order, Aux_Smth(i_s_s))

  if ((Symmetry == 1) .AND. (MOD(N,4) /= 0)) then
    write (6,110) N ; write (99,110) N ; STOP
  endif

  if (Symmetry == 0) then
    i_fft = FFT_FWAVE_START_FULL(N1)
    i_bft = FFT_BWAVE_START_FULL(N1)
  else
    i_cos = FFT_CWAVE_START_HALF(N1)
    i_sin = FFT_SWAVE_START_HALF(N1)
  endif

#if defined(ESSL)
  if (Symmetry == 0) then
    call ESSL_FFT_Init (N0, ldy, BlockSize, Index, Aux_Smth(i_fft),  1)
    call ESSL_FFT_Init (N0, ldy, BlockSize, Index, Aux_Smth(i_bft), -1)
  else
    call ESSL_COS_Init (N0, ldy, BlockSize, Index, Aux_Smth(i_cos))
    call ESSL_SIN_Init (N0, ldy, BlockSize, Index, Aux_Smth(i_sin))
  endif
#endif

#if defined(VFFT)
  if (Symmetry == 0) then
    call VRFFTI        (N0, Aux_Smth(i_fft))
  else
    call VFFT_COS_Init (N0, Aux_Smth(i_cos))
    call VFFT_SIN_Init (N0, Aux_Smth(i_sin))
  endif
#endif

#if defined(CRAY)
  if (Symmetry == 0) then
    call Cray_FFT_Init (N0, Aux_Smth(i_fft))
  else
    call Cray_COS_Init (N0, Aux_Smth(i_cos))
    call Cray_SIN_Init (N0, Aux_Smth(i_sin))
  endif
#endif

 100 format (1x/1x,'FATAL Error : F_Setup_Aux_Smth_FFT'/  &
                1x,'              The Fourier Gauss Radau point type is ',&
                                  'not supported' / &
                1x,'              for the Transform Algorithm'/)
 110 format (1x/1x,'FATAL Error : F_Setup_Aux_Smth_FFT'/  &
                1x,'              Invalid Transform Length, N=',i5,' !'/ &
                1x,'              For Symmetry case, Mod(N,2) = 0'/)

  END Subroutine F_Setup_Aux_Smth_FFT 

END MODULE Fourier_Smoothing_FFT
