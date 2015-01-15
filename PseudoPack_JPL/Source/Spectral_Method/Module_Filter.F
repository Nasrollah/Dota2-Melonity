#include "pseudopack.h"



! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! FUNCTION NAME: Filter 
! DESCRIPTION  : Setup Filtering Coefficient Sigma.
!
!  Filter_Choice = Choice of Filtering function
!                    1. Exponential Function (M_CutOffutOff, alpha, Order)
!                    2. Lanczos filter (2nd order)
!                    3. Raised CoSINe  (2nd order)
!                    4. 8th order Sharpened Rasied CoSINe
!
!  Mode_N+1  = Total number of Modes (a_k, k=0,...,Mode_N)
!  M_CutOff  = Cut-Off Frequency for the Exponential Filter
!  alpha     = constant in which the value of the Nth mode = Exp(-alpha)
!  Order     = Order of EXPonential filter = Exp(-a(k/Mode_N)**Order)
!  Sigma     = Filter function of size Mode_N+1
!
! ----------------------------------------------------------------------

MODULE Filtering

implicit NONE

PRIVATE

PUBLIC  :: Filter
PUBLIC  :: Change_Filter
PUBLIC  :: Smoothing_ON_OFF

CONTAINS

  Subroutine Filter (Filter_Choice, Mode_N, M_CutOff, alpha, Order, Sigma)

  integer  :: Filter_Choice, Mode_N, M_CutOff
  REALTYPE :: alpha, Order

  REALTYPE, dimension(0:Mode_N) :: Sigma

  integer  :: i
  REALTYPE :: PiN, Arg

#ifdef DEBUG
  ASSERT(Mode_N>0,"Filter:Mode_N<=0")
  ASSERT(Filter_Choice>=0,"Filter:Choice<0")
  ASSERT(Filter_Choice<=4,"Filter:Choice>4")
  ASSERT(M_CutOff>=0,"Filter:M_CutOff<0")
  ASSERT(M_CutOff<Mode_N,"Filter:M_CutOff>=Mode_N")
#endif

  PiN = PI/Mode_N

  Sigma(0) = ONE

  SELECT CASE (Filter_Choice)
    CASE DEFAULT
      Sigma = ONE

    CASE (1)
      if (alpha == ZERO) alpha = -LOG(EPSILON(ONE))

      Sigma(1:M_CutOff) = ONE

      if (Order == ZERO) then ; Sigma(M_CutOff+1:Mode_N) = ZERO ; RETURN ; endif

      do i = M_CutOff+1,Mode_N
             Arg = (i-M_CutOff)*ONE/(Mode_N-M_CutOff)
        Sigma(i) = EXP(-alpha*Arg**Order)
      enddo
 
    CASE (2)
      do i = 1,Mode_N
             Arg = i*PiN
        Sigma(i) = SIN(Arg)/Arg
      enddo

    CASE (3)
      do i = 1,Mode_N
             Arg = i*PiN
        Sigma(i) = HALF*(ONE+COS(Arg))
      enddo

    CASE (4)
      do i = 1,Mode_N
             Arg = HALF*(ONE+COS(i*PiN))
        Sigma(i) = (35+Arg*(-84+Arg*(70-20*Arg)))*Arg**4
      enddo
 
    END SELECT 

  END Subroutine Filter 

! ----------------------------------------------------------------------
! FUNCTION NAME: Change_Filter 
! DESCRIPTION  : Replace the built-in filter inside Aux_Diff or Aux_Smth 
!                with a user specified filter function Sigma.
!                Fourier and Chebyshev Transform Algorithm ONLY.
!
!  Mode_N     = The last resolvable mode by the methods
!  Smooth_1   = Turn on Filtering before differentiation for Chebyshev Method
!               Turn on Filtering for Fourier Differentiation
!               Turn on Filtering for Smoothing
!  Smooth_2   = Turn on Filtering aftter differentiation for Chebyshev Method
!               No used otherwise.
!  Operator   = Arrays holding filtering parameters (Aux_Diff or Aux_Smth)
!  Sigma      = New Filter, array of length (Mode_N+1)
!
!  Note : 
!         Aux_???? must be defined by called to Setup_Operator_???? piror to
!           call to this routine.
!
!         Sigma can be computed by call to Filter
!
!         This is essentially the same as the Setup_Operator_???? except that
!         the user will supply the filtering function Sigma.  
!         Make sure that parameters other than those related to filtering
!         are exactly the same as the one used to define Aux_????.
!         ???? can either be Diff or Smth
! ----------------------------------------------------------------------

  Subroutine Change_Filter (Mode_N, Smooth_1,Smooth_2, Sigma, Operator)

  integer  :: N_Pnts, N, N1, Mode_N, i_s_s
  integer  :: Method, Symmetry, Smooth_1, Smooth_2

  REALTYPE, dimension(0:Mode_N) :: Sigma
  REALTYPE, dimension(*)        :: Operator

#if defined (DEBUG)
  integer  :: N0, Algorithm, IO_Error

  Algorithm  = INT(Operator(_I_ALGORITHM))
#endif
  Method     = INT(Operator(_I_METHOD))
  Symmetry   = INT(Operator(_I_SYMMETRY))
  N_Pnts     = INT(Operator(_I_N_PNTS))

#if defined (DEBUG)
  IO_Error = 0

  if (Method == 0) N0 = 2*Mode_N
  if (Method == 1) N0 = Mode_N+1

  if (N_Pnts /= N0) then
    write (6,*) 'Method = ',Method
    write (6,*) 'Mode_N = ',Mode_N
    write (6,*) 'N      = ',N0    ,' computed  from Mode_N'
    write (6,*) 'N_Pnts = ',N_Pnts,' retrieved from Operator'
    write (6,*) 'Inconsistent Transfrom Length : N, N_Pnts'
    IO_Error = 1
  endif

  if ((Method == 2) .OR. (Algorithm /= 2)) then
    write (6,*) 'Filter can only be changed for Transform Algorithm'
    IO_Error = 1
  endif
    
  if (IO_Error == 1) STOP
#endif

  call Smoothing_ON_OFF (Smooth_1, Smooth_2, Operator)

  SELECT CASE (Method)
    CASE (0)
      N  = N_Pnts

      if (Symmetry == 0) then
        N1 = N
        i_s_s = FFT_SIGMA_START_FULL(N1)
      else
        N1 = 2*(N/2)+1
        i_s_s = FFT_SIGMA_START_HALF(N1)
      endif

    CASE (1)
      N  = N_Pnts-1
      N1 = N+1

      if (Symmetry == 0) then
        i_s_s = CFT_SIGMA_START_FULL(N1)
      else
        i_s_s = CFT_SIGMA_START_HALF(N1)
      endif

  END SELECT

  Operator(i_s_s:i_s_s+Mode_N) = Sigma

  END Subroutine Change_Filter 

! ----------------------------------------------------------------------
! FUNCTION NAME: Smoothing_ON_OFF 
! DESCRIPTION  : Routines to turn on/off the filtering for Operator
! ----------------------------------------------------------------------
  Subroutine Smoothing_ON_OFF (Smooth_1, Smooth_2, Operator)

  integer  :: Smooth_1, Smooth_2

  REALTYPE, dimension(*) :: Operator

  Operator(_I_SMOOTH_1) = Smooth_1
  Operator(_I_SMOOTH_2) = Smooth_2

  END Subroutine Smoothing_ON_OFF 

END MODULE Filtering
