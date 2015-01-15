#if defined (LOCAL_GAMMA)

#include "pseudopack.h"

#if defined (SINGLE_PRECISION)
# define GAMMA_CONSTANT		0.9375E0
#else
# define GAMMA_CONSTANT		0.9375D0
#endif

! ----------------------------------------------------------------
! Author      : Wai Sun Don (10-18-98)
! FUNCTION    : Gamma(x)
! Description : Gamma(x) compute the complete Gamma function with argument x 
! 
!  This is an adaptation of the code DGAMMA from SLATEC Library by Fullerton, W.
!  Please see the original code for details (see below at the bottom).
!
!  This version cleans up the original code and 
!  was rewritten with C-Preprocessor. 
!  The gamma function is now valid for both single and double precision.
!
!  Warning and Fatal error messages are written to the unit 6 (standard Output)
!  and unit 99.
!
!  Must be declared as   External 
!
!  Series for Gamma on the interval [0,1] <Double Precision>
!                   with weighted error   5.79E-32
!                     log weighted error  31.24
!           significant figures required  30.00
!                decimal places required  32.05
! ------------------------------------------------------------------

MODULE Gamma_Function

implicit NONE

integer , SAVE :: N_Gamma_CS
logical , SAVE :: Initialize_Gamma = .TRUE.
REALTYPE, SAVE :: X_Min, X_Max, dx_Relative

integer , SAVE :: N_L_Gamma_CS
logical , SAVE :: Initialize_L_Gamma = .TRUE.
REALTYPE, SAVE :: C_Big, C_Max

#if defined (SINGLE_PRECISION)
REALTYPE, dimension(42), parameter ::                     &
  Gamma_CS = (/ +0.8571195590989331421920062399942E-2   , &
                +0.4415381324841006757191315771652E-2   , &
                +0.5685043681599363378632664588789E-1   , &
                -0.4219835396418560501012500186624E-2   , &
                +0.1326808181212460220584006796352E-2   , &
                -0.1893024529798880432523947023886E-3   , &
                +0.3606925327441245256578082217225E-4   , &
                -0.6056761904460864218485548290365E-5   , &
                +0.1055829546302283344731823509093E-5   , &
                -0.1811967365542384048291855891166E-6   , &
                +0.3117724964715322277790254593169E-7   , &
                -0.5354219639019687140874081024347E-8   , &
                +0.9193275519859588946887786825940E-9   , &
                -0.1577941280288339761767423273953E-9   , &
                +0.2707980622934954543266540433089E-10  , &
                -0.4646818653825730144081661058933E-11  , &
                +0.7973350192007419656460767175359E-12  , &
                -0.1368078209830916025799499172309E-12  , &
                +0.2347319486563800657233471771688E-13  , &
                -0.4027432614949066932766570534699E-14  , &
                +0.6910051747372100912138336975257E-15  , &
                -0.1185584500221992907052387126192E-15  , &
                +0.2034148542496373955201026051932E-16  , &
                -0.3490054341717405849274012949108E-17  , &
                +0.5987993856485305567135051066026E-18  , &
                -0.1027378057872228074490069778431E-18  , &
                +0.1762702816060529824942759660748E-19  , &
                -0.3024320653735306260958772112042E-20  , &
                +0.5188914660218397839717833550506E-21  , &
                -0.8902770842456576692449251601066E-22  , &
                +0.1527474068493342602274596891306E-22  , &
                -0.2620731256187362900257328332799E-23  , &
                +0.4496464047830538670331046570666E-24  , &
                -0.7714712731336877911703901525333E-25  , &
                +0.1323635453126044036486572714666E-25  , &
                -0.2270999412942928816702313813333E-26  , &
                +0.3896418998003991449320816639999E-27  , &
                -0.6685198115125953327792127999999E-28  , &
                +0.1146998663140024384347613866666E-28  , &
                -0.1967938586345134677295103999999E-29  , &
                +0.3376448816585338090334890666666E-30  , &
                -0.5793070335782135784625493333333E-31  /)

REALTYPE, dimension(15), parameter ::                      &
  L_Gamma_CS = (/ +0.1666389480451863247205729650822E+0  , &
                  -0.1384948176067563840732986059135E-4  , &
                  +0.9810825646924729426157171547487E-8  , &
                  -0.1809129475572494194263306266719E-10 , &
                  +0.6221098041892605227126015543416E-13 , &
                  -0.3399615005417721944303330599666E-15 , &
                  +0.2683181998482698748957538846666E-17 , &
                  -0.2868042435334643284144622399999E-19 , &
                  +0.3962837061046434803679306666666E-21 , &
                  -0.6831888753985766870111999999999E-23 , &
                  +0.1429227355942498147573333333333E-24 , &
                  -0.3547598158101070547199999999999E-26 , &
                  +0.1025680058010470912000000000000E-27 , &
                  -0.3401102254316748799999999999999E-29 , &
                  +0.1276642195630062933333333333333E-30 /)

REALTYPE, parameter ::                             &
  Pi      = 3.14159265358979323846264338327950E0 , &
  L_Sq2Pi = 0.91893853320467274178032973640562E0  

#else

REALTYPE, dimension(42), parameter ::                     &
  Gamma_CS = (/ +0.8571195590989331421920062399942D-2   , &
                +0.4415381324841006757191315771652D-2   , &
                +0.5685043681599363378632664588789D-1   , &
                -0.4219835396418560501012500186624D-2   , &
                +0.1326808181212460220584006796352D-2   , &
                -0.1893024529798880432523947023886D-3   , &
                +0.3606925327441245256578082217225D-4   , &
                -0.6056761904460864218485548290365D-5   , &
                +0.1055829546302283344731823509093D-5   , &
                -0.1811967365542384048291855891166D-6   , &
                +0.3117724964715322277790254593169D-7   , &
                -0.5354219639019687140874081024347D-8   , &
                +0.9193275519859588946887786825940D-9   , &
                -0.1577941280288339761767423273953D-9   , &
                +0.2707980622934954543266540433089D-10  , &
                -0.4646818653825730144081661058933D-11  , &
                +0.7973350192007419656460767175359D-12  , &
                -0.1368078209830916025799499172309D-12  , &
                +0.2347319486563800657233471771688D-13  , &
                -0.4027432614949066932766570534699D-14  , &
                +0.6910051747372100912138336975257D-15  , &
                -0.1185584500221992907052387126192D-15  , &
                +0.2034148542496373955201026051932D-16  , &
                -0.3490054341717405849274012949108D-17  , &
                +0.5987993856485305567135051066026D-18  , &
                -0.1027378057872228074490069778431D-18  , &
                +0.1762702816060529824942759660748D-19  , &
                -0.3024320653735306260958772112042D-20  , &
                +0.5188914660218397839717833550506D-21  , &
                -0.8902770842456576692449251601066D-22  , &
                +0.1527474068493342602274596891306D-22  , &
                -0.2620731256187362900257328332799D-23  , &
                +0.4496464047830538670331046570666D-24  , &
                -0.7714712731336877911703901525333D-25  , &
                +0.1323635453126044036486572714666D-25  , &
                -0.2270999412942928816702313813333D-26  , &
                +0.3896418998003991449320816639999D-27  , &
                -0.6685198115125953327792127999999D-28  , &
                +0.1146998663140024384347613866666D-28  , &
                -0.1967938586345134677295103999999D-29  , &
                +0.3376448816585338090334890666666D-30  , &
                -0.5793070335782135784625493333333D-31  /)

REALTYPE, dimension(15), parameter ::                      &
  L_Gamma_CS = (/ +0.1666389480451863247205729650822D+0  , &
                  -0.1384948176067563840732986059135D-4  , &
                  +0.9810825646924729426157171547487D-8  , &
                  -0.1809129475572494194263306266719D-10 , &
                  +0.6221098041892605227126015543416D-13 , &
                  -0.3399615005417721944303330599666D-15 , &
                  +0.2683181998482698748957538846666D-17 , &
                  -0.2868042435334643284144622399999D-19 , &
                  +0.3962837061046434803679306666666D-21 , &
                  -0.6831888753985766870111999999999D-23 , &
                  +0.1429227355942498147573333333333D-24 , &
                  -0.3547598158101070547199999999999D-26 , &
                  +0.1025680058010470912000000000000D-27 , &
                  -0.3401102254316748799999999999999D-29 , &
                  +0.1276642195630062933333333333333D-30 /)

REALTYPE, parameter ::                             &
  Pi      = 3.14159265358979323846264338327950D0 , &
  L_Sq2Pi = 0.91893853320467274178032973640562D0  
#endif

PRIVATE

PUBLIC  :: GAMMA

CONTAINS

  REALTYPE FUNCTION GAMMA (X)

  integer  :: i, N
  REALTYPE :: X, Y, SinPiY
 
  if (X == ZERO) then ; write (6 ,100) ; write (99,100) ; STOP ; endif

  if (Initialize_Gamma) then
    N_Gamma_CS  = Gamma_Accuracy (Gamma_CS, 42, EPSILON(ONE)/RADIX(ONE)/TEN)
    dx_Relative = SQRT(EPSILON(ONE))

    call Gamma_Limit (X_Min, X_Max)

    Initialize_Gamma = .FALSE.
  endif

! Compute Gamma(X) for ABS(X) > 10, recall Y = ABS(X).
  if (ABS(X) >  TEN) then
    if (X  >  X_Max) then ; write (6 ,101) ; write (99,101) ; STOP ; endif
 
    GAMMA = ZERO
    if (X <  X_Min) then ; write (6 ,102) ; write (99,102) ; RETURN ; endif
 
    Y = ABS(X)

    GAMMA = EXP((Y-HALF)*LOG(Y) - Y + L_Sq2Pi + L_Gamma(Y))

    if (X >  ZERO) RETURN
 
    SinPiY = SIN(Pi*Y)
    if (SinPiY == ZERO) then ; write (6 ,103) ; write (99,103) ; STOP ; endif
 
    if (ABS((X-INT(X-HALF))/X) <  dx_Relative) then
      write (6 ,104) ; write (99,104)
    endif
 
    GAMMA = -Pi/(Y*SinPiY*GAMMA)

    RETURN
  endif
 
! Compute Gamma(X) for |X|<=XBND.  
! Reduce Interval and Find Gamma(1+Y) for 0 <= Y < 1 first of all.
 
  N = X
  if (X <  ZERO) N = N - 1
  Y = X - N
  N = N - 1

  GAMMA = GAMMA_CONSTANT + CS_at_X (TWO*Y-ONE, Gamma_CS, N_Gamma_CS)

  if (N == 0) RETURN
 
  if (N >  0) then
! Compute Gamma(X) for X>=2 AND X<=10

    do i = 1,N
      GAMMA = (Y+i)*GAMMA
    enddo

  else
! Compute Gamma(X) for X < 1

    N = -N

    if ((X <  ZERO) .AND. (X+N-2 == ZERO)) then
      write (6 ,103) ; write (99,103) ; STOP
    endif
 
    if ((X <  -HALF) .AND. (ABS((X-INT(X-HALF))/X) <  dx_Relative)) then
      write (6 ,104) ; write (99,104) 
    endif
 
    do i = 1,N
      GAMMA = GAMMA/(X+i-1)
    enddo

  endif

 100 format (1x,'FATAL ERROR : <Gamma> X=0')
 101 format (1x,'FATAL ERROR : <Gamma> X is too big, Gamma overflowed')
 102 format (1x,'WARNING : <Gamma> X is too small, Gamma underflowed')
 103 format (1x,'FATAL ERROR : <Gamma> X is a negative integer')
 104 format (1x,'WARNING : <Gamma> Answer<epsilon/2, X is too close'/ &
             1x,'          to the negative integer')

  CONTAINS

! ----------------------------------------------------------------
! FUNCTION    : Gamma_Limit 
! Description : Calculate the minimum and maximum legal bounds for X 
!                 in Gamma(X).
!              
!               XMIN and XMAX are not the only bounds, but they are the 
!               only non-trivial ones to calculate.
!
! X_MIN   Minimum legal value of X in gamma(X).  
!         Any smaller value of X might result in underflow.
! X_MAX   Maximum legal value of X in gamma(X).  
!         Any larger value of X might cause overflow.
! ----------------------------------------------------------------

    Subroutine Gamma_Limit (X_Min, X_Max)

    integer  :: i
    logical  :: Error
    REALTYPE :: X_Old, X_Min, X_Max, Log_X
    REALTYPE :: X_Min_c, X_Max_c

#if defined (SINGLE_PRECISION)
# define X_MIN_C	0.2258E0
# define X_MAX_C	0.9189E0
#else
# define X_MIN_C	0.2258D0
# define X_MAX_C	0.9189D0
#endif

    X_Min = -LOG(TINY(ONE)) ; X_Min_c = X_MIN_C+X_Min
    do i = 1,10
      X_Old = X_Min ; Log_X = LOG(X_Min)

      X_Min =  X_Min*(-HALF*(Log_X-ONE)+X_Min+X_Min_c)/(X_Min*Log_X+HALF)

      if (ABS(X_Min-X_Old) <  HALF*HUNDREDTH) then 
        Error = .TRUE. ; EXIT 
      endif
    enddo

    if (Error) then ; write (6 ,100) ; write (99,100) ; STOP ; endif
 
    X_Min = -X_Min + HUNDREDTH
 
    X_Max =  LOG(HUGE(ONE)) ; X_Max_c = X_MAX_C-X_Max
    do i = 1,10
      X_Old = X_Max ; Log_X = LOG(X_Max)

      X_Max = -X_Max*(-HALF*(Log_X-ONE)-X_Max+X_Max_c)/(X_Max*Log_X-HALF)

      if (ABS(X_Max-X_Old) <  HALF*HUNDREDTH) then 
        Error = .TRUE. ; EXIT 
      endif
    enddo

   if (Error) then ; write (6 ,101) ; write (99,101) ; STOP ; endif

   X_Max = X_Max - HUNDREDTH

   X_Min = MAX(X_Min, -X_Max+ONE)
 
 100 format (1x,'FATAL ERROR : <Gamma_Limit> Unable to find X_Min')
 101 format (1x,'FATAL ERROR : <Gamma_Limit> Unable to find X_Max')

    END Subroutine Gamma_Limit 

  END FUNCTION GAMMA 
! ----------------------------------------------------------------
! FUNCTION    : L_Gamma (x)
! Description : Compute the log Gamma correction factor so that
!               LOG(Gamma(X)) = LOG(SQRT(2*PI))+(X-5.)*LOG(X)-X+L_Gamma(X).
!               for x <= 10.
!
! Series for L_Gamma on the interval  0 to 0.01 <Double Precision>
!                    with weighted error   1.28E-31
!                     log weighted error  30.89
!           significant figures required  29.81
!                decimal places required  31.48
! ----------------------------------------------------------------

  REALTYPE FUNCTION L_Gamma (X)

  REALTYPE :: X

  if (X <  TEN) then ; write (6 ,100) ; write (99,100) ; STOP ; endif

  if (Initialize_L_Gamma) then
    N_L_Gamma_CS = Gamma_Accuracy (L_Gamma_CS, 15, EPSILON(ONE)/RADIX(ONE))
    C_Big        = ONE/SQRT(EPSILON(ONE)/RADIX(ONE))
    C_Max        = EXP(MIN(LOG(HUGE(ONE)/TWELVE),-LOG(TWELVE*TINY(ONE))))

    Initialize_L_Gamma = .FALSE.
  endif
 
  if (X >= C_Max) then
    L_Gamma = ZERO

    write (6 ,101) ; write (99,101)
  else
    if (X <  C_Big) then
      L_Gamma = CS_at_X (TWO*(TEN/X)**2-ONE, L_Gamma_CS, N_L_Gamma_CS) /X
    else
      L_Gamma = ONE/(TWELVE*X)
    endif
  endif

 100 format (1x,'FATAL ERROR : <L_Gamma> X must be >= 10')
 101 format (1x,'WARNING : <L_Gamma> X is too big, L_Gamma underflowed')

  END FUNCTION L_Gamma
! ----------------------------------------------------------------
! FUNCTION    : CS_at_X (x)
! Description : Evaluate the N-term Chebyshev series CS at X.  
!               Adapted from a method presented in the paper by Broucke.
!
!  X    value at which the series is to be evaluated.
!  CS   array of N terms of a Chebyshev series.  
!       In evaluating CS, only half the first coefficient is summed.
!  N    number of terms in array CS.
!
! ***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
!         Chebyshev series, Algorithm 446, Communications of
!         the A.C.M. 16, (1973) pp. 254-256.
!       L. Fox and I. B. Parker, Chebyshev Polynomials in
!         Numerical Analysis, Oxford University Press, 1968, page 56.
! ----------------------------------------------------------------

  REALTYPE FUNCTION CS_at_X (X, CS, N)

  integer  :: i, N
  REALTYPE :: X, B0, B1, B2

  REALTYPE, dimension(1:N) :: CS

  if (N <  1   ) then ; write (6 ,100) ; write (99,100) ; STOP ; endif
  if (N >  1000) then ; write (6 ,101) ; write (99,101) ; STOP ; endif

  if (ABS(X) >  ONE+EPSILON(ONE)) then ; write (6 ,102) ; write (99,102) ; endif
 
  B0 = ZERO ; B1 = ZERO
  do i = N,1,-1
    B2 = B1
    B1 = B0
    B0 = TWO*X*B1 - B2 + CS(i)
  enddo
 
  CS_at_X = HALF*(B0-B2)
 
 100 format (1x,'FATAL ERROR : <CS_at_X> Number of terms <= 0')
 101 format (1x,'FATAL ERROR : <CS_at_X> Number of terms >  1000')
 102 format (1x,'WARNING : <CS_at_X> |X| > 1')

  END FUNCTION CS_at_X 
! ----------------------------------------------------------------
! Author      : Wai Sun Don (10-18-98)
! FUNCTION    : Gamma_Accuracy
! Description : Determine the number of terms needed in an orthogonal
!               polynomial series so that it meets a specified accuracy.
!
! OS               Array of N term coefficients in an orthogonal series.
! N                Number of coefficients in OS.
! Error_Tolerance  Scalar containing requested accuracy of series.
! ----------------------------------------------------------------

  Integer FUNCTION Gamma_Accuracy (OS, N, Error_Tolerance)

  integer  :: N, i
  REALTYPE :: Error, Error_Tolerance

  REALTYPE, dimension(*) :: OS

  if (N <  1) then ; write (6 ,100) ; write (99,100) ; endif
 
  Error = ZERO
  do i = N,1,-1
    Error = Error + ABS(OS(i))
    if (Error >  Error_Tolerance) EXIT
  enddo
 
  if (i == N) then ; write (6 ,101) i ; write (99,101) i ; endif

  Gamma_Accuracy = i
 
 100 format (1x,'WARNING : <Gamma_Accuracy> Number of Coefficients < 1')
 101 format (1x,'WARNING : <Gamma_Accuracy> Chebyshev Series is too short'/ &
             1x,'          for the specified accuracy! N terms used = ',i4)

  END FUNCTION Gamma_Accuracy 

END MODULE Gamma_Function

#else

MODULE Gamma_Function

implicit NONE

END MODULE Gamma_Function

#endif
