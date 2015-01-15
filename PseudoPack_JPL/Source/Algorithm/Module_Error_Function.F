#include "pseudopack.h"

! This subprogram computes approximate values for erf(x).
! This subprogram computes approximate values for erfc(x).
! This subprogram computes approximate values for exp(x*x) * erfc(x).
!------------------------------------------------------------------
!
! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!   for a real argument  x.  It contains three FUNCTION type
!   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
!   and one SUBROUTINE type subprogram, CALERF.  The calling
!   statements for the primary entries are:
!
!                   Y=ERF(X)     (or   Y=DERF(X)),
!
!                   Y=ERFC(X)    (or   Y=DERFC(X)),
!   and
!                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
!
!   The routine  CALERF  is intended for internal packet use only,
!   all computations within the packet being concentrated in this
!   routine.  The function subprograms invoke  CALERF  with the
!   statement
!
!          CALL CALERF(ARG,RESULT,JINT)
!
!   where the parameter usage is as follows
!
!      Function                     Parameters for CALERF
!       call              ARG                  Result          JINT
!
!     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
!     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
!     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
!
!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XMIN   = the smallest positive floating-point number.
!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!            the negative of the solution to the equation
!            2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!            A conservative value is the largest machine number X
!            such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!            W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!            machine precision.  A conservative value is
!            1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!            of XINF and 1/[sqrt(pi)*XMIN].
!
!   Approximate values for some important machines are:
!
!                          XMIN       XINF        XNEG     XSMALL
!
!  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
!
!
!                          XBIG       XHUGE       XMAX
!
!  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ERFC = 0      for  ARG .GE. XBIG;
!
!                       ERFCX = XINF  for  ARG .LT. XNEG;
!      and
!                       ERFCX = 0     for  ARG .GE. XMAX.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, EXP
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: March 19, 1990
! --------------------------------------------------------------
!
! This subroutine and associated functions are a modification of 
! the error function erf found in NetLib/Specfun.
! 
! Modification : March 7, 2008
! By           : Wai Sun Don


MODULE Error_Function

IMPLICIT NONE

!------------------------------------------------------------------
! Coefficients for approximation to  erf  in first interval
!------------------------------------------------------------------
REALTYPE, dimension(5), PARAMETER ::                              &
           A = (/ 3.16112374387056560d00, 1.13864154151050156d02, &
                  3.77485237685302021d02, 3.20937758913846947d03, &
                  1.85777706184603153d-1                          /)

REALTYPE, dimension(4), PARAMETER ::                              &
           B = (/ 2.36012909523441209d01, 2.44024637934444173d02, &
                  1.28261652607737228d03, 2.84423683343917062d03  /)

!------------------------------------------------------------------
! Coefficients for approximation to  erfc  in second interval
!------------------------------------------------------------------
REALTYPE, dimension(9), PARAMETER ::                              &
           C = (/ 5.64188496988670089d-1, 8.88314979438837594d00, &
                  6.61191906371416295d01, 2.98635138197400131d02, &
                  8.81952221241769090d02, 1.71204761263407058d03, &
                  2.05107837782607147d03, 1.23033935479799725d03, &
                  2.15311535474403846d-8                          /)

REALTYPE, dimension(8), PARAMETER ::                              &
           D = (/ 1.57449261107098347d01, 1.17693950891312499d02, &
                  5.37181101862009858d02, 1.62138957456669019d03, &
                  3.29079923573345963d03, 4.36261909014324716d03, &
                  3.43936767414372164d03, 1.23033935480374942d03  /)

!------------------------------------------------------------------
! Coefficients for approximation to  erfc  in third interval
!------------------------------------------------------------------
REALTYPE, dimension(6), PARAMETER ::                              &
           P = (/ 3.05326634961232344d-1, 3.60344899949804439d-1, &
                  1.25781726111229246d-1, 1.60837851487422766d-2, &
                  6.58749161529837803d-4, 1.63153871373020978d-2  /)

REALTYPE, dimension(5), PARAMETER ::                              &
           Q = (/ 2.56852019228982242d00, 1.87295284992346047d00, &
                  5.27905102951428412d-1, 6.05183413124413191d-2, &
                  2.33520497626869185d-3                          /)

REALTYPE, PARAMETER :: Threshold = 0.468750d0
REALTYPE, PARAMETER :: SQRT_Pi   = 5.64189583547756286950d-1
REALTYPE, PARAMETER :: Sixteen   = 16.0d0

PRIVATE
PUBLIC  :: Erf, Erfc, Erfcx

CONTAINS

  REALTYPE FUNCTION Erf (x)

  REALTYPE :: x, Answer

  CALL CalErf (x, Answer, 0)

  Erf = Answer

  END FUNCTION Erf

!--------------------------------------------------------------------
  REALTYPE FUNCTION Erfc (x)

  REALTYPE :: x, Answer

  CALL CalErf (x, Answer, 1)

  Erfc = Answer

  END FUNCTION Erfc

!------------------------------------------------------------------
  REALTYPE FUNCTION Erfcx (x)

  REALTYPE :: x, Answer

  CALL CalErf (x, Answer, 2)

  Erfcx = Answer

  END FUNCTION Erfcx

!--------------------------------------------------------------------

  SUBROUTINE CalErf (ARG, Answer, Case)

  INTEGER  :: i, Case
  REALTYPE :: Arg, Answer, X, Y, Z, U
  REALTYPE :: X_Inf, X_Neg, X_Small, X_Big, X_Min, X_Max, X_Huge
  REALTYPE :: X_Num, X_Den

  X_Inf   = HUGE(ONE)                             !   1.79d+308
  X_Neg   = -SQRT(LOG(HALF*X_Inf))                ! -26.628
  X_Small = HALF*EPSILON(ONE)                     !   1.11d-16
  X_Min   = HALF*TINY(ONE)                        !   1.18d-308
!  X_Max   = HALF*MIN(X_Inf, ONE/(SQRT(PI)*X_Min)) !   2.53d+307
  X_Max   = HALF*MIN(X_Inf, ONE/(SQRT(PI)*HALF*TINY(ONE))) !2.53d+307
  X_Huge  = ONE/SQRT(TWO*X_Small)                 !   6.71d07
  X_Big   = 26.5430d0                             !  IEEE
#if defined (SINGLE_PRECISION)
  X_Big   =  9.1940e0
#endif

  X = ARG ; Y = ABS(X)

  if (Y <= Threshold) then ! |X| <= 0.46875
    Z = ZERO ; if (Y >  X_Small) Z = Y*Y

      X_Num = Z*       A(5)  ; X_Den = Z
    do i = 1,3
      X_Num = Z*(X_Num+A(i)) ; X_Den = Z*(X_Den + B(i))
    enddo
      X_Num = X*(X_Num+A(4)) ; X_Den =    X_Den + B(4)

    Answer = X_Num/X_Den

    if (Case /= 0) Answer = ONE-Answer
    if (Case == 2) Answer = EXP(Z) * Answer

    RETURN
  else if (Y <= FOUR) then ! erfc  for 0.46875 <= |X| <= 4.0
      X_Num = Y*       C(9)  ; X_Den = Y
    do i = 1,7
      X_Num = Y*(X_Num+C(i)) ; X_Den = Y*(X_Den + D(i))
    enddo
      X_Num =    X_Num+C(8)  ; X_Den =    X_Den + D(8)

    Answer = X_Num/X_Den

    if (Case /= 2) then
      Z = INT(Y*Sixteen)/Sixteen ; U = (Y-Z)*(Y+Z)

      Answer = Answer*EXP(-Z*Z)*EXP(-U)
    endif
  else ! erfc  for |X| > 4.0
    Answer = ZERO

    if (Y >= X_Big) then
      if ((Case /= 2) .OR. (Y >= X_Max)) GOTO 100

      if (Y >= X_Huge) then ; Answer = SQRT_Pi/Y ; GOTO 100 ; endif
    endif

    Z = ONE/(Y*Y)

      X_Num = Z*       P(6)  ; X_Den = Z
    do i = 1,4
      X_Num = Z*(X_Num+P(i)) ; X_Den = Z*(X_Den + Q(i))
    enddo
      X_Num = Z*(X_Num+P(5)) ; X_Den =    X_Den + Q(5)

    Answer = X_Num/X_Den ; Answer = (SQRT_Pi-Answer)/Y

    if (Case /= 2) then
      Z = INT(Y*Sixteen)/Sixteen ; U = (Y-Z)*(Y+Z)

      Answer = Answer*EXP(-Z*Z)*EXP(-U)
    endif
  endif

100 CONTINUE

  if (Case == 0) Answer = (HALF-Answer) + HALF

  if (X >= ZERO) RETURN

  SELECT CASE (Case)
    CASE (0)
      Answer = -Answer

    CASE (1)
      Answer = TWO-Answer

    CASE (2)
      if (X <  x_Neg) then
        Answer = X_Inf
      else
        Z = INT(X*Sixteen)/Sixteen ; U = (X-Z)*(X+Z)

        Answer = TWO*EXP(Z*Z)*EXP(U)-Answer
      endif

  END SELECT

  END SUBROUTINE CalErf 

END MODULE Error_Function
