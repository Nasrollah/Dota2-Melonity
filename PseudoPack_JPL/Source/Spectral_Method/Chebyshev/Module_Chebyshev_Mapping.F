#include "pseudopack.h"


! ----------------------------------------------------------------------
!  Date         : 3-18-2000
!  Bug          : In Semi_Inf_Log Mapping, a1 = -a0/E2b, not a1 = -a0*E2b
!  Modification : In Semi_Infinite Mapping, allows mapping from
!                 [-1,1] <=> [x_Left,infinity] or
!                 [-1,1] <=> [x_Left, beta] for the Semi_Inf_Log Mapping.
!                 instead of mapping [-1,1] <=> [0,....] only.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------
 
! ----------------------------------------------------------------------
! FUNCTION NAME: C_Mapping
! DESCRIPTION  : For Chebyshev and Legendre Methods.
!                Computes various things associated with Mapping including
!                the mapped grid Points x and its Transformation Metrics.
!
!  N       = Total   number of Points regardless Symmetry
!  N0+1    = Acutual number of Points taking into account of symmetry
!  x_Left  = Left  Physical Domain Points
!  x_Right = Right Physical Domain Points
!  x       = Computational Grid Points (Input)
!            Mapped (if Mapping used) Phsyical Grid Points (Output)
!  A       = coefficients of the Transformation Metrics
!            that make up the Mapped k-th Derivative
!
!  if Map  <  0  use user defined Mapping functions. (see Note below)
!  if Map  =  0  No mapping is used
!  if Map  >  0  use built in Mapping functions.
!
!  if  1 < ABS(Map) <= 10, It is mapping for FINITE        domain mapping
!          =  1  Tal-Ezer Map (alpha, Manual)
!          =  2  TAN Map  (alpha, beta)
!
!  if 10 < ABS(Map) <= 20, It is mapping for SEMI-INFINITE domain mapping
!          = 11  ALGEBRAIC Map   (alpha)
!          = 12  EXP Map (alpha)
!          = 13  LOG Map (alpha, beta)
!
!  if 20 < ABS(Map) <= 30, It is mapping for INFINITE      domain mapping
!          = 21  TAN Map   (alpha)
!          = 22  ARCTANH Map   (alpha)
!          = 23  ALGEBRAIC Map (alpha)
!
!          = 99  Linear Mapping is used since (x_Right-x_Left) .ne. 2 (Output)
!
!  For Chebyshev and Legendre Method (Method=1,2) and Tal-Ezer Map (Map=1) ONLY
!    Manual = 0, the value of alpha is computed from machine zero
!    Manual = 1, the value of alpha is input
!
!  Note :
!
!  5.  Here is the meaning for N, N0, N2 for the Chebyshev/Legendre methods.
!      N0=N   ; if (Symmetry=1) N0=N/2 
!      N2=N0  ; if (x(0) = -99) N2=N   (If full length of x computed regardness
!                                         of symmetry, use only for computing x)
!
!      xi(N)=-1 , xi(N-1) ,..., xi(N/2) ,..., xi(1) , xi(0)=1
!       x(N)    ,  x(N-1) ,...,  x(N/2) ,...,  x(1) , x(0)     x=x(xi)
!      
!      In the other words, if (x(N2)=-1) N2=N 
!                          if (N = N0) No symmetry else symmetry
!                          if (No symmetry and the derivative of the mapping 
!                              is written in terms of x) then one must set the 
!                            derivatives of the mapping at N0=N to be zero for
!                            infinite domain.
! 
!  6.  Here is the meaning for N, N0, N2 for the Fourier method.
!      N0=N-1 ; if (Symmetry=1) N0=N/2 
!      N2=N0  ; if (x(0) = -99) N2=N-1 (If full length of x computed regardness
!                                         of symmetry, use only for computing x)
!
!      x(0)=0 , x(1) ,..., x(N/2) ,..., x(N-1) , x(N)=2 Pi
!      
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! FUNCTION NAME: Tal_Ezer_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE Tal_Ezer_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: Tal_Ezer_Map
PUBLIC :: Map_Tal_Ezer
PUBLIC :: Inverse_Map_Tal_Ezer

CONTAINS

  Subroutine Tal_Ezer_Map                                            &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i
  REALTYPE :: kappa, gamma, delta

#ifdef DEBUG
  if (Manual /= 0) then
    ASSERT(alpha> ZERO, "Tal_Ezer_Map:alpha<=0")
    ASSERT(alpha< ONE,  "Tal_Ezer_Map:alpha> 1")
  endif
#endif
 
  if (Manual == 0) alpha = MAX(ONE/COSH(-LOG(EPSILON(ONE))/N),ZERO)

  if (alpha == ZERO) Map = 0

  if (Map == 0) RETURN

  kappa = ASIN(alpha)
  gamma = -kappa*kappa
  delta =  kappa/alpha

! Computed D_Map_D_x_? based on the Chebyshev Grid Points.

  do i = 0,N0
    D_Map_D_x_1(i) = delta*SQRT(ONE-(alpha*x(i))**2)
    D_Map_D_x_2(i) = gamma*x(i)
    D_Map_D_x_3(i) = gamma*D_Map_D_x_1(i)
    D_Map_D_x_4(i) = gamma*D_Map_D_x_2(i)
  enddo

! ------------------- Mapped Grid Points -------------------------
  x = ASIN(alpha*x)/kappa
! ----------------------------------------------------------------

#if 0
! Computed D_Map_D_x_? based on the Tal-Ezer Mapped Chebyshev Grid Points.

  d1 = kappa
  d2 = d1*d1
  d3 = d2*d1

  do i = 0,N0
    dc = delta*COS(kappa*x(i))
    ds = delta*SIN(kappa*x(i))

    D_Map_D_x_1(i) =  dc
    D_Map_D_x_2(i) = -ds*d1
    D_Map_D_x_3(i) = -dc*d2
    D_Map_D_x_4(i) =  ds*d3
  enddo
#endif

  END Subroutine Tal_Ezer_Map 
! ----------------------------------------------------------------------
! FUNCTION NAME: Tal_Ezer
! DESCRIPTION  : Map the Chebyshev Point to the Tal_Ezer Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Tal_Ezer (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, x_Left, x_Right, xx

#if defined (DEBUG)
  if ((x <  -ONE) .OR. (x >  ONE)) then
    write (6 ,100) x ; write (99,100) x ; STOP
  endif

 100 format (1x,'Map_Tal_Ezer : x is out of range [-1,1]'/1x,'   x = ',f14.7)
#endif

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)
 
  if ((Map == 0) .OR. (alpha == ZERO)) then
    Map_Tal_Ezer = x ; RETURN
  endif

  xx = ASIN(alpha*x)/ASIN(alpha)

  call Map_Computational_2_Physical (1,   -ONE,    ONE, xx, &
                                        x_Left,x_Right, xx)

  if (x == -ONE) xx = x_Left
  if (x ==  ONE) xx = x_Right

  Map_Tal_Ezer = xx

  END FUNCTION Map_Tal_Ezer 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Tal_Ezer
! DESCRIPTION  : Map the Tal_Ezer Mapped Physical Point back to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Tal_Ezer (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, x_Left, x_Right, xx

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)
 
#if defined (DEBUG)
  if (((x_Left <  x_Right) .AND. ((x <  x_Left) .OR. (x >  x_Right))) .OR. &
      ((x_Left >  x_Right) .AND. ((x >  x_Left) .OR. (x <  x_Right)))) then
    write (6 ,100) x, x_Left, x_Right ; write (99,100) x, x_Left, x_Right ; STOP
  endif

 100 format (1x,'Inverse_Map_Tal_Ezer : x is out of range [x_Left,x_Right]'/ &
             1x,'   x, x_Left, x_Right = ',3(1x,f14.7))
#endif

  if ((Map /= 0) .AND. (alpha /= ZERO)) then
    call Map_Physical_2_Computational (1, x_Left, x_Right,  x, &
                                            -ONE,     ONE, xx)

    xx = SIN(ASIN(alpha)*xx)/alpha
  else
    xx = x
  endif

  if (x == x_Left ) xx = -ONE
  if (x == x_Right) xx =  ONE

  Inverse_Map_Tal_Ezer = xx

  END FUNCTION Inverse_Map_Tal_Ezer 

END MODULE Tal_Ezer_Mapping



! ----------------------------------------------------------------------
! FUNCTION NAME: TAN_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE TAN_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: TAN_Map
PUBLIC :: Map_TAN
PUBLIC :: Inverse_Map_TAN

CONTAINS

  Subroutine TAN_Map                                                 &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i
  REALTYPE :: kappa, gamma, delta, omega, x_0, s_0,  &
              a0, a1, d2, d4, y, c, s2, c2, f0, f1, f2, f3

#ifdef DEBUG
  ASSERT(alpha> 0, "TAN_Map:alpha<=0")
  ASSERT(beta>=x_Left, "TAN_Map:beta<x_Left")
  ASSERT(beta<=x_Right, "TAN_Map:beta>x_Right")
#endif
 
  call Map_Physical_2_Computational (1, x_Left, x_Right, beta, &
                                          -ONE,     ONE,  x_0)

  if ((alpha == ONE) .AND. (x_0 == ZERO)) Map = 0

  if (Map == 0) RETURN

  kappa = ATAN(alpha*(ONE-x_0))
  gamma = ATAN(alpha*(ONE+x_0))
  delta = HALF*(kappa+gamma)
  omega = HALF*(kappa-gamma)
  s_0   = -omega/delta

! Computed D_Map_D_x_? based on the Chebyshev Grid Points.
    
  a0 = alpha/delta
  d2 = TWO*delta
  d4 = d2*d2

  do i = 0,N0
    y   = delta*(x(i)-s_0)
    c   = COS(y)
    s2  = alpha*SIN(TWO*y)
    c2  = alpha*COS(TWO*y)

    f0 =  a0*c*c
    f1 = -   s2
    f2 = -d2*c2
    f3 =  d4*s2

    D_Map_D_x_1(i) = f0
    D_Map_D_x_2(i) = f0*f1
    D_Map_D_x_3(i) = f0*f1**2+f2*f0**2
    D_Map_D_x_4(i) = f0*f1**3+FOUR*f1*f2*f0**2+f3*f0**3
  enddo

! ------------------- Mapped Grid Points -------------------------
  x = x_0+TAN(delta*x+omega)/alpha
! ----------------------------------------------------------------

#if 0
! Computed D_Map_D_x_? based on the TAN Mapped Chebyshev Grid Points.

  a0  = alpha/delta
  a1  = TWO*alpha*alpha*a0

  do i = 0,N0
    y   = alpha*(x(i)-x_0)
    y2  = y*y
    z   = ONE/(ONE+y2)

    z1  = a0*z
    z2  = a1*z*z
    z3  = z2*z
    z4  = z3*z
     
    D_Map_D_x_1(i) =  z1
    D_Map_D_x_2(i) = -z2*y/alpha
    D_Map_D_x_3(i) =  z3*  (THREE*y2-ONE)
    D_Map_D_x_4(i) = -z4*y*(  y2-ONE)*TWELVE
  enddo
#endif

  END Subroutine TAN_Map
! ----------------------------------------------------------------------
! FUNCTION NAME: Map_TAN
! DESCRIPTION  : Map the Chebyshev Point to the TAN Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_TAN (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: kappa, gamma, delta, omega, x_0, xx

#if defined (DEBUG)
  if ((x <  -ONE) .OR. (x >  ONE)) then
    write (6 ,100) x ; write (99,100) x ; STOP
  endif

 100 format (1x,'Map_TAN : x is out of range [-1,1]'/1x,'   x = ',f14.7)
#endif

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

  if (Map == 0) then
    Map_TAN = x ; RETURN
  endif

  call Map_Physical_2_Computational (1, x_Left, x_Right, beta, &
                                          -ONE,     ONE,  x_0)

  if ((alpha == ONE) .AND. (x_0 == ZERO)) then
    Map_TAN = x ; RETURN
  endif

  kappa = ATAN(alpha*(ONE-x_0))
  gamma = ATAN(alpha*(ONE+x_0))
  delta = HALF*(kappa+gamma)
  omega = HALF*(kappa-gamma)

  xx = x_0+TAN(delta*x+omega)/alpha

  call Map_Computational_2_Physical (1,   -ONE,    ONE, xx, &
                                        x_Left,x_Right, xx)

  Map_TAN = xx

  END FUNCTION Map_TAN 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_TAN
! DESCRIPTION  : Map the TAN Mapped Physical Point back to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_TAN (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: kappa, gamma, delta, omega, x_0, s_0, xx

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

#if defined (DEBUG)
  if (((x_Left <  x_Right) .AND. ((x <  x_Left) .OR. (x >  x_Right))) .OR. &
      ((x_Left >  x_Right) .AND. ((x >  x_Left) .OR. (x <  x_Right)))) then
    write (6 ,100) x, x_Left, x_Right ; write (99,100) x, x_Left, x_Right ; STOP
  endif

 100 format (1x,'Inverse_Map_TAN : x is out of range [x_Left,x_Right]'/ &
             1x,'   x, x_Left, x_Right = ',3(1x,f14.7))
#endif

  if (Map == 0) then
    Inverse_Map_TAN = x ; RETURN
  endif

  call Map_Physical_2_Computational (1, x_Left, x_Right, beta, &
                                          -ONE,     ONE,  x_0)

  if ((alpha == ONE) .AND. (x_0 == ZERO)) then
    Inverse_Map_TAN = x ; RETURN
  endif

  call Map_Physical_2_Computational (1, x_Left, x_Right,  x, &
                                          -ONE,     ONE, xx)

  kappa = ATAN(alpha*(ONE-x_0))
  gamma = ATAN(alpha*(ONE+x_0))
  delta = HALF*(kappa+gamma)
  omega = HALF*(kappa-gamma)
  s_0   = -omega/delta

  xx = s_0+ATAN(alpha*(xx-x_0))/delta

  Inverse_Map_TAN = xx

  END FUNCTION Inverse_Map_TAN 

END MODULE TAN_Mapping

! ----------------------------------------------------------------------
! FUNCTION NAME: Semi_Infinite_ALGEBRAIC_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE Semi_INF_ALGEBRAIC_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: Semi_INF_ALGEBRAIC_Map
PUBLIC :: Map_Semi_INF_ALGEBRAIC
PUBLIC :: Inverse_Map_Semi_INF_ALGEBRAIC

CONTAINS

  Subroutine Semi_INF_ALGEBRAIC_Map                                  &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i
  REALTYPE :: a0, z, f0, f1, f2, f3, xi

#ifdef DEBUG
  ASSERT(ABS(alpha)> 0, "Semi_INF_ALGEBRAIC_Map:alpha=0")
#endif
 
  if (alpha == ZERO) Map = 0

  if (Map == 0) RETURN

  alpha = ABS(alpha)

! Computed D_Map_D_x_? based on the Chebyshev Grid Points.

  a0 = ONE/alpha

  do i = 0,N0
    z = (x(i)-ONE)

    f0 = HALF*a0*z*z
    f1 =  a0*z
    f2 =  a0
    f3 = ZERO

    D_Map_D_x_1(i) = f0
    D_Map_D_x_2(i) = f0*f1
    D_Map_D_x_3(i) = f0*f1**2+f2*f0**2
    D_Map_D_x_4(i) = f0*f1**3+FOUR*f1*f2*f0**2+f3*f0**3
  enddo

! ------------------- Mapped Grid Points -------------------------
  do i = 0,N2
    if (x(i) == ONE) then
      xi = x(i)*(ONE-EPSILON(ONE))
    else
      xi = x(i)
    endif

    x(i) = x_Left + alpha*(ONE+xi)/(ONE-xi)
  enddo

! ----------------------------------------------------------------

#if 0
! Computed D_Map_D_x_? based on the Mapped Chebyshev Grid Point.

    D_Map_D_x_1(0) = ZERO
    D_Map_D_x_2(0) = ZERO
    D_Map_D_x_3(0) = ZERO
    D_Map_D_x_4(0) = ZERO
 
  do i = 1,N0
    z  = ONE/(x(i)+alpha)

    z1 =  TWO  *z *z *alpha
    z2 = -TWO  *z1*z
    z3 = -THREE*z2*z
    z4 = -FOUR *z3*z

    D_Map_D_x_1(i) = z1
    D_Map_D_x_2(i) = z2
    D_Map_D_x_3(i) = z3
    D_Map_D_x_4(i) = z4
  enddo
#endif

  END Subroutine Semi_INF_ALGEBRAIC_Map 
! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Semi_INF_ALGEBRAIC
! DESCRIPTION  : Map the Semi_INF_ALGEBRAIC Mapped Physical Point from
!                Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Semi_INF_ALGEBRAIC (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx, xi

#if defined (DEBUG)
  if ((x <  -ONE) .OR. (x >  ONE)) then
    write (6 ,100) x ; write (99,100) x ; STOP
  endif

 100 format (1x,'Map_Semi_INF_ALGEBRAIC : x is out of range [-1,1]'/ &
             1x,'   x = ',f14.7)
#endif

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

  if (Map == 0) then
    xx = x 
  else
    if (x == ONE) then 
      xi = x*(ONE-EPSILON(ONE)) 
    else 
      xi = x 
    endif

    xx = x_Left + alpha*(ONE+xi)/(ONE-xi)
  endif

  Map_Semi_INF_ALGEBRAIC = xx

  END FUNCTION Map_Semi_INF_ALGEBRAIC 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Semi_INF_ALGEBRAIC
! DESCRIPTION  : Map the Semi_INF_ALGEBRAIC Mapped Physical Point back 
!                to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Semi_INF_ALGEBRAIC (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx, xi

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

#if defined (DEBUG)
  if (x <  x_Left) then
    write (6 ,100) x, x_Left          ; write (99,100) x, x_Left          ; STOP
  endif

 100 format (1x,'Inverse_Map_Semi_INF_ALGEBRAIC : x is out of ', &
                'range [x_Left,Infinite]'/  &
             1x,'   x, x_Left          = ',2(1x,f14.7))
#endif

  if (Map == 0) then
    xx = x 
  else
    xi = x - x_Left
    xx = (xi-alpha)/(xi+alpha)
  endif

  Inverse_Map_Semi_INF_ALGEBRAIC = xx

  END FUNCTION Inverse_Map_Semi_INF_ALGEBRAIC 

END MODULE Semi_INF_ALGEBRAIC_Mapping




! ----------------------------------------------------------------------
! FUNCTION NAME: Semi_Infinite_EXP_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! Note : CHQZ page 73.
! ----------------------------------------------------------------------

MODULE Semi_INF_EXP_Mapping  
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: Semi_INF_EXP_Map  
PUBLIC :: Map_Semi_INF_EXP
PUBLIC :: Inverse_Map_Semi_INF_EXP

CONTAINS

  Subroutine Semi_INF_EXP_Map                                        &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i
  REALTYPE :: a0, f0, f1, f2, f3, xi

#ifdef DEBUG
  ASSERT(ABS(alpha)> 0, "Semi_INF_EXP_Map:alpha=0")
#endif
 
  if (alpha == ZERO) Map = 0

  if (Map == 0) RETURN

! Computed D_Map_D_x_? based on the Chebyshev Grid Points.

  a0 = ONE/alpha

  do i = 0,N0
    f0 = -a0*(x(i)-ONE)
    f1 = -a0
    f2 = ZERO
    f3 = ZERO

    D_Map_D_x_1(i) = f0
    D_Map_D_x_2(i) = f0*f1
    D_Map_D_x_3(i) = f0*f1**2+f2*f0**2
    D_Map_D_x_4(i) = f0*f1**3+FOUR*f1*f2*f0**2+f3*f0**3
  enddo

! ------------------- Mapped Grid Points -------------------------
  do i = 0,N2
    if (x(i) == ONE) then
      xi = x(i)*(ONE-EPSILON(ONE))
    else
      xi = x(i)
    endif

    x(i) = x_Left - alpha*LOG(HALF*(ONE-xi))
  enddo

! ----------------------------------------------------------------

#if 0
! Computed D_Map_D_x_? based on the Mapped Chebyshev Grid Point.

  a0 = ONE/alpha

    D_Map_D_x_1(0) = ZERO
    D_Map_D_x_2(0) = ZERO
    D_Map_D_x_3(0) = ZERO
    D_Map_D_x_4(0) = ZERO

  do i = 1,N0
    z  = TWO*EXP(-a0*x(i))

    z1 =  a0*z
    z2 = -a0*z1
    z3 = -a0*z2
    z4 = -a0*z3
  
    D_Map_D_x_1(i) = z1
    D_Map_D_x_2(i) = z2
    D_Map_D_x_3(i) = z3
    D_Map_D_x_4(i) = z4
  enddo
#endif

  END Subroutine Semi_INF_EXP_Map 
! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Semi_INF_EXP
! DESCRIPTION  : Map the Chebyshev Point to the  Semi_INF_EXP 
!      Mapped Physical Point 
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Semi_INF_EXP (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx, xi

#if defined (DEBUG)
  if ((x <  -ONE) .OR. (x >  ONE)) then
    write (6 ,100) x ; write (99,100) x ; STOP
  endif

 100 format (1x,'Map_Semi_INF_EXP : x is out of range [-1,1]'/ &
             1x,'   x = ',f14.7)
#endif

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

  if (Map == 0) then
    xx = x
  else
    if (x == ONE) then ; xi = x*(ONE-EPSILON(ONE)) ; else ; xi = x ; endif

    xx = x_Left - alpha*LOG(HALF*(ONE-xi))
  endif

  Map_Semi_INF_EXP = xx

  END FUNCTION Map_Semi_INF_EXP 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Semi_INF_EXP
! DESCRIPTION  : Map the Semi_INF_EXP Mapped Physical Point back 
!         to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Semi_INF_EXP (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx, xi

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

#if defined (DEBUG)
  if (x <  x_Left) then
    write (6 ,100) x, x_Left          ; write (99,100) x, x_Left          ; STOP
  endif

 100 format (1x,'Inverse_Map_Semi_INF_EXP : x is out of ', &
                'range [x_Left,Infinite]'/  &
             1x,'   x, x_Left          = ',2(1x,f14.7))
#endif

  if (Map == 0) then
    xx = x
  else
    xi = x - x_Left
    xx = ONE-TWO*EXP(-xi/alpha)
  endif

  Inverse_Map_Semi_INF_EXP = xx

  END FUNCTION Inverse_Map_Semi_INF_EXP 

END MODULE Semi_INF_EXP_Mapping  



! ----------------------------------------------------------------------
! FUNCTION NAME: Semi_Infinite_LOG_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! Note : CHQZ page 73.
! ----------------------------------------------------------------------

MODULE Semi_INF_LOG_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: Semi_INF_LOG_Map
PUBLIC :: Map_Semi_INF_LOG
PUBLIC :: Inverse_Map_Semi_INF_LOG

CONTAINS

  Subroutine Semi_INF_LOG_Map                                        &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i
  REALTYPE :: E2a, E2b, a0, a1, f0, f1, f2, f3

#ifdef DEBUG
  ASSERT(ABS(alpha)> 0, "Semi_INF_LOG_Map:alpha=0")
  ASSERT(beta-x_Left> 0, "Semi_INF_LOG_Map:beta<x_Left")
#endif
 
  if (alpha == ZERO) Map = 0

  if (Map == 0) RETURN

  E2a = EXP(TWO*alpha)
  E2b = (beta-x_Left)/(ONE-E2a)
  a0  = ONE/alpha

! Computed D_Map_D_x_? based on the Chebyshev Grid Points.

! BUG : 3-18-2000
! a1 = -a0*E2b
  a1 = -a0/E2b

  do i = 0,N0
    f0 = a1*EXP(alpha*(x(i)-ONE)) 
    f1 = alpha*f0
    f2 = alpha*f1
    f3 = alpha*f2

    D_Map_D_x_1(i) = f0
    D_Map_D_x_2(i) = f0*f1
    D_Map_D_x_3(i) = f0*f1**2+f2*f0**2
    D_Map_D_x_4(i) = f0*f1**3+FOUR*f1*f2*f0**2+f3*f0**3
  enddo

! ------------------- Mapped Grid Points -------------------------
  x = x_Left + E2b*(EXP(alpha*(ONE-x))-E2a)
! ----------------------------------------------------------------

#if 0
! Computed D_Map_D_x_? based on the Mapped Chebyshev Grid Point.

  a1 = E2a*E2b

  do i = 0,N0
    z  = -ONE/(x(i)+a1)

    z1 =   a0*z
    z2 =   z1*z
    z3 = TWO  *z2*z
    z4 = THREE*z3*z

    D_Map_D_x_1(i) = z1
    D_Map_D_x_2(i) = z2
    D_Map_D_x_3(i) = z3
    D_Map_D_x_4(i) = z4
  enddo
#endif

  END Subroutine Semi_INF_LOG_Map 
! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Semi_INF_LOG
! DESCRIPTION  : Map the Chebyshev Point to the  Semi_INF_LOG 
!      Mapped Physical Point 
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Semi_INF_LOG (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx, E2a, E2b

#if defined (DEBUG)
  if ((x <  -ONE) .OR. (x >  ONE)) then
    write (6 ,100) x ; write (99,100) x ; STOP
  endif

 100 format (1x,'Map_Semi_INF_LOG : x is out of range [-1,1]'/ &
             1x,'   x = ',f14.7)
#endif

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

#if defined (DEBUG)
  ASSERT(beta-x_Left> 0, "Map_Semi_INF_LOG:beta<x_Left")
#endif

  if (Map == 0) then
    xx = x
  else
    E2a = EXP(TWO*alpha)
    E2b = (beta-x_Left)/(ONE-E2a)

    xx = x_Left + E2b*(EXP(alpha*(ONE-x))-E2a)
  endif

  Map_Semi_INF_LOG = xx

  END FUNCTION Map_Semi_INF_LOG 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Semi_INF_LOG
! DESCRIPTION  : Map the Semi_INF_LOG Mapped Physical Point back 
!                to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Semi_INF_LOG (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx, E2a, E2b

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

#if defined (DEBUG)
  if ((x <  x_Left) .OR. (x >  beta)) then
    write (6 ,100) x, x_Left, beta    ; write (99,100) x, x_Left, beta    ; STOP
  endif

 100 format (1x,'Inverse_Map_Semi_INF_LOG : x is out of ', &
                'range [x_Left,beta]'/  &
             1x,'   x, x_Left, beta    = ',3(1x,f14.7))
#endif

  if (Map == 0) then
    xx = x
  else
    E2a = EXP(TWO*alpha)
    E2b = (beta-x_Left)/(ONE-E2a)

    xx = ONE-LOG((x-x_Left)/E2b+E2a)/alpha
  endif

  Inverse_Map_Semi_INF_LOG = xx

  END FUNCTION Inverse_Map_Semi_INF_LOG 

END MODULE Semi_INF_LOG_Mapping


! ----------------------------------------------------------------------
! FUNCTION NAME: Infinite_TAN_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE Infinite_TAN_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: Infinite_TAN_Map
PUBLIC :: Map_Infinite_TAN
PUBLIC :: Inverse_Map_Infinite_TAN

CONTAINS

  Subroutine Infinite_TAN_Map                                        &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i
  REALTYPE :: Pi2, a0, d0, d1, d2, y, c, s2, c2, f0, f1, f2, f3, xi

#ifdef DEBUG
  ASSERT(ABS(alpha)> 0, "Infinite TAN_Map:alpha=0")
#endif
 
  if (alpha == ZERO) Map = 0

  if (Map == 0) Return

  alpha = ABS(alpha)
  Pi2   = HALF*PI

! Computed D_Map_D_x_? based on the Chebyshev Grid Points.
    
  a0 = ONE/alpha
  d0 = ONE/Pi2
  d1 = -PI
  d2 = d1*d1

  do i = 0,N0
    y  = Pi2*x(i)

    c  = COS(y)
    s2 = SIN(TWO*y)
    c2 = COS(TWO*y)

    f0 =  d0*a0*c*c
    f1 = -   a0*s2
    f2 =  d1*a0*c2
    f3 =  d2*a0*s2

    D_Map_D_x_1(i) = f0
    D_Map_D_x_2(i) = f0*f1
    D_Map_D_x_3(i) = f0*f1**2+f2*f0**2
    D_Map_D_x_4(i) = f0*f1**3+FOUR*f1*f2*f0**2+f3*f0**3
  enddo

! ------------------- Mapped Grid Points -------------------------
  do i = 0,N2
    if (ABS(x(i)) == ONE) then
      xi = x(i)*(ONE-EPSILON(ONE))
    else
      xi = x(i)
    endif

    x(i) = alpha*TAN(Pi2*xi)
  enddo

! ----------------------------------------------------------------

#if 0
! Computed D_Map_D_x_? based on the Mapped Chebyshev Grid Point.

  a0 = ONE/Pi2

    D_Map_D_x_1(0) = ZERO
    D_Map_D_x_2(0) = ZERO
    D_Map_D_x_3(0) = ZERO
    D_Map_D_x_4(0) = ZERO

  if (N == N0) then
    D_Map_D_x_1(N) = ZERO
    D_Map_D_x_2(N) = ZERO
    D_Map_D_x_3(N) = ZERO
    D_Map_D_x_4(N) = ZERO

    NN = N0-1
  else
    NN = N0
  endif

  do i = 1,NN
    y   = x(i)/alpha
    y2  = y*y
    z   = ONE/(alpha*(ONE+y2))

    z1  =     a0*z
    z2  = TWO*z1*z
    z3  =     z2*z
    z4  =     z3*z
     
    D_Map_D_x_1(i) =  z1
    D_Map_D_x_2(i) = -z2*y
    D_Map_D_x_3(i) =  z3*  (THREE*y2-ONE)
    D_Map_D_x_4(i) = -z4*y*(  y2-ONE)*TWELVE
  enddo
#endif

  END Subroutine Infinite_TAN_Map 
! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Infinite_TAN
! DESCRIPTION  : Map the Chebyshev Point to the TAN Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Infinite_TAN (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx, xi

#if defined (DEBUG)
  if ((x <  -ONE) .OR. (x >  ONE)) then
    write (6 ,100) x ; write (99,100) x ; STOP
  endif

 100 format (1x,'Map_Infinite_TAN : x is out of range [-1,1]'/ &
             1x,'   x = ',f14.7)
#endif

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)

  if (Map == 0) then
    xx = x
  else
    if (ABS(x) == ONE) then ; xi = x*(ONE-EPSILON(ONE)) ; else ; xi = x ; endif

    xx = alpha*TAN(HALF*PI*xi)
  endif

  Map_Infinite_TAN = xx

  END FUNCTION Map_Infinite_TAN 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Infinite_TAN
! DESCRIPTION  : Map the TAN Mapped Physical Point back to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Infinite_TAN (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

  if (Map == 0) then
    xx = x
  else
    xx = (TWO/PI)*ATAN(x/alpha)
  endif

  Inverse_Map_Infinite_TAN = xx

  END FUNCTION Inverse_Map_Infinite_TAN 

END MODULE Infinite_TAN_Mapping


! ----------------------------------------------------------------------
! FUNCTION NAME: Infinite_ARCTANH_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE Infinite_ARCTANH_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: Infinite_ARCTANH_Map
PUBLIC :: Map_Infinite_ARCTANH
PUBLIC :: Inverse_Map_Infinite_ARCTANH

CONTAINS

  Subroutine Infinite_ARCTANH_Map                                    &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i
  REALTYPE :: a0, f0, f1, f2, f3, xi, ArcTanh

#ifdef DEBUG
  ASSERT(ABS(alpha)> 0, "Finite_ARCTANH_Map:alpha=0")
#endif
 
  if (alpha == ZERO) Map = 0

  if (Map == 0) Return

! Computed D_Map_D_x_? based on the Chebyshev Grid Points.
    
  a0 = ONE/alpha

  do i = 0,N0
    f0 = -a0*(x(i)*x(i)-ONE)
    f1 = -a0*  TWO*x(i)
    f2 = -a0*  TWO
    f3 = ZERO

    D_Map_D_x_1(i) = f0
    D_Map_D_x_2(i) = f0*f1
    D_Map_D_x_3(i) = f0*f1**2+f2*f0**2
    D_Map_D_x_4(i) = f0*f1**3+FOUR*f1*f2*f0**2+f3*f0**3
  enddo

! ------------------- Mapped Grid Points -------------------------
! ArcTanh(x) = ArcCosh[1/sqrt(1-x**2)]

  do i = 0,N2
    if (ABS(x(i)) == ONE) then
      xi = x(i)*(ONE-EPSILON(ONE))
    else
      xi = x(i)
    endif

!   ArcTanh = alpha*ATANH(xi)
!   ArcTanh = alpha*ACOSH( ONE/SQRT(ONE-xi*xi) )
    ArcTanh = alpha*LOG( SQRT( (ONE+ABS(xi))/(ONE-ABS(xi)) ) )

    if (x(i) >  ZERO) x(i) =  ArcTanh
    if (x(i) == ZERO) x(i) =  ZERO
    if (x(i) <  ZERO) x(i) = -ArcTanh
  enddo

! ----------------------------------------------------------------

#if 0
! Computed D_Map_D_x_? based on the Mapped Chebyshev Grid Point.

    D_Map_D_x_1(0) = ZERO
    D_Map_D_x_2(0) = ZERO
    D_Map_D_x_3(0) = ZERO
    D_Map_D_x_4(0) = ZERO

  if (N == N0) then
    D_Map_D_x_1(N) = ZERO
    D_Map_D_x_2(N) = ZERO
    D_Map_D_x_3(N) = ZERO
    D_Map_D_x_4(N) = ZERO

    NN = N0-1
  else
    NN = N0
  endif

  a0 = ONE/alpha
  
  do i = 1,NN
    y  = a0*x(i)
    z  = a0/COSH(y)

    z1 =    z
    z2 = TWO*z1*z
    z3 =     z2*z*z
    z4 =     z3*z
 
    D_Map_D_x_1(i) =  z1
    D_Map_D_x_2(i) = -z2*   TANH(y)
    D_Map_D_x_3(i) =  z3*alpha*(COSH(  TWO*y)-TWO)
    D_Map_D_x_4(i) = -z4*alpha*(SINH(THREE*y)-ELEVEN*SINH(y)))
  enddo
#endif

  END Subroutine Infinite_ARCTANH_Map 
! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Infinite_ARCTANH
! DESCRIPTION  : Map the Chebyshev Point to the ARCTANH Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Infinite_ARCTANH (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha

  REALTYPE :: xx, xi, ArcTanh

#if defined (DEBUG)
  if ((x <  -ONE) .OR. (x >  ONE)) then
    write (6 ,100) x ; write (99,100) x ; STOP
  endif

 100 format (1x,'Map_Infinite_ARCTANH : x is out of range [-1,1]'/ &
             1x,'   x = ',f14.7)
#endif

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)

  if (Map == 0) then
    xx = x
  else
    if (ABS(x) == ONE) then ; xi = x*(ONE-EPSILON(ONE)) ; else ; xi = x ; endif

!    ArcTanh = alpha*ATANH(xi)
!    ArcTanh = alpha*ACOSH( ONE/SQRT(ONE-xi*xi) )
     ArcTanh = alpha*LOG( SQRT( (ONE+ABS(xi))/(ONE-ABS(xi)) ) )

     if (x >  ZERO) xx =  ArcTanh
     if (x == ZERO) xx =  ZERO
     if (x <  ZERO) xx = -ArcTanh
  endif

  Map_Infinite_ARCTANH = xx

  END FUNCTION Map_Infinite_ARCTANH 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Infinite_ARCTANH
! DESCRIPTION  : Map the ARCTANH Mapped Physical Point back to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Infinite_ARCTANH (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha

  REALTYPE :: xx

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)

  if (Map == 0) then
    xx = x
  else
    xx = TANH(x/alpha)
  endif

  Inverse_Map_Infinite_ARCTANH = xx

  END FUNCTION Inverse_Map_Infinite_ARCTANH 

END MODULE Infinite_ARCTANH_Mapping


! ----------------------------------------------------------------------
! FUNCTION NAME: Infinite_ALGEBRAIC_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE Infinite_ALGEBRAIC_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: Infinite_ALGEBRAIC_Map
PUBLIC :: Map_Infinite_ALGEBRAIC
PUBLIC :: Inverse_Map_Infinite_ALGEBRAIC

CONTAINS

  Subroutine Infinite_ALGEBRAIC_Map                                  &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i
  REALTYPE :: a0, y, y2, z, z3, z4, z5, z6, xi

#ifdef DEBUG
  ASSERT(ABS(alpha)> 0, "Finite_ALGEBRAIC_Map:alpha=0")
#endif
 
  if (alpha == ZERO) Map = 0

  if (Map == 0) Return

! Computed D_Map_D_x_? based on the Chebyshev Grid Points.
    
  a0 = ONE/alpha

  do i = 0,N0
    y  = x(i)
    y2 = y*y
    z  = a0*SQRT(ONE-y2)

    z3 = alpha*alpha*z*z*z
    z4 = THREE*z3*z
    z5 =   z4*z
    z6 =   z5*z

    D_Map_D_x_1(i) =  z3
    D_Map_D_x_2(i) = -z4*y
    D_Map_D_x_3(i) =  z5  *( FIVE*y2-ONE)
    D_Map_D_x_4(i) = -z6*y*(SEVEN*y2-THREE)*FIVE
  enddo

! ------------------- Mapped Grid Points -------------------------
  do i = 0,N2
    if (ABS(x(i)) == ONE) then
      xi = x(i)*(ONE-EPSILON(ONE))
    else
      xi = x(i)
    endif

    x(i) = alpha*xi/SQRT(ONE-xi*xi)
  enddo

! ----------------------------------------------------------------

#if 0
! Computed D_Map_D_x_? based on the Mapped Chebyshev Grid Point.

    D_Map_D_x_1(0) = ZERO
    D_Map_D_x_2(0) = ZERO
    D_Map_D_x_3(0) = ZERO
    D_Map_D_x_4(0) = ZERO

  if (N == N0) then
    D_Map_D_x_1(N) = ZERO
    D_Map_D_x_2(N) = ZERO
    D_Map_D_x_3(N) = ZERO
    D_Map_D_x_4(N) = ZERO

    NN = N0-1
  else
    NN = N0
  endif

  a0 = ONE/alpha
  
  do i = 1,NN
    y   = x(i)*a0
    y2  = y*y
    z   = ONE/SQRT(ONE+y2)

    z2  = a0*z*z
    z3  = z2*z
    z5  = THREE*z3*z2
    z7  =   z5*z2
    z9  =   z7*z2
     
    D_Map_D_x_1(i) =  z3
    D_Map_D_x_2(i) = -z5*y
    D_Map_D_x_3(i) =  z7*  (FOUR*y2-ONE)
    D_Map_D_x_4(i) = -z9*y*(FOUR*y2-THREE)*FIVE
  enddo
#endif

  END Subroutine Infinite_ALGEBRAIC_Map 
! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Infinite_ALGEBRAIC
! DESCRIPTION  : Map the Chebyshev Point to the ALGEBRAIC Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Infinite_ALGEBRAIC (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx, xi

#if defined (DEBUG)
  if ((x <  -ONE) .OR. (x >  ONE)) then
    write (6 ,100) x ; write (99,100) x ; STOP
  endif

 100 format (1x,'Map_Infinite_ALGEBRAIC : x is out of range [-1,1]'/ &
             1x,'   x = ',f14.7)
#endif

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)

  if (Map == 0) then
    xx = x
  else
    if (ABS(x) == ONE) then ; xi = x*(ONE-EPSILON(ONE)) ; else ; xi = x ; endif

    xx = alpha*xi/SQRT(ONE-xi*xi)
  endif

  Map_Infinite_ALGEBRAIC = xx

  END FUNCTION Map_Infinite_ALGEBRAIC 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Infinite_ALGEBRAIC
! DESCRIPTION  : Map the ALGEBRAIC Mapped Physical Point back to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Infinite_ALGEBRAIC (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: xx, xi

  Map     = OPERATOR_MAP(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

  if (Map == 0) then
    xx = x
  else
    xi = x/alpha
    xx = xi/SQRT(ONE+xi*xi)
  endif

  Inverse_Map_Infinite_ALGEBRAIC = xx

  END FUNCTION Inverse_Map_Infinite_ALGEBRAIC 

END MODULE Infinite_ALGEBRAIC_Mapping

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Linear_Mapping
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE C_Linear_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: C_Linear_Map
PUBLIC ::         Map_Linear_C
PUBLIC :: Inverse_Map_Linear_C

CONTAINS

  Subroutine C_Linear_Map                                              &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right, Slope

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i

  if (Map == 0) RETURN

! Computed D_Map_D_x_? based on the Chebyshev Grid Points.

  Slope = (x_Right-x_Left)/TWO

  do i = 0,N0
    D_Map_D_x_1(i) = Slope
    D_Map_D_x_2(i) = ZERO
    D_Map_D_x_3(i) = ZERO
    D_Map_D_x_4(i) = ZERO
  enddo

! ------------------- Mapped Grid Points -------------------------
  x = x_Right+Slope*(x-ONE)
! ----------------------------------------------------------------

  END Subroutine C_Linear_Map 
! ----------------------------------------------------------------------
! FUNCTION NAME: Linear
! DESCRIPTION  : Map the Chebyshev Point to the Linear Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Linear_C (Aux_Diff, x)

  REALTYPE               :: x, y
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: x_Left, x_Right

  Map     = OPERATOR_MAP(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)
 
#if defined (DEBUG)
  if ((x <  -ONE) .OR. (x >  ONE)) then
    write (6 ,100) x ; write (99,100) x ; STOP
  endif

 100 format (1x,'Map_Linear_C : x is out of range [-1,1]'/1x,'   x = ',f14.7)
#endif

  call Map_Computational_2_Physical (1,   -ONE,    ONE, x, &
                                        x_Left,x_Right, y)

  Map_Linear_C = y

  END FUNCTION Map_Linear_C 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Linear_C
! DESCRIPTION  : Map the Linear Mapped Physical Point back to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Linear_C (Aux_Diff, x)

  REALTYPE               :: x, y
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map
  REALTYPE :: x_Left, x_Right

  Map     = OPERATOR_MAP(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)
 
#if defined (DEBUG)
  if (((x_Left <  x_Right) .AND. ((x <  x_Left) .OR. (x >  x_Right))) .OR. &
      ((x_Left >  x_Right) .AND. ((x >  x_Left) .OR. (x <  x_Right)))) then
    write (6 ,100) x, x_Left, x_Right ; write (99,100) x, x_Left, x_Right ; STOP
  endif

 100 format (1x,'Inverse_Map_Linear_C : x is out of range [x_Left,x_Right]'/ &
             1x,'   x, x_Left, x_Right = ',3(1x,f14.7))
#endif

  call Map_Physical_2_Computational (1, x_Left, x_Right, x, &
                                          -ONE,     ONE, y)

  Inverse_Map_Linear_C = y

  END FUNCTION Inverse_Map_Linear_C 

END MODULE C_Linear_Mapping
!
! ----------------------------------------------------------------------
!
MODULE Chebyshev_Mapping
  USE Grid_Auxillary_Subroutine

  USE Tal_Ezer_Mapping
  USE TAN_Mapping 
  USE Semi_INF_ALGEBRAIC_Mapping 
  USE Semi_INF_EXP_Mapping 
  USE Semi_INF_LOG_Mapping 
  USE Infinite_TAN_Mapping 
  USE Infinite_ARCTANH_Mapping
  USE Infinite_ALGEBRAIC_Mapping
  USE C_Linear_Mapping

  USE User_Defined_Mapping 

implicit NONE

PUBLIC  ::       C_Mapping
PUBLIC  :: C_Forward_Mapping
PUBLIC  :: C_Inverse_Mapping

CONTAINS

  Subroutine C_Mapping                                               &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  D_Map_D_x_1 = ONE
  D_Map_D_x_2 = ZERO
  D_Map_D_x_3 = ZERO
  D_Map_D_x_4 = ZERO

  SELECT CASE (Map)
    CASE (0)
      RETURN

    CASE (1)
      call Tal_Ezer_Map                                            &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (2)
      call TAN_Map                                                 &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (11)
      call Semi_INF_ALGEBRAIC_Map                                  &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (12)
      call Semi_INF_EXP_Map                                        &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (13)
      call Semi_INF_LOG_Map                                        &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (21)
      call Infinite_TAN_Map                                        &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (22)
      call Infinite_ARCTANH_Map                                    &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (23)
      call Infinite_ALGEBRAIC_Map                                  &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (99)
      call C_Linear_Map                                            &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (:-1)
      call User_Defined_Map                                        &
           (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

#ifdef DEBUG
    CASE DEFAULT
      write ( 6,*) 'Unsupported Chebyshev/Legendre mapping Map=',Map
      write (99,*) 'Unsupported Chebyshev/Legendre mapping Map=',Map

      STOP
#endif

  END SELECT

  END Subroutine C_Mapping
!
! ---------------------------------------------------------------------
!
  REALTYPE FUNCTION C_Forward_Mapping (Operator, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Operator

  SELECT CASE (INT(Operator(_I_MAP)))
    CASE (0)
      C_Forward_Mapping = Map_Linear_C           (Operator, x)
!     C_Forward_Mapping = x

    CASE (1)
      C_Forward_Mapping = Map_Tal_Ezer           (Operator, x)

    CASE (2)
      C_Forward_Mapping = Map_TAN                (Operator, x)

    CASE (11)
      C_Forward_Mapping = Map_Semi_INF_ALGEBRAIC (Operator, x)

    CASE (12)
      C_Forward_Mapping = Map_Semi_INF_EXP       (Operator, x)

    CASE (13)
      C_Forward_Mapping = Map_Semi_INF_LOG       (Operator, x)

    CASE (21)
      C_Forward_Mapping = Map_Infinite_TAN       (Operator, x)

    CASE (22)
      C_Forward_Mapping = Map_Infinite_ARCTANH   (Operator, x)

    CASE (23)
      C_Forward_Mapping = Map_Infinite_ALGEBRAIC (Operator, x)

    CASE (99)
      C_Forward_Mapping = Map_Linear_C           (Operator, x)

    CASE (:-1)
      C_Forward_Mapping = Map_User_Defined       (Operator, x)

#ifdef DEBUG
    CASE DEFAULT
      write (6 ,*) 'Unsupported Chebyshev/Legendre mapping =',Operator(_I_MAP)
      write (99,*) 'Unsupported Chebyshev/Legendre mapping =',Operator(_I_MAP)

      STOP
#endif

  END SELECT

  END FUNCTION C_Forward_Mapping 
!
! ---------------------------------------------------------------------
!
  REALTYPE FUNCTION C_Inverse_Mapping (Operator, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Operator

  SELECT CASE (INT(Operator(_I_MAP)))
    CASE (0)
      C_Inverse_Mapping = Inverse_Map_Linear_C           (Operator, x)
!     C_Inverse_Mapping = x

    CASE (1)
      C_Inverse_Mapping = Inverse_Map_Tal_Ezer           (Operator, x)

    CASE (2)
      C_Inverse_Mapping = Inverse_Map_TAN                (Operator, x)

    CASE (11)
      C_Inverse_Mapping = Inverse_Map_Semi_INF_ALGEBRAIC (Operator, x)

    CASE (12)
      C_Inverse_Mapping = Inverse_Map_Semi_INF_EXP       (Operator, x)

    CASE (13)
      C_Inverse_Mapping = Inverse_Map_Semi_INF_LOG       (Operator, x)

    CASE (21)
      C_Inverse_Mapping = Inverse_Map_Infinite_TAN       (Operator, x)

    CASE (22)
      C_Inverse_Mapping = Inverse_Map_Infinite_ARCTANH   (Operator, x)

    CASE (23)
      C_Inverse_Mapping = Inverse_Map_Infinite_ALGEBRAIC (Operator, x)

    CASE (99)
      C_Inverse_Mapping = Inverse_Map_Linear_C           (Operator, x)

    CASE (:-1)
      C_Inverse_Mapping = Inverse_Map_User_Defined       (Operator, x)

#ifdef DEBUG
    CASE DEFAULT
      write (6 ,*) 'Unsupported Chebyshev/Legendre mapping =',Operator(_I_MAP)
      write (99,*) 'Unsupported Chebyshev/Legendre mapping =',Operator(_I_MAP)

      STOP
#endif

  END SELECT

  END FUNCTION C_Inverse_Mapping 

END MODULE Chebyshev_Mapping
