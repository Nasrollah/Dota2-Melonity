#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: F_Mapping
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : For Fourier Methods 
!                Computes various things associated with Mapping including
!                the mapped grid Points x and its Transformation Metrics.
!
!  if the first elements of x is coded with -99, then regardless of the
!    SYMMETRY setting, the collocations Points x are returned in
!    full domain rather than half domain even if SYMMETRY=1.
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
!  Map     < -1, user defined mapping
!          =  0, No mapping is used
!          =  1, Angular Mapping is used with parameter alpha and beta=Angle.
!
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! FUNCTION NAME: Fourier_Angular_Map
! DESCRIPTION  : Compute the Mapped grid Points AND its four transfORmation
!                 metrics D_Map_D_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE Fourier_Angular_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: Fourier_Angular_Map
PUBLIC :: Map_Angular
PUBLIC :: Inverse_Map_Angular

CONTAINS

  Subroutine Fourier_Angular_Map                                     &
             (N, N0, N2, Map, Map_F,  alpha, beta, x_Left, x_Right,  &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Map_F
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2,   &
                               D_Map_D_x_3, D_Map_D_x_4

  integer  :: i
  REALTYPE :: gamma, kappa, delta, Omega, x_Omega, Angle
  REALTYPE :: s, c, d, f0, f1, f2, f3

#ifdef DEBUG
  ASSERT(ABS(alpha)< ONE, "Fourier_Angular_Map:|alpha|>=1")

  if (((x_Left <  x_Right) .AND. ((beta <  x_Left) .OR. (beta >  x_Right)))   &
                                .OR.                                          &
      ((x_Left >  x_Right) .AND. ((beta >  x_Left) .OR. (beta <  x_Right))))  &
    then
    write (6 ,100) Map_F, beta, x_Left, x_Right
    write (99,100) Map_F, beta, x_Left, x_Right

    STOP
  endif

 100 format (1x,'Fourier_Angular_Map : beta is out of range [x_Left, x_Right]'/&
             1x,'   Map_F, beta, x_Left, x_Right = ',i3, 1x, 3(1x,f14.7))
#endif
 
  if (alpha == ZERO) Map = 0

  if (Map == 0) RETURN

  call Map_Physical_2_Computational (1, x_Left, x_Right, beta, &
                                          ZERO,   TWOPI, Angle)

  gamma = ONE-alpha*alpha
  kappa = ONE+alpha*alpha
  delta = TWO*alpha

  do i = 0,N0
    Omega = x(i)-Angle
    s     = SIN(Omega)
    c     = COS(Omega)
    d     = delta/gamma

    f0 = (kappa-delta*c)/gamma
    f1 =  d*s
    f2 =  d*c
    f3 = -d*s

    D_Map_D_x_1(i) = f0
    D_Map_D_x_2(i) = f0*f1
    D_Map_D_x_3(i) = f0*f1**2+f2*f0**2
    D_Map_D_x_4(i) = f0*f1**3+FOUR*f1*f2*f0**2+f3*f0**3
  enddo

  do i = 0,N2
      Omega = x(i)-Angle
    x_Omega = ATAN2(gamma*SIN(Omega),kappa*COS(Omega)-delta)
       x(i) = D_2_R*MOD(x_Omega*R_2_D+CIRCLE, CIRCLE)
  enddo

#if 0
  do i = 0,N0
    Omega = x(i)
    s     = SIN(Omega)
    c     = COS(Omega)
    d     = delta/gamma

    z1 = gamma/(kappa+delta*c)
    z2 = d*  s*z1*z1
    z3 = d*( c*z1*z1+ TWO*s*z1*z2)
    z4 = d*(-s*z1*z1+FOUR*c*z1*z2+TWO*s*(z1*z3+z2*z2))

    D_Map_D_x_1(i) = z1
    D_Map_D_x_2(i) = z2
    D_Map_D_x_3(i) = z3
    D_Map_D_x_4(i) = z4
  enddo
#endif

  END Subroutine Fourier_Angular_Map
! ----------------------------------------------------------------------
! FUNCTION NAME: Map_Angular
! DESCRIPTION  : Map the Fourier Point to the Angular Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Angular (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map, Map_F
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: gamma, kappa, delta, Omega, x_Omega, Angle, xx

  Map     = OPERATOR_MAP(Aux_Diff)
  Map_F   = OPERATOR_MANUAL(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

#if defined (DEBUG)
  if ((x <  ZERO) .OR. (x >  TWOPI)) then
    write (6 ,100) Map_F, x ; write (99,100) Map_F, x ; STOP
  endif

 100 format (1x,'Map_Angular: x is out of range [0, 2*Pi]'/ &
             1x,'   Map_F, x = ',i3, 1x, f14.7)
#endif

  if ((Map == 0) .OR. (alpha == ZERO)) then
    Map_Angular = x ; RETURN
  endif

  call Map_Physical_2_Computational (1, x_Left, x_Right, beta, &
                                          ZERO,   TWOPI, Angle)

  gamma = ONE-alpha*alpha
  kappa = ONE+alpha*alpha
  delta = TWO*alpha
  Omega = x-Angle

  x_Omega = ATAN2(gamma*SIN(Omega),kappa*COS(Omega)-delta)
       xx = D_2_R*MOD(x_Omega*R_2_D+CIRCLE, CIRCLE)

  call Map_Computational_2_Physical (1,   ZERO,  TWOPI, xx, &
                                        x_Left,x_Right, xx)
  
  Map_Angular = xx

  END FUNCTION Map_Angular 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Angular
! DESCRIPTION  : Map the Angular Mapped Physical Point back to Fourier Point.
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Angular (Aux_Diff, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map, Map_F
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE :: gamma, kappa, delta, Omega, x_Omega, Angle, xx

  Map     = OPERATOR_MAP(Aux_Diff)
  Map_F   = OPERATOR_MANUAL(Aux_Diff)
  alpha   = OPERATOR_ALPHA(Aux_Diff)
  beta    = OPERATOR_BETA(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)

#if defined (DEBUG)
  if (((x_Left <  x_Right) .AND. ((x <  x_Left) .OR. (x >  x_Right))) .OR.  &
      ((x_Left >  x_Right) .AND. ((x >  x_Left) .OR. (x <  x_Right)))) then
    write (6 ,100) Map_F, x, x_Left, x_Right
    write (99,100) Map_F, x, x_Left, x_Right

    STOP
  endif

 100 format (1x,'Inverse_Map_Angular : x is out of range [x_Left, x_Right]'/ &
             1x,'   Map_F, x, x_Left, x_Right = ',i3, 1x, 3(1x,f14.7))
#endif
 
  if ((Map == 0) .OR. (alpha == ZERO)) then
    Inverse_Map_Angular = x ; RETURN
  endif

  call Map_Physical_2_Computational (1, x_Left, x_Right, beta, &
                                          ZERO,   TWOPI, Angle)
  
  call Map_Physical_2_Computational (1, x_Left, x_Right,  x, &
                                          ZERO,   TWOPI, xx)

  gamma = ONE-alpha*alpha
  kappa = ONE+alpha*alpha
  delta = TWO*alpha
  Omega = xx

  x_Omega = ATAN2(gamma*SIN(Omega),kappa*COS(Omega)+delta)
       xx = D_2_R*MOD((x_Omega+Angle)*R_2_D+CIRCLE,CIRCLE)

  Inverse_Map_Angular = xx

  END FUNCTION Inverse_Map_Angular 

END MODULE Fourier_Angular_Mapping

! ----------------------------------------------------------------------
! FUNCTION NAME: F_Linear_Mapping
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_Map_D_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE F_Linear_Mapping
  USE Grid_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC :: F_Linear_Map
PUBLIC ::         Map_Linear_F
PUBLIC :: Inverse_Map_Linear_F

CONTAINS

  Subroutine F_Linear_Map                                              &
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

  Slope = (x_Right-x_Left)/TWOPI

  do i = 0,N0
    D_Map_D_x_1(i) = Slope
    D_Map_D_x_2(i) = ZERO
    D_Map_D_x_3(i) = ZERO
    D_Map_D_x_4(i) = ZERO
  enddo

! ------------------- Mapped Grid Points -------------------------
  x = x_Right+Slope*(x-TWOPI)
! ----------------------------------------------------------------

  END Subroutine F_Linear_Map 
! ----------------------------------------------------------------------
! FUNCTION NAME: Linear
! DESCRIPTION  : Map the Chebyshev Point to the Linear Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_Linear_F (Aux_Diff, x)

  REALTYPE               :: x, y
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map, Map_F
  REALTYPE :: x_Left, x_Right

  Map     = OPERATOR_MAP(Aux_Diff)
  Map_F   = OPERATOR_MANUAL(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)
 
#if defined (DEBUG)
  if ((x <  ZERO) .OR. (x >  TWOPI)) then
    write (6 ,100) Map_F, x ; write (99,100) Map_F, x ; STOP
  endif

 100 format (1x,'Map_Linear_F: x is out of range [0, 2*Pi]'/ &
             1x,'   Map_F, x = ',i3, 1x, f14.7)
#endif

  call Map_Computational_2_Physical (1,   ZERO,  TWOPI, x, &
                                        x_Left,x_Right, y)

  Map_Linear_F = y

  END FUNCTION Map_Linear_F 
! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_Linear_F
! DESCRIPTION  : Map the Linear Mapped Physical Point back to Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_Linear_F (Aux_Diff, x)

  REALTYPE               :: x, y
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: Map, Map_F
  REALTYPE :: x_Left, x_Right

  Map     = OPERATOR_MAP(Aux_Diff)
  Map_F   = OPERATOR_MANUAL(Aux_Diff)
  x_Left  = OPERATOR_X_LEFT(Aux_Diff)
  x_Right = OPERATOR_X_RIGHT(Aux_Diff)
 
#if defined (DEBUG)
  if (((x_Left <  x_Right) .AND. ((x <  x_Left) .OR. (x >  x_Right))) .OR.  &
      ((x_Left >  x_Right) .AND. ((x >  x_Left) .OR. (x <  x_Right)))) then
    write (6 ,100) Map_F, x, x_Left, x_Right
    write (99,100) Map_F, x, x_Left, x_Right

    STOP
  endif

 100 format (1x,'Inverse_Map_Linear_F : x is out of range [x_Left, x_Right]'/ &
             1x,'   Map_F, x, x_Left, x_Right = ',i3, 1x, 3(1x,f14.7))
#endif

  call Map_Physical_2_Computational (1, x_Left, x_Right, x, &
                                          ZERO,   TWOPI, y)

  Inverse_Map_Linear_F = y

  END FUNCTION Inverse_Map_Linear_F 

END MODULE F_Linear_Mapping

!
! ----------------------------------------------------------------------
!
MODULE Fourier_Mapping
  USE Grid_Auxillary_Subroutine

  USE Fourier_Angular_Mapping
  USE F_Linear_Mapping
  USE User_Defined_Mapping

implicit NONE

PUBLIC  ::       F_Mapping
PUBLIC  :: F_Forward_Mapping
PUBLIC  :: F_Inverse_Mapping

CONTAINS

  Subroutine F_Mapping                                                &
             (N, N0, N2, Map, Map_F,  alpha, beta, x_Left, x_Right,   &
              x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

  integer  :: N, N0, N2, Map, Map_F
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
      call Fourier_Angular_Map                                     &
           (N, N0, N2, Map, Map_F,  alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (99)
      call F_Linear_Map                                            &
           (N, N0, N2, Map, Map_F,  alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)

    CASE (:-1)
      call User_Defined_Map                                        &
           (N, N0, N2, Map, Map_F,  alpha, beta, x_Left, x_Right,  &
            x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)
  
#ifdef DEBUG
    CASE DEFAULT
      write (6 ,*) 'Unsupported Fourier mapping Map=',Map
      write (99,*) 'Unsupported Fourier mapping Map=',Map

      STOP
#endif

  END SELECT

  END Subroutine F_Mapping
!
! ---------------------------------------------------------------------
!
  REALTYPE FUNCTION F_Forward_Mapping (Operator, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Operator

  SELECT CASE (INT(Operator(_I_MAP)))
    CASE (0)
      F_Forward_Mapping = Map_Linear_F             (Operator, x)

    CASE (1)
      F_Forward_Mapping = Map_Angular            (Operator, x)

    CASE (99)
      F_Forward_Mapping = Map_Linear_F             (Operator, x)

    CASE (:-1)
      F_Forward_Mapping = Map_User_Defined       (Operator, x)

#ifdef DEBUG
    CASE DEFAULT
      write (6 ,*) 'Unsupported Fourier mapping =',Operator(_I_MAP)
      write (99,*) 'Unsupported Fourier mapping =',Operator(_I_MAP)

      STOP
#endif

  END SELECT

  END FUNCTION F_Forward_Mapping 
!
! ---------------------------------------------------------------------
!
  REALTYPE FUNCTION F_Inverse_Mapping (Operator, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: Operator

  SELECT CASE (INT(Operator(_I_MAP)))
    CASE (0)
      F_Inverse_Mapping = Inverse_Map_Linear_F           (Operator, x)

    CASE (1)
      F_Inverse_Mapping = Inverse_Map_Angular            (Operator, x)

    CASE (99)
      F_Inverse_Mapping = Inverse_Map_Linear_F           (Operator, x)

    CASE (:-1)
      F_Inverse_Mapping = Inverse_Map_User_Defined       (Operator, x)

#ifdef DEBUG
    CASE DEFAULT
      write (6 ,*) 'Unsupported Fourier mapping =',Operator(_I_MAP)
      write (99,*) 'Unsupported Fourier mapping =',Operator(_I_MAP)

      STOP
#endif

  END SELECT

  END FUNCTION F_Inverse_Mapping 

END MODULE Fourier_Mapping
