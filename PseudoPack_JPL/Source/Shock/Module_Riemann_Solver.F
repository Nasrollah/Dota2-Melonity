#include "pseudopack.h"

! ----------------------------------------------------------------------
! FUNCTION NAME: PS_Riemann_Solver
! AUTHOR       : Wai-Sun Don
! DESCRIPTION  : Exact solution of the Riemann problem of 1D Euler at s
!                See Gottlieb & Groth, JCP 78, p. 437-458
!
! s = (x-x_s)/t        Speed at location x relative to the x_s, the location of
!                        the known discontinuity and time interval (Input).
! Q_Left               Array of size 3 (Input)
!                        containing (Density, Velocity, Pressure) or 
!                                   (Density, Moemntum, Total Energy) at Left.
! Q_Right              Array of size 3 (Input)
!                        containing (Density, Velocity, Pressure) or 
!                                   (Density, Moemntum, Total Energy) at Right.
! Q_at_s               Array of size 3 (Output)
!                        containing (Density, Velocity, Pressure) or 
!                                   (Density, Moemntum, Total Energy) at s.
! Variable_Type        Character of size 1 (Input, OPTIONAL)
!                        Conservative ('C') or Primitive ('P') Variable 
!                        Default : 'C'
! Specific_Heat_Ratio  specific heat ratio (gamma) (Input, OPTIONAL)
!                        Default : 1.4
! Tolerance            Epsilon_Tol Level for the Newton's method. (Input, OPTIONAL)
! Max_Iteration        Maximum number of iteration for the Newton's method.
!                                                             (Input, OPTIONAL) 
! Shock_Case           Character of size 3 (Output, OPTIONAL)
!                        Which Case is the being used?
!                        SCR, SCS, RCR, RCS
! 
! Small_Zero           Small number added to avoid division by zero at
!                        the contact discontinuity where the shock speed is
!                        computed as (d1*u1-d2*u2)/(d1-d2+Small_Zero) and 
!                        d1=d2, u1=u2=0 (Input, OPTIONAL)
!
!*********************************************************************

MODULE Riemann_Solver

implicit NONE

PRIVATE
PUBLIC  :: PS_Riemann_Solver

CONTAINS

!-----------------------------------------------------------------

  Subroutine PS_Riemann_Solver (s, Q_Left, Q_Right, Q_at_s,                    &
                                Variable_Type, Specific_Heat_Ratio, Tolerance, &
                                Max_Iteration, Shock_Case, Small_Zero)

  REALTYPE                   :: s
  REALTYPE, dimension(3)     :: Q_Left, Q_Right, Q_at_s

  character(LEN=1), OPTIONAL :: Variable_Type
  REALTYPE        , OPTIONAL :: Specific_Heat_Ratio
  REALTYPE        , OPTIONAL :: Tolerance
  integer         , OPTIONAL :: Max_Iteration
  character(LEN=3), OPTIONAL :: Shock_Case
  REALTYPE        , OPTIONAL :: Small_Zero

  REALTYPE         :: Epsilon_Tol
  integer          :: Number_Of_Iteration
  character(LEN=1) :: Var_type
  character(LEN=3) :: Case

  REALTYPE :: gamma, gp1, gm1, g1, g2, g3, g4, g5
  REALTYPE :: Vacuum
  REALTYPE :: d1, u1, p1, c1
  REALTYPE :: d2, u2, p2, c2, s2
  REALTYPE :: d3, u3, p3, c3, s3
  REALTYPE :: d4, u4, p4, c4
  REALTYPE :: y, z
  REALTYPE :: Density, Velocity, Pressure
  REALTYPE :: Contact_Zero

  Var_Type            = 'C'
  gamma               = 1.40d0
  Epsilon_Tol         = 1.0d-3 
  Number_Of_Iteration = 100
  Contact_Zero        = ZERO                                 

  if (PRESENT(Variable_Type      ))            Var_Type = Variable_Type
  if (PRESENT(Tolerance          ))         Epsilon_Tol = ABS(Tolerance)
  if (PRESENT(Max_Iteration      )) Number_Of_Iteration = Max_Iteration
  if (PRESENT(Specific_Heat_Ratio))               gamma = Specific_Heat_Ratio
  if (PRESENT(Small_Zero         ))        Contact_Zero = Small_Zero

  gp1     = gamma+ONE
  gm1     = gamma-ONE
  g1      = gp1/gm1         !  (g+1)/(g-1)
  g2      = TWO/gp1         !      2/(g+1)
  g3      = TWO/gm1         !      2/(g-1)
  g4      = ONE/(gamma*g3)  !  (g-1)/(2*g)
  g5      = g3*SQRT(g4)     !  2*sqrt((g-1)/(2*g))/(g-1)

  SELECT CASE (Var_Type)
    CASE DEFAULT
      d1 = Q_Left(1)                     ; d4 = Q_Right(1)
      u1 = Q_Left(2)/d1                  ; u4 = Q_Right(2)/d4
      p1 = gm1*(Q_Left(3)-HALF*d1*u1*u1) ; p4 = gm1*(Q_Right(3)-HALF*d4*u4*u4)
      c1 = SQRT(gamma*p1/d1)             ; c4 = SQRT(gamma*p4/d4)

    CASE ('P','p')
      d1 = Q_Left(1)                     ; d4 = Q_Right(1)
      u1 = Q_Left(2)                     ; u4 = Q_Right(2)
      p1 = Q_Left(3)                     ; p4 = Q_Right(3)
      c1 = SQRT(gamma*p1/d1)             ; c4 = SQRT(gamma*p4/d4)

  END SELECT

!  Test If a non-vacuum solution exists.
  Vacuum = g3*(c4+c1)-(u4-u1)
  if (Vacuum <  Epsilon_Tol) Stop 'VACUUM!  Stopped....'

!  Determine which case we are dealing with:
!  S-Shock, C-Contact, R-Rarefaction

  call Riemann_Shock_Case 

  call Newton_Solver 

  SELECT CASE (Case)
    CASE ('RCR')              !    1-Rarefaction   2-Contact  3-Rarefaction
      call Rarefaction_Contact_Rarefaction

    CASE ('SCR')              !    1-Shock         2-Contact  3-Rarefaction
      call Shock_Contact_Rarefaction

    CASE ('RCS')              !    1-Rarefaction   2-Contact  3-Shock
      call Rarefaction_Contact_Shock

    CASE ('SCS')              !    1-Shock         2-Contact  3-Shock
      call Shock_Contact_Shock

  END SELECT

  SELECT CASE (Var_Type)
    CASE DEFAULT
      Q_at_s(1) = Density
      Q_at_s(2) = Velocity*Density
      Q_at_s(3) = Pressure/gm1+HALF*Density*Velocity*Velocity

    CASE ('P','p')
      Q_at_s(1) = Density
      Q_at_s(2) = Velocity
      Q_at_s(3) = Pressure

  END SELECT

  if (PRESENT(Shock_Case)) Shock_Case = Case

  CONTAINS
 
!---------------------------------------------------
    Subroutine Riemann_Shock_Case 

    REALTYPE :: P_41, P_14
    REALTYPE :: u_SCN, u_NCR, u_NCS, u_RCN

    P_41 = p4/p1
    P_14 = p1/p4

    if (p4 >  p1) then
      u_SCN = u1 + c1*g5*(ONE-P_41    )/SQRT(ONE+g1*P_41)
      u_NCR = u1 + c4*g3*(ONE-P_14**g4)

      Case = 'SCR'
      if (u4 <  u_SCN) Case = 'SCS'
      if (u4 >  u_NCR) Case = 'RCR'
    else
      u_NCS = u1 + c4*g5*(ONE-P_14    )/SQRT(ONE+g1*P_14)
      u_RCN = u1 + c1*g3*(ONE-P_41**g4)

      Case = 'RCS'
      if (u4 <  u_NCS) Case = 'SCS'
      if (u4 >  u_RCN) Case = 'RCR'
    endif

    END Subroutine Riemann_Shock_Case 
 
!---------------------------------------------------
! Using Newton methods to solve for x of the nonlinear equation
!        H(x)+c4/c1*H(x+B)-C = 0
!---------------------------------------------------
    Subroutine Newton_Solver 

    integer  :: Iteration
    REALTYPE :: B, C, D, H, dH, dy
    REALTYPE :: Arg, Arg1, Arg2, Arg3

    B = LOG(p4/p1)
    C = (u4-u1)/c1
    D = c4/c1

    y = ZERO
    do Iteration = 1,Number_Of_Iteration

      if (y >  ZERO) then
        Arg  = EXP(-g4*y)

        H    = g3*(ONE-Arg)
        dH   = Arg/gamma
      else
        Arg1 = EXP(-y)
        Arg2 = ONE+g1*Arg1
        Arg3 = SQRT(Arg2)

        H    =      g5*(   ONE-Arg1)/          Arg3
        dH   = Arg1*g5*(g1+ONE+Arg2)/(TWO*Arg2*Arg3)
      endif

      z = y + B
      if (z >  ZERO) then
        Arg  = EXP(-g4*z)

        H    = H  + D*g3*(ONE-Arg)                         - C
        dH   = dH + D*Arg/gamma
      else
        Arg1 = EXP(-z)
        Arg2 = ONE+g1*Arg1
        Arg3 = SQRT(Arg2)

        H    = H  + D*     g5*(   ONE-Arg1)/          Arg3 - C
        dH   = dH + D*Arg1*g5*(g1+ONE+Arg2)/(TWO*Arg2*Arg3)
      endif

      dy = H/dH

      if (ABS(dy) <  Epsilon_Tol) RETURN

      y = y - dy
    enddo

    write (6,*) '              s : ',s
    write (6,*) 'Iteration Error : ',dy

    Stop 'No Convergent !  Stopped ....'

    END Subroutine Newton_Solver

!---------------------------------------------------
    Subroutine Rarefaction_Contact_Rarefaction

    REALTYPE :: Arg, Arg3, Arg4, Arg5

    Arg3 = EXP(-z/gamma)
    Arg4 = EXP(-y/gamma)
    Arg5 = (d4/d1)*Arg3/Arg4

    d2 = d1*Arg4                   
    u2 = u1 + c1*g3*(ONE-EXP(-g4*y)) 
    p2 = p1*EXP(-y)                
    c2 = SQRT(gamma*p2/d2)         

    d3 = d2*Arg5 
    u3 = u2
    p3 = p2
    c3 = SQRT(gamma*p3/d3)

    if (s <  (u1-c1)) then
      Density = d1 ; Velocity = u1 ; Pressure = p1

      else if (s <  (u2-c2)) then
             Arg = LOG(g2+(u1-s)/c1/g1)/g4

        Density  = d1*EXP(Arg/gamma)
        Velocity = u1 + c1*g3*(ONE-EXP(g4*Arg))
        Pressure = p1*EXP(Arg)

        else if (s <  u2) then
          Density = d2 ; Velocity = u2 ; Pressure = p2

          else if (s <  (u3+c3)) then
            Density = d3 ; Velocity = u3 ; Pressure = p3

            else if (s <  (u4+c4)) then
                   Arg = LOG(g2+(s-u3)/c3/g1)/g4

              Density  = d3*EXP(Arg/gamma)
              Velocity = u3 + c3*g3*(EXP(g4*Arg)-ONE)
              Pressure = p3*EXP(Arg)

              else
                Density = d4 ; Velocity = u4 ; Pressure = p4
    endif

    END Subroutine Rarefaction_Contact_Rarefaction

!---------------------------------------------------
    Subroutine Shock_Contact_Rarefaction

    REALTYPE :: Arg, Arg0, Arg1, Arg3, Arg4, Arg5

    Arg0 = EXP(y)
    Arg1 = ONE/Arg0

    Arg3 = EXP(-z/gamma)
    Arg4 = (g1+Arg0)/(ONE+g1*Arg0)
    Arg5 = (d4/d1)*Arg3/Arg4

    d2 = d1*Arg4                               
    u2 = u1 + c1*g5*(ONE-Arg1)/SQRT(ONE+g1*Arg1) 
    p2 = p1*Arg1                               
    c2 = SQRT(gamma*p2/d2)                     
!   s2 = (d2*u2-d1*u1)/(d2-d1)
    s2 = (d2*u2-d1*u1)/(d2-d1 + Contact_Zero)

    d3 = d2*Arg5
    u3 = u2
    p3 = p2
    c3 = SQRT(gamma*p3/d3)

    if (s <  s2) then
      Density = d1 ; Velocity = u1 ; Pressure = p1

      else if (s <  u2) then
        Density = d2 ; Velocity = u2 ; Pressure = p2

        else if (s <  (u3+c3)) then
          Density = d3 ; Velocity = u3 ; Pressure = p3

          else if (s <  (u4+c4)) then
                 Arg = LOG(g2+(s-u3)/c3/g1)/g4

            Density  = d3*EXP(Arg/gamma)
            Velocity = u3 + c3*g3*(EXP(g4*Arg)-ONE)
            Pressure = p3*EXP(Arg)

            else
              Density = d4 ; Velocity = u4 ; Pressure = p4
    endif

    END Subroutine Shock_Contact_Rarefaction

!---------------------------------------------------
    Subroutine Rarefaction_Contact_Shock

    REALTYPE :: Arg, Arg0, Arg3, Arg4, Arg5

    Arg0 = EXP(z)

    Arg3 = (g1+Arg0)/(ONE+g1*Arg0)
    Arg4 = EXP(-y/gamma)
    Arg5 = (d4/d1)*Arg3/Arg4

    d2 = d1*Arg4
    u2 = u1 + c1*g3*(ONE-EXP(-g4*y))
    p2 = p1*EXP(-y)
    c2 = SQRT(gamma*p2/d2)

    d3 = d2*Arg5
    u3 = u2
    p3 = p2
    c3 = SQRT(gamma*p3/d3)
!   s3 = (d4*u4-d3*u3)/(d4-d3)
    s3 = (d4*u4-d3*u3)/(d4-d3 + Contact_Zero)

    if (s <  (u1-c1)) then
      Density = d1 ; Velocity = u1 ; Pressure = p1

      else if (s <  (u2-c2)) then
             Arg = LOG(g2+(u1-s)/c1/g1)/g4

        Density  = d1*EXP(Arg/gamma)
        Velocity = u1 + c1*g3*(ONE-EXP(g4*Arg))
        Pressure = p1*EXP(Arg)

        else if (s <  u2) then
          Density = d2 ; Velocity = u2 ; Pressure = p2

          else if (s <  s3) then
            Density = d3 ; Velocity = u3 ; Pressure = p3

            else
              Density = d4 ; Velocity = u4 ; Pressure = p4
    endif

    END Subroutine Rarefaction_Contact_Shock

!---------------------------------------------------
    Subroutine Shock_Contact_Shock 

    REALTYPE :: Arg0, Arg1, Arg2, Arg3, Arg4, Arg5

    Arg0 = EXP(z)
    Arg1 = EXP(y)
    Arg2 = ONE/Arg1

    Arg3 = (g1+Arg0)/(ONE+g1*Arg0)
    Arg4 = (g1+Arg1)/(ONE+g1*Arg1)
    Arg5 = (d4/d1)*Arg3/Arg4

    d2 = d1*Arg1
    u2 = u1 + c1*g5*(ONE-Arg2)/SQRT(ONE+g1*Arg2)
    p2 = p1*Arg2
    c2 = SQRT(gamma*p2/d2)
!   s2 = (d2*u2-d1*u1)/(d2-d1)
    s2 = (d2*u2-d1*u1)/(d2-d1 + Contact_Zero)

    d3 = d2*Arg5
    u3 = u2
    p3 = p2
    c3 = SQRT(gamma*p3/d3)
!   s3 = (d4*u4-d3*u3)/(d4-d3)
    s3 = (d4*u4-d3*u3)/(d4-d3 + Contact_Zero)

    if (s <  s2) then
      Density = d1 ; Velocity = u1 ; Pressure = p1

      else if (s <  u2) then
        Density = d2 ; Velocity = u2 ; Pressure = p2

        else if (s <  s3) then
          Density = d3 ; Velocity = u3 ; Pressure = p3

          else
            Density = d4 ; Velocity = u4 ; Pressure = p4
    endif

    END Subroutine Shock_Contact_Shock

  END Subroutine PS_Riemann_Solver 

END MODULE Riemann_Solver
