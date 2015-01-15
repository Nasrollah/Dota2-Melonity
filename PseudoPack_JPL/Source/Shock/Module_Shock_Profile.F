#include "pseudopack.h"

!==================================================================
!
! Function    : Setup_Shock_Profile
! Author      : Wai Sun Don
! Description : Specify Initial Conditions AND Update Shock Location
!                for One Dimensional Problems.
!                  Q(x,t) = q_Left     if x <  Shock_Location
!                         = q_Right    if x => Shock_Location
!               with additional perturbations in density for 
!               Shock_Profile_Type=1,2
!
!  Shock_Profile_Type    Type of Initial Conditions
!               0   Step Functions  (123, SOD, LAX, LeBlanc)  (DEFAULT)
!
!               1   Moving Shock-Density Wave Interaction
!                       rho = rho+Amplitude*sin(Wave_Number*pi*x)
!                     where
!                       rho is the density to the right of the shock.
!
!               2   Moving Shock-Entropy Wave Interaction
!                       rho = rho*exp(-Amplitude/P*COS(Wave_Number*pi*x)
!                     where
!                       rho is the density to the right of the shock.
!
!               3   Moving Shock-Entropy Wave Interaction
!                       rho = rho*exp(-Amplitude*SIN(Wave_Number*x))
!                     where
!                       rho is the density to the right of the shock.
!
!               4   Moving Shock-Entropy Wave Interaction
!                       rho = rho*exp(-Amplitude*SIN(Wave_Number*x)/Gamma)
!                     where
!                       rho is the density to the right of the shock.
!
!  x           Collocation grid of size (:)          (Input) 
!  q           Conservative variables of size (:,3)  (Output)
!
!  Shock_Location Shock Location      (Input/Output)
!
!  q_Left      Array of size 3 (Input)
!                Conservative Variables at the Inflow  (Left of the shock)
!  q_Right     Array of size 3 (Input)
!                Conservative Variables at the Outflow (Right of the shock)
!
!  Gamma       Specific Heat ratio (Input)
!  Wave_Number Wave Number of the perturbation for Shock_Profile_Type=1,2
!  Amplitude   Amplitude   of the perturbation for Shock_Profile_Type=1,2
!
! -----------------------------------------------------------------

MODULE Setup_Shock_Profile

implicit NONE

INTERFACE PS_Shock_Profile_Setup
  MODULE PROCEDURE Shock_Profile_1D_A
  MODULE PROCEDURE Shock_Profile_1D_B
END INTERFACE

PRIVATE

PUBLIC  :: PS_Shock_Profile_Setup

CONTAINS
!
! =========================================================================
!
  Subroutine Shock_Profile_1D_B (Shock_Profile_Type, Shock_Location,  &
                                 q_Left, q_Right, x, q)
                                 

  integer            :: Shock_Profile_Type
  REALTYPE           :: Shock_Location, Gamma, Wave_Number, Amplitude

  REALTYPE, dimension(3)   :: q_Left, q_Right
  REALTYPE, dimension(:)   :: x
  REALTYPE, dimension(:,:) :: q

  Gamma = 1.40d0 ; Wave_Number = 0 ; Amplitude = ZERO
  
  call Shock_Profile_1D_A (Shock_Profile_Type, Shock_Location,        &
                           q_Left, q_Right,                           &
                           Gamma, Wave_Number, Amplitude, x, q)

  END Subroutine Shock_Profile_1D_B 
!
! =========================================================================
!
  Subroutine Shock_Profile_1D_A (Shock_Profile_Type, Shock_Location,  &
                                 q_Left, q_Right,                     &
                                 Gamma, Wave_Number, Amplitude, x, q)
                                 

  integer            :: Shock_Profile_Type
  REALTYPE           :: Shock_Location, Gamma, Wave_Number, Amplitude

  REALTYPE, dimension(3)   :: q_Left, q_Right
  REALTYPE, dimension(:)   :: x 
  REALTYPE, dimension(:,:) :: q

  integer :: i, N, N_q

  N = SIZE(x) ; N_q = SIZE(q, DIM=1)

#if defined (DEBUG)
  if (Shock_Profile_Type >  5) STOP 'ERROR Shock_Profile_1D : Type >  5'

  if (N /= N_q) STOP 'ERROR Shock_Profile_1D : Size mismatch between x and q'
#endif

  SELECT CASE (Shock_Profile_Type)
    CASE DEFAULT
      call Shock_Profile_1D_0 

    CASE (1)
      call Shock_Profile_1D_1 

    CASE (2)
      call Shock_Profile_1D_2 

    CASE (3)
      call Shock_Profile_1D_3 

    CASE (4)
      call Shock_Profile_1D_4 

    CASE (5)
      call Shock_Profile_1D_5 

  END SELECT

! call Update_Shock_Location_1D 

  CONTAINS
!
!==================================================================
!
    Subroutine Shock_Profile_1D_0 

    do i = 1,N
      if (x(i) <  Shock_Location) then
        q(i,:) = q_Left
      else
        q(i,:) = q_Right
      endif
    enddo
   
    END Subroutine Shock_Profile_1D_0 
!
!==================================================================
!
    Subroutine Shock_Profile_1D_1 

    REALTYPE :: r, u, P, Perturbation

    do i = 1,N
      if (x(i) <  Shock_Location) then
        q(i,:) = q_Left
      else
        r = q_Right(1)
        u = q_Right(2)/r
        P = (Gamma-ONE)*(q_Right(3)-HALF*r*u*u)

#if 0
        Perturbation = Amplitude*SIN(Wave_Number*PI*x(i))
#else
        Perturbation = Amplitude*SIN(Wave_Number*   x(i))
#endif
  
        q(i,1) = r+Perturbation
        q(i,2) = u*q(i,1)
        q(i,3) = P/(Gamma-ONE)+HALF*q(i,2)*q(i,2)/q(i,1)
      endif
    enddo
 
    END Subroutine Shock_Profile_1D_1 
!
!==================================================================
!
    Subroutine Shock_Profile_1D_2 

    REALTYPE :: r, u, P, Perturbation

    do i = 1,N
      if (x(i) <  Shock_Location) then
        q(i,:) = q_Left
      else
        r = q_Right(1)
        u = q_Right(2)/r
        P = (Gamma-ONE)*(q_Right(3)-HALF*r*u*u)

        Perturbation = EXP(-Amplitude/P*SIN(Wave_Number*PI*x(i))) 
  
        q(i,1) = r*Perturbation
        q(i,2) = u*q(i,1)
        q(i,3) = P/(Gamma-ONE)+HALF*q(i,2)*q(i,2)/q(i,1)
      endif
    enddo
 
    END Subroutine Shock_Profile_1D_2 
!
!==================================================================
!
    Subroutine Shock_Profile_1D_3 

    REALTYPE :: r, u, P, Perturbation

    do i = 1,N
      if (x(i) <  Shock_Location) then
        q(i,:) = q_Left
      else
        r = q_Right(1)
        u = q_Right(2)/r
        P = (Gamma-ONE)*(q_Right(3)-HALF*r*u*u)

        Perturbation = EXP(-Amplitude*SIN(Wave_Number*x(i)))
  
        q(i,1) = r*Perturbation
        q(i,2) = u*q(i,1)
        q(i,3) = P/(Gamma-ONE)+HALF*q(i,2)*q(i,2)/q(i,1)
      endif
    enddo
 
    END Subroutine Shock_Profile_1D_3 
!
!==================================================================
!
    Subroutine Shock_Profile_1D_4

    REALTYPE :: r, u, P, Perturbation

    do i = 1,N
      if (x(i) <  Shock_Location) then
        q(i,:) = q_Left
      else
        r = q_Right(1)
        u = q_Right(2)/r
        P = (Gamma-ONE)*(q_Right(3)-HALF*r*u*u)

        Perturbation = EXP(-Amplitude*SIN(Wave_Number*x(i))/Gamma)
 
        q(i,1) = r*Perturbation
        q(i,2) = u*q(i,1)
        q(i,3) = P/(Gamma-ONE)+HALF*q(i,2)*q(i,2)/q(i,1)
      endif
    enddo

    END Subroutine Shock_Profile_1D_4
!
!==================================================================
!
    Subroutine Shock_Profile_1D_5 

    REALTYPE :: r, u, P

    REALTYPE, parameter :: Shock_Location_1 = 0.10d0
    REALTYPE, parameter :: Shock_Location_2 = 0.90d0

    do i = 1,N
      if  (x(i) <  Shock_Location_1) q(i,:) = q_Left

      if ((x(i) >= Shock_Location_1) .AND. (x(i) <= Shock_Location_2)) then
        r = ONE
        u = ZERO
        P = ONE/HUNDRED

        q(i,1) = r
        q(i,2) = r*u
        q(i,3) = P/(Gamma-ONE)+HALF*q(i,2)*q(i,2)/q(i,1)
      endif

      if  (x(i) >  Shock_Location_2) q(i,:) = q_Right
    enddo
   
    END Subroutine Shock_Profile_1D_5 

#if 0
!
!==================================================================
!
    Subroutine Update_Shock_Location_1D 
  
#if defined (DEBUG)
    if (x(N) <  x(1)) then
      if (Shock_Location <  x(N)) STOP 'ERROR USL1D: Shock_Location < x(N)'
      if (Shock_Location >  x(1)) STOP 'ERROR USL1D: Shock_Location > x(1)'
    endif
    if (x(N) >  x(1)) then
      if (Shock_Location <  x(1)) STOP 'ERROR USL1D: Shock_Location < x(1)'
      if (Shock_Location >  x(N)) STOP 'ERROR USL1D: Shock_Location > x(N)'
    endif
#endif

    do i = 1,N
      if ((Shock_Location >  x(i+1)) .AND. (Shock_Location <= x(i))) then
        Shock_Location = HALF*(x(i)+x(i+1))
        EXIT
      endif
    enddo

    END Subroutine Update_Shock_Location_1D 
#endif

  END Subroutine Shock_Profile_1D_A

END MODULE Setup_Shock_Profile
