#include "pseudopack.h"
 
!==================================================================
!
! Function    : Shock_Flow_Specification
! Author      : Wai-Sun Don
! Description : Specify the type of 1-D Shock Problems and 
!               Post and Pre Shock Variables
!
!  Shock_Type  Type of Shock definitions (Input)
!
!               -1   Custom Shock Tube Proglem
!
!                1   Shock Tube Problem - 123
!                      (Final_Time=1.0, x_Left=-5, x_Right=5, Gamma=1.4,
!                                                             Shock_Location=0)
!
!                2   Shock Tube Problem - SOD
!                      (Final_Time=2.0, x_Left=-5, x_Right=5, Gamma=1.4,
!                                                             Shock_Location=0)
!
!                3   Shock Tube Problem - LAX
!                      (Final_Time=1.3, x_Left=-5, x_Right=5, Gamma=1.4,
!                                                             Shock_Location=0)
!
!                4   Shock Tube Problem - LeBlanc
!                      (Final_Time=1.0d-4, x_Left=-10, x_Right=10, Gamma=1.4,
!                                                             Shock_Location=0)
!
!                5   Interacting Blast Wave (Gamma)
!
!                6   Moving Shock (Gamma, Mach)
!               -6   Moving Shock (Gamma, Mach, Rho_Right, P_Right)
!               -7   Moving Shock (Gamma, Mach, Rho_Right, T_Right, 
!                                         Molecular_Weight)
!
!                7   Fixed  Shock (Gamma, Mach)
!
!                8   Moving Shock (Gamma, P_Ratio)
!
! Variable_Type  Character of size 1 (Input)
!                  Conservative ('C') or Primitive ('P') Variable
!                  Default : 'C'
!
!  Mach        Shock Mach Number   (Input/Output)
!  P_Ratio     Pressure ratio      (Input)
!  Gamma       Specific Heat ratio (Input)
!  Shock_Speed Shock Speed         (Output)
!
!  Rho_Right / Rho_Left Density     to the Right / Left of the Shock 
!    U_Right / U_Left   Velocity    to the Right / Left of the Shock 
!    P_Right / P_Left   Velocity    to the Right / Left of the Shock 
!    C_Right / C_Left   Sound Speed to the Right / Left of the Shock 
!              ( Input)  if Shock_Type < 0
!              (Output)  if Shock_Type > 0
!
!  q_Left    Array of size 3 (Output)
!                  Left  stage of the shock 
!  q_Right   Array of size 3 (Output)
!                  Right stage of the shock 
!
! -----------------------------------------------------------------
 
MODULE Shock_Flow_Specification

implicit NONE

INTERFACE PS_Shock_Specification
  MODULE PROCEDURE Shock_Flow_1D_A
  MODULE PROCEDURE Shock_Flow_1D_B
  MODULE PROCEDURE Shock_Flow_1D_C
  MODULE PROCEDURE Shock_Flow_1D_D
END INTERFACE

REALTYPE, parameter    :: R_Universal = 8.314410d7        ! erg/K

PRIVATE

PUBLIC  :: PS_Shock_Specification
PUBLIC  :: PS_Convert_Variable_Type
PUBLIC  :: PS_Conservative_To_Primitive 
PUBLIC  :: PS_Primitive_To_Conservative

CONTAINS
!
!===================================================================
!
  Subroutine Shock_Flow_1D_B (Shock_Type, Variable_Type, q_Left, q_Right, &
                              Mach, P_Ratio, Gamma)

  integer                :: Shock_Type
  character(LEN=1)       :: Variable_Type
  REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
  REALTYPE :: Final_Time, x_Left, x_Right, Shock_Location
  REALTYPE :: Rho_Left , U_Left , P_Left , T_Left
  REALTYPE :: Rho_Right, U_Right, P_Right, T_Right
  REALTYPE :: C_Left, C_Right
  REALTYPE, dimension(3) :: q_Left, q_Right

  REALTYPE :: MW = ONE

  call Shock_Flow_1D_A (Shock_Type, Variable_Type,                   &
                        Mach, P_Ratio, MW, Gamma, Shock_Speed,       &
                        Rho_Left , U_Left , P_Left , T_Left ,        &
                        Rho_Right, U_Right, P_Right, T_Right,        &
                        Final_Time, x_Left, x_Right, Shock_Location, &
                        C_Left, C_Right, q_Left, q_Right)

  END Subroutine Shock_Flow_1D_B
!
!===================================================================
!
  Subroutine Shock_Flow_1D_C (Shock_Type, Variable_Type, q_Left, q_Right,  &
                              Final_Time, x_Left, x_Right, Shock_Location, &
                              Mach, P_Ratio, Gamma)

  integer                :: Shock_Type
  character(LEN=1)       :: Variable_Type
  REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
  REALTYPE :: Final_Time, x_Left, x_Right, Shock_Location
  REALTYPE :: Rho_Left , U_Left , P_Left , T_Left
  REALTYPE :: Rho_Right, U_Right, P_Right, T_Right
  REALTYPE :: C_Left, C_Right
  REALTYPE, dimension(3) :: q_Left, q_Right

  REALTYPE :: MW = ONE

  call Shock_Flow_1D_A (Shock_Type, Variable_Type,                   &
                        Mach, P_Ratio, MW, Gamma, Shock_Speed,       &
                        Rho_Left , U_Left , P_Left , T_Left ,        &
                        Rho_Right, U_Right, P_Right, T_Right,        &
                        Final_Time, x_Left, x_Right, Shock_Location, &
                        C_Left, C_Right, q_Left, q_Right)

  END Subroutine Shock_Flow_1D_C
!
!===================================================================
!
  Subroutine Shock_Flow_1D_D (Shock_Type, Variable_Type,                   &
                              Mach, P_Ratio, Gamma, Shock_Speed,           &
                              Rho_Left , U_Left , P_Left ,                 &
                              Rho_Right, U_Right, P_Right,                 &
                              Final_Time, x_Left, x_Right, Shock_Location, &
                              C_Left, C_Right, q_Left, q_Right)

  integer                :: Shock_Type
  character(LEN=1)       :: Variable_Type
  REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
  REALTYPE :: Final_Time, x_Left, x_Right, Shock_Location
  REALTYPE :: Rho_Left , U_Left , P_Left , T_Left
  REALTYPE :: Rho_Right, U_Right, P_Right, T_Right
  REALTYPE :: C_Left, C_Right
  REALTYPE, dimension(3) :: q_Left, q_Right

  REALTYPE :: MW = ONE

  call Shock_Flow_1D_A (Shock_Type, Variable_Type,                   &
                        Mach, P_Ratio, MW, Gamma, Shock_Speed,       &
                        Rho_Left , U_Left , P_Left , T_Left ,        &
                        Rho_Right, U_Right, P_Right, T_Right,        &
                        Final_Time, x_Left, x_Right, Shock_Location, &
                        C_Left, C_Right, q_Left, q_Right)

  END Subroutine Shock_Flow_1D_D
!
!===================================================================
!
  Subroutine Shock_Flow_1D_A (Shock_Type, Variable_Type,                   &
                              Mach, P_Ratio, Molecular_Weight, Gamma,      &
                              Shock_Speed,                                 &
                              Rho_Left , U_Left , P_Left , T_Left ,        &
                              Rho_Right, U_Right, P_Right, T_Right,        &
                              Final_Time, x_Left, x_Right, Shock_Location, &
                              C_Left, C_Right, q_Left, q_Right) 

  integer                :: Shock_Type
  character(LEN=1)       :: Variable_Type
  REALTYPE :: Mach, P_Ratio, Molecular_Weight, Gamma, Shock_Speed
  REALTYPE :: Final_Time, x_Left, x_Right, Shock_Location
  REALTYPE :: Rho_Left , U_Left , P_Left , T_Left
  REALTYPE :: Rho_Right, U_Right, P_Right, T_Right
  REALTYPE :: C_Left, C_Right
  REALTYPE, dimension(3) :: q_Left, q_Right
 
  REALTYPE :: gm1, gp1

#if defined (DEBUG)
  if (Shock_Type < -7) Stop 'ERROR Type_Flow : Shock_Type < -7'
  if (Shock_Type >  8) Stop 'ERROR Type_Flow : Shock_Type >  8'
#endif
 
#if 0
  write (6,*)                 Shock_Type, Variable_Type,                   &
                              Mach, P_Ratio, Molecular_Weight, Gamma,      &
                              Shock_Speed,                                 &
                              Rho_Left , U_Left , P_Left , T_Left ,        &
                              Rho_Right, U_Right, P_Right, T_Right,        &
                              Final_Time, x_Left, x_Right, Shock_Location, &
                              C_Left, C_Right, q_Left, q_Right
#endif

#if 0
 if (Shock_Type == -7) call Type_Custom_Moving_Shock_2
 if (Shock_Type == -6) call Type_Custom_Moving_Shock_1
 if (Shock_Type == -5) call Type_Custom_Shock_Tube
 if (Shock_Type == -4) call Type_Custom_Shock_Tube
 if (Shock_Type == -3) call Type_Custom_Shock_Tube
 if (Shock_Type == -2) call Type_Custom_Shock_Tube
 if (Shock_Type == -1) call Type_Custom_Shock_Tube
 if (Shock_Type ==  1) call Type_123
 if (Shock_Type ==  2) call Type_SOD
 if (Shock_Type ==  3) call Type_LAX
 if (Shock_Type ==  4) call Type_LeBlanc
 if (Shock_Type ==  5) call Type_Interacting_Blast_Wave
 if (Shock_Type ==  6) call Type_Moving
 if (Shock_Type ==  7) call Type_Fixed
 if (Shock_Type ==  8) call Type_P_Ratio

 if (Shock_Type ==  6) call Type_Moving
#endif

#if 1
 SELECT CASE (Shock_Type)
   CASE (-7)
     call Type_Custom_Moving_Shock_2 (Molecular_Weight,                      &
                              Mach, P_Ratio, Gamma, Shock_Speed,             &
                              Rho_Right, U_Right, P_Right, C_Right, T_Right, &
                              Rho_Left , U_Left , P_Left , C_Left , T_Left )

   CASE (-6)
     call Type_Custom_Moving_Shock_1 (                                       &
                              Mach, P_Ratio, Gamma, Shock_Speed,             &
                              Rho_Right, U_Right, P_Right, C_Right, T_Right, &
                              Rho_Left , U_Left , P_Left , C_Left , T_Left )

   CASE (-5:-1)
     call Type_Custom_Shock_Tube     (                                       &
                                       Mach, P_Ratio, Gamma, Shock_Speed,    &
                                       Rho_Right, U_Right, P_Right, C_Right, &
                                       Rho_Left , U_Left , P_Left , C_Left )
   CASE (1)
     call Type_123 (Mach, P_Ratio, Gamma,                          &
                    Shock_Speed, Shock_Location, Final_Time,       &
                    x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                    x_Left , Rho_Left , U_Left , P_Left , C_Left )

   CASE (2)
     call Type_SOD (Mach, P_Ratio, Gamma,                          &
                    Shock_Speed, Shock_Location, Final_Time,       &
                    x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                    x_Left , Rho_Left , U_Left , P_Left , C_Left )

   CASE (3)
     call Type_LAX (Mach, P_Ratio, Gamma,                          &
                    Shock_Speed, Shock_Location, Final_Time,       &
                    x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                    x_Left , Rho_Left , U_Left , P_Left , C_Left )

   CASE (4)
     call Type_LeBlanc (Mach, P_Ratio, Gamma,                          &
                        Shock_Speed, Shock_Location, Final_Time,       &
                         x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                         x_Left , Rho_Left , U_Left , P_Left , C_Left )

   CASE (5)
     call Type_Interacting_Blast_Wave (                                     &
                             Mach, P_Ratio, Gamma,                          &
                             Shock_Speed, Shock_Location, Final_Time,       &
                             x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                             x_Left , Rho_Left , U_Left , P_Left , C_Left )

   CASE (6)
     call Type_Moving (Mach, P_Ratio, Gamma, Shock_Speed,             &
                       x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                       x_Left , Rho_Left , U_Left , P_Left , C_Left )

   CASE (7)
     call Type_Fixed  (Mach, P_Ratio, Gamma, Shock_Speed,             &
                       x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                       x_Left , Rho_Left , U_Left , P_Left , C_Left )
 

   CASE (8)
     call Type_P_Ratio (Mach, P_Ratio, Gamma, Shock_Speed,             &
                        x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                        x_Left , Rho_Left , U_Left , P_Left , C_Left )

 END SELECT
#endif
!
!  ** Define Conservative Variables **
!
  SELECT CASE (Variable_Type)
    CASE DEFAULT
      q_Left(1)  = Rho_Left
      q_Left(2)  = Rho_Left *U_Left
      q_Left(3)  = P_Left /(Gamma-ONE) + HALF*Rho_Left *U_Left **2

      q_Right(1) = Rho_Right
      q_Right(2) = Rho_Right*U_Right
      q_Right(3) = P_Right/(Gamma-ONE) + HALF*Rho_Right*U_Right**2

    CASE ('P','p')
      q_Left(1)  = Rho_Left
      q_Left(2)  =   U_Left
      q_Left(3)  =   P_Left

      q_Right(1) = Rho_Right
      q_Right(2) =   U_Right
      q_Right(3) =   P_Right

  END SELECT
  
  END Subroutine Shock_Flow_1D_A 
!
!==================================================================
! Description : Custom Shock Tube Problem
!
!  Gamma    specific heat ratio (Input)
! -----------------------------------------------------------------
    Subroutine Type_Custom_Shock_Tube     (                                  &
                                       Mach, P_Ratio, Gamma, Shock_Speed,    &
                                       Rho_Right, U_Right, P_Right, C_Right, &
                                       Rho_Left , U_Left , P_Left , C_Left )

    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: Rho_Left , U_Left , P_Left , C_Left 
    REALTYPE :: Rho_Right, U_Right, P_Right, C_Right
 
      C_Right = SQRT(Gamma*P_Right/Rho_Right)
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)
 
        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    if (Rho_Left == Rho_Right) then
      Shock_Speed = ZERO
    else
      Shock_Speed = (Rho_Left*U_Left-Rho_Right*U_Right)/(Rho_Left-Rho_Right)
    endif
           Mach = Shock_Speed/C_Right
  
    END Subroutine Type_Custom_Shock_Tube
!
!==================================================================
!  Description : Moving Normal shock condition  (R-H relation) 
!
!  Mach      Shock Mach number   (Input)
!  Gamma     specific heat ratio (Input)
!  Rho_Right Pre Shock Density   (Input)
!    P_Right Pre Shock Pressure  (Input)
! -----------------------------------------------------------------
    Subroutine Type_Custom_Moving_Shock_1 (                                  &
                              Mach, P_Ratio, Gamma, Shock_Speed,             &
                              Rho_Right, U_Right, P_Right, C_Right, T_Right, &
                              Rho_Left , U_Left , P_Left , C_Left , T_Left )

    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: Rho_Left , U_Left , P_Left , C_Left , T_Left
    REALTYPE :: Rho_Right, U_Right, P_Right, C_Right, T_Right
 
    REALTYPE :: gm1, gp1

    gm1 = Gamma-ONE ; gp1 = Gamma+ONE
     
      C_Right = SQRT(Gamma*P_Right/Rho_Right)
      U_Right = Mach*C_Right
  
      P_Left  = P_Right*(TWO*Gamma*Mach*Mach-gm1)/gp1
    Rho_Left  = Rho_Right*(gp1*Mach*Mach)/(gm1*Mach*Mach+TWO)
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)
      U_Left  = U_Right*Rho_Right/Rho_Left
 
        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    Shock_Speed = Mach*C_Right
        U_Right = Shock_Speed-U_Right
        U_Left  = Shock_Speed-U_Left
 
    END Subroutine Type_Custom_Moving_Shock_1
!
!==================================================================
!  Description : Moving Normal shock condition  (R-H relation)
!
!  Mach      Shock Mach number     (Input)
!  Gamma     specific heat ratio   (Input)
!  Rho_Right Pre Shock Density     (Input)
!    T_Right Pre Shock Temperature (Input)
!            Molecular Weight      (Input)
! -----------------------------------------------------------------
    Subroutine Type_Custom_Moving_Shock_2 (Molecular_Weight,                 &
                              Mach, P_Ratio, Gamma, Shock_Speed,             &
                              Rho_Right, U_Right, P_Right, C_Right, T_Right, &
                              Rho_Left , U_Left , P_Left , C_Left , T_Left )

    REALTYPE :: Molecular_Weight
    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: Rho_Left , U_Left , P_Left , C_Left , T_Left
    REALTYPE :: Rho_Right, U_Right, P_Right, C_Right, T_Right

    REALTYPE :: R
    REALTYPE :: gm1, gp1

    gm1 = Gamma-ONE ; gp1 = Gamma+ONE

      R       = R_Universal/Molecular_Weight

      T_Right = T_Right                            ! K
    Rho_Right = Rho_Right                          ! g/cm^3
      P_Right = R*Rho_Right*T_Right                ! g/(cm*s^2) = dyne/cm^2
      C_Right = SQRT(Gamma*P_Right/Rho_Right)      ! cm/s
      U_Right = Mach*C_Right                       ! cm/s

      P_Left  = P_Right*(TWO*Gamma*Mach*Mach-gm1)/gp1          ! g/(cm*s^2)
    Rho_Left  = Rho_Right*(gp1*Mach*Mach)/(gm1*Mach*Mach+TWO)  ! g/cm^3
      T_Left  = P_Left/(R*Rho_Left)                            ! K
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)                    ! cm/s
      U_Left  = U_Right*Rho_Right/Rho_Left                     ! cm/s

        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    Shock_Speed = Mach*C_Right
        U_Right = Shock_Speed-U_Right
        U_Left  = Shock_Speed-U_Left

    END Subroutine Type_Custom_Moving_Shock_2
!
!==================================================================
! Description : 123 Shock Tube Problem
!
!  Gamma    specific heat ratio (Input)
! -----------------------------------------------------------------
    Subroutine Type_123 (Mach, P_Ratio, Gamma,                          &
                         Shock_Speed, Shock_Location, Final_Time,       &
                         x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                         x_Left , Rho_Left , U_Left , P_Left , C_Left )
 
    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: Shock_Location, Final_Time
    REALTYPE :: x_Left , Rho_Left , U_Left , P_Left , C_Left
    REALTYPE :: x_Right, Rho_Right, U_Right, P_Right, C_Right
 
    Rho_Right =  ONE
      U_Right =  TWO
      P_Right =  0.4d0
      C_Right = SQRT(Gamma*P_Right/Rho_Right)
 
    Rho_Left  =  ONE
      U_Left  = -TWO
      P_Left  =  P_Right
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)
 
        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    Shock_Speed = ZERO
           Mach = Shock_Speed/C_Right
  
        x_Left     = -FIVE 
        x_Right    =  FIVE 
    Shock_Location =  ZERO
    Final_Time     =  ONE

    END Subroutine Type_123     
!
!==================================================================
! Description : SOD Shock Tube Problem
!
!  Gamma    specific heat ratio (Input)
! -----------------------------------------------------------------
    Subroutine Type_SOD (Mach, P_Ratio, Gamma,                          &
                         Shock_Speed, Shock_Location, Final_Time,       &
                         x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                         x_Left , Rho_Left , U_Left , P_Left , C_Left )
 
    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: Shock_Location, Final_Time
    REALTYPE :: x_Left , Rho_Left , U_Left , P_Left , C_Left
    REALTYPE :: x_Right, Rho_Right, U_Right, P_Right, C_Right

    Rho_Right = ONE
      U_Right = ZERO
      P_Right = ONE
      C_Right = SQRT(Gamma*P_Right/Rho_Right)

    Rho_Left  = 0.125d0
      U_Left  = ZERO
      P_Left  = 0.1d0
      C_Left  = SQRT(Gamma*P_Left /Rho_Left )

        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    Shock_Speed = (Rho_Left*U_Left-Rho_Right*U_Right)/(Rho_Left-Rho_Right)
           Mach = Shock_Speed/C_Right
  
        x_Left     = -FIVE 
        x_Right    =  FIVE 
    Shock_Location =  ZERO
    Final_Time     =  TWO

    END Subroutine Type_SOD     
!
!==================================================================
! Description : LAX Shock Tube Problem
!
!  Gamma    specific heat ratio (Input)
! -----------------------------------------------------------------
    Subroutine Type_LAX (Mach, P_Ratio, Gamma,                          &
                         Shock_Speed, Shock_Location, Final_Time,       &
                         x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                         x_Left , Rho_Left , U_Left , P_Left , C_Left )
 
    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: Shock_Location, Final_Time
    REALTYPE :: x_Left , Rho_Left , U_Left , P_Left , C_Left
    REALTYPE :: x_Right, Rho_Right, U_Right, P_Right, C_Right

    Rho_Right = HALF
      U_Right = ZERO
      P_Right = 0.571d0
      C_Right = SQRT(Gamma*P_Right/Rho_Right)
 
    Rho_Left  = 0.445d0
      U_Left  = 0.698d0
      P_Left  = 3.528d0
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)
 
        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    Shock_Speed = (Rho_Left*U_Left-Rho_Right*U_Right)/(Rho_Left-Rho_Right)
           Mach = Shock_Speed/C_Right
 
        x_Left     = -FIVE 
        x_Right    =  FIVE 
    Shock_Location =  ZERO
    Final_Time     = 1.30d0 

    END Subroutine Type_LAX     
!
!==================================================================
! Description : LeBlanc Shock Tube Problem
!
!  Gamma    specific heat ratio (Input)
! -----------------------------------------------------------------
    Subroutine Type_LeBlanc (Mach, P_Ratio, Gamma,                          &
                             Shock_Speed, Shock_Location, Final_Time,       &
                             x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                             x_Left , Rho_Left , U_Left , P_Left , C_Left )
 
    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: Shock_Location, Final_Time
    REALTYPE :: x_Left , Rho_Left , U_Left , P_Left , C_Left
    REALTYPE :: x_Right, Rho_Right, U_Right, P_Right, C_Right

    Rho_Right = TWO
      U_Right = ZERO
      P_Right = 1.0d9
      C_Right = SQRT(Gamma*P_Right/Rho_Right)

    Rho_Left  = 0.001d0
      U_Left  = ZERO
      P_Left  = ONE
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)

        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    Shock_Speed = (Rho_Left*U_Left-Rho_Right*U_Right)/(Rho_Left-Rho_Right)
           Mach = Shock_Speed/C_Right

        x_Left     = -TEN
        x_Right    =  TEN
    Shock_Location =  ZERO
    Final_Time     =  1.0d-4

    END Subroutine Type_LeBlanc
!
!==================================================================
! Description : 1D Interacting Blast Wave Problem
!
!  Gamma    specific heat ratio (Input)
! -----------------------------------------------------------------
    Subroutine Type_Interacting_Blast_Wave (                                &
                             Mach, P_Ratio, Gamma,                          &
                             Shock_Speed, Shock_Location, Final_Time,       &
                             x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                             x_Left , Rho_Left , U_Left , P_Left , C_Left )
 
    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: Shock_Location, Final_Time
    REALTYPE :: x_Left , Rho_Left , U_Left , P_Left , C_Left
    REALTYPE :: x_Right, Rho_Right, U_Right, P_Right, C_Right
 
    Rho_Right = ONE
      U_Right = ZERO
      P_Right = HUNDRED
      C_Right = SQRT(Gamma*P_Right/Rho_Right)

    Rho_Left  = ONE
      U_Left  = ZERO
      P_Left  = THOUSAND
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)

        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    Shock_Speed = ZERO
           Mach = Shock_Speed/C_Right

        x_Left     =  ZERO
        x_Right    =  ONE
    Shock_Location =  0.10d0
    Final_Time     =  0.0380d0

    END Subroutine Type_Interacting_Blast_Wave
!
!==================================================================
!  Description : Moving Normal shock condition  (R-H relation) 
!
!  Mach     Shock Mach number  (Input)
!  Gamma    specific heat ratio (Input)
!
!     if Mach = 3 these are the Left and Right states
!         P_Right = 1          P_Left  = 10.3333333
!       Rho_Right = 1        Rho_Left  =  3.857143
!         U_Right = 0          U_Left  =  2.629369
! -----------------------------------------------------------------
    Subroutine Type_Moving (Mach, P_Ratio, Gamma, Shock_Speed,             &
                            x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                            x_Left , Rho_Left , U_Left , P_Left , C_Left )
 
    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: x_Left , Rho_Left , U_Left , P_Left , C_Left
    REALTYPE :: x_Right, Rho_Right, U_Right, P_Right, C_Right

    REALTYPE :: gm1, gp1
 
    gm1 = Gamma-ONE ; gp1 = Gamma+ONE
     
      P_Right = ONE
    Rho_Right = ONE
      C_Right = SQRT(Gamma*P_Right/Rho_Right)
      U_Right = Mach*C_Right
  
      P_Left  = P_Right*(TWO*Gamma*Mach*Mach-gm1)/gp1
    Rho_Left  = Rho_Right*(gp1*Mach*Mach)/(gm1*Mach*Mach+TWO)
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)
      U_Left  = U_Right*Rho_Right/Rho_Left
 
        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    Shock_Speed = Mach*C_Right
        U_Right = Shock_Speed-U_Right
        U_Left  = Shock_Speed-U_Left
 
    END Subroutine Type_Moving 
!
!==================================================================
!  Description : Fixed Normal shock condition  (R-H relation) 
!                Stationary Shock with u at the Right not equal zero.
!
!  Mach     Shock Mach number (Input)
!  Gamma    specific heat ratio (Input)
! -----------------------------------------------------------------
    Subroutine Type_Fixed  (Mach, P_Ratio, Gamma, Shock_Speed,             &
                            x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                            x_Left , Rho_Left , U_Left , P_Left , C_Left )
 
    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: x_Left , Rho_Left , U_Left , P_Left , C_Left
    REALTYPE :: x_Right, Rho_Right, U_Right, P_Right, C_Right

    REALTYPE :: gm1, gp1
 
    gm1 = Gamma-ONE ; gp1 = Gamma+ONE
     
      P_Right = ONE
    Rho_Right = ONE
      C_Right = SQRT(Gamma*P_Right/Rho_Right)
      U_Right = Mach*C_Right
  
      P_Left  =   P_Right*(TWO*Gamma*Mach*Mach-gm1)/gp1
    Rho_Left  = Rho_Right*(gp1*Mach*Mach)/(gm1*Mach*Mach+TWO)
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)
      U_Left  = U_Right*Rho_Right/Rho_Left
 
        P_Ratio = MAX(P_Left/P_Right, P_Right/P_Left)
    Shock_Speed = Mach*C_Left
        U_Right = Shock_Speed-U_Right
        U_Left  = Shock_Speed-U_Left
 
    END Subroutine Type_Fixed  
!
!==================================================================
!  Description : Moving Normal shock condition  (R-H relation) 
!
!  P_Ratio  Pressure ratio between the Left and Right of the shock (Input)
!  Gamma    specific heat ratio (Input)
! -----------------------------------------------------------------
    Subroutine Type_P_Ratio (Mach, P_Ratio, Gamma, Shock_Speed,             &
                             x_Right, Rho_Right, U_Right, P_Right, C_Right, &
                             x_Left , Rho_Left , U_Left , P_Left , C_Left )
 
    REALTYPE :: Mach, P_Ratio, Gamma, Shock_Speed
    REALTYPE :: x_Left , Rho_Left , U_Left , P_Left , C_Left
    REALTYPE :: x_Right, Rho_Right, U_Right, P_Right, C_Right

    REALTYPE :: gm1, gp1
    REALTYPE :: alpha, beta, tau

    gm1 = Gamma-ONE ; gp1 = Gamma+ONE
     
    alpha = gm1/TWO ; beta = gp1/gm1 ; tau = alpha/Gamma
 
    Rho_Right = ONE
      U_Right = ONE
      P_Right = ONE
      C_Right = SQRT(Gamma*P_Right/Rho_Right)
   
    Rho_Left  = Rho_Right*(ONE+beta*P_Ratio)/(beta+P_Ratio)
      U_Left  =   U_Right+C_Right*SQRT(tau)/alpha*(P_Ratio-ONE)/ &
                                  SQRT(ONE+beta*P_Ratio)
      P_Left  =   P_Right*P_Ratio
      C_Left  = SQRT(Gamma*P_Left/Rho_Left)
 
    Shock_Speed = (Rho_Left*U_Left-Rho_Right*U_Right)/(Rho_Left-Rho_Right)
           Mach = Shock_Speed/C_Right
  
    END Subroutine Type_P_Ratio 
!
! ==================================================================
!
  Subroutine PS_Convert_Variable_Type (Operation, Gamma, C, P)

  integer                :: Operation
  REALTYPE               :: Gamma
  REALTYPE, dimension(3) :: C, P

  SELECT CASE (Operation)
    CASE (2)
      P = C

    CASE (1)
      P(1) = C(1)
      P(2) = C(2)/C(1)
      P(3) = (Gamma-ONE)*(C(3)-HALF*P(1)*P(2)**2)

    CASE (-1)
      C(1) = P(1)
      C(2) = P(2)*P(1)
      C(3) = P(3)/(Gamma-ONE)+HALF*P(1)*P(2)**2

    CASE (-2)
      C = P

  END SELECT

  END Subroutine PS_Convert_Variable_Type
!
! ==================================================================
!
  Subroutine PS_Primitive_To_Conservative (Gamma, P, C)

  REALTYPE               :: Gamma
  REALTYPE, dimension(3) :: P, C

  C(1) = P(1)
  C(2) = P(2)*P(1)
  C(3) = P(3)/(Gamma-ONE)+HALF*P(1)*P(2)**2

  END Subroutine PS_Primitive_To_Conservative 
!
! ==================================================================
!
  Subroutine PS_Conservative_To_Primitive (Gamma, C, P)

  REALTYPE               :: Gamma
  REALTYPE, dimension(3) :: C, P

  P(1) = C(1)
  P(2) = C(2)/C(1)
  P(3) = (Gamma-ONE)*(C(3)-HALF*P(1)*P(2)**2)

  END Subroutine PS_Conservative_To_Primitive 

END MODULE Shock_Flow_Specification
