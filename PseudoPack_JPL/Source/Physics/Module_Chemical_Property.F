#include "pseudopack.h"


!------------------------------------------------------------------------------
! Units : Joule = NT-m ! NT = kg-m/s**2 ! Joule/kg = 9.869233e-03 cm**3-atm/g
!         Rankine (R) = 1.8 Kelvin (K)
!         cal = 4.1868 joule
! Universal Gas Constant (Joule/kg-mole K= 82.056156 cm**3-atm/g-mole-K)
! use UNIX utility  units   to do unit conversion, see man units for detail.
!
! MW, Molecular Weight.  (kg/kg-mole)
! hto, Sum of sensible enthalpy and Chemical potential at 0 K (Joule/kg).
!
! Coeffficients for specific heat of species Cp and Gibbs Energy of Species g.
! Cp = (a+T*(b+T*(c+T*(d+e*T))))*R/MW
! g  = (a*(1-ln(T))-T*(b+T*(c/3+T*(d/6+e/10*T)))/2+f)*T*R+hto*MW
!
! R  with units (J/kg-K) = R_Universal_Gas_Constant*f_i
! Cp with units (J/kg-K)
! g  with units (J/kg-mole)
! acp, bcp, ccp, dcp, ecp, fcp with units kg/kg-mole
!
! Coeffficients for laminar viscosity from Sutherland's Law.
! (mu/mu0) = (T/Tmu)**(3/2)*(Tmu+Smu)/(T+Smu)   (kg/m-sec)
!    Tmu, Smu are dimensioned in Kelvin (K)
!
! Coeffficients for thermal conductivity from Sutherland's Law.
! (Tc/kp0) = (T/Tkp)**(3/2)*(Tkp+Skp)/(T+Skp)   (J/m-K-sec)
!    Tkp, Skp are dimensioned in Kelvin (K)
!
! Coeffficients for Binary Diffusion Coefficients.
! abdc, Collision Diameter (Angstroms).
! bbdc, Molecular Temperature Parameter (K).
!
! Reaction rate variables for forward reaction rates, kf, of each reaction.
! kf=akf*T**pkf*exp(-ekf/(1.987*T))
! Pre-exponential akf; Temperature power pkf; Activation Energy ekf, cal/gm-mole
!
!------------------------------------------------------------------------------

MODULE Chemical_Property

implicit NONE

REALTYPE, parameter :: R_0                      = 82.056156d0
REALTYPE, parameter :: R_Universal_Gas_Constant = 8314.340d0 
REALTYPE, parameter :: Kelvin_to_Rankine        = 1.80d0
REALTYPE, parameter :: Atm_to_Metric            = 101325.0d0
      
integer             ::             Model = 0
integer             ::               NCS = 2
integer             ::               NCR = 0
character(LEN=100)  ::         Fuel_Type = 'He'
integer             ::        Fuel_Index = 2
integer             ::            f_Fuel = ONE

REALTYPE            ::  Frequency_Factor = 5.5410d14
REALTYPE            ::             Power =    0.0d0
REALTYPE            ::        Activation = 7200.0d0

REALTYPE            ::    Schmidt_Number = 0.220d0

TYPE Chemical
  integer           :: Model
  integer           :: NCS, NCR
  integer           :: Fuel_Index
  Character(LEN=100):: Fuel_Type
  integer           :: f_Fuel

  REALTYPE          :: Frequency_Factor, Power, Activation
  REALTYPE          :: Schmidt_Number

  REALTYPE, dimension(:)  , POINTER :: R_n
  REALTYPE, dimension(:)  , POINTER :: MW, hto
  REALTYPE, dimension(:)  , POINTER :: acp ,bcp ,ccp ,dcp ,ecp ,fcp
  REALTYPE, dimension(:)  , POINTER :: acpb,bcpb,ccpb,dcpb,ecpb
  REALTYPE, dimension(:)  , POINTER :: abdc, bbdc
  REALTYPE, dimension(:)  , POINTER :: mu0, Tmu, Smu, Cmu
  REALTYPE, dimension(:)  , POINTER :: kp0, Tkp, Skp, Ckp
  REALTYPE, dimension(:)  , POINTER :: akf, pkf, ekf
  REALTYPE, dimension(:,:), POINTER :: ac
  integer , dimension(:,:), POINTER :: LHS, RHS
  REALTYPE, dimension(:,:), POINTER :: MW_1, MW_2
END TYPE Chemical

INTERFACE PS_Chemical_Property
  MODULE PROCEDURE PS_Chemical_Species
END INTERFACE

PRIVATE
PUBLIC  :: Chemical
PUBLIC  :: R_0, R_Universal_Gas_Constant, Kelvin_to_Rankine, Atm_to_Metric
PUBLIC  :: Model, NCS, NCR, Fuel_Type, Fuel_Index, f_Fuel
PUBLIC  :: Frequency_Factor, Power, Activation, Schmidt_Number

PUBLIC  :: PS_Chemical_Property
PUBLIC  :: PS_Chemical_Species

CONTAINS
  
  Subroutine PS_Chemical_Species (Index, Name, MW,  hto,               &
                                  acp , bcp , ccp , dcp , ecp , fcp , &
                                  acpb, bcpb, ccpb, dcpb, ecpb,       &
                                   mu0,  Tmu,  Smu,  Cmu,             &
                                   kp0,  Tkp,  Skp,  Ckp,             &
                                  abdc, bbdc,  R_n, lid50)

  integer  :: Index
  REALTYPE ::  R_n
  REALTYPE ::   MW, hto
  REALTYPE :: acp , bcp , ccp , dcp , ecp , fcp
  REALTYPE :: acpb, bcpb, ccpb, dcpb, ecpb
  REALTYPE ::  mu0,  Tmu,  Smu,  Cmu
  REALTYPE ::  kp0,  Tkp,  Skp,  Ckp
  REALTYPE :: abdc, bbdc

  Character(LEN=100) :: Name

  integer , OPTIONAL :: lid50

  SELECT CASE (Index)
    CASE (1)                   
        Name = 'O_2'

        MW = 32.000       ;  hto = -272918.21             ! O2

       acp =  3.0809      ;  bcp =  0.16962d-2  ;  ccp = -0.76334d-6
       dcp =  0.17140d-9  ;  ecp = -0.14116d-13 ;  fcp =  6.6490

       mu0 = 0.1919d-4    ;  Tmu = 273.111      ;  Smu = 138.889
       kp0 = 0.0245781    ;  Tkp = 273.111      ;  Skp = 222.222

      abdc = 3.467        ; bbdc = 106.7

    CASE (2)                        
        Name = 'H_2'

        MW =  2.016       ;  hto = -4280070.46             ! H2

       acp = 3.4990       ;  bcp = -0.18651d-3  ;  ccp =  0.46064d-6
       dcp = -0.13157d-9  ;  ecp =  0.11679d-13 ;  fcp = -4.1500

       mu0 = 0.08411d-4   ;  Tmu = 273.111      ;  Smu =  96.6667
       kp0 = 0.162815     ;  Tkp = 273.111      ;  Skp =  166.667

      abdc = 2.827        ; bbdc =  59.7

    CASE (3)
        Name = 'H2_O'

        MW = 18.016       ;  hto =  -13973684.55           ! H2-O

       acp = 3.4990       ;  bcp =  0.14878d-2  ;  ccp =  0.87544d-7
       dcp = -0.11499d-9  ;  ecp =  0.13495d-13 ;  fcp =  2.3296

       mu0 = 0.1703d-4    ;  Tmu = 416.667      ;  Smu = 861.111
       kp0 = 0.0179442    ;  Tkp = 273.111      ;  Skp = 1277.78

      abdc = 2.641        ; bbdc = 809.1

    CASE (4)
        Name = 'O-H'

        MW = 17.008       ;  hto =    1776329.80           ! O-H

       acp = 3.5811       ;  bcp = -0.28556d-3  ;  ccp =  0.56412d-6
       dcp = -0.17127d-9  ;  ecp =  0.15762d-13 ;  fcp =  1.7590

       mu0 = 0.1826d-4    ;  Tmu = 273.111      ;  Smu = 138.70
       kp0 = 0.0453456    ;  Tkp = 273.111      ;  Skp = 230.783

      abdc = 3.147        ; bbdc =  79.8

    CASE (5)
        Name = 'H'

        MW =  1.008       ;  hto =  210118548.80           ! H

       acp = 2.5000       ;  bcp =  0.63743d-14 ;  ccp = -0.51755d-17
       dcp =  0.15299d-20 ;  ecp = -0.14782d-24 ;  fcp = -0.4510

       mu0 = 0.07006d-4   ;  Tmu = 273.111      ;  Smu = 117.667
       kp0 = 0.217028     ;  Tkp = 273.111      ;  Skp = 117.817

      abdc = 2.708        ; bbdc =  37.0

    CASE (6)
        Name = 'O'

        MW = 16.000       ;  hto =   15179906.26           ! O

       acp = 2.6746       ;  bcp = -0.28666d-3  ;  ccp =  0.16002d-6
       dcp = -0.36096d-10 ;  ecp =  0.30438d-14 ;  fcp =  4.2098

       mu0 = 0.1743d-4    ;  Tmu = 273.111      ;  Smu = 156.972
       kp0 = 0.0358538    ;  Tkp = 273.111      ;  Skp = 130.561

      abdc = 3.050        ; bbdc = 106.7

    CASE (7)
        Name = 'H-O_2'

        MW = 33.008       ;  hto =     316548.45           ! H-O2

       acp = 3.2263       ;  bcp =  0.38264d-2  ;  ccp = -0.16480d-5
       dcp =  0.32550d-9  ;  ecp = -0.23968d-13 ;  fcp =  7.8090

       mu0 = 0.0893656d-4 ;  Tmu = 273.111      ;  Smu = 324.506
       kp0 = 0.0152561    ;  Tkp = 273.111      ;  Skp = 117.189

      abdc = 4.196        ; bbdc = 289.3

    CASE (8)
        Name = 'H_2-O_2'

        MW = 34.016       ;  hto =   -4393772.81           ! H2-O2

       acp = 3.8344       ;  bcp =  0.57626d-2  ;  ccp = -0.23656d-5
       dcp =  0.45571d-9  ;  ecp = -0.33155d-13 ;  fcp =  4.5980

       mu0 = 0.0893656d-4 ;  Tmu = 273.111      ;  Smu = 324.506
       kp0 = 0.0152561    ;  Tkp = 273.111      ;  Skp = 117.189

      abdc = 4.196        ; bbdc = 289.3

    CASE (9)
        Name = 'N_2'

        MW = 28.016       ;  hto =    -302736.23           ! N2

       acp = 3.1459       ;  bcp =  0.99154d-3  ;  ccp = -0.22912d-6
       dcp =  0.12181d-10 ;  ecp =  0.11024d-14 ;  fcp =  4.8370

       mu0 = 0.1663d-4    ;  Tmu = 273.111      ;  Smu = 106.667
       kp0 = 0.0245781    ;  Tkp = 273.111      ;  Skp = 166.667

      abdc = 3.798        ; bbdc =  71.4

    CASE (10)
        Name = 'He'

        MW =  4.0026      ;  hto =   -1548318.90           ! He

       acp = 2.5000       ;  bcp =  0.0d0       ;  ccp =  0.0d0
       dcp =  0.0d0       ;  ecp =  0.0d0       ;  fcp =  0.2870d-1

       mu0 = 0.18699d-4   ;  Tmu = 273.111      ;  Smu =   79.4444
       kp0 = 0.21702837   ;  Tkp = 273.111      ;  Skp =  117.81667

      abdc = 0.0d0        ; bbdc =   0.0d0

    CASE DEFAULT
      call Read_Chemical_Property

  END SELECT 

  Cmu = mu0*(Tmu+Smu)/(SQRT(Tmu)*Tmu)
  Ckp = kp0*(Tkp+Skp)/(SQRT(Tkp)*Tkp)

  acpb = acp ; bcpb = bcp/TWO  ; ccpb = ccp/THREE
               dcpb = dcp/FOUR ; ecpb = ecp/FIVE

  R_n = R_Universal_Gas_Constant/MW

  CONTAINS
    
    Subroutine Read_Chemical_Property

    integer  :: lid

    lid = 50 ; if (PRESENT(lid50)) lid = lid50

    read (lid,*) Index 
    read (lid,*) Name
    read (lid,*)   MW, hto
    read (lid,*)  acp, bcp, ccp, dcp, ecp, fcp 
    read (lid,*)  mu0, Tmu, Smu
    read (lid,*)  kp0, Tkp, Skp
    read (lid,*) abdc, bbdc

    END Subroutine Read_Chemical_Property

  END Subroutine PS_Chemical_Species

END MODULE Chemical_Property
