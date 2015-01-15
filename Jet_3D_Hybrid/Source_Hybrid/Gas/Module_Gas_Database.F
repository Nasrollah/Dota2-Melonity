#include "pseudopack.h"


MODULE Gas_Properties

  USE PseudoPack, ONLY: I_Am, First_Processor

implicit NONE

REALTYPE, parameter :: Celius_2_Kelvin  =    273.150d0
REALTYPE, parameter :: Kelvin_2_Rankine =      1.80d0
REALTYPE, parameter :: Atm_2_Pascal     = 101325.0d0         ! Pascal=kg m/s^2
REALTYPE, parameter :: R_Universal      =      8.314410      ! J/K-mole
REALTYPE, parameter :: KJKg_2_Jmole     = 29.23950d0/1.010d0  

integer  :: Gas_Setting
logical  :: Unit_Conversion

REALTYPE :: Reference_Temperature      
REALTYPE :: Reference_Pressure         

REALTYPE :: R_Specific = 287                   !  = cp-cv = R_Universal/MW_Gas
REALTYPE :: Gamma      = 1.40d0                !  Specific Gas Constant
REALTYPE :: Pr         = 0.71                  !  Pantle mumber
REALTYPE :: Re         = ZERO                  !  Reynold Number
REALTYPE :: Sc         = ONE                   !  Schmidt Number = mu/zeta

REALTYPE ::    mu_default = 1.51d-5    ! m^2 s
REALTYPE :: kappa_default = 2.60d-2    ! W/(m K)
REALTYPE ::  zeta_default = HUGE(ONE)  ! m^-2 s^-1

REALTYPE ::    mu_0       = 1.51d-5    ! m^2 s
REALTYPE :: kappa_0       = 2.60d-2    ! W/(m K)
REALTYPE ::  zeta_0       = HUGE(ONE)  ! m^-2 s^-1

integer  ::   Gas_Composition = 1
integer  :: Gamma_Formulation = 1

REALTYPE :: MW_Gas_1, Rho_Gas_1, Gamma_Gas_1, Cp_Gas_1, Cv_Gas_1, P_Gas_1
REALTYPE :: MW_Gas_2, Rho_Gas_2, Gamma_Gas_2, Cp_Gas_2, Cv_Gas_2, P_Gas_2
REALTYPE :: MW_Gas_3, Rho_Gas_3, Gamma_Gas_3, Cp_Gas_3, Cv_Gas_3, P_Gas_3

! Gamma_Gas(:,n) stores the Gamma, Cp, Cv, MW, Rho of gas n
REALTYPE, dimension(257,5) :: Gamma_Gas

#include "Gas_Properties.i"

TYPE Fluid_Properties
  Character(LEN=8) :: Name

  integer  :: Index

  REALTYPE ::     f   ! Mass fraction

  REALTYPE :: Gamma   

  REALTYPE ::    Cp   ! J/(mol-K) 
  REALTYPE ::    Cv   ! J/(mol-K)

  REALTYPE ::    MW   ! g/(mol)
  REALTYPE ::   Rho   ! g/cm^3

  REALTYPE ::   nu    ! (cm s)/g       ! Kinematic viscosity = mu/rho

  REALTYPE :: T_mu                     ! Reference temperatures in Kelvin
  REALTYPE ::   mu    ! (cm^2 s^-1)

  REALTYPE :: T_kappa
  REALTYPE ::   kappa ! W/(cm K)

  REALTYPE :: T_zeta
  REALTYPE ::   zeta  ! cm^2 s^-1

  REALTYPE ::    Pr 
  REALTYPE ::    Re

  REALTYPE ::    R
END TYPE

TYPE (Fluid_Properties)               :: Xenon, Argon
TYPE (Fluid_Properties)               :: Hydrogen, Nitrogen, Oxygen
TYPE (Fluid_Properties)               :: Air, Acetone, SF6 
TYPE (Fluid_Properties)               :: Air_Acetone, AA
TYPE (Fluid_Properties)               :: CHBr, CRF
TYPE (Fluid_Properties), dimension(9) :: Gas
 
PUBLIC

PUBLIC  :: Fluid_Properties

PUBLIC  :: Xenon, Argon
PUBLIC  :: Hydrogen, Nitrogen, Oxygen
PUBLIC  :: Air, Acetone, SF6 
PUBLIC  :: Air_Acetone, AA
PUBLIC  :: CHBr, CRF
PUBLIC  :: Gas

PUBLIC  ::   Gas_Composition
PUBLIC  :: Gamma_Formulation

PUBLIC  :: Read_Gas_Properties
PUBLIC  :: Gas_Unit_Conversion

CONTAINS

  Subroutine Read_Gas_Properties (Filename)

  integer            :: lid, IOS
  character(LEN=256) :: Filename

  if (I_Am == First_Processor) &
    write (6,*) 'Reading .... ', TRIM(Filename)

  lid = 66 

  OPEN (lid, File = TRIM(Filename), STATUS='OLD', IOSTAT=IOS)

  read (lid,'(//)')
  read (lid,*) Gas_Setting, Unit_Conversion
  read (lid,*) Reference_Temperature
  read (lid,*) Reference_Pressure
  read (lid,'(///)')
  read (lid,*) Air
  read (lid,*) SF6
  read (lid,*) Acetone
  read (lid,*) AA
  read (lid,*) Xenon
  read (lid,*) Argon
  read (lid,*) Hydrogen
  read (lid,*) Oxygen
  read (lid,*) Nitrogen
  read (lid,*) CHBr
  read (lid,*) CRF

  CLOSE (lid)

  if (I_Am == First_Processor) &
    write (6,*) '....... Done ', TRIM(Filename)

#if 0
  call Debug
#endif

  Reference_Temperature = Reference_Temperature + Celius_2_Kelvin
  Reference_Pressure    = Reference_Pressure    * Atm_2_Pascal 

  Air_Acetone = AA

CONTAINS
!
! ====================================================
!
  SUBROUTINE Debug

  if (I_Am == First_Processor) then
     write (6,100) Air
     write (6,100) SF6
     write (6,100) CRF
   endif

  if (I_Am == First_Processor) write (6,*) 'Check the Cp for Acetone!'
  Acetone%Cp = Acetone%Cp*KJKg_2_Jmole

100 format (1x,a10,i5,100(g12.5,:))

  END SUBROUTINE Debug

  END Subroutine Read_Gas_Properties
! 
! =========================================================================
!
  Subroutine Gas_Unit_Conversion (Gas)

  TYPE (Fluid_Properties) :: Gas
  
if (I_Am == First_Processor) then
  write (6,*) ' *********  Before Unit Conversion *******'
  write (6,100) Gas
endif

  Gas%MW    = Gas%MW   /1000      ! kg/mole
  Gas%Rho   = Gas%Rho  *1000      ! kg/m^3

  Gas%nu    = Gas%nu   /10        ! kg/(m-s)
  Gas%mu    = Gas%mu   /10000     ! m^2/s
  Gas%kappa = Gas%kappa*100       ! W/(m-K)
  Gas%zeta  = Gas%zeta /10000     ! m^2/s

  Gas%R     = R_Universal/Gas%MW  ! J/(kg-K)  ! J = N m, N = kg-m/s^2
  Gas%Cp    = Gas%Cp     /Gas%MW  ! J/(kg-K)
  Gas%Cv    = Gas%Cv     /Gas%MW  ! J/(kg-K)

if (I_Am == First_Processor) then
  write (6,*) ' *********  After Unit Conversion *******'
  write (6,100) Gas
  write (6,*) '*****************************************'
endif

100 format (1x,a10,i5,100(g12.5,:))

  END Subroutine Gas_Unit_Conversion

END MODULE Gas_Properties
