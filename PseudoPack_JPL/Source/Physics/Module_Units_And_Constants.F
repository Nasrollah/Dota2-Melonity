#include "pseudopack.h"

#define PARATYPE	REALTYPE, parameter :: 

MODULE Units_And_Constants 

implicit NONE

PARATYPE PS_Kelvin_to_Rankine = 1.80d0
PARATYPE PS_Celius_to_Kelvin  = 273.15
PARATYPE PS_Atm_to_Metric     = 101325.0d0

PARATYPE PS_YOTTA          = 1.0d24 
PARATYPE PS_ZETTA          = 1.0d21 
PARATYPE PS_EXA            = 1.0d18 
PARATYPE PS_PETA           = 1.0d15 
PARATYPE PS_TERA           = 1.0d12 
PARATYPE PS_GIGA           = 1.0d9 
PARATYPE PS_MEGA           = 1.0d6 
PARATYPE PS_KILO           = 1.0d3 
PARATYPE PS_MILLI          = 1.0d-3 
PARATYPE PS_MICRO          = 1.0d-6 
PARATYPE PS_NANO           = 1.0d-9 
PARATYPE PS_PICO           = 1.0d-12 
PARATYPE PS_FEMTO          = 1.0d-15 
PARATYPE PS_ATTO           = 1.0d-18 
PARATYPE PS_ZEPTO          = 1.0d-21 
PARATYPE PS_YOCTO          = 1.0d-24 

PARATYPE PS_FINE_STRUCTURE = 7.2973525330d-3   
PARATYPE PS_AVOGADRO       = 6.0221419900d23                !  1/mol 

!  CGS Units

PARATYPE PS_CGS_SPEED_OF_LIGHT           = 2.99792458d10    ! cm / s  
PARATYPE PS_CGS_GRAVITATIONAL_CONSTANT   = 6.673d-8         ! cm^3 / g s^2 
PARATYPE PS_CGS_PLANCKS_CONSTANT_H       = 6.62606876d-27   ! g cm^2 / s 
PARATYPE PS_CGS_PLANCKS_CONSTANT_HBAR    = 1.05457159642d-27! g cm^2 / s 
PARATYPE PS_CGS_ASTRONOMICAL_UNIT        = 1.49597870691d13 ! cm 
PARATYPE PS_CGS_LIGHT_YEAR               = 9.46053620707d17 ! cm 
PARATYPE PS_CGS_PARSEC                   = 3.08567758135d18 ! cm  
PARATYPE PS_CGS_GRAV_ACCEL               = 9.80665d2        ! cm / s^2 
PARATYPE PS_CGS_ELECTRON_VOLT            = 1.602176462d-12  ! g cm^2 / s^2 
PARATYPE PS_CGS_MASS_ELECTRON            = 9.10938188d-28   ! g 
PARATYPE PS_CGS_MASS_MUON                = 1.88353109d-25   ! g 
PARATYPE PS_CGS_MASS_PROTON              = 1.67262158d-24   ! g 
PARATYPE PS_CGS_MASS_NEUTRON             = 1.67492716d-24   ! g  
PARATYPE PS_CGS_RYDBERG                  = 2.17987190389d-11! g cm^2 / s^2 
PARATYPE PS_CGS_BOLTZMANN                = 1.3806503d-16    ! g cm^2 / K s^2 
PARATYPE PS_CGS_BOHR_MAGNETON            = 9.27400899d-20   ! A cm^2  
PARATYPE PS_CGS_NUCLEAR_MAGNETON         = 5.05078317d-23   ! A cm^2  
PARATYPE PS_CGS_ELECTRON_MAGNETIC_MOMENT = 9.28476362d-20   ! A cm^2 
PARATYPE PS_CGS_PROTON_MAGNETIC_MOMENT   = 1.410606633d-22  ! A cm^2 
PARATYPE PS_CGS_MOLAR_GAS                = 8.314472d7       ! g cm^2/K mol s^2
PARATYPE PS_CGS_STANDARD_GAS_VOLUME      = 2.2710981d4      ! cm^3 / mol 
PARATYPE PS_CGS_MINUTE                   = 6.0d1            ! s 
PARATYPE PS_CGS_HOUR                     = 3.6d3            ! s 
PARATYPE PS_CGS_DAY                      = 8.64d4           ! s 
PARATYPE PS_CGS_WEEK                     = 6.048d5          ! s 
PARATYPE PS_CGS_INCH                     = 2.54d0           ! cm 
PARATYPE PS_CGS_FOOT                     = 3.048d1          ! cm 
PARATYPE PS_CGS_YARD                     = 9.144d1          ! cm 
PARATYPE PS_CGS_MILE                     = 1.609344d5       ! cm 
PARATYPE PS_CGS_NAUTICAL_MILE            = 1.852d5          ! cm 
PARATYPE PS_CGS_FATHOM                   = 1.8288d2         ! cm 
PARATYPE PS_CGS_MIL                      = 2.54d-3          ! cm 
PARATYPE PS_CGS_POINT                    = 3.52777777778d-2 ! cm 
PARATYPE PS_CGS_TEXPOINT                 = 3.51459803515d-2 ! cm 
PARATYPE PS_CGS_MICRON                   = 1.0d-4           ! cm 
PARATYPE PS_CGS_ANGSTROM                 = 1.0d-8           ! cm 
PARATYPE PS_CGS_HECTARE                  = 1.0d8            ! cm^2 
PARATYPE PS_CGS_ACRE                     = 4.04685642241d7  ! cm^2 
PARATYPE PS_CGS_BARN                     = 1.0d-24          ! cm^2 
PARATYPE PS_CGS_LITER                    = 1.0d3            ! cm^3 
PARATYPE PS_CGS_US_GALLON                = 3.78541178402d3  ! cm^3 
PARATYPE PS_CGS_QUART                    = 9.46352946004d2  ! cm^3 
PARATYPE PS_CGS_PINT                     = 4.73176473002d2  ! cm^3 
PARATYPE PS_CGS_CUP                      = 2.36588236501d2  ! cm^3 
PARATYPE PS_CGS_FLUID_OUNCE              = 2.95735295626d1  ! cm^3 
PARATYPE PS_CGS_TABLESPOON               = 1.47867647813d1  ! cm^3 
PARATYPE PS_CGS_TEASPOON                 = 4.92892159375d0  ! cm^3 
PARATYPE PS_CGS_CANADIAN_GALLON          = 4.54609d3        ! cm^3 
PARATYPE PS_CGS_UK_GALLON                = 4.546092d3       ! cm^3 
PARATYPE PS_CGS_MILES_PER_HOUR           = 4.4704d1         ! cm / s 
PARATYPE PS_CGS_KILOMETERS_PER_HOUR      = 2.77777777778d1  ! cm / s 
PARATYPE PS_CGS_KNOT                     = 5.14444444444d1  ! cm / s 
PARATYPE PS_CGS_POUND_MASS               = 4.5359237d2      ! g 
PARATYPE PS_CGS_OUNCE_MASS               = 2.8349523125d1   ! g 
PARATYPE PS_CGS_TON                      = 9.0718474d5      ! g 
PARATYPE PS_CGS_METRIC_TON               = 1.0d6            ! g 
PARATYPE PS_CGS_UK_TON                   = 1.0160469088d6   ! g 
PARATYPE PS_CGS_TROY_OUNCE               = 3.1103475d1      ! g 
PARATYPE PS_CGS_CARAT                    = 2.0d-1           ! g 
PARATYPE PS_CGS_UNIFIED_ATOMIC_MASS      = 1.66053873d-24   ! g 
PARATYPE PS_CGS_GRAM_FORCE               = 9.80665d2        ! cm g / s^2 
PARATYPE PS_CGS_POUND_FORCE              = 4.44822161526d5  ! cm g / s^2 
PARATYPE PS_CGS_KILOPOUND_FORCE          = 4.44822161526d8  ! cm g / s^2 
PARATYPE PS_CGS_POUNDAL                  = 1.38255d4        ! cm g / s^2 
PARATYPE PS_CGS_CALORIE                  = 4.1868d7         ! g cm^2 / s^2 
PARATYPE PS_CGS_BTU                      = 1.05505585262d10 ! g cm^2 / s^2 
PARATYPE PS_CGS_THERM                    = 1.05506d15       ! g cm^2 / s^2  
PARATYPE PS_CGS_HORSEPOWER               = 7.457d9          ! g cm^2 / s^3 
PARATYPE PS_CGS_BAR                      = 1.0d6            ! g / cm s^2  
PARATYPE PS_CGS_STD_ATMOSPHERE           = 1.01325d6        ! g / cm s^2 
PARATYPE PS_CGS_TORR                     = 1.33322368421d3  ! g / cm s^2 
PARATYPE PS_CGS_METER_OF_MERCURY         = 1.33322368421d6  ! g / cm s^2 
PARATYPE PS_CGS_INCH_OF_MERCURY          = 3.38638815789d4  ! g / cm s^2 
PARATYPE PS_CGS_INCH_OF_WATER            = 2.490889d3       ! g / cm s^2 
PARATYPE PS_CGS_PSI                      = 6.89475729317d4  ! g / cm s^2 
PARATYPE PS_CGS_POISE                    = 1.0d0            ! g / cm s 
PARATYPE PS_CGS_STOKES                   = 1.0d0            ! cm^2 / s 
PARATYPE PS_CGS_FARADAY                  = 9.6485341472d4   ! A s / mol 
PARATYPE PS_CGS_ELECTRON_CHARGE          = 1.602176462d-19  ! A s 
PARATYPE PS_CGS_GAUSS                    = 1.0d-1           ! g / A s^2 
PARATYPE PS_CGS_STILB                    = 1.0d0            ! cd / cm^2 
PARATYPE PS_CGS_LUMEN                    = 1.0d0            ! cd sr 
PARATYPE PS_CGS_LUX                      = 1.0d-4           ! cd sr / cm^2 
PARATYPE PS_CGS_PHOT                     = 1.0d0            ! cd sr / cm^2 
PARATYPE PS_CGS_FOOTCANDLE               = 1.076d-3         ! cd sr / cm^2 
PARATYPE PS_CGS_LAMBERT                  = 1.0d0            ! cd sr / cm^2 
PARATYPE PS_CGS_FOOTLAMBERT              = 1.07639104d-3    ! cd sr / cm^2 
PARATYPE PS_CGS_CURIE                    = 3.7d10           ! 1 / s 
PARATYPE PS_CGS_ROENTGEN                 = 2.58d-7          ! A s / g 
PARATYPE PS_CGS_RAD                      = 1.0d2            ! cm^2 / s^2 
PARATYPE PS_CGS_SOLAR_MASS               = 1.98892d33       ! g 
PARATYPE PS_CGS_BOHR_RADIUS              = 5.291772083d-9   ! cm 
PARATYPE PS_CGS_NEWTON                   = 1.0d5            ! g cm   / s^2 
PARATYPE PS_CGS_DYNE                     = 1.0d0            ! g cm   / s^2 
PARATYPE PS_CGS_JOULE                    = 1.0d7            ! g cm^2 / s^2 
PARATYPE PS_CGS_ERG                      = 1.0d0            ! g cm^2 / s^2 

PARATYPE PS_CGS_STEFAN_BOLTZMANN_CONST    = 5.67039934436d-5  ! g / K^4 s^3 
PARATYPE PS_CGS_THOMSON_CROSS_SECTION     = 6.65245853542d-25 ! cm^2 

#if 0
PARATYPE PS_CGSM_BOHR_MAGNETON            = 9.27400899d-21   ! Abamp cm^2  
PARATYPE PS_CGSM_NUCLEAR_MAGNETON         = 5.05078317d-24   ! Abamp cm^2  
PARATYPE PS_CGSM_ELECTRON_MAGNETIC_MOMENT = 9.28476362d-21   ! Abamp cm^2 
PARATYPE PS_CGSM_PROTON_MAGNETIC_MOMENT   = 1.410606633d-23  ! Abamp cm^2 
PARATYPE PS_CGSM_FARADAY                  = 9.6485341472d3   ! Abamp s / mol 
PARATYPE PS_CGSM_ELECTRON_CHARGE          = 1.602176462d-20  ! Abamp s 
PARATYPE PS_CGSM_GAUSS                    = 1.0d0            ! g / Abamp s^2 
PARATYPE PS_CGSM_ROENTGEN                 = 2.58d-8          ! Abamp s / g 
#endif

!  MKS Units

PARATYPE PS_MKS_SPEED_OF_LIGHT           = 2.99792458d8      ! m / s 
PARATYPE PS_MKS_GRAVITATIONAL_CONSTANT   = 6.673d-11         ! m^3 / kg s^2 
PARATYPE PS_MKS_PLANCKS_CONSTANT_H       = 6.62606876d-34    ! kg m^2 / s 
PARATYPE PS_MKS_PLANCKS_CONSTANT_HBAR    = 1.05457159642d-34 ! kg m^2 / s 
PARATYPE PS_MKS_ASTRONOMICAL_UNIT        = 1.49597870691d11  ! m 
PARATYPE PS_MKS_LIGHT_YEAR               = 9.46053620707d15  ! m 
PARATYPE PS_MKS_PARSEC                   = 3.08567758135d16  ! m 
PARATYPE PS_MKS_GRAV_ACCEL               = 9.80665d0         ! m / s^2 
PARATYPE PS_MKS_ELECTRON_VOLT            = 1.602176462d-19   ! kg m^2 / s^2 
PARATYPE PS_MKS_MASS_ELECTRON            = 9.10938188d-31    ! kg 
PARATYPE PS_MKS_MASS_MUON                = 1.88353109d-28    ! kg 
PARATYPE PS_MKS_MASS_PROTON              = 1.67262158d-27    ! kg 
PARATYPE PS_MKS_MASS_NEUTRON             = 1.67492716d-27    ! kg 
PARATYPE PS_MKS_RYDBERG                  = 2.17987190389d-18 ! kg m^2 / s^2 
PARATYPE PS_MKS_BOLTZMANN                = 1.3806503d-23     ! kg m^2 / K s^2 
PARATYPE PS_MKS_BOHR_MAGNETON            = 9.27400899d-24    ! A m^2 
PARATYPE PS_MKS_NUCLEAR_MAGNETON         = 5.05078317d-27    ! A m^2 
PARATYPE PS_MKS_ELECTRON_MAGNETIC_MOMENT = 9.28476362d-24    ! A m^2 
PARATYPE PS_MKS_PROTON_MAGNETIC_MOMENT   = 1.410606633d-26   ! A m^2 
PARATYPE PS_MKS_MOLAR_GAS                = 8.314472d0        ! kg m^2/K mol s^2
PARATYPE PS_MKS_STANDARD_GAS_VOLUME      = 2.2710981d-2      ! m^3 / mol 
PARATYPE PS_MKS_MINUTE                   = 6.0d1             ! s 
PARATYPE PS_MKS_HOUR                     = 3.6d3             ! s 
PARATYPE PS_MKS_DAY                      = 8.64d4            ! s 
PARATYPE PS_MKS_WEEK                     = 6.048d5           ! s 
PARATYPE PS_MKS_INCH                     = 2.54d-2           ! m 
PARATYPE PS_MKS_FOOT                     = 3.048d-1          ! m 
PARATYPE PS_MKS_YARD                     = 9.144d-1          ! m 
PARATYPE PS_MKS_MILE                     = 1.609344d3        ! m 
PARATYPE PS_MKS_NAUTICAL_MILE            = 1.852d3           ! m 
PARATYPE PS_MKS_FATHOM                   = 1.8288d0          ! m 
PARATYPE PS_MKS_MIL                      = 2.54d-5           ! m 
PARATYPE PS_MKS_POINT                    = 3.52777777778d-4  ! m 
PARATYPE PS_MKS_TEXPOINT                 = 3.51459803515d-4  ! m 
PARATYPE PS_MKS_MICRON                   = 1.0d-6            ! m 
PARATYPE PS_MKS_ANGSTROM                 = 1.0d-10           ! m 
PARATYPE PS_MKS_HECTARE                  = 1.0d4             ! m^2 
PARATYPE PS_MKS_ACRE                     = 4.04685642241d3   ! m^2 
PARATYPE PS_MKS_BARN                     = 1.0d-28           ! m^2 
PARATYPE PS_MKS_LITER                    = 1.0d-3            ! m^3 
PARATYPE PS_MKS_US_GALLON                = 3.78541178402d-3  ! m^3 
PARATYPE PS_MKS_QUART                    = 9.46352946004d-4  ! m^3 
PARATYPE PS_MKS_PINT                     = 4.73176473002d-4  ! m^3 
PARATYPE PS_MKS_CUP                      = 2.36588236501d-4  ! m^3 
PARATYPE PS_MKS_FLUID_OUNCE              = 2.95735295626d-5  ! m^3 
PARATYPE PS_MKS_TABLESPOON               = 1.47867647813d-5  ! m^3 
PARATYPE PS_MKS_TEASPOON                 = 4.92892159375d-6  ! m^3 
PARATYPE PS_MKS_CANADIAN_GALLON          = 4.54609d-3        ! m^3 
PARATYPE PS_MKS_UK_GALLON                = 4.546092d-3       ! m^3 
PARATYPE PS_MKS_MILES_PER_HOUR           = 4.4704d-1         ! m / s 
PARATYPE PS_MKS_KILOMETERS_PER_HOUR      = 2.77777777778d-1  ! m / s 
PARATYPE PS_MKS_KNOT                     = 5.14444444444d-1  ! m / s 
PARATYPE PS_MKS_POUND_MASS               = 4.5359237d-1      ! kg 
PARATYPE PS_MKS_OUNCE_MASS               = 2.8349523125d-2   ! kg 
PARATYPE PS_MKS_TON                      = 9.0718474d2       ! kg 
PARATYPE PS_MKS_METRIC_TON               = 1.0d3             ! kg 
PARATYPE PS_MKS_UK_TON                   = 1.0160469088d3    ! kg 
PARATYPE PS_MKS_TROY_OUNCE               = 3.1103475d-2      ! kg 
PARATYPE PS_MKS_CARAT                    = 2.0d-4            ! kg 
PARATYPE PS_MKS_UNIFIED_ATOMIC_MASS      = 1.66053873d-27    ! kg 
PARATYPE PS_MKS_GRAM_FORCE               = 9.80665d-3        ! kg m / s^2 
PARATYPE PS_MKS_POUND_FORCE              = 4.44822161526d0   ! kg m / s^2 
PARATYPE PS_MKS_KILOPOUND_FORCE          = 4.44822161526d3   ! kg m / s^2 
PARATYPE PS_MKS_POUNDAL                  = 1.38255d-1        ! kg m / s^2 
PARATYPE PS_MKS_CALORIE                  = 4.1868d0          ! kg m^2 / s^2 
PARATYPE PS_MKS_BTU                      = 1.05505585262d3   ! kg m^2 / s^2 
PARATYPE PS_MKS_THERM                    = 1.05506d8         ! kg m^2 / s^2 
PARATYPE PS_MKS_HORSEPOWER               = 7.457d2           ! kg m^2 / s^3 
PARATYPE PS_MKS_BAR                      = 1.0d5             ! kg / m s^2 
PARATYPE PS_MKS_STD_ATMOSPHERE           = 1.01325d5         ! kg / m s^2 
PARATYPE PS_MKS_TORR                     = 1.33322368421d2   ! kg / m s^2 
PARATYPE PS_MKS_METER_OF_MERCURY         = 1.33322368421d5   ! kg / m s^2 
PARATYPE PS_MKS_INCH_OF_MERCURY          = 3.38638815789d3   ! kg / m s^2 
PARATYPE PS_MKS_INCH_OF_WATER            = 2.490889d2        ! kg / m s^2 
PARATYPE PS_MKS_PSI                      = 6.89475729317d3   ! kg / m s^2 
PARATYPE PS_MKS_POISE                    = 1.0d-1            ! kg m^-1 s^-1 
PARATYPE PS_MKS_STOKES                   = 1.0d-4            ! m^2 / s 
PARATYPE PS_MKS_FARADAY                  = 9.6485341472d4    ! A s / mol 
PARATYPE PS_MKS_ELECTRON_CHARGE          = 1.602176462d-19   ! A s 
PARATYPE PS_MKS_GAUSS                    = 1.0d-4            ! kg / A s^2 
PARATYPE PS_MKS_STILB                    = 1.0d4             ! cd / m^2 
PARATYPE PS_MKS_LUMEN                    = 1.0d0             ! cd sr 
PARATYPE PS_MKS_LUX                      = 1.0d0             ! cd sr / m^2 
PARATYPE PS_MKS_PHOT                     = 1.0d4             ! cd sr / m^2 
PARATYPE PS_MKS_FOOTCANDLE               = 1.076d1           ! cd sr / m^2 
PARATYPE PS_MKS_LAMBERT                  = 1.0d4             ! cd sr / m^2 
PARATYPE PS_MKS_FOOTLAMBERT              = 1.07639104d1      ! cd sr / m^2 
PARATYPE PS_MKS_CURIE                    = 3.7d10            ! 1 / s 
PARATYPE PS_MKS_ROENTGEN                 = 2.58d-4           ! A s / kg 
PARATYPE PS_MKS_RAD                      = 1.0d-2            ! m^2 / s^2 
PARATYPE PS_MKS_SOLAR_MASS               = 1.98892d30        ! kg 
PARATYPE PS_MKS_BOHR_RADIUS              = 5.291772083d-11   ! m 
PARATYPE PS_MKS_NEWTON                   = 1.0d0             ! kg m / s^2 
PARATYPE PS_MKS_DYNE                     = 1.0d-5            ! kg m / s^2 
PARATYPE PS_MKS_JOULE                    = 1.0d0             ! kg m^2 / s^2 
PARATYPE PS_MKS_ERG                      = 1.0d-7            ! kg m^2 / s^2 
PARATYPE PS_MKS_STEFAN_BOLTZMANN_CONST   = 5.67039934436d-8  ! kg/K^4 s^3 
PARATYPE PS_MKS_THOMSON_CROSS_SECTION    = 6.65245853542d-29 ! m^2 
PARATYPE PS_MKS_VACUUM_PERMITTIVITY      = 8.854187817d-12   ! A^2s^4/kg m^3 
PARATYPE PS_MKS_VACUUM_PERMEABILITY      = 1.25663706144d-6  ! kg m/A^2 s^2 

PUBLIC

#undef PARATYPE

END MODULE Units_And_Constants 
