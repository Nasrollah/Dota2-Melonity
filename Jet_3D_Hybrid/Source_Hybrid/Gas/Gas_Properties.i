REALTYPE, parameter :: Translation_Factor = 29.23950d0/1.010d0 ! KJ/kg to J/mol

! ============= Air_Acetone ===================

REALTYPE ::    Pr_AA  = HUGE(ONE)
REALTYPE ::    MW_AA  =  35.9250d0             ! g/g-mole
REALTYPE ::   Rho_AA  = 1.3510d-3              ! g/cm^3 = 1000 g/l
REALTYPE ::    mu_AA  = HUGE(ONE)     ! m^2 s^-1
REALTYPE :: kappa_AA  = HUGE(ONE)
REALTYPE ::  zeta_AA  = HUGE(ONE)     ! m^-1 s^-1
REALTYPE :: Gamma_AA  = 1.2760d0
REALTYPE ::     f_AA  = 1.0d0

REALTYPE ::    Cp_AA  = HUGE(ONE) ! = f_AA_Acetone*Cp_Acetone + f_AA_Air*Cp_Air
REALTYPE ::    Cv_AA  = HUGE(ONE) ! = Cp_AA/Gamma_AA

REALTYPE ::  Cp_Acetone    = 1.470d0*Translation_Factor

#if 1
REALTYPE ::   f_AA_Acetone    = 0.250d0
REALTYPE ::   f_AA_Air        = 0.750d0
#else
REALTYPE ::   f_AA_Acetone    = 0.500d0
REALTYPE ::   f_AA_Air        = 0.500d0
#endif

! ============= Air ===================

REALTYPE ::    Pr_Air =   0.710d0
REALTYPE ::    MW_Air =  29.040d0
REALTYPE ::   Rho_Air =   1.180d0     ! kg kmol^-1
REALTYPE ::    mu_Air =  15.70d-6     ! m^2 s^-1
REALTYPE :: kappa_Air = HUGE(ONE)
REALTYPE ::  zeta_Air =  20.40d0      ! m^-1 s^-1
REALTYPE :: Gamma_Air =  1.400d0      
REALTYPE ::     f_Air = 1.0d0

REALTYPE ::    Cp_Air     = 29.2395                ! J/(mol.K)
REALTYPE ::    Cv_Air     = 20.844                 ! J/(mol.K)

#if 0
REALTYPE ::    MW_Air     =  39.9480d0             ! g/g-mole
REALTYPE ::   Rho_Air     = 0.890d-3               ! g/cm^3
REALTYPE :: Gamma_Air     = 1.40280d0
#endif

! ============= SF6 ===================

REALTYPE ::    Pr_SF6 =   0.900d0
REALTYPE ::    MW_SF6 = 146.070d0
REALTYPE ::   Rho_SF6 =   5.970d0     ! kg kmol^-1
REALTYPE ::    mu_SF6 =   2.470d-6    ! m^2 s^-1
REALTYPE :: kappa_SF6 = HUGE(ONE)
REALTYPE ::  zeta_SF6 =   9.70d0      ! m^-2 s^-1
REALTYPE :: Gamma_SF6 =  1.090d0      
REALTYPE ::     f_SF6     = 1.0d0

REALTYPE ::    Cp_SF6     = 97.85350d0             ! J/mol.K
REALTYPE ::    Cv_SF6     = 89.52750d0             ! J/mol.K

! ============= Xenon ===================

REALTYPE ::    MW_Xe      = 131.2900d0             ! g/g-mole
REALTYPE ::   Rho_Xe      = 2.950d-3               ! g/cm^3
!REALTYPE :: Gamma_Xe      = 5.0d0/3.0d0
REALTYPE :: Gamma_Xe      = 1.6670d0
REALTYPE ::     f_Xe      = 1.0d0

REALTYPE ::   Cp_Xe       = 35.31970d0             ! J/(mol.K)
!REALTYPE ::   Cv_Xe       = 12.60480d0             ! J/(mol.K)
REALTYPE ::   Cv_Xe       = 21.18760d0             ! J/(mol.K)

! ============= Argon ===================

REALTYPE ::    MW_Ar      =  39.9480d0             ! g/g-mole
REALTYPE ::   Rho_Ar      = 0.890d-3               ! g/cm^3
!REALTYPE :: Gamma_Ar      = 5.0d0/3.0d0
REALTYPE :: Gamma_Ar      = 1.6640d0
REALTYPE ::     f_Ar      = 1.0d0

REALTYPE ::    Cp_Ar      = 20.80d0                ! J/mol.K
REALTYPE ::    Cv_Ar      = 12.50d0                ! J/mol.K

! ============= CHBr ===================

REALTYPE ::    MW_CHBr    = 12.0d0                 ! g/g-mole
REALTYPE ::   Rho_CHBr    = 1.420d-3               ! g/cm^3 = 1000 g/l
REALTYPE :: Gamma_CHBr    = 5.0d0/3.0d0
REALTYPE ::     f_CHBr    = 1.0d0

REALTYPE ::    Cv_CHBr    = 2.590d0*KJKg_2_Jmole    ! J/mol.K = 74.98050d0     
REALTYPE ::    Cp_CHBr    = 124.96750d0             ! J/mol.K

! ============= CRF ===================

REALTYPE ::    MW_CRF     = 12.0d0                 ! g/g-mole
REALTYPE ::   Rho_CRF     = 0.050d-3               ! g/cm^3 = 1000 g/l
REALTYPE :: Gamma_CRF     = 1.46d0
REALTYPE ::     f_CRF     = 1.0d0

REALTYPE ::    Cv_CRF     = 2.590d0*KJKg_2_Jmole    ! J/mol.K = 74.98050d0
REALTYPE ::    Cp_CRF     = 109.471530d0            ! J/mol.K

! ============= O2 ===================

REALTYPE ::    MW_O2      = 32.00000d0             ! g/g-mole
REALTYPE ::   Rho_O2      = 4.4750d-3              ! g/cm^3
!REALTYPE :: Gamma_O2      = 5.0d0/3.0d0
REALTYPE :: Gamma_O2      = 1.3933650d0
REALTYPE ::     f_O2      = 0.23140d0

REALTYPE ::    Cp_O2      = 29.40d0                ! J/(mol.K)
REALTYPE ::    Cv_O2      = 21.10d0                ! J/(mol.K)

! ============= H2 ===================

REALTYPE ::    MW_H2      =  2.01600d0             ! g/g-mole
REALTYPE ::   Rho_H2      =  0.0d0                 ! g/cm^3
!REALTYPE :: Gamma_H2      = 5.0d0/3.0d0
REALTYPE :: Gamma_H2      = 1.384259
REALTYPE ::     f_H2      = 1.0d0

REALTYPE ::    Cp_H2      = 29.90d0                ! J/(mol.K)
REALTYPE ::    Cv_H2      = 21.60d0                ! J/(mol.K)

! ============= N2 ===================

REALTYPE ::    MW_N2      = 28.01600d0             ! g/g-mole
REALTYPE ::   Rho_N2      = 4.6140d-3              ! g/cm^3
!REALTYPE :: Gamma_N2      = 5.0d0/3.0d0
REALTYPE :: Gamma_N2      = 1.4038460d0
REALTYPE ::     f_N2      = 0.76860d0

REALTYPE ::    Cp_N2      = 29.20d0                ! J/(mol.K)
REALTYPE ::    Cv_N2      = 20.80d0                ! J/(mol.K)

! ***************************************************************************
