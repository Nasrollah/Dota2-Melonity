#include "pseudopack.h"


MODULE Gas_Mixture_Formulation

  USE Gas_Properties

implicit NONE

INTERFACE PS_Viscous_Dissipativity
  MODULE PROCEDURE Viscous_Dissipativity_0D
  MODULE PROCEDURE Viscous_Dissipativity_1D
  MODULE PROCEDURE Viscous_Dissipativity_2D
  MODULE PROCEDURE Viscous_Dissipativity_3D
END INTERFACE

INTERFACE PS_Thermal_Conductivity
  MODULE PROCEDURE Thermal_Conductivity_0D
  MODULE PROCEDURE Thermal_Conductivity_1D
  MODULE PROCEDURE Thermal_Conductivity_2D
  MODULE PROCEDURE Thermal_Conductivity_3D
END INTERFACE

INTERFACE PS_Molecular_Diffusivity
  MODULE PROCEDURE Molecular_Diffusivity_0D
  MODULE PROCEDURE Molecular_Diffusivity_1D
  MODULE PROCEDURE Molecular_Diffusivity_2D
  MODULE PROCEDURE Molecular_Diffusivity_3D
END INTERFACE

INTERFACE PS_Species_Coefficient
  MODULE PROCEDURE Species_Coefficient_0D
  MODULE PROCEDURE Species_Coefficient_1D
  MODULE PROCEDURE Species_Coefficient_2D
  MODULE PROCEDURE Species_Coefficient_3D
END INTERFACE

INTERFACE PS_Mixture_Molecular_Weight
  MODULE PROCEDURE Mixture_Molecular_Weight_0D
  MODULE PROCEDURE Mixture_Molecular_Weight_1D
  MODULE PROCEDURE Mixture_Molecular_Weight_2D
  MODULE PROCEDURE Mixture_Molecular_Weight_3D
END INTERFACE

INTERFACE PS_Mixture_Coefficient
  MODULE PROCEDURE Mixture_Coefficient_0D
  MODULE PROCEDURE Mixture_Coefficient_1D
  MODULE PROCEDURE Mixture_Coefficient_2D
  MODULE PROCEDURE Mixture_Coefficient_3D
END INTERFACE

REALTYPE, parameter :: Power = 0.7860d0

PUBLIC

PUBLIC  :: PS_Viscous_Dissipativity
PUBLIC  :: PS_Thermal_Conductivity
PUBLIC  :: PS_Molecular_Diffusivity
PUBLIC  :: PS_Species_Coefficient
PUBLIC  :: PS_Mixture_Molecular_Weight
PUBLIC  :: PS_Mixture_Coefficient

PRIVATE :: Power

CONTAINS
!
!===================================================================
!
  Subroutine Mixture_Molecular_Weight_0D (f_1, MW_1, f_2, MW_2,  MW)

  REALTYPE                   :: MW_1, MW_2
  REALTYPE                   :: f_1, f_2, MW

  MW = ONE/(f_1/MW_1 + f_2/MW_2)

  END Subroutine Mixture_Molecular_Weight_0D
!
!===================================================================
!
  Subroutine Mixture_Molecular_Weight_1D (f_1, MW_1, f_2, MW_2,  MW)

  REALTYPE                   :: MW_1, MW_2
  REALTYPE, dimension(:)     :: f_1, f_2, MW

  MW = ONE/(f_1/MW_1 + f_2/MW_2)

  END Subroutine Mixture_Molecular_Weight_1D
!
!===================================================================
!
  Subroutine Mixture_Molecular_Weight_2D (f_1, MW_1, f_2, MW_2,  MW)

  REALTYPE                   :: MW_1, MW_2
  REALTYPE, dimension(:,:)   :: f_1, f_2, MW

  MW = ONE/(f_1/MW_1 + f_2/MW_2)

  END Subroutine Mixture_Molecular_Weight_2D
!
!===================================================================
!
  Subroutine Mixture_Molecular_Weight_3D (f_1, MW_1, f_2, MW_2,  MW)

  REALTYPE                   :: MW_1, MW_2
  REALTYPE, dimension(:,:,:) :: f_1, f_2, MW

  MW = ONE/(f_1/MW_1 + f_2/MW_2)

  END Subroutine Mixture_Molecular_Weight_3D

!
!===================================================================
!
  Subroutine Viscous_Dissipativity_0D (T_0, mu_0, T, mu)
 
  REALTYPE                   :: T_0, mu_0
  REALTYPE                   :: T  , mu

  mu = mu_0*(T/T_0)**Power

  END Subroutine Viscous_Dissipativity_0D
!
!===================================================================
!
  Subroutine Viscous_Dissipativity_1D (T_0, mu_0, T, mu)
 
  REALTYPE                   :: T_0, mu_0
  REALTYPE, dimension(:)     :: T  , mu

  mu = mu_0*(T/T_0)**Power

  END Subroutine Viscous_Dissipativity_1D
!
!===================================================================
!
  Subroutine Viscous_Dissipativity_2D (T_0, mu_0, T, mu)
 
  REALTYPE                   :: T_0, mu_0
  REALTYPE, dimension(:,:)   :: T  , mu

  mu = mu_0*(T/T_0)**Power

  END Subroutine Viscous_Dissipativity_2D
!
!===================================================================
!
  Subroutine Viscous_Dissipativity_3D (T_0, mu_0, T, mu)
 
  REALTYPE                   :: T_0, mu_0
  REALTYPE, dimension(:,:,:) :: T  , mu

  mu = mu_0*(T/T_0)**Power

  END Subroutine Viscous_Dissipativity_3D

!
!===================================================================
!
  Subroutine Thermal_Conductivity_0D (T_0, kappa_0, T, kappa)

  REALTYPE                   :: T_0, kappa_0
  REALTYPE                   :: T  , kappa

  kappa = kappa_0*(T/T_0)**Power

  END Subroutine Thermal_Conductivity_0D
!
!===================================================================
!
  Subroutine Thermal_Conductivity_1D (T_0, kappa_0, T, kappa)

  REALTYPE                   :: T_0, kappa_0
  REALTYPE, dimension(:)     :: T  , kappa

  kappa = kappa_0*(T/T_0)**Power

  END Subroutine Thermal_Conductivity_1D
!
!===================================================================
!
  Subroutine Thermal_Conductivity_2D (T_0, kappa_0, T, kappa)

  REALTYPE                   :: T_0, kappa_0
  REALTYPE, dimension(:,:)   :: T  , kappa

  kappa = kappa_0*(T/T_0)**Power

  END Subroutine Thermal_Conductivity_2D
!
!===================================================================
!
  Subroutine Thermal_Conductivity_3D (T_0, kappa_0, T, kappa)

  REALTYPE                   :: T_0, kappa_0
  REALTYPE, dimension(:,:,:) :: T  , kappa

  kappa = kappa_0*(T/T_0)**Power

  END Subroutine Thermal_Conductivity_3D

!
!===================================================================
!
  Subroutine Molecular_Diffusivity_0D (T_0, zeta_0, T, zeta)

  REALTYPE                   :: T_0, zeta_0
  REALTYPE                   :: T  , zeta

  zeta = zeta_0*(T/T_0)**Power

  END Subroutine Molecular_Diffusivity_0D
!
!===================================================================
!
  Subroutine Molecular_Diffusivity_1D (T_0, zeta_0, T, zeta)

  REALTYPE                   :: T_0, zeta_0
  REALTYPE, dimension(:)     :: T  , zeta

  zeta = zeta_0*(T/T_0)**Power

  END Subroutine Molecular_Diffusivity_1D
!
!===================================================================
!
  Subroutine Molecular_Diffusivity_2D (T_0, zeta_0, T, zeta)

  REALTYPE                   :: T_0, zeta_0
  REALTYPE, dimension(:,:)   :: T  , zeta

  zeta = zeta_0*(T/T_0)**Power

  END Subroutine Molecular_Diffusivity_2D
!
!===================================================================
!
  Subroutine Molecular_Diffusivity_3D (T_0, zeta_0, T, zeta)

  REALTYPE                   :: T_0, zeta_0
  REALTYPE, dimension(:,:,:) :: T  , zeta

  zeta = zeta_0*(T/T_0)**Power

  END Subroutine Molecular_Diffusivity_3D

!
!===================================================================
!
  Subroutine Species_Coefficient_0D (T_0, mu_0, T, mu)

  REALTYPE                   :: T_0, mu_0
  REALTYPE                   :: T  , mu

  mu = mu_0*(T/T_0)**Power

  END Subroutine Species_Coefficient_0D
!
!===================================================================
!
  Subroutine Species_Coefficient_1D (T_0, mu_0, T, mu)

  REALTYPE                   :: T_0, mu_0
  REALTYPE, dimension(:)     :: T  , mu

  mu = mu_0*(T/T_0)**Power

  END Subroutine Species_Coefficient_1D
!
!===================================================================
!
  Subroutine Species_Coefficient_2D (T_0, mu_0, T, mu)

  REALTYPE                   :: T_0, mu_0
  REALTYPE, dimension(:,:)   :: T  , mu

  mu = mu_0*(T/T_0)**Power

  END Subroutine Species_Coefficient_2D
!
!===================================================================
!
  Subroutine Species_Coefficient_3D (T_0, mu_0, T, mu)

  REALTYPE                   :: T_0, mu_0
  REALTYPE, dimension(:,:,:) :: T  , mu

  mu = mu_0*(T/T_0)**Power

  END Subroutine Species_Coefficient_3D

!
!===================================================================
!
  Subroutine Mixture_Coefficient_0D (mu_1, f_1, MW_1,     &
                                     mu_2, f_2, MW_2,  mu)

  REALTYPE               :: MW_1, MW_2
  REALTYPE               :: f_1, f_2, mu_1, mu_2, mu

  integer                             :: N
  REALTYPE                            :: A, B, C

  A = f_1*SQRT(MW_1) ; B = f_2*SQRT(MW_2) ; C = A+B ; A = A/C ; B = B/C

  mu = A*mu_1 + B*mu_2

  END Subroutine Mixture_Coefficient_0D 
!
!===================================================================
!
  Subroutine Mixture_Coefficient_1D (mu_1, f_1, MW_1,     &
                                     mu_2, f_2, MW_2,  mu)

  REALTYPE               :: MW_1, MW_2
  REALTYPE, dimension(:) :: f_1, f_2, mu_1, mu_2, mu

  REALTYPE :: SQRT_MW_1, SQRT_MW_2
  REALTYPE, dimension(SIZE(mu,DIM=1)) :: A, B, C

  SQRT_MW_1 = SQRT(MW_1) ; SQRT_MW_2 = SQRT(MW_2)

  A = f_1*SQRT_MW_1 ; B = f_2*SQRT_MW_2 ; C = A+B ; A = A/C ; B = B/C

  mu = A*mu_1 + B*mu_2

  END Subroutine Mixture_Coefficient_1D 
!
!===================================================================
!
  Subroutine Mixture_Coefficient_2D (mu_1, f_1, MW_1,     &
                                     mu_2, f_2, MW_2,  mu)

  REALTYPE                 :: MW_1, MW_2
  REALTYPE, dimension(:,:) :: f_1, f_2, mu_1, mu_2, mu

  REALTYPE :: SQRT_MW_1, SQRT_MW_2
  REALTYPE, dimension(SIZE(mu,DIM=1),SIZE(mu,DIM=2)) :: A, B, C

  SQRT_MW_1 = SQRT(MW_1) ; SQRT_MW_2 = SQRT(MW_2)

  A = f_1*SQRT_MW_1 ; B = f_2*SQRT_MW_2 ; C = A+B ; A = A/C ; B = B/C

  mu = A*mu_1 + B*mu_2

  END Subroutine Mixture_Coefficient_2D 
!
!===================================================================
!
  Subroutine Mixture_Coefficient_3D (mu_1, f_1, MW_1,     &
                                     mu_2, f_2, MW_2,  mu)

  REALTYPE                   :: MW_1, MW_2
  REALTYPE, dimension(:,:,:) :: f_1, f_2, mu_1, mu_2, mu

  integer                                 :: N, M, K

  REALTYPE :: SQRT_MW_1, SQRT_MW_2
  REALTYPE, dimension(SIZE(mu,DIM=1),SIZE(mu,DIM=2),SIZE(mu,DIM=3)) :: A, B, C

  SQRT_MW_1 = SQRT(MW_1) ; SQRT_MW_2 = SQRT(MW_2)

  A = f_1*SQRT_MW_1 ; B = f_2*SQRT_MW_2 ; C = A+B ; A = A/C ; B = B/C

  mu = A*mu_1 + B*mu_2

  END Subroutine Mixture_Coefficient_3D 

END MODULE Gas_Mixture_Formulation
