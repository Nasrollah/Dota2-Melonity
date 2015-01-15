#include "pseudopack.h"


MODULE Navier_Stokes_Package

  USE Navier_Stokes_Coefficients

#if defined (NAVIER_STOKES)
  USE Thermal_Gradient_Boundary_Condition
  USE Thermal_Gradient_module

  USE Viscous_Dissipation_Boundary_Condition
  USE Viscous_Dissipation_module

  USE Molecular_Diffusion_Boundary_Condition
  USE Molecular_Diffusion_module

  USE Navier_Stokes_Flux

  USE Navier_Stokes_CFL_Conditions
#endif

END MODULE Navier_Stokes_Package
