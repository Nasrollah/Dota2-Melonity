#include "macros.h"
#include "dimension.h"

#include "Module_HOPE_File_System.F"

#include "Module_Boundary_Condition_z.F"
#include "Module_Injector_3D.F"

#include "Module_Gas_Database.F"
#include "Module_Gas_Mixture_Formulation.F"
#include "Module_Gas_Mixture_Coefficient.F"

!----
#include "Module_Temperature.F"             

#include "Module_Navier_Stokes_Coefficients.F"

#if defined (NAVIER_STOKES)
#include "Module_Navier_Stokes_Flux_Boundary_Condition.F"

#include "Module_Thermal_Gradient_Boundary_Condition.F"
#include "Module_Molecular_Diffusion_Boundary_Condition.F"
#include "Module_Viscous_Dissipation_Boundary_Condition.F"

#include "Module_Stress_Tensor.F"
#include "Module_Thermal_Gradient.F"
#include "Module_Molecular_Diffusion.F"
#include "Module_Viscous_Dissipation.F"

#include "Module_Navier_Stokes_Flux.F"
#include "Module_Navier_Stokes_CFL_Conditions.F"
#endif

#include "Module_Navier_Stokes_Package.F"

#include "Module_Adjust_D_F_Epsilon.F"
