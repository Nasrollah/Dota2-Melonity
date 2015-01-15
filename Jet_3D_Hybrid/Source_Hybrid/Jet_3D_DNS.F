#include "Jet_3D_DNS.h"

Program Jet_3D_DNS

  USE PSEUDOpack

  USE HOPE_File_System
  USE Gas_Properties
  USE Temperature
  USE Navier_Stokes_Package
  USE Adjust_D_F_Epsilon
  USE Boundary_Condition_z
  USE Injector_3D

IMPLICIT NONE

#include "Jet_3D_DNS.d"

! --------------------------------------------------------------

#if defined (PARALLEL_MPI)
call PS_Initialization
#endif

call File_System_Setup (Dimensions)

!------
! Input Data
!------
call Input (7)

#if defined (WITH_SPECIES_EQUATION)
  Number_Of_Species = 2
#else
  Number_Of_Species = 1
#endif

ND = Dimensions ; NF = ND+2 ; NS = Number_Of_Species 

NV = NF+NS-1 

#include "Processor_Topology.i"
#include "Hybrid_Grid.i"
#include "WENO_BC_Setup.i"  
#include "Variables_Allocation.i"
#include "Case_Study.i"

call Geometry

#if defined (HYBRID)
#include "Central_Difference_Operator_Setup.i"
#endif

#if defined (NAVIER_STOKES)
#include "Navier_Stokes_Difference_Operator_Setup.i"
#endif


        ! KAUSHIK ADDED THIS
       Time = ZERO
       dt = 0.0d0

call Initial_Condition(Time,dt/3.0d0) 

Step = 0 ; Time = ZERO ; call Time_Step (dt, Euler_dt, NS_dt)

if (Restart /= 0) then
#include "Restart_File_Read.i"
endif

dt_Original = dt

PS_Counter(1) = Time/MAX( Restart_Resolution, EPSILON(ONE))
PS_Counter(2) = Time/MAX(   Movie_Resolution, EPSILON(ONE))

if (Time /= Final_Time) then
#include "Initial_Output.i"
endif

#include "CPU_Timing_Output_Header.i"

CPU_Start = TIMER

Do While (Time <  Final_Time)

#include "Time_Step.i"
#include "Adjust_Time_Step.i"
#include "Adjust_MR_Activation.i"

  CPU_RK_Begin = TIMER

  call Runge_Kutta

  CPU_RK_End   = TIMER

#include "Intermediate_Output.i"
enddo

#include "Finalize_Output.i"
#include "Done.i"

CONTAINS
!
! ==================================================
!
#include "Input.F"
#include "Geometry.F"

#include "Initial_Condition.F"
#include "Uniform_Flow.F"

#include "Time_Step.F"
#include "Runge_Kutta.F"

#include "Flux.F"
#include "WENO_Flux.F"

#if defined (HYBRID)
#include "CD_Euler_Flux.F"
#include "Hybrid_Flux.F"
#include "Central_Flux.F"
#include "Rusanov_Flux.F"
#include "Gottlieb_Groth.F"
#include "HLLC.F"
#include "WENO_kaushik.F"
#include "Viscous_Flux_kaushik.F"
#include "WENO_cons_Riemann.F"
#include "les_flux.F"
#include "sponge_bc.F"
#endif

#include "Central_Smoothing.F"

#include "Boundary_Condition.F"
#include "Adjust_Mass_Fraction.F"

#include "Graphics.F"

#if defined (MULTI_RESOLUTION) 
#include "MR_Flag_Graphics.F"
#endif

END Program Jet_3D_DNS
