#include "pseudopack.h"


MODULE Thermal_Gradient_Boundary_Condition

  USE PseudoPack

! USE Boundary_Condition_Parameter
  USE Navier_Stokes_Flux_Boundary_Condition

implicit NONE

INTERFACE PS_Thermal_Gradient_BC
  MODULE PROCEDURE Thermal_Gradient_1D
  MODULE PROCEDURE Thermal_Gradient_2D
  MODULE PROCEDURE Thermal_Gradient_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Thermal_Gradient_BC

CONTAINS
!
!===================================================================
!
  Subroutine Thermal_Gradient_1D (N0, N5, N2, N3,      &
                                  Grad_T, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,1)                      :: Grad_T
  REALTYPE, dimension(3,1)                          :: BC_IEO

  BC_IEO =  ONE 

  call PS_Flux_Boundary_Condition (N0, N5, N2, N3,   1, Grad_T,                &
                                   BC_On_Off, BC_Distributed, BC_Order,        &
                                   BC_Method, BC_Symmetry   , BC_Type, BC_IEO, &
                                   MPI_Comm_Type)

  END Subroutine Thermal_Gradient_1D 
!
!===================================================================
!
  Subroutine Thermal_Gradient_2D (N0, N5, N2, N3,      &
                                  M0, M5, M2, M3,      &
                                  Grad_T, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,2)                :: Grad_T
  REALTYPE, dimension(3,2)                          :: BC_IEO

  BC_IEO =  ONE 

  call PS_Flux_Boundary_Condition (N0, N5, N2, N3,                             &
                                   M0, M5, M2, M3,   2, Grad_T,                &
                                   BC_On_Off, BC_Distributed, BC_Order,        &
                                   BC_Method, BC_Symmetry   , BC_Type, BC_IEO, &
                                   MPI_Comm_Type)

  END Subroutine Thermal_Gradient_2D 
!
!===================================================================
!
  Subroutine Thermal_Gradient_3D (N0, N5, N2, N3,      &
                                  M0, M5, M2, M3,      &
                                  K0, K5, K2, K3,      &
                                  Grad_T, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,3)          :: Grad_T
  REALTYPE, dimension(3,3)                          :: BC_IEO

  BC_IEO =  ONE 

  call PS_Flux_Boundary_Condition (N0, N5, N2, N3,                             &
                                   M0, M5, M2, M3,                             &
                                   K0, K5, K2, K3,   3, Grad_T,                &
                                   BC_On_Off, BC_Distributed, BC_Order,        &
                                   BC_Method, BC_Symmetry   , BC_Type, BC_IEO, &
                                   MPI_Comm_Type)

  END Subroutine Thermal_Gradient_3D

END MODULE Thermal_Gradient_Boundary_Condition
