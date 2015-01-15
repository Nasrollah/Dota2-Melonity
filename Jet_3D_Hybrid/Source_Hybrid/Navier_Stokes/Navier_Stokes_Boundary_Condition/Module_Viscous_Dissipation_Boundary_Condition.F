#include "pseudopack.h"


MODULE Viscous_Dissipation_Boundary_Condition

  USE PseudoPack

! USE Boundary_Condition_Parameter
  USE Navier_Stokes_Flux_Boundary_Condition

implicit NONE

INTERFACE PS_Viscous_Dissipation_BC
  MODULE PROCEDURE Stress_Tensor_1D
  MODULE PROCEDURE Stress_Tensor_2D
  MODULE PROCEDURE Stress_Tensor_3D
END INTERFACE

INTERFACE PS_Stress_Tensor_BC
  MODULE PROCEDURE Stress_Tensor_1D
  MODULE PROCEDURE Stress_Tensor_2D
  MODULE PROCEDURE Stress_Tensor_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Stress_Tensor_BC
PUBLIC :: PS_Viscous_Dissipation_BC

CONTAINS
!
!===================================================================
!
  Subroutine Stress_Tensor_1D (N0, N5, N2, N3,         &
                               Tau, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,1)                      :: Tau
  REALTYPE, dimension(3,1)                          :: BC_IEO

  BC_IEO =  ONE 

  call PS_Flux_Boundary_Condition (N0, N5, N2, N3,   1, Tau,                   &
                                   BC_On_Off, BC_Distributed, BC_Order,        &
                                   BC_Method, BC_Symmetry   , BC_Type, BC_IEO, &
                                   MPI_Comm_Type)

  END Subroutine Stress_Tensor_1D 
!
!===================================================================
!
  Subroutine Stress_Tensor_2D (N0, N5, N2, N3,         &
                               M0, M5, M2, M3,         &
                               Tau, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,3)                :: Tau
  REALTYPE, dimension(3,3)                          :: BC_IEO

  BC_IEO =  ONE ; BC_IEO(1,3) = -ONE  ; BC_IEO(2,3) = -ONE

  call PS_Flux_Boundary_Condition (N0, N5, N2, N3,                             &
                                   M0, M5, M2, M3,   3, Tau,                   &
                                   BC_On_Off, BC_Distributed, BC_Order,        &
                                   BC_Method, BC_Symmetry   , BC_Type, BC_IEO, &
                                   MPI_Comm_Type)

  END Subroutine Stress_Tensor_2D 
!
!===================================================================
!
  Subroutine Stress_Tensor_3D (N0, N5, N2, N3,         &
                               M0, M5, M2, M3,         &
                               K0, K5, K2, K3,         &
                               Tau, MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,6)          :: Tau
  REALTYPE, dimension(3,6)                          :: BC_IEO

  BC_IEO =  ONE ; BC_IEO(1,4) = -ONE ; BC_IEO(1,5) = -ONE
                  BC_IEO(2,4) = -ONE ; BC_IEO(2,6) = -ONE
                  BC_IEO(3,5) = -ONE ; BC_IEO(3,6) = -ONE

  call PS_Flux_Boundary_Condition (N0, N5, N2, N3,                             &
                                   M0, M5, M2, M3,                             &
                                   K0, K5, K2, K3,   6, Tau,                   &
                                   BC_On_Off, BC_Distributed, BC_Order,        &
                                   BC_Method, BC_Symmetry   , BC_Type, BC_IEO, &
                                   MPI_Comm_Type)

  END Subroutine Stress_Tensor_3D

END MODULE Viscous_Dissipation_Boundary_Condition
