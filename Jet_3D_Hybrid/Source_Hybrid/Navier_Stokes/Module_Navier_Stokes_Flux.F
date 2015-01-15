#include "pseudopack.h"


MODULE Navier_Stokes_Flux

  USE PseudoPack

  USE Navier_Stokes_Coefficients, ONLY:   Thermal_Conduction 
  USE Navier_Stokes_Coefficients, ONLY:   Viscous_Dissipation
  USE Navier_Stokes_Coefficients, ONLY: Molecular_Diffusion  

  USE Navier_Stokes_Coefficients, ONLY: CD_TC_x, CD_TC_y, CD_TC_z
  USE Navier_Stokes_Coefficients, ONLY: CD_VD_x, CD_VD_y, CD_VD_z
  USE Navier_Stokes_Coefficients, ONLY: CD_MD_x, CD_MD_y, CD_MD_z

  USE Gas_Properties            , ONLY: mu_0, kappa_0, zeta_0
  USE Gas_Mixture_Coefficient

  USE   Thermal_Gradient_module
  USE   Viscous_Dissipation_module
  USE Molecular_Diffusion_module

implicit NONE

INTERFACE PS_CD_NS_Flux
  MODULE PROCEDURE Flux_1D
  MODULE PROCEDURE Flux_2D
  MODULE PROCEDURE Flux_3D
END INTERFACE

PRIVATE

PUBLIC :: PS_CD_NS_Flux

CONTAINS
!
!===================================================================
!
  Subroutine Flux_1D (N0, N5, N2, N3,                &
                      NV, NF, NS,                    &
                      Q, U, T, D_Flux, D_F_Epsilon,  &
                      MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: NV, NF, NS

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,NV)             :: Q, D_Flux
  REALTYPE, dimension(N0:N5)                :: U, T
  REALTYPE, dimension(N0:N5)                :: mu, kappa, zeta
  REALTYPE, dimension(N0:N5,1)              :: Grad_T

  REALTYPE, OPTIONAL :: D_F_Epsilon

  if (.NOT. (Viscous_Dissipation .OR. Thermal_Conduction &
                                 .OR. Molecular_Diffusion)) RETURN

  D_Flux = ZERO 

  call PS_Gas_Mixture_Coefficient (Q, T,                                   &
                                     Viscous_Dissipation,    mu_0, mu   ,  &
                                     Thermal_Conduction , kappa_0, kappa,  &
                                   Molecular_Diffusion  ,  zeta_0, zeta    )

  if (Thermal_Conduction) then
      call PS_Thermal_Gradient    (N0, N5, N2, N3, CD_TC_x,     &
                                   kappa, T, Grad_T,            &
                                   MPI_Comm_Type)

    if (Viscous_Dissipation) &
      call PS_Viscous_Dissipation (N0, N5, N2, N3, CD_VD_x,     &
                                   NV, mu, U, D_Flux, Grad_T,   &
                                   MPI_Comm_Type=MPI_Comm_Type)
  else
    if (Viscous_Dissipation) &
      call PS_Viscous_Dissipation (N0, N5, N2, N3, CD_VD_x,     &
                                   NV, mu, U, D_Flux,           &
                                   MPI_Comm_Type=MPI_Comm_Type)
  endif

    if (Molecular_Diffusion) &
      call PS_Molecular_Diffusion (N0, N5, N2, N3, CD_MD_x,     &
                                   NV, NF, NS, zeta, Q, D_Flux, &
                                   MPI_Comm_Type)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO

  END Subroutine Flux_1D
!
!===================================================================
!
  Subroutine Flux_2D (N0, N5, N2, N3,                  &
                      M0, M5, M2, M3,                  &
                      NV, NF, NS,                      &
                      Q, U, V, T, D_Flux, D_F_Epsilon, &
                      MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: NV, NF, NS

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)       :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5)          :: U, V, T
  REALTYPE, dimension(N0:N5,M0:M5)          :: mu, kappa, zeta
  REALTYPE, dimension(N0:N5,M0:M5,2)        :: Grad_T

  REALTYPE, OPTIONAL :: D_F_Epsilon

  REALTYPE :: CPU0, CPU1, CPU2, CPU3

  if (.NOT. (Viscous_Dissipation .OR. Thermal_Conduction &
                                 .OR. Molecular_Diffusion)) RETURN

  D_Flux = ZERO 

!CPU0 = TIMER

  call PS_Gas_Mixture_Coefficient (Q, T,                                   &
                                     Viscous_Dissipation,    mu_0, mu   ,  &
                                     Thermal_Conduction , kappa_0, kappa,  &
                                   Molecular_Diffusion  ,  zeta_0, zeta    )

!CPU1 = TIMER
!write (6,*) 'Gas Mixture :', CPU1-CPU0

  if (Thermal_Conduction) then
!CPU0 = TIMER
      call PS_Thermal_Gradient    (N0, N5, N2, N3, CD_TC_x,      &
                                   M0, M5, M2, M3, CD_TC_y,      &
                                   kappa, T, Grad_T,             &
                                   MPI_Comm_Type)
!CPU1 = TIMER
!write (6,*) 'Thermal_Gradient :', CPU1-CPU0
!CPU0 = TIMER
    if (Viscous_Dissipation) &
      call PS_Viscous_Dissipation (N0, N5, N2, N3, CD_VD_x,      &
                                   M0, M5, M2, M3, CD_VD_y,      &
                                   NV, mu, U, V, D_Flux, Grad_T, &
                                   MPI_Comm_Type=MPI_Comm_Type)
  else
!CPU0 = TIMER
    if (Viscous_Dissipation) &
      call PS_Viscous_Dissipation (N0, N5, N2, N3, CD_VD_x,      &
                                   M0, M5, M2, M3, CD_VD_y,      &
                                   NV, mu, U, V, D_Flux,         &
                                   MPI_Comm_Type=MPI_Comm_Type)
  endif
!CPU1 = TIMER
!write (6,*) 'Viscous Dissipation :', CPU1-CPU0
!CPU0 = TIMER  
  if (Molecular_Diffusion) &
      call PS_Molecular_Diffusion (N0, N5, N2, N3, CD_MD_x,      &
                                   M0, M5, M2, M3, CD_MD_y,      &
                                   NV, NF, NS, zeta, Q, D_Flux,  &
                                   MPI_Comm_Type)
!CPU1 = TIMER
!write (6,*) 'Molecular Diffusion :', CPU1-CPU0
!CPU0 = TIMER
  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO
!CPU1 = TIMER
!write (6,*) 'D_F_Epsilon : ', CPU1-CPU0

  END Subroutine Flux_2D
!
!===================================================================
!
  Subroutine Flux_3D (N0, N5, N2, N3,                     &
                      M0, M5, M2, M3,                     &
                      K0, K5, K2, K3,                     &
                      NV, NF, NS,                         &
                      Q, U, V, W, T, D_Flux, D_F_Epsilon, &
                      MPI_Comm_Type)

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3
  integer  :: K0, K5, K2, K3
  integer  :: NV, NF, NS

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, T
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: mu, kappa, zeta
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,3)  :: Grad_T

  REALTYPE, OPTIONAL :: D_F_Epsilon

  if (.NOT. (Viscous_Dissipation .OR. Thermal_Conduction &
                                 .OR. Molecular_Diffusion)) RETURN

  D_Flux = ZERO 

  call PS_Gas_Mixture_Coefficient (Q, T,                                   &
                                     Viscous_Dissipation,    mu_0, mu   ,  &
                                     Thermal_Conduction , kappa_0, kappa,  &
                                   Molecular_Diffusion  ,  zeta_0, zeta    )

  if (Thermal_Conduction) then
      call PS_Thermal_Gradient    (N0, N5, N2, N3, CD_TC_x,         &
                                   M0, M5, M2, M3, CD_TC_y,         &
                                   K0, K5, K2, K3, CD_TC_z,         &
                                   kappa, T, Grad_T,                &
                                   MPI_Comm_Type)

    if (Viscous_Dissipation) &
      call PS_Viscous_Dissipation (N0, N5, N2, N3, CD_VD_x,         &
                                   M0, M5, M2, M3, CD_VD_y,         &
                                   K0, K5, K2, K3, CD_VD_z,         &
                                   NV, mu, U, V, W, D_Flux, Grad_T, &
                                   MPI_Comm_Type=MPI_Comm_Type)
  else
    if (Viscous_Dissipation) &
      call PS_Viscous_Dissipation (N0, N5, N2, N3, CD_VD_x,         &
                                   M0, M5, M2, M3, CD_VD_y,         &
                                   K0, K5, K2, K3, CD_VD_z,         &
                                   NV, mu, U, V, W, D_Flux,         &
                                   MPI_Comm_Type=MPI_Comm_Type)
  endif

    if (Molecular_Diffusion) &
      call PS_Molecular_Diffusion (N0, N5, N2, N3, CD_MD_x,         &
                                   M0, M5, M2, M3, CD_MD_y,         &
                                   K0, K5, K2, K3, CD_MD_z,         &
                                   NV, NF, NS, zeta, Q, D_Flux,     &
                                   MPI_Comm_Type)

  if (PRESENT(D_F_Epsilon)) WHERE (ABS(D_Flux) <  D_F_Epsilon) D_Flux = ZERO

  END Subroutine Flux_3D

END MODULE Navier_Stokes_Flux
