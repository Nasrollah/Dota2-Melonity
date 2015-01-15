#include "pseudopack.h"


MODULE Navier_Stokes_CFL_Conditions

  USE PseudoPack

  USE Temperature
! USE Primitive_Variable

  USE Navier_Stokes_Coefficients

  USE Gas_Properties
  USE Gas_Mixture_Coefficient
  USE Gas_Mixture_Formulation

implicit NONE

INTERFACE PS_Navier_Stokes_Time_Step
  MODULE PROCEDURE Time_Step_1D
  MODULE PROCEDURE Time_Step_2D
  MODULE PROCEDURE Time_Step_3D

  MODULE PROCEDURE Time_Step_1D_G
  MODULE PROCEDURE Time_Step_2D_G
  MODULE PROCEDURE Time_Step_3D_G
END INTERFACE

PRIVATE
PUBLIC  :: PS_Navier_Stokes_Time_Step

CONTAINS
!
!==========================================================================
!
  Subroutine Time_Step_1D         (N0, N2, N3, N5, NV,                   &
                                   x, Q, Gamma,                          &
                                   CFL, dt                               )

  integer  :: N0, N2, N3, N5
  integer  :: NV

  REALTYPE                                  :: CFL, dt
  REALTYPE, dimension(N0:N5)                :: x
  REALTYPE, dimension(N0:N5,NV)             :: Q
  REALTYPE                                  :: Gamma
  REALTYPE, dimension(N0:N5)                :: mu, kappa, zeta
  REALTYPE, dimension(N0:N5)                :: U, P, T

  integer  :: i
  integer  :: N6
  REALTYPE :: dx, dL, dt_NS, dt_mu, dt_kappa, dt_zeta

  if (.NOT. (Viscous_Dissipation .OR. Thermal_Conduction &
                                 .OR. Molecular_Diffusion)) then
    dt = HUGE(ONE) ; RETURN
  endif

  call PS_Primitive_Variable (N0, N5, NV, Gamma, Q, U, P)
  call PS_Temperature        (N0, N5, NV, Q, P, T)

  call PS_Gas_Mixture_Coefficient (Q, T,                                   &
                                     Viscous_Dissipation,    mu_0, mu   ,  &
                                     Thermal_Conduction , kappa_0, kappa,  &
                                   Molecular_Diffusion  ,  zeta_0, zeta    )

  mu = mu/Q(:,1)    ! Kinematic Viscosity

                             dt       = -ONE
                             dt_NS    = -ONE
  if (  Viscous_Dissipation) dt_mu    = -ONE
  if (  Thermal_Conduction ) dt_kappa = -ONE
  if (Molecular_Diffusion  ) dt_zeta  = -ONE

    N6 = N3-1 ; if (N2 == N3) N6 = N2

!$OMP PARALLEL DO REDUCTION (MAX : dt_NS)
      do i = N2,N6
        dx = ABS(x(i+1)-x(i))

        dL = ONE/dx**2 

        if (  Viscous_Dissipation) dt_mu    = MAX(dt_mu   ,    mu(i)*dL)
        if (  Thermal_Conduction ) dt_kappa = MAX(dt_kappa, kappa(i)*dL)
        if (Molecular_Diffusion  ) dt_zeta  = MAX(dt_zeta ,  zeta(i)*dL)
      enddo
!$OMP END PARALLEL DO

  if (  Viscous_Dissipation) dt_NS = MAX(dt_NS, dt_mu   )
  if (  Thermal_Conduction ) dt_NS = MAX(dt_NS, dt_kappa)
  if (Molecular_Diffusion  ) dt_NS = MAX(dt_NS, dt_zeta )

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_NS
  else
    call MPI_ALLREDUCE (dt_NS, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                      MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
    dt = CFL/dt_NS
#endif

  dt = HALF*dt

  END Subroutine Time_Step_1D
!
!==========================================================================
!
  Subroutine Time_Step_2D         (N0, N2, N3, N5,                       &
                                   M0, M2, M3, M5, NV,                   &
                                   x, y, Q, Gamma,                       &
                                   CFL, dt                               )

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5
  integer  :: NV

  REALTYPE                                  :: CFL, dt
  REALTYPE, dimension(N0:N5,M0:M5)          :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,NV)       :: Q
  REALTYPE                                  :: Gamma
  REALTYPE, dimension(N0:N5,M0:M5)          :: mu, kappa, zeta
  REALTYPE, dimension(N0:N5,M0:M5)          :: U, V, P, T

  integer  :: i, j
  integer  :: N6, M6
  REALTYPE :: dx, dy, dL, dt_NS, dt_mu, dt_kappa, dt_zeta

  if (.NOT. (Viscous_Dissipation .OR. Thermal_Conduction &
                                 .OR. Molecular_Diffusion)) then
    dt = HUGE(ONE) ; RETURN
  endif

  call PS_Primitive_Variable (N0, N5, M0, M5, NV, Gamma, Q, U, V, P)
  call PS_Temperature        (N0, N5, M0, M5, NV, Q, P, T)

  call PS_Gas_Mixture_Coefficient (Q, T,                                   &
                                     Viscous_Dissipation,    mu_0, mu   ,  &
                                     Thermal_Conduction , kappa_0, kappa,  &
                                   Molecular_Diffusion  ,  zeta_0, zeta    )

  mu = mu/Q(:,:,1)    ! Kinematic Viscosity

                             dt       = -ONE
                             dt_NS    = -ONE
  if (  Viscous_Dissipation) dt_mu    = -ONE
  if (  Thermal_Conduction ) dt_kappa = -ONE
  if (Molecular_Diffusion  ) dt_zeta  = -ONE

    N6 = N3-1 ; if (N2 == N3) N6 = N2
    M6 = M3-1 ; if (M2 == M3) M6 = M2

!$OMP PARALLEL DO REDUCTION (MAX : dt_NS)
    do j = M2,M6
      do i = N2,N6
        dx = ABS(x(i+1,j  )-x(i,j))
        dy = ABS(y(i  ,j+1)-y(i,j))

        dL = ONE/dx**2 + ONE/dy**2

        if (  Viscous_Dissipation) dt_mu    = MAX(dt_mu   ,    mu(i,j)*dL)
        if (  Thermal_Conduction ) dt_kappa = MAX(dt_kappa, kappa(i,j)*dL)
        if (Molecular_Diffusion  ) dt_zeta  = MAX(dt_zeta ,  zeta(i,j)*dL)
      enddo
    enddo
!$OMP END PARALLEL DO

  if (  Viscous_Dissipation) dt_NS = MAX(dt_NS, dt_mu   )
  if (  Thermal_Conduction ) dt_NS = MAX(dt_NS, dt_kappa)
  if (Molecular_Diffusion  ) dt_NS = MAX(dt_NS, dt_zeta )

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_NS
  else
    call MPI_ALLREDUCE (dt_NS, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                      MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
    dt = CFL/dt_NS
#endif

  dt = HALF*dt

  END Subroutine Time_Step_2D
!
!==========================================================================
!
  Subroutine Time_Step_3D         (N0, N2, N3, N5,                       &
                                   M0, M2, M3, M5,                       &
                                   K0, K2, K3, K5, NV,                   &
                                   x, y, z, Q, Gamma,                    &
                                   CFL, dt                               )

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5
  integer  :: K0, K2, K3, K5
  integer  :: NV

  REALTYPE                                  :: CFL, dt
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE                                  :: Gamma
  REALTYPE, dimension(:,:,:), ALLOCATABLE   :: mu, kappa, zeta
  REALTYPE, dimension(:,:,:), ALLOCATABLE   :: U, V, W, P, T

  integer  :: i, j, k
  integer  :: N6, M6, K6
  REALTYPE :: dx, dy, dz, dL, dt_NS, dt_mu, dt_kappa, dt_zeta

  if (.NOT. (Viscous_Dissipation .OR. Thermal_Conduction &
                                 .OR. Molecular_Diffusion)) then
    dt = HUGE(ONE) ; RETURN
  endif

    ALLOCATE (U(N0:N5,M0:M5,K0:K5), V(N0:N5,M0:M5,K0:K5))
    ALLOCATE (W(N0:N5,M0:M5,K0:K5), P(N0:N5,M0:M5,K0:K5))

  call PS_Primitive_Variable (N0, N5, M0, M5, K0, K5, NV, Gamma, Q, U, V, W, P)

  DEALLOCATE (U, V, W)
    ALLOCATE (T(N0:N5,M0:M5,K0:K5))

  call PS_Temperature        (N0, N5, M0, M5, K0, K5, NV, Q, P, T)

  DEALLOCATE (P)

  if (  Viscous_Dissipation)   ALLOCATE (   mu(N0:N5,M0:M5,K0:K5))
  if (  Thermal_Conduction )   ALLOCATE (kappa(N0:N5,M0:M5,K0:K5))
  if (Molecular_Diffusion  )   ALLOCATE ( zeta(N0:N5,M0:M5,K0:K5))

  call PS_Gas_Mixture_Coefficient (Q, T,                                   &
                                     Viscous_Dissipation,    mu_0, mu   ,  &
                                     Thermal_Conduction , kappa_0, kappa,  &
                                   Molecular_Diffusion  ,  zeta_0, zeta    )

  mu = mu/Q(:,:,:,1)    ! Kinematic Viscosity

  DEALLOCATE (T)

                             dt       = -ONE
                             dt_NS    = -ONE
  if (  Viscous_Dissipation) dt_mu    = -ONE
  if (  Thermal_Conduction ) dt_kappa = -ONE
  if (Molecular_Diffusion  ) dt_zeta  = -ONE

    N6 = N3-1 ; if (N2 == N3) N6 = N2
    M6 = M3-1 ; if (M2 == M3) M6 = M2
    K6 = K3-1 ; if (K2 == K3) K6 = K2

!$OMP PARALLEL DO REDUCTION (MAX : dt_NS)
  do k = K2,K6
    do j = M2,M6
      do i = N2,N6
        dx = ABS(x(i+1,j  ,k  )-x(i,j,k))
        dy = ABS(y(i  ,j+1,k  )-y(i,j,k))
        dz = ABS(z(i  ,j  ,k+1)-z(i,j,k))

        dL = ONE/dx**2 + ONE/dy**2 + ONE/dz**2

        if (  Viscous_Dissipation) dt_mu    = MAX(dt_mu   ,    mu(i,j,k)*dL)
        if (  Thermal_Conduction ) dt_kappa = MAX(dt_kappa, kappa(i,j,k)*dL)
        if (Molecular_Diffusion  ) dt_zeta  = MAX(dt_zeta ,  zeta(i,j,k)*dL)
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

  if (  Viscous_Dissipation) DEALLOCATE (mu)
  if (  Thermal_Conduction ) DEALLOCATE (kappa)
  if (Molecular_Diffusion  ) DEALLOCATE (zeta)

  if (  Viscous_Dissipation) dt_NS = MAX(dt_NS, dt_mu   )

    ! DANGER - KAUSHIK COMMENTED THIS; DO NOT WANT THERMAL TIME SCALE AS IT IS TOO LOW
!  if (  Thermal_Conduction ) dt_NS = MAX(dt_NS, dt_kappa)

  if (Molecular_Diffusion  ) dt_NS = MAX(dt_NS, dt_zeta )

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_NS
  else
    call MPI_ALLREDUCE (dt_NS, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                      MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
    dt = CFL/dt_NS
#endif

  dt = HALF*dt

  END Subroutine Time_Step_3D
!
!==========================================================================
!
  Subroutine Time_Step_1D_G       (N0, N2, N3, N5, NV,                   &
                                   x, Q, Gamma,                          &
                                   CFL, dt                               )

  integer  :: N0, N2, N3, N5
  integer  :: NV

  REALTYPE                                  :: CFL, dt
  REALTYPE, dimension(N0:N5)                :: x
  REALTYPE, dimension(N0:N5,NV)             :: Q
  REALTYPE, dimension(N0:N5)                :: Gamma
  REALTYPE, dimension(N0:N5)                :: mu, kappa, zeta
  REALTYPE, dimension(N0:N5)                :: U, P, T

  integer  :: i
  integer  :: N6
  REALTYPE :: dx, dL, dt_NS, dt_mu, dt_kappa, dt_zeta

  if (.NOT. (Viscous_Dissipation .OR. Thermal_Conduction &
                                 .OR. Molecular_Diffusion)) then
    dt = HUGE(ONE) ; RETURN
  endif

  call PS_Primitive_Variable (N0, N5, NV, Gamma, Q, U, P)
  call PS_Temperature        (N0, N5, NV, Q, P, T)

  call PS_Gas_Mixture_Coefficient (Q, T,                                   &
                                     Viscous_Dissipation,    mu_0, mu   ,  &
                                     Thermal_Conduction , kappa_0, kappa,  &
                                   Molecular_Diffusion  ,  zeta_0, zeta    )

  mu = mu/Q(:,1) ! Kinematic Viscosity

                             dt       = -ONE
                             dt_NS    = -ONE
  if (  Viscous_Dissipation) dt_mu    = -ONE
  if (  Thermal_Conduction ) dt_kappa = -ONE
  if (Molecular_Diffusion  ) dt_zeta  = -ONE

    N6 = N3-1 ; if (N2 == N3) N6 = N2

!$OMP PARALLEL DO REDUCTION (MAX : dt_NS)
      do i = N2,N6
        dx = ABS(x(i+1)-x(i))

        dL = ONE/dx**2 

        if (  Viscous_Dissipation) dt_mu    = MAX(dt_mu   ,    mu(i)*dL)
        if (  Thermal_Conduction ) dt_kappa = MAX(dt_kappa, kappa(i)*dL)
        if (Molecular_Diffusion  ) dt_zeta  = MAX(dt_zeta ,  zeta(i)*dL)
      enddo
!$OMP END PARALLEL DO

  if (  Viscous_Dissipation) dt_NS = MAX(dt_NS, dt_mu   )
  if (  Thermal_Conduction ) dt_NS = MAX(dt_NS, dt_kappa)
  if (Molecular_Diffusion  ) dt_NS = MAX(dt_NS, dt_zeta )

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_NS
  else
    call MPI_ALLREDUCE (dt_NS, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                      MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
    dt = CFL/dt_NS
#endif

  dt = HALF*dt

  END Subroutine Time_Step_1D_G
!
!==========================================================================
!
  Subroutine Time_Step_2D_G       (N0, N2, N3, N5,                       &
                                   M0, M2, M3, M5, NV,                   &
                                   x, y, Q, Gamma,                       &
                                   CFL, dt                               )

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5
  integer  :: NV

  REALTYPE                                  :: CFL, dt
  REALTYPE, dimension(N0:N5,M0:M5)          :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,NV)       :: Q
  REALTYPE, dimension(N0:N5,M0:M5)          :: Gamma
  REALTYPE, dimension(N0:N5,M0:M5)          :: mu, kappa, zeta
  REALTYPE, dimension(N0:N5,M0:M5)          :: U, V, P, T

  integer  :: i, j
  integer  :: N6, M6
  REALTYPE :: dx, dy, dL, dt_NS, dt_mu, dt_kappa, dt_zeta

  if (.NOT. (Viscous_Dissipation .OR. Thermal_Conduction &
                                 .OR. Molecular_Diffusion)) then
    dt = HUGE(ONE) ; RETURN
  endif

  call PS_Primitive_Variable (N0, N5, M0, M5, NV, Gamma, Q, U, V, P)
  call PS_Temperature        (N0, N5, M0, M5, NV, Q, P, T)

  call PS_Gas_Mixture_Coefficient (Q, T,                                   &
                                     Viscous_Dissipation,    mu_0, mu   ,  &
                                     Thermal_Conduction , kappa_0, kappa,  &
                                   Molecular_Diffusion  ,  zeta_0, zeta    )

  mu = mu/Q(:,:,1) ! Kinematic Viscosity

#if defined (DEBUG_HOPE_NS_NU_KAPPA_ZETA)
#include "NS_Debug.i"
#endif

                             dt       = -ONE
                             dt_NS    = -ONE
  if (  Viscous_Dissipation) dt_mu    = -ONE
  if (  Thermal_Conduction ) dt_kappa = -ONE
  if (Molecular_Diffusion  ) dt_zeta  = -ONE

    N6 = N3-1 ; if (N2 == N3) N6 = N2
    M6 = M3-1 ; if (M2 == M3) M6 = M2

!$OMP PARALLEL DO REDUCTION (MAX : dt_NS)
    do j = M2,M6
      do i = N2,N6
        dx = ABS(x(i+1,j  )-x(i,j))
        dy = ABS(y(i  ,j+1)-y(i,j))

        dL = ONE/dx**2 + ONE/dy**2

        if (  Viscous_Dissipation) dt_mu    = MAX(dt_mu   ,    mu(i,j)*dL)
        if (  Thermal_Conduction ) dt_kappa = MAX(dt_kappa, kappa(i,j)*dL)
        if (Molecular_Diffusion  ) dt_zeta  = MAX(dt_zeta ,  zeta(i,j)*dL)
      enddo
    enddo
!$OMP END PARALLEL DO

  if (  Viscous_Dissipation) dt_NS = MAX(dt_NS, dt_mu   )
  if (  Thermal_Conduction ) dt_NS = MAX(dt_NS, dt_kappa)
  if (Molecular_Diffusion  ) dt_NS = MAX(dt_NS, dt_zeta )

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_NS
  else
    call MPI_ALLREDUCE (dt_NS, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                      MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
    dt = CFL/dt_NS
#endif

  dt = HALF*dt

  END Subroutine Time_Step_2D_G
!
!==========================================================================
!
  Subroutine Time_Step_3D_G       (N0, N2, N3, N5,                       &
                                   M0, M2, M3, M5,                       &
                                   K0, K2, K3, K5, NV,                   &
                                   x, y, z, Q, Gamma,                    &
                                   CFL, dt                               )

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5
  integer  :: K0, K2, K3, K5
  integer  :: NV

  REALTYPE                                  :: CFL, dt
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: Gamma
  REALTYPE, dimension(:,:,:), ALLOCATABLE   :: mu, kappa, zeta
  REALTYPE, dimension(:,:,:), ALLOCATABLE   :: U, V, W, P, T

  integer  :: i, j, k
  integer  :: N6, M6, K6
  REALTYPE :: dx, dy, dz, dL, dt_NS, dt_mu, dt_kappa, dt_zeta

  if (.NOT. (Viscous_Dissipation .OR. Thermal_Conduction &
                                 .OR. Molecular_Diffusion)) then
    dt = HUGE(ONE) ; RETURN
  endif

    ALLOCATE (U(N0:N5,M0:M5,K0:K5), V(N0:N5,M0:M5,K0:K5))
    ALLOCATE (W(N0:N5,M0:M5,K0:K5), P(N0:N5,M0:M5,K0:K5))

  call PS_Primitive_Variable (N0, N5, M0, M5, K0, K5, NV, Gamma, Q, U, V, W, P)

  DEALLOCATE (U, V, W)
    ALLOCATE (T(N0:N5,M0:M5,K0:K5))

  call PS_Temperature        (N0, N5, M0, M5, K0, K5, NV, Q, P, T)

  DEALLOCATE (P)

  if (  Viscous_Dissipation)   ALLOCATE (   mu(N0:N5,M0:M5,K0:K5))
  if (  Thermal_Conduction )   ALLOCATE (kappa(N0:N5,M0:M5,K0:K5))
  if (Molecular_Diffusion  )   ALLOCATE ( zeta(N0:N5,M0:M5,K0:K5))

  call PS_Gas_Mixture_Coefficient (Q, T,                                   &
                                     Viscous_Dissipation,    mu_0, mu   ,  &
                                     Thermal_Conduction , kappa_0, kappa,  &
                                   Molecular_Diffusion  ,  zeta_0, zeta    )

  mu = mu/Q(:,:,:,1) ! Kinematic Viscosity

  DEALLOCATE (T)

                             dt       = -ONE
                             dt_NS    = -ONE
  if (  Viscous_Dissipation) dt_mu    = -ONE
  if (  Thermal_Conduction ) dt_kappa = -ONE
  if (Molecular_Diffusion  ) dt_zeta  = -ONE

    N6 = N3-1 ; if (N2 == N3) N6 = N2
    M6 = M3-1 ; if (M2 == M3) M6 = M2
    K6 = K3-1 ; if (K2 == K3) K6 = K2

!$OMP PARALLEL DO REDUCTION (MAX : dt_NS)
  do k = K2,K6
    do j = M2,M6
      do i = N2,N6
        dx = ABS(x(i+1,j  ,k  )-x(i,j,k))
        dy = ABS(y(i  ,j+1,k  )-y(i,j,k))
        dz = ABS(z(i  ,j  ,k+1)-z(i,j,k))

        dL = ONE/dx**2 + ONE/dy**2 + ONE/dz**2

        if (  Viscous_Dissipation) dt_mu    = MAX(dt_mu   ,    mu(i,j,k)*dL)
        if (  Thermal_Conduction ) dt_kappa = MAX(dt_kappa, kappa(i,j,k)*dL)
        if (Molecular_Diffusion  ) dt_zeta  = MAX(dt_zeta ,  zeta(i,j,k)*dL)
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

  if (  Viscous_Dissipation) DEALLOCATE (mu)
  if (  Thermal_Conduction ) DEALLOCATE (kappa)
  if (Molecular_Diffusion  ) DEALLOCATE (zeta)

  if (  Viscous_Dissipation) dt_NS = MAX(dt_NS, dt_mu   )
  if (  Thermal_Conduction ) dt_NS = MAX(dt_NS, dt_kappa)
  if (Molecular_Diffusion  ) dt_NS = MAX(dt_NS, dt_zeta )

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_NS
  else
    call MPI_ALLREDUCE (dt_NS, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                      MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
    dt = CFL/dt_NS
#endif

  dt = HALF*dt

  END Subroutine Time_Step_3D_G

END MODULE Navier_Stokes_CFL_Conditions
