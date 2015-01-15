#include "pseudopack.h"

!*******************************************************************
! Low Storage 4 stages 3th order Runge Kutta Scheme by 
!   Carpenter and Kennedy
!
!   do m = 1,4
!     Q1 = A(m)*Q1 + dt*DFLUX(Q)
!     Q  =      Q  + B(m)*Q1
!   enddo
!
!   Q1 is a temporary storage  (Don't use it to compute the DFlux)
!   Q  is the intermediate solution (use this to compute the DFlux) at stage m
!*******************************************************************


MODULE Runge_Kutta_LS4_3

implicit NONE

REALTYPE, dimension(1:4), parameter ::              &
     A = (/   0.0d0       , &
            -11.0d0/15.0d0, &
             -5.0d0/ 3.0d0, &
             -1.0d0        /)

REALTYPE, dimension(1:4), parameter ::              &
     B = (/   1.0d0/ 3.0d0, &
              5.0d0/ 6.0d0, &
              3.0d0/ 5.0d0, &
             -1.0d0/ 4.0d0 /)
                                                 
REALTYPE, dimension(0:4), parameter ::              &
     C = (/   0.0d0       , &
              3.0d0/ 9.0d0, &
              5.0d0/ 9.0d0, &
              8.0d0/ 9.0d0, &
              1.0d0        /)
                                                 
INTERFACE PS_Runge_Kutta_LS4_3
  MODULE PROCEDURE Runge_Kutta_0D_LS4_3G
  MODULE PROCEDURE Runge_Kutta_1D_LS4_3G
  MODULE PROCEDURE Runge_Kutta_2D_LS4_3G
  MODULE PROCEDURE Runge_Kutta_3D_LS4_3G

  MODULE PROCEDURE Runge_Kutta_0D_LS4_3
  MODULE PROCEDURE Runge_Kutta_1D_LS4_3
  MODULE PROCEDURE Runge_Kutta_2D_LS4_3
  MODULE PROCEDURE Runge_Kutta_3D_LS4_3
END INTERFACE

PRIVATE

PUBLIC :: PS_Runge_Kutta_LS4_3

CONTAINS
!
!======================================================================
!
  Subroutine Runge_Kutta_0D_LS4_3G (Stage, N0, N5, N2, N3,         &
                                           NV, dt, Q, Q1, D_Flux,  &
                                           Time, Time_Now)

  integer  :: Stage, N0, N5, N2, N3,  NV
  REALTYPE :: dt
  REALTYPE, dimension(N0:N5) :: Q, Q1, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: i, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do i = N2,N3
        Q1(i) =                dt*D_Flux(i)
        Q (i) =      Q (i) +   B(m)*Q1(i)
      enddo
!$OMP END PARALLEL DO

    CASE (2,3,4)
!$OMP PARALLEL DO
      do i = N2,N3
        Q1(i) = A(m)*Q1(i) + dt*D_Flux(i)
        Q (i) =      Q (i) +   B(m)*Q1(i)
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(Time_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_0D_LS4_3G
!
!======================================================================
!
  Subroutine Runge_Kutta_1D_LS4_3G (Stage, N0, N5, N2, N3,         &
                                           M0, M5, M2, M3,         &
                                           NV, dt, Q, Q1, D_Flux,  &
                                           Time, Time_Now)
      
  integer  :: Stage, N0, N5, N2, N3
  integer  ::        M0, M5, M2, M3,  NV
  REALTYPE :: dt
  REALTYPE, dimension(N0:N5,NV) :: Q, Q1, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: i, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do i = N2,N3
        Q1(i,:) =                dt*D_Flux(i,:)
        Q (i,:) =      Q (i,:) +   B(m)*Q1(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2,3,4)
!$OMP PARALLEL DO
      do i = N2,N3
        Q1(i,:) = A(m)*Q1(i,:) + dt*D_Flux(i,:)
        Q (i,:) =      Q (i,:) +   B(m)*Q1(i,:)
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(Time_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_1D_LS4_3G
!
!======================================================================
!
  Subroutine Runge_Kutta_2D_LS4_3G (Stage, N0, N5, N2, N3,         &
                                           M0, M5, M2, M3,         &
                                           NV, dt, Q, Q1, D_Flux,  &
                                           Time, Time_Now)
      
  integer  :: Stage, N0, N5, N2, N3
  integer  ::        M0, M5, M2, M3,  NV
  REALTYPE :: dt
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, Q1, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: j, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do j = M2,M3
        Q1(N2:N3,j,:) =                      dt*D_Flux(N2:N3,j,:)
        Q (N2:N3,j,:) =      Q (N2:N3,j,:) +   B(m)*Q1(N2:N3,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2,3,4)
!$OMP PARALLEL DO
      do j = M2,M3
        Q1(N2:N3,j,:) = A(m)*Q1(N2:N3,j,:) + dt*D_Flux(N2:N3,j,:)
        Q (N2:N3,j,:) =      Q (N2:N3,j,:) +   B(m)*Q1(N2:N3,j,:)
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(Time_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_2D_LS4_3G
!
!======================================================================
!
  Subroutine Runge_Kutta_3D_LS4_3G (Stage, N0, N5, N2, N3,         &
                                           M0, M5, M2, M3,         &
                                           K0, K5, K2, K3,         &
                                           NV, dt, Q, Q1, D_Flux,  &
                                           Time, Time_Now)
      
  integer  :: Stage, N0, N5, N2, N3
  integer  ::        M0, M5, M2, M3
  integer  ::        K0, K5, K2, K3,  NV
  REALTYPE :: dt
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, Q1, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: k, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do k = K2,K3
        Q1(N2:N3,M2:M3,k,:) =                    dt*D_Flux(N2:N3,M2:M3,k,:)
        Q (N2:N3,M2:M3,k,:) = Q(N2:N3,M2:M3,k,:) + B(m)*Q1(N2:N3,M2:M3,k,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2,3,4)
!$OMP PARALLEL DO
      do k = K2,K3
        Q1(N2:N3,M2:M3,k,:) = A(m)*Q1(N2:N3,M2:M3,k,:) + &
                                                 dt*D_Flux(N2:N3,M2:M3,k,:)
        Q (N2:N3,M2:M3,k,:) =      Q (N2:N3,M2:M3,k,:) + &
                                                   B(m)*Q1(N2:N3,M2:M3,k,:)
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(Time_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_3D_LS4_3G
!
!======================================================================
!
  Subroutine Runge_Kutta_0D_LS4_3 (Stage, dt, Q, Q1, D_Flux, Time, Time_Now)

  integer  :: Stage
  REALTYPE :: dt
  REALTYPE, dimension(:) :: Q, Q1, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: i, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q1(i) =                dt*D_Flux(i)
        Q (i) =      Q (i) +   B(m)*Q1(i)
      enddo
!$OMP END PARALLEL DO

    CASE (2,3,4)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q1(i) = A(m)*Q1(i) + dt*D_Flux(i)
        Q (i) =      Q (i) +   B(m)*Q1(i)
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(Time_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_0D_LS4_3
!
!======================================================================
!
  Subroutine Runge_Kutta_1D_LS4_3 (Stage, dt, Q, Q1, D_Flux, Time, Time_Now)

  integer  :: Stage
  REALTYPE :: dt
  REALTYPE, dimension(:,:) :: Q, Q1, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: i, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q1(i,:) =                dt*D_Flux(i,:)
        Q (i,:) =      Q (i,:) +   B(m)*Q1(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2,3,4)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q1(i,:) = A(m)*Q1(i,:) + dt*D_Flux(i,:)
        Q (i,:) =      Q (i,:) +   B(m)*Q1(i,:)
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(Time_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_1D_LS4_3
!
!======================================================================
!
  Subroutine Runge_Kutta_2D_LS4_3 (Stage, dt, Q, Q1, D_Flux, Time, Time_Now)

  integer  :: Stage
  REALTYPE :: dt
  REALTYPE, dimension(:,:,:) :: Q, Q1, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: j, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do j = LBOUND(Q,DIM=2),UBOUND(Q,DIM=2)
        Q1(:,j,:) =                  dt*D_Flux(:,j,:)
        Q (:,j,:) =      Q (:,j,:) +   B(m)*Q1(:,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2,3,4)
!$OMP PARALLEL DO
      do j = LBOUND(Q,DIM=2),UBOUND(Q,DIM=2)
        Q1(:,j,:) = A(m)*Q1(:,j,:) + dt*D_Flux(:,j,:)
        Q (:,j,:) =      Q (:,j,:) +   B(m)*Q1(:,j,:)
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(Time_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_2D_LS4_3
!
!======================================================================
!
  Subroutine Runge_Kutta_3D_LS4_3 (Stage, dt, Q, Q1, D_Flux, Time, Time_Now)

  integer  :: Stage
  REALTYPE :: dt
  REALTYPE, dimension(:,:,:,:) :: Q, Q1, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: k, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do k = LBOUND(Q,DIM=3),UBOUND(Q,DIM=3)
        Q1(:,:,k,:) =                    dt*D_Flux(:,:,k,:)
        Q (:,:,k,:) =      Q (:,:,k,:) +   B(m)*Q1(:,:,k,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2,3,4)
!$OMP PARALLEL DO
      do k = LBOUND(Q,DIM=3),UBOUND(Q,DIM=3)
        Q1(:,:,k,:) = A(m)*Q1(:,:,k,:) + dt*D_Flux(:,:,k,:)
        Q (:,:,k,:) =      Q (:,:,k,:) +   B(m)*Q1(:,:,k,:)
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(Time_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_3D_LS4_3

END MODULE Runge_Kutta_LS4_3
