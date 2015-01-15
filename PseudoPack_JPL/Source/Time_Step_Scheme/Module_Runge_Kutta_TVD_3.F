#include "pseudopack.h"

!*******************************************************************
! Low Storage 3 stages 3th order Runge Kutta TVD Scheme by 
!   Osher and Shu
!
!  Q1 =    Q +          dt*DFLUX(Q)
!  Q1 = (3*Q +   Q1 +   dt*DFLUX(Q1))/4
!  Q  = (  Q + 2*Q1 + 2*dt*DFLUX(Q1))/3
!
!   Q1 is a temporary storage  
!   Q  is the intermediate solution 
!*******************************************************************

MODULE Runge_Kutta_TVD_3

implicit NONE

INTERFACE PS_Runge_Kutta_TVD_3
  MODULE PROCEDURE Runge_Kutta_0D_TVD_3G
  MODULE PROCEDURE Runge_Kutta_1D_TVD_3G
  MODULE PROCEDURE Runge_Kutta_2D_TVD_3G
  MODULE PROCEDURE Runge_Kutta_3D_TVD_3G

  MODULE PROCEDURE Runge_Kutta_0D_TVD_3
  MODULE PROCEDURE Runge_Kutta_1D_TVD_3
  MODULE PROCEDURE Runge_Kutta_2D_TVD_3
  MODULE PROCEDURE Runge_Kutta_3D_TVD_3
END INTERFACE

REALTYPE, dimension(0:3) :: C = (/ ZERO, ONE, HALF, ONE /)

PRIVATE

PUBLIC :: PS_Runge_Kutta_TVD_3

CONTAINS
!
!======================================================================
!
  Subroutine Runge_Kutta_0D_TVD_3G (Stage, N0, N5, N2, N3,        &
                                           NV, dt, Q, Q1, D_Flux, &
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
        Q1(i) = Q(i) + dt*D_Flux(i)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do i = N2,N3
        Q1(i) = (THREE*Q(i) + Q1(i) + dt*D_Flux(i))/FOUR
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = N2,N3
        Q (i) = (Q(i) + TWO*Q1(i) + TWO*dt*D_Flux(i))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_0D_TVD_3G
!
!======================================================================
!
  Subroutine Runge_Kutta_1D_TVD_3G (Stage, N0, N5, N2, N3,        &
                                           NV, dt, Q, Q1, D_Flux, &
                                           Time, Time_Now)

  integer  :: Stage, N0, N5, N2, N3,  NV
  REALTYPE :: dt
  REALTYPE, dimension(N0:N5,NV) :: Q, Q1, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: i, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do i = N2,N3
        Q1(i,:) = Q(i,:) + dt*D_Flux(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do i = N2,N3
        Q1(i,:) = (THREE*Q(i,:) + Q1(i,:) + dt*D_Flux(i,:))/FOUR
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = N2,N3
        Q (i,:) = (Q(i,:) + TWO*Q1(i,:) + TWO*dt*D_Flux(i,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_1D_TVD_3G
!
!======================================================================
!
  Subroutine Runge_Kutta_2D_TVD_3G (Stage, N0, N5, N2, N3,        &
                                           M0, M5, M2, M3,        &
                                           NV, dt, Q, Q1, D_Flux, &
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
        Q1(N2:N3,j,:) = Q(N2:N3,j,:) + dt*D_Flux(N2:N3,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do j = M2,M3
        Q1(N2:N3,j,:) = (THREE*Q(N2:N3,j,:) + Q1(N2:N3,j,:) + &
                                      dt*D_Flux(N2:N3,j,:))/FOUR
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do j = M2,M3
        Q (N2:N3,j,:) = (Q(N2:N3,j,:) + TWO*Q1(N2:N3,j,:) + &
                                  TWO*dt*D_Flux(N2:N3,j,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_2D_TVD_3G
!
!======================================================================
!
  Subroutine Runge_Kutta_3D_TVD_3G (Stage, N0, N5, N2, N3,        &
                                           M0, M5, M2, M3,        &
                                           K0, K5, K2, K3,        &
                                           NV, dt, Q, Q1, D_Flux, &
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
        Q1(N2:N3,M2:M3,k,:) = Q(N2:N3,M2:M3,k,:) + dt*D_Flux(N2:N3,M2:M3,k,:)
      enddo
!$OMP END PARALLEL DO
      
    CASE (2)
!$OMP PARALLEL DO
      do k = K2,K3
        Q1(N2:N3,M2:M3,k,:) = (THREE*Q(N2:N3,M2:M3,k,:) + Q1(N2:N3,M2:M3,k,:) +&
                                      dt*D_Flux(N2:N3,M2:M3,k,:))/FOUR
      enddo
!$OMP END PARALLEL DO
      
    CASE (3)
!$OMP PARALLEL DO
      do k = K2,K3
        Q (N2:N3,M2:M3,k,:) = (Q(N2:N3,M2:M3,k,:) + TWO*Q1(N2:N3,M2:M3,k,:) + &
                                  TWO*dt*D_Flux(N2:N3,M2:M3,k,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_3D_TVD_3G
!
!======================================================================
!
  Subroutine Runge_Kutta_0D_TVD_3 (Stage, dt, Q, Q1, D_Flux, Time, Time_Now)

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
        Q1(i) = Q(i) + dt*D_Flux(i)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q1(i) = (THREE*Q(i) + Q1(i) + dt*D_Flux(i))/FOUR
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q (i) = (Q(i) + TWO*Q1(i) + TWO*dt*D_Flux(i))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_0D_TVD_3
!
!======================================================================
!
  Subroutine Runge_Kutta_1D_TVD_3 (Stage, dt, Q, Q1, D_Flux, Time, Time_Now)

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
        Q1(i,:) = Q(i,:) + dt*D_Flux(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q1(i,:) = (THREE*Q(i,:) + Q1(i,:) + dt*D_Flux(i,:))/FOUR
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q (i,:) = (Q(i,:) + TWO*Q1(i,:) + TWO*dt*D_Flux(i,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_1D_TVD_3
!
!======================================================================
!
  Subroutine Runge_Kutta_2D_TVD_3 (Stage, dt, Q, Q1, D_Flux, Time, Time_Now)

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
        Q1(:,j,:) = Q(:,j,:) + dt*D_Flux(:,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do j = LBOUND(Q,DIM=2),UBOUND(Q,DIM=2)
        Q1(:,j,:) = (THREE*Q(:,j,:) + Q1(:,j,:) + dt*D_Flux(:,j,:))/FOUR
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do j = LBOUND(Q,DIM=2),UBOUND(Q,DIM=2)
        Q (:,j,:) = (Q(:,j,:) + TWO*Q1(:,j,:) + TWO*dt*D_Flux(:,j,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_2D_TVD_3
!
!======================================================================
!
  Subroutine Runge_Kutta_3D_TVD_3 (Stage, dt, Q, Q1, D_Flux, Time, Time_Now)

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
        Q1(:,:,k,:) = Q(:,:,k,:) + dt*D_Flux(:,:,k,:)
      enddo
!$OMP END PARALLEL DO
      
    CASE (2)
!$OMP PARALLEL DO
      do k = LBOUND(Q,DIM=3),UBOUND(Q,DIM=3)
        Q1(:,:,k,:) = (THREE*Q(:,:,k,:) + Q1(:,:,k,:) + &
                                      dt*D_Flux(:,:,k,:))/FOUR
      enddo
!$OMP END PARALLEL DO
      
    CASE (3)
!$OMP PARALLEL DO
      do k = LBOUND(Q,DIM=3),UBOUND(Q,DIM=3)
        Q (:,:,k,:) = (Q(:,:,k,:) + TWO*Q1(:,:,k,:) + &
                                  TWO*dt*D_Flux(:,:,k,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_3D_TVD_3

END MODULE Runge_Kutta_TVD_3
