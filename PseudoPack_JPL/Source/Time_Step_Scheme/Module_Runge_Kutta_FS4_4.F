#include "pseudopack.h"

!*******************************************************************
! Full Storage 4 stages 4th order Runge Kutta Scheme
!
!  Q1 =    Q +                  dt/2*DFLUX(Q )
!  Q2 =    Q +                  dt/2*DFLUX(Q1)
!  Q3 =    Q +                  dt  *DFLUX(Q2)
!  Q  = ( -Q + Q1 + 2*Q2 + Q3 + dt/2*DFLUX(Q3))/3
!
!   Q1, Q2, Q3 is a temporary storage at some intermediate stage.
!   Q  is the Initial and Final solution 
!*******************************************************************

MODULE Runge_Kutta_FS4_4

implicit NONE

INTERFACE PS_Runge_Kutta_FS4_4
  MODULE PROCEDURE Runge_Kutta_0D_FS4_4G
  MODULE PROCEDURE Runge_Kutta_1D_FS4_4G
  MODULE PROCEDURE Runge_Kutta_2D_FS4_4G
  MODULE PROCEDURE Runge_Kutta_3D_FS4_4G

  MODULE PROCEDURE Runge_Kutta_0D_FS4_4
  MODULE PROCEDURE Runge_Kutta_1D_FS4_4
  MODULE PROCEDURE Runge_Kutta_2D_FS4_4
  MODULE PROCEDURE Runge_Kutta_3D_FS4_4
END INTERFACE

REALTYPE, dimension(0:4) :: C = (/ ZERO, HALF, HALF, ONE, ONE /)

PRIVATE

PUBLIC :: PS_Runge_Kutta_FS4_4

CONTAINS
!
!======================================================================
!
  Subroutine Runge_Kutta_0D_FS4_4G (Stage, N0, N5, N2, N3,                &
                                           NV, dt, Q, Q1, Q2, Q3, D_Flux, &
                                           Time, Time_Now)

  integer  :: Stage, N0, N5, N2, N3,  NV
  REALTYPE :: dt
  REALTYPE, dimension(N0:N5) :: Q, Q1, Q2, Q3, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: i, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do i = N2,N3
        Q1(i) = Q(i) + HALF*dt*D_Flux(i)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do i = N2,N3
        Q2(i) = Q(i) + HALF*dt*D_Flux(i)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = N2,N3
        Q3(i) = Q(i) +      dt*D_Flux(i)
      enddo
!$OMP END PARALLEL DO

    CASE (4)
!$OMP PARALLEL DO
      do i = N2,N3
        Q (i) = (-Q (i) + Q1(i) + TWO*Q2(i) + Q3(i) + &
                   HALF*dt*D_Flux(i))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_0D_FS4_4G
!
!======================================================================
!
  Subroutine Runge_Kutta_1D_FS4_4G (Stage, N0, N5, N2, N3,                &
                                           NV, dt, Q, Q1, Q2, Q3, D_Flux, &
                                           Time, Time_Now)

  integer  :: Stage, N0, N5, N2, N3,  NV
  REALTYPE :: dt
  REALTYPE, dimension(N0:N5,NV) :: Q, Q1, Q2, Q3, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: i, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do i = N2,N3
        Q1(i,:) = Q(i,:) + HALF*dt*D_Flux(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do i = N2,N3
        Q2(i,:) = Q(i,:) + HALF*dt*D_Flux(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = N2,N3
        Q3(i,:) = Q(i,:) +      dt*D_Flux(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (4)
!$OMP PARALLEL DO
      do i = N2,N3
        Q (i,:) = (-Q (i,:) + Q1(i,:) + TWO*Q2(i,:) + Q3(i,:) + &
                   HALF*dt*D_Flux(i,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_1D_FS4_4G
!
!======================================================================
!
  Subroutine Runge_Kutta_2D_FS4_4G (Stage, N0, N5, N2, N3,                &
                                           M0, M5, M2, M3,                &
                                           NV, dt, Q, Q1, Q2, Q3, D_Flux, &
                                           Time, Time_Now)

  integer  :: Stage, N0, N5, N2, N3
  integer  ::        M0, M5, M2, M3,  NV
  REALTYPE :: dt
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q, Q1, Q2, Q3, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: j, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do j = M2,M3
        Q1(N2:N3,j,:) = Q(N2:N3,j,:) + HALF*dt*D_Flux(N2:N3,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do j = M2,M3
        Q2(N2:N3,j,:) = Q(N2:N3,j,:) + HALF*dt*D_Flux(N2:N3,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do j = M2,M3
        Q3(N2:N3,j,:) = Q(N2:N3,j,:) +      dt*D_Flux(N2:N3,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (4)
!$OMP PARALLEL DO
      do j = M2,M3
        Q (N2:N3,j,:) = (   -Q (N2:N3,j,:) + Q1(N2:N3,j,:) + &
                         TWO*Q2(N2:N3,j,:) + Q3(N2:N3,j,:) + &
                         HALF*dt*D_Flux(N2:N3,j,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_2D_FS4_4G
!
!======================================================================
!
  Subroutine Runge_Kutta_3D_FS4_4G (Stage, N0, N5, N2, N3,                &
                                           M0, M5, M2, M3,                &
                                           K0, K5, K2, K3,                &
                                           NV, dt, Q, Q1, Q2, Q3, D_Flux, &
                                           Time, Time_Now)
      
  integer  :: Stage, N0, N5, N2, N3
  integer  ::        M0, M5, M2, M3
  integer  ::        K0, K5, K2, K3,  NV
  REALTYPE :: dt
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q, Q1, Q2, Q3, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: k, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do k = K2,K3
        Q1(N2:N3,M2:M3,k,:) =              Q(N2:N3,M2:M3,k,:) +  &
                              HALF*dt*D_Flux(N2:N3,M2:M3,k,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do k = K2,K3
        Q2(N2:N3,M2:M3,k,:) =              Q(N2:N3,M2:M3,k,:) +  &
                              HALF*dt*D_Flux(N2:N3,M2:M3,k,:)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do k = K2,K3
        Q3(N2:N3,M2:M3,k,:) = Q(N2:N3,M2:M3,k,:) + dt*D_Flux(N2:N3,M2:M3,k,:)
      enddo
!$OMP END PARALLEL DO

    CASE (4)
!$OMP PARALLEL DO
      do k = K2,K3
        Q (N2:N3,M2:M3,k,:) = (   -Q (N2:N3,M2:M3,k,:) + Q1(N2:N3,M2:M3,k,:) + &
                               TWO*Q2(N2:N3,M2:M3,k,:) + Q3(N2:N3,M2:M3,k,:) + &
                               HALF*dt*D_Flux(N2:N3,M2:M3,k,:))/THREE
                    
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_3D_FS4_4G
!
!======================================================================
!
  Subroutine Runge_Kutta_0D_FS4_4 (Stage, dt, Q, Q1, Q2, Q3, D_Flux, &
                                          Time, Time_Now)

  integer  :: Stage
  REALTYPE :: dt
  REALTYPE, dimension(:) :: Q, Q1, Q2, Q3, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: i, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q1(i) = Q(i) + HALF*dt*D_Flux(i)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q2(i) = Q(i) + HALF*dt*D_Flux(i)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q3(i) = Q(i) +      dt*D_Flux(i)
      enddo
!$OMP END PARALLEL DO

    CASE (4)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q (i) = (-Q (i) + Q1(i) + TWO*Q2(i) + &
                  Q3(i) + HALF*dt*D_Flux(i))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_0D_FS4_4
!
!======================================================================
!
  Subroutine Runge_Kutta_1D_FS4_4 (Stage, dt, Q, Q1, Q2, Q3, D_Flux, &
                                          Time, Time_Now)

  integer  :: Stage
  REALTYPE :: dt
  REALTYPE, dimension(:,:) :: Q, Q1, Q2, Q3, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: i, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q1(i,:) = Q(i,:) + HALF*dt*D_Flux(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q2(i,:) = Q(i,:) + HALF*dt*D_Flux(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q3(i,:) = Q(i,:) +      dt*D_Flux(i,:)
      enddo
!$OMP END PARALLEL DO

    CASE (4)
!$OMP PARALLEL DO
      do i = LBOUND(Q,DIM=1),UBOUND(Q,DIM=1)
        Q (i,:) = (-Q (i,:) + Q1(i,:) + TWO*Q2(i,:) + &
                    Q3(i,:) + HALF*dt*D_Flux(i,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_1D_FS4_4
!
!======================================================================
!
  Subroutine Runge_Kutta_2D_FS4_4 (Stage, dt, Q, Q1, Q2, Q3, D_Flux, &
                                          Time, Time_Now)

  integer  :: Stage
  REALTYPE :: dt
  REALTYPE, dimension(:,:,:) :: Q, Q1, Q2, Q3, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: j, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do j = LBOUND(Q,DIM=2),UBOUND(Q,DIM=2)
        Q1(:,j,:) = Q(:,j,:) + HALF*dt*D_Flux(:,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do j = LBOUND(Q,DIM=2),UBOUND(Q,DIM=2)
        Q2(:,j,:) = Q(:,j,:) + HALF*dt*D_Flux(:,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do j = LBOUND(Q,DIM=2),UBOUND(Q,DIM=2)
        Q3(:,j,:) = Q(:,j,:) +      dt*D_Flux(:,j,:)
      enddo
!$OMP END PARALLEL DO

    CASE (4)
!$OMP PARALLEL DO
      do j = LBOUND(Q,DIM=2),UBOUND(Q,DIM=2)
        Q (:,j,:) = (-Q (:,j,:) + Q1(:,j,:) + TWO*Q2(:,j,:) + &
                      Q3(:,j,:) + HALF*dt*D_Flux(:,j,:))/THREE
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_2D_FS4_4
!
!======================================================================
!
  Subroutine Runge_Kutta_3D_FS4_4 (Stage, dt, Q, Q1, Q2, Q3, D_Flux, &
                                          Time, Time_Now)

  integer  :: Stage
  REALTYPE :: dt
  REALTYPE, dimension(:,:,:,:) :: Q, Q1, Q2, Q3, D_Flux

  REALTYPE, OPTIONAL :: Time, Time_Now

  integer  :: k, m

  m = Stage

  SELECT CASE (Stage)
    CASE (1)
!$OMP PARALLEL DO
      do k = LBOUND(Q,DIM=3),UBOUND(Q,DIM=3)
        Q1(:,:,k,:) = Q(:,:,k,:) + HALF*dt*D_Flux(:,:,k,:)
      enddo
!$OMP END PARALLEL DO

    CASE (2)
!$OMP PARALLEL DO
      do k = LBOUND(Q,DIM=3),UBOUND(Q,DIM=3)
        Q2(:,:,k,:) = Q(:,:,k,:) + HALF*dt*D_Flux(:,:,k,:)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do k = LBOUND(Q,DIM=3),UBOUND(Q,DIM=3)
        Q3(:,:,k,:) = Q(:,:,k,:) +      dt*D_Flux(:,:,k,:)
      enddo
!$OMP END PARALLEL DO

    CASE (4)
!$OMP PARALLEL DO
      do k = LBOUND(Q,DIM=3),UBOUND(Q,DIM=3)
        Q (:,:,k,:) = (-Q (:,:,k,:) + Q1(:,:,k,:) + TWO*Q2(:,:,k,:) +  &
                        Q3(:,:,k,:) + HALF*dt*D_Flux(:,:,k,:))/THREE
                    
      enddo
!$OMP END PARALLEL DO

  END SELECT

  if (PRESENT(Time) .AND. PRESENT(TIME_Now)) Time_Now = Time + C(m)*dt

  END Subroutine Runge_Kutta_3D_FS4_4

END MODULE Runge_Kutta_FS4_4
