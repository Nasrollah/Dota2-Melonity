#include "pseudopack.h"

! ----------------------------------------------------------------------
! FUNCTION NAME: Euler_Time_Step
! AUTHOR       : Wai-Sun Don
! DESCRIPTION  : Compute the time step based on the inviscid part of the
!                Navier-Stokes Equations with Perfect gas law P = R*rho*T
!               
! Gamma      The specific heat ratio (Input)
! x, y, z    Arrays contain the collocation grid (Input)
!              of size (0:N)         for 1D problem (y and z are ignored)
!              of size (0:N,0:M)     for 2D problem (z is ignored)
!              of size (0:N,0:M,0:L) for 3D problem
! Q          Array contains the Conservative variables (Input)
!              of size (0:N,       ,3) for a 2D problem
!              of size (0:N,0:M,   ,4) for a 2D problem
!              of size (0:N,0:M,0:L,5) for a 3D problem
! CFL        CFL number of the Euler equations (Input)
! dt         Stable Euler time step  (Output)
!
!*********************************************************************


MODULE Euler_CFL_Conditions

#if defined (PARALLEL_MPI)
  USE Processor
#endif

implicit NONE

INTERFACE PS_Euler_Time_Step
  MODULE PROCEDURE Euler_Time_Step_1D
  MODULE PROCEDURE Euler_Time_Step_2D
  MODULE PROCEDURE Euler_Time_Step_3D

  MODULE PROCEDURE Euler_Time_Step_1D_G
  MODULE PROCEDURE Euler_Time_Step_2D_G
  MODULE PROCEDURE Euler_Time_Step_3D_G
END INTERFACE

PRIVATE
PUBLIC  :: PS_Euler_Time_Step

CONTAINS

!------
#if 1
!------
!
!==========================================================================
!
  Subroutine Euler_Time_Step_1D   (N0, N2, N3, N5, NV,                   &
                                   x,       Q, Gamma,                    &
                                   CFL, dt,                              &
                                   V_Moving_Frame)

  integer  :: N0, N2, N3, N5

  integer  :: NV

  REALTYPE                      :: CFL, dt
  REALTYPE, dimension(N0:N5)    :: x
  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE                      :: Gamma
  REALTYPE, dimension(N0:N5)    :: U, c2

  REALTYPE, dimension(1), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(1)           :: S

  integer  :: i, N6
  REALTYPE :: dx, dL, dU, c, dt_Euler

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  dt = -ONE ; dt_Euler = -ONE

      U = Q(:,2)/Q(:,1)
     c2 = Gamma*(Gamma-ONE)*(Q(:,3)/Q(:,1) - HALF*(U*U))

    N6 = N3-1 ; if (N2 == N3) N6 = N2 
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
      do i = N2,N6
        dx = ABS(x(i+1)-x(i))
 
        dL = SQRT(ONE/dx**2)

        dU = ABS(U(i)-S(1))/dx 

        if (c2(i) >  ZERO) then
          c = SQRT(c2(i))

          dt_Euler = MAX(dt_Euler, dU+c*dL, dU-c*dL) 
        endif
      enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_1D
!
!==========================================================================
!
  Subroutine Euler_Time_Step_2D   (N0, N2, N3, N5,                       &
                                   M0, M2, M3, M5, NV,                   &
                                   x, y,    Q, Gamma,                    &
                                   CFL, dt,                              &
                                   V_Moving_Frame)

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5

  integer  :: NV

  REALTYPE                            :: CFL, dt
  REALTYPE, dimension(N0:N5,M0:M5)    :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE                            :: Gamma
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, c2

  REALTYPE, dimension(2), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(2)           :: S

  integer  :: i, j, N6, M6
  REALTYPE :: dx, dy, dL, dU, c, dt_Euler

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  dt = -ONE ; dt_Euler = -ONE

      U = Q(:,:,2)/Q(:,:,1)
      V = Q(:,:,3)/Q(:,:,1)
     c2 = Gamma*(Gamma-ONE)*(Q(:,:,4)/Q(:,:,1) - HALF*(U*U+V*V))

    N6 = N3-1 ; if (N2 == N3) N6 = N2 
    M6 = M3-1 ; if (M2 == M3) M6 = M2 
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
    do j = M2,M6
      do i = N2,N6
        dx = ABS(x(i+1,j  )-x(i,j))
        dy = ABS(y(i  ,j+1)-y(i,j))
 
        dL = SQRT(ONE/dx**2 + ONE/dy**2)

        dU = ABS(U(i,j)-S(1))/dx + ABS(V(i,j)-S(2))/dy 

        if (c2(i,j) >  ZERO) then
          c = SQRT(c2(i,j))

          dt_Euler = MAX(dt_Euler, dU+c*dL, dU-c*dL) 
        endif
      enddo
    enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_2D
!
!==========================================================================
!
  Subroutine Euler_Time_Step_3D   (N0, N2, N3, N5,                       &
                                   M0, M2, M3, M5,                       &
                                   K0, K2, K3, K5, NV,                   &
                                   x, y, z, Q, Gamma,                    &
                                   CFL, dt,                              &
                                   V_Moving_Frame)

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5
  integer  :: K0, K2, K3, K5

  integer  :: NV

  REALTYPE                                  :: CFL, dt
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z 
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE                                  :: Gamma
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, c2

  REALTYPE, dimension(3), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(3)           :: S

  integer  :: N6, M6, K6
  integer  :: i, j, k
  REALTYPE :: dx, dy, dz, dL, dU, c, dt_Euler

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  dt = -ONE ; dt_Euler = -ONE

      U = Q(:,:,:,2)/Q(:,:,:,1)
      V = Q(:,:,:,3)/Q(:,:,:,1)
      W = Q(:,:,:,4)/Q(:,:,:,1)
     c2 = Gamma*(Gamma-ONE)*(Q(:,:,:,5)/Q(:,:,:,1) - HALF*(U*U+V*V+W*W))

    N6 = N3-1 ; if (N2 == N3) N6 = N2 
    M6 = M3-1 ; if (M2 == M3) M6 = M2 
    K6 = K3-1 ; if (K2 == K3) K6 = K2 
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
  do k = K2,K6
    do j = M2,M6
      do i = N2,N6
        dx = ABS(x(i+1,j  ,k  )-x(i,j,k))
        dy = ABS(y(i  ,j+1,k  )-y(i,j,k))
        dz = ABS(z(i  ,j  ,k+1)-z(i,j,k))
 
        dL = SQRT(ONE/dx**2 + ONE/dy**2 + ONE/dz**2)

        dU = ABS(U(i,j,k)-S(1))/dx + ABS(V(i,j,k)-S(2))/dy + &
                                     ABS(W(i,j,k)-S(3))/dz

        if (c2(i,j,k) >  ZERO) then
          c = SQRT(c2(i,j,k))

          dt_Euler = MAX(dt_Euler, dU+c*dL, dU-c*dL) 
        endif
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_3D
!
!==========================================================================
!
  Subroutine Euler_Time_Step_1D_G (N0, N2, N3, N5, NV,                   &
                                   x,       Q, Gamma,                    &
                                   CFL, dt,                              &
                                   V_Moving_Frame)

  integer  :: N0, N2, N3, N5

  integer  :: NV

  REALTYPE                      :: CFL, dt
  REALTYPE, dimension(N0:N5)    :: x
  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: Gamma
  REALTYPE, dimension(N0:N5)    :: U, c2

  REALTYPE, dimension(1), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(1)           :: S

  integer  :: N6
  integer  :: i
  REALTYPE :: dx, dL, dU, c, dt_Euler

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  dt = -ONE ; dt_Euler = -ONE

      U = Q(:,2)/Q(:,1)
     c2 = Gamma*(Gamma-ONE)*(Q(:,3)/Q(:,1) - HALF*(U*U))

    N6 = N3-1 ; if (N2 == N3) N6 = N2 
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
      do i = N2,N6
        dx = ABS(x(i+1)-x(i))
 
        dL = SQRT(ONE/dx**2)

        dU = ABS(U(i)-S(1))/dx 

        if (c2(i) >  ZERO) then
          c = SQRT(c2(i))

          dt_Euler = MAX(dt_Euler, dU+c*dL, dU-c*dL) 
        endif
      enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_1D_G
!
!==========================================================================
!
  Subroutine Euler_Time_Step_2D_G (N0, N2, N3, N5,                       &
                                   M0, M2, M3, M5, NV,                   &
                                   x, y,    Q, Gamma,                    &
                                   CFL, dt,                              &
                                   V_Moving_Frame)

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5

  integer  :: NV

  REALTYPE                            :: CFL, dt
  REALTYPE, dimension(N0:N5,M0:M5)    :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: Gamma
  REALTYPE, dimension(N0:N5,M0:M5)    :: U, V, c2

  REALTYPE, dimension(2), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(2)           :: S

  integer  :: N6, M6
  integer  :: i, j
  REALTYPE :: dx, dy, dL, dU, c, dt_Euler

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  dt = -ONE ; dt_Euler = -ONE

      U = Q(:,:,2)/Q(:,:,1)
      V = Q(:,:,3)/Q(:,:,1)
     c2 = Gamma*(Gamma-ONE)*(Q(:,:,4)/Q(:,:,1) - HALF*(U*U+V*V))

    N6 = N3-1 ; if (N2 == N3) N6 = N2 
    M6 = M3-1 ; if (M2 == M3) M6 = M2 
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
    do j = M2,M6
      do i = N2,N6
        dx = ABS(x(i+1,j  )-x(i,j))
        dy = ABS(y(i  ,j+1)-y(i,j))
 
        dL = SQRT(ONE/dx**2 + ONE/dy**2)

        dU = ABS(U(i,j)-S(1))/dx + ABS(V(i,j)-S(2))/dy 

        if (c2(i,j) >  ZERO) then
          c = SQRT(c2(i,j))

          dt_Euler = MAX(dt_Euler, dU+c*dL, dU-c*dL) 
        endif
      enddo
    enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_2D_G
!
!==========================================================================
!
  Subroutine Euler_Time_Step_3D_G (N0, N2, N3, N5,                       &
                                   M0, M2, M3, M5,                       &
                                   K0, K2, K3, K5, NV,                   &
                                   x, y, z, Q, Gamma,                    &
                                   CFL, dt,                              &
                                   V_Moving_Frame)

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5
  integer  :: K0, K2, K3, K5

  integer  :: NV

  REALTYPE                                  :: CFL, dt
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z 
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: Gamma
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, c2

  REALTYPE, dimension(3), OPTIONAL :: V_Moving_Frame
  REALTYPE, dimension(3)           :: S

  integer  :: N6, M6, K6
  integer  :: i, j, k
  REALTYPE :: dx, dy, dz, dL, dU, c, dt_Euler

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  dt = -ONE ; dt_Euler = -ONE

      U = Q(:,:,:,2)/Q(:,:,:,1)
      V = Q(:,:,:,3)/Q(:,:,:,1)
      W = Q(:,:,:,4)/Q(:,:,:,1)
     c2 = Gamma*(Gamma-ONE)*(Q(:,:,:,5)/Q(:,:,:,1) - HALF*(U*U+V*V+W*W))

    N6 = N3-1 ; if (N2 == N3) N6 = N2 
    M6 = M3-1 ; if (M2 == M3) M6 = M2 
    K6 = K3-1 ; if (K2 == K3) K6 = K2 
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
  do k = K2,K6
    do j = M2,M6
      do i = N2,N6
        dx = ABS(x(i+1,j  ,k  )-x(i,j,k))
        dy = ABS(y(i  ,j+1,k  )-y(i,j,k))
        dz = ABS(z(i  ,j  ,k+1)-z(i,j,k))
 
        dL = SQRT(ONE/dx**2 + ONE/dy**2 + ONE/dz**2)

        dU = ABS(U(i,j,k)-S(1))/dx + ABS(V(i,j,k)-S(2))/dy + &
                                     ABS(W(i,j,k)-S(3))/dz

        if (c2(i,j,k) >  ZERO) then
          c = SQRT(c2(i,j,k))

          dt_Euler = MAX(dt_Euler, dU+c*dL, dU-c*dL) 
        endif
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_3D_G

!------
#else
!------

  Subroutine Euler_Time_Step_1D (x, Q, Gamma, CFL, dt, V_Moving_Frame, n_0)

  REALTYPE                        :: Gamma, CFL, dt
  REALTYPE, dimension(0:)         :: x
  REALTYPE, dimension(0:,:)       :: Q

  integer ,               OPTIONAL :: n_0
  REALTYPE, dimension(1), OPTIONAL :: V_Moving_Frame

  integer  :: N6
  integer  :: N, i, ii, i_Offset
  REALTYPE :: gm1, dx, dL, u, T, c2, dt_Euler
  REALTYPE :: ABS_U

  REALTYPE, dimension(1) :: S

  i_Offset = 0 ; if (PRESENT(n_0)) i_Offset = n_0

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  N = UBOUND(Q,DIM=1) 

  gm1 = Gamma - ONE

  dt = -ONE ; dt_Euler = -ONE

    N6 = N-1 ; if (N == 0) N6 = 0
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
  do i = 0,N6
    ii = i_Offset+i 

    dx = ABS(x(ii+1)-x(ii))
 
    dL = ONE/dx

    u  = Q(i,2)/Q(i,1)
    T  = Q(i,3)/Q(i,1) - HALF*(u**2)
    c2 = Gamma*gm1*T

    ABS_U = ABS(u-S(1))

    if (c2 >  ZERO) then
      dt_Euler = MAX(dt_Euler,               &
                     ABS_U/dx + SQRT(c2)*dL, &
                     ABS_U/dx - SQRT(c2)*dL)
    endif
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_1D
!
!==========================================================================
!
  Subroutine Euler_Time_Step_2D (x, y, Q, Gamma, CFL, dt, &
                                       V_Moving_Frame, n_0, m_0)

  REALTYPE                     :: Gamma, CFL, dt
  REALTYPE, dimension(0:,0:)   :: x, y
  REALTYPE, dimension(0:,0:,:) :: Q

  integer ,               OPTIONAL :: n_0, m_0
  REALTYPE, dimension(2), OPTIONAL :: V_Moving_Frame

  integer  :: N6, M6
  integer  :: N, M, i, j, ii, jj, i_Offset, j_Offset
  REALTYPE :: gm1, dx, dy, dL, u, v, T, c2, dt_Euler
  REALTYPE :: ABS_U, ABS_V

  REALTYPE, dimension(2) :: S

  i_Offset = 0 ; if (PRESENT(n_0)) i_Offset = n_0
  j_Offset = 0 ; if (PRESENT(m_0)) j_Offset = m_0

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  N = UBOUND(Q,DIM=1) ; M = UBOUND(Q,DIM=2)

  gm1 = Gamma - ONE

  dt = -ONE ; dt_Euler = -ONE

    N6 = N-1 ; if (N == 0) N6 = 0
    M6 = M-1 ; if (M == 0) M6 = 0
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
  do j = 0,M6
    do i = 0,N6
      ii = i_Offset+i ; jj = j_Offset+j

      dx = ABS(x(ii+1,jj  )-x(ii,jj))
      dy = ABS(y(ii  ,jj+1)-y(ii,jj))
      dL = SQRT(ONE/dx**2 + ONE/dy**2)

      u  = Q(i,j,2)/Q(i,j,1)
      v  = Q(i,j,3)/Q(i,j,1)
      T  = Q(i,j,4)/Q(i,j,1) - HALF*(u**2 + v**2)
      c2 = Gamma*gm1*T

      ABS_U = ABS(u-S(1)) ; ABS_V = ABS(v-S(2))

      if (c2 >  ZERO) then
        dt_Euler = MAX(dt_Euler,                          &
                       ABS_U/dx + ABS_V/dy + SQRT(c2)*dL, &
                       ABS_U/dx + ABS_V/dy - SQRT(c2)*dL)
      endif
    enddo
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_2D
!
!==========================================================================
!
  Subroutine Euler_Time_Step_3D (x, y, z, Q, Gamma, CFL, dt, &
                                          V_Moving_Frame, n_0, m_0, k_0)

  REALTYPE                        :: Gamma, CFL, dt
  REALTYPE, dimension(0:,0:,0:)   :: x, y, z
  REALTYPE, dimension(0:,0:,0:,:) :: Q

  integer ,               OPTIONAL :: n_0, m_0, k_0
  REALTYPE, dimension(3), OPTIONAL :: V_Moving_Frame

  integer  :: N6, M6, K6
  integer  :: N, M, L, i, j, k, ii, jj, kk, i_Offset, j_Offset, k_Offset
  REALTYPE :: gm1, dx, dy, dz, dL, u, v, w, T, c2, dt_Euler
  REALTYPE :: ABS_U, ABS_V, ABS_W

  REALTYPE, dimension(3) :: S

  i_Offset = 0 ; if (PRESENT(n_0)) i_Offset = n_0
  j_Offset = 0 ; if (PRESENT(m_0)) j_Offset = m_0
  k_Offset = 0 ; if (PRESENT(k_0)) k_Offset = k_0

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  N = UBOUND(Q,DIM=1) ; M = UBOUND(Q,DIM=2) ; L = UBOUND(Q,DIM=3)

  dt = -ONE ; dt_Euler = -ONE

    N6 = N-1 ; if (N == 0) N6 = 0
    M6 = M-1 ; if (M == 0) M6 = 0
    K6 = L-1 ; if (L == 0) K6 = 0
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
  do k = 0,K6
    do j = 0,M6
      do i = 0,N6
        ii = i_Offset+i ; jj = j_Offset+j ; kk = k_Offset+k

        dx = ABS(x(ii+1,jj  ,kk  )-x(ii,jj,kk))
        dy = ABS(y(ii  ,jj+1,kk  )-y(ii,jj,kk))
        dz = ABS(z(ii  ,jj  ,kk+1)-z(ii,jj,kk))
 
        dL = SQRT(ONE/dx**2 + ONE/dy**2 + ONE/dz**2)

        u  = Q(i,j,k,2)/Q(i,j,k,1)
        v  = Q(i,j,k,3)/Q(i,j,k,1)
        w  = Q(i,j,k,4)/Q(i,j,k,1)
        T  = Q(i,j,k,5)/Q(i,j,k,1) - HALF*(u**2 + v**2 + w**2)
        c2 = Gamma*(Gamma-ONE)*T

        ABS_U = ABS(u-S(1)) ; ABS_V = ABS(v-S(2)) ; ABS_W = ABS(w-S(3))

        if (c2 >  ZERO) then
          dt_Euler = MAX(dt_Euler,                                     &
                         ABS_U/dx + ABS_V/dy + ABS_W/dz + SQRT(c2)*dL, &
                         ABS_U/dx + ABS_V/dy + ABS_W/dz - SQRT(c2)*dL)
        endif
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_3D
!
!==========================================================================
!
  Subroutine Euler_Time_Step_1D_G (x, Q, Specific_Heat_Ratio, CFL, dt, &
                                      V_Moving_Frame, n_0)

  REALTYPE                        :: CFL, dt
  REALTYPE, dimension(0:)         :: x
  REALTYPE, dimension(0:,:)       :: Q
  REALTYPE, dimension(0:)         :: Specific_Heat_Ratio 

  integer ,               OPTIONAL :: n_0
  REALTYPE, dimension(1), OPTIONAL :: V_Moving_Frame

  integer  :: N6
  integer  :: N, i, ii, i_Offset
  REALTYPE :: Gamma, gm1, dx, dL, u, T, c2, dt_Euler
  REALTYPE :: ABS_U

  REALTYPE, dimension(1) :: S

  i_Offset = 0 ; if (PRESENT(n_0)) i_Offset = n_0

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  N = UBOUND(Q,DIM=1) 

  dt = -ONE ; dt_Euler = -ONE

    N6 = N-1 ; if (N == 0) N6 = 0
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
  do i = 0,N6
    ii = i_Offset+i 

    dx = ABS(x(ii+1)-x(ii))
 
    dL = ONE/dx

    u  = Q(i,2)/Q(i,1)
    T  = Q(i,3)/Q(i,1) - HALF*(u**2)

    Gamma = Specific_Heat_Ratio(i)
      gm1 = Gamma - ONE

    c2 = Gamma*gm1*T

    ABS_U = ABS(u-S(1)) 

    if (c2 >  ZERO) then
      dt_Euler = MAX(dt_Euler,               &
                     ABS_U/dx + SQRT(c2)*dL, &
                     ABS_U/dx - SQRT(c2)*dL)
    endif
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_1D_G
!
!==========================================================================
!
  Subroutine Euler_Time_Step_2D_G (x, y, Q, Specific_Heat_Ratio, CFL, dt, &
                                         V_Moving_Frame, n_0, m_0)

  REALTYPE                     :: CFL, dt
  REALTYPE, dimension(0:,0:)   :: x, y
  REALTYPE, dimension(0:,0:,:) :: Q
  REALTYPE, dimension(0:,0:)   :: Specific_Heat_Ratio 

  integer ,               OPTIONAL :: n_0, m_0
  REALTYPE, dimension(2), OPTIONAL :: V_Moving_Frame

  integer  :: N6, M6
  integer  :: N, M, i, j, ii, jj, i_Offset, j_Offset
  REALTYPE :: Gamma, gm1, dx, dy, dL, u, v, T, c2, dt_Euler
  REALTYPE :: ABS_U, ABS_V

  REALTYPE, dimension(2) :: S

  i_Offset = 0 ; if (PRESENT(n_0)) i_Offset = n_0
  j_Offset = 0 ; if (PRESENT(m_0)) j_Offset = m_0

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  N = UBOUND(Q,DIM=1) ; M = UBOUND(Q,DIM=2)

  dt = -ONE ; dt_Euler = -ONE

    N6 = N-1 ; if (N == 0) N6 = 0
    M6 = M-1 ; if (M == 0) M6 = 0
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
  do j = 0,M6
    do i = 0,N6
      ii = i_Offset+i ; jj = j_Offset+j

      dx = ABS(x(ii+1,jj  )-x(ii,jj))
      dy = ABS(y(ii  ,jj+1)-y(ii,jj))
      dL = SQRT(ONE/dx**2 + ONE/dy**2)

      u  = Q(i,j,2)/Q(i,j,1)
      v  = Q(i,j,3)/Q(i,j,1)
      T  = Q(i,j,4)/Q(i,j,1) - HALF*(u**2 + v**2)

      Gamma = Specific_Heat_Ratio(i,j)
        gm1 = Gamma - ONE

      c2 = Gamma*gm1*T

      ABS_U = ABS(u-S(1)) ; ABS_V = ABS(v-S(2)) 

      if (c2 >  ZERO) then
        dt_Euler = MAX(dt_Euler,                          &
                       ABS_U/dx + ABS_V/dy + SQRT(c2)*dL, &
                       ABS_U/dx + ABS_V/dy - SQRT(c2)*dL)
      endif
    enddo
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_2D_G
!
!==========================================================================
!
  Subroutine Euler_Time_Step_3D_G (x, y, z, Q, Specific_Heat_Ratio, CFL, dt, &
                                            V_Moving_Frame, n_0, m_0, k_0)

  REALTYPE                        :: CFL, dt
  REALTYPE, dimension(0:,0:,0:)   :: x, y, z 
  REALTYPE, dimension(0:,0:,0:,:) :: Q
  REALTYPE, dimension(0:,0:,0:)   :: Specific_Heat_Ratio

  integer ,               OPTIONAL :: n_0, m_0, k_0
  REALTYPE, dimension(3), OPTIONAL :: V_Moving_Frame

  integer  :: N6, M6, K6
  integer  :: N, M, L, i, j, k, ii, jj, kk, i_Offset, j_Offset, k_Offset
  REALTYPE :: Gamma, gm1, dx, dy, dz, dL, u, v, w, T, c2, dt_Euler
  REALTYPE :: ABS_U, ABS_V, ABS_W

  REALTYPE, dimension(3) :: S

  i_Offset = 0 ; if (PRESENT(n_0)) i_Offset = n_0
  j_Offset = 0 ; if (PRESENT(m_0)) j_Offset = m_0
  k_Offset = 0 ; if (PRESENT(k_0)) k_Offset = k_0

  S = ZERO     ; if (PRESENT(V_Moving_Frame)) S = V_Moving_Frame

  N = UBOUND(Q,DIM=1) ; M = UBOUND(Q,DIM=2) ; L = UBOUND(Q,DIM=3)

  dt = -ONE ; dt_Euler = -ONE

    N6 = N-1 ; if (N == 0) N6 = 0
    M6 = M-1 ; if (M == 0) M6 = 0
    K6 = L-1 ; if (L == 0) K6 = 0
!$OMP PARALLEL DO REDUCTION (MAX : dt_Euler)
  do k = 0,K6
    do j = 0,M6
      do i = 0,N6
        ii = i_Offset+i ; jj = j_Offset+j ; kk = k_Offset+k

        dx = ABS(x(ii+1,jj  ,kk  )-x(ii,jj,kk))
        dy = ABS(y(ii  ,jj+1,kk  )-y(ii,jj,kk))
        dz = ABS(z(ii  ,jj  ,kk+1)-z(ii,jj,kk))
 
        dL = SQRT(ONE/dx**2 + ONE/dy**2 + ONE/dz**2)

        u  = Q(i,j,k,2)/Q(i,j,k,1)
        v  = Q(i,j,k,3)/Q(i,j,k,1)
        w  = Q(i,j,k,4)/Q(i,j,k,1)
        T  = Q(i,j,k,5)/Q(i,j,k,1) - HALF*(u**2 + v**2 + w**2)

        Gamma = Specific_Heat_Ratio(i,j,k)   
          gm1 = Gamma - ONE

        c2 = Gamma*gm1*T

        ABS_U = ABS(u-S(1)) ; ABS_V = ABS(v-S(2)) ; ABS_W = ABS(w-S(3))

        if (c2 >  ZERO) then
          dt_Euler = MAX(dt_Euler,                                     &
                         ABS_U/dx + ABS_V/dy + ABS_W/dz + SQRT(c2)*dL, &
                         ABS_U/dx + ABS_V/dy + ABS_W/dz - SQRT(c2)*dL)
        endif
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    dt = CFL/dt_Euler
  else
    call MPI_ALLREDUCE (dt_Euler, dt, 1, MPI_REALTYPE, MPI_MAX, &
                                         MPI_Comm_World, MPI_Error_Status)

    dt = CFL/dt
  endif
#else
   dt = CFL/dt_Euler
#endif

  END Subroutine Euler_Time_Step_3D_G
!------
#endif
!------

END MODULE Euler_CFL_Conditions
