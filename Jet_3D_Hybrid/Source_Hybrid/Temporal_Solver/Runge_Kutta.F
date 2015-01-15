#if defined (MULTI_RESOLUTION)
#define _D_FLUX_FLAG    D_Flux, MR_Flag, Flag_WENO
#else
#define _D_FLUX_FLAG    D_Flux
#endif

  Subroutine Runge_Kutta

  integer  :: Runge_Kutta_Stage, i, j, k
  REALTYPE :: Time_n, Time_Now

  integer  :: M_s, M_s_0

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q1, D_Flux
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: f

  Time_n = Time-dt

  Q1 = Q ; D_Flux = ZERO 

! --------------------------
#if defined (MULTI_RESOLUTION)
#include "Multi_Resolution.i"

  if (MR_Activation == 1) then
    if (MR_Output_Screen == 1) then
      write ( 6,100) Step, MAXVAL(ABS(MR_Scale)), MR_Tolerance
    endif
      write (79,100) Step, MAXVAL(ABS(MR_Scale)), MR_Tolerance
  endif
100 format (1x,'Step, Max MR Scale/Tolerance : ',1x,i10, 1x, 2(g20.10,:))
#endif

! --------------------------
! Stage 1 :

  Runge_Kutta_Stage = 1 ; Time_Now = Time_n + dt

  call Flux (Q , _D_FLUX_FLAG)

  call PS_Runge_Kutta_TVD_3 (1, N0, N5, N2, N3, M0, M5, M2, M3,         &
                                K0, K5, K2, K3,  NV, dt, Q, Q1, D_Flux )

#if defined (WITH_SPECIES_EQUATION)
  call PS_Adjust_Mass_Fraction (N0, N5, N2, N3, M0, M5, M2, M3,            &
                                K0, K5, K2, K3,                            &
                                NV, NF, NS, Q1,                            &
                                Adjust_Mass_Fraction, Check_Mass_Fraction, &
                                F_Zero) 
#endif

  call BC   (Q1,Time_n+dt/3.0d0,dt/3.0d0)


! --------------------------
! Stage 2 :

  Runge_Kutta_Stage = 2 ; Time_Now = Time_n + dt/2

  call Flux (Q1, _D_FLUX_FLAG)

  call PS_Runge_Kutta_TVD_3 (2, N0, N5, N2, N3, M0, M5, M2, M3,         &
                                K0, K5, K2, K3,  NV, dt, Q, Q1, D_Flux )

#if defined (WITH_SPECIES_EQUATION)
  call PS_Adjust_Mass_Fraction (N0, N5, N2, N3, M0, M5, M2, M3,            &
                                K0, K5, K2, K3,                            &
                                NV, NF, NS, Q1,                            &
                                Adjust_Mass_Fraction, Check_Mass_Fraction, &
                                F_Zero) 
#endif

  call BC   (Q1,Time_n+2.0d0*dt/3.0d0,dt/3.0d0)


! --------------------------
! Stage 3:

  Runge_Kutta_Stage = 3 ; Time_Now = Time_n + dt

  call Flux (Q1, _D_FLUX_FLAG)

  call PS_Runge_Kutta_TVD_3 (3, N0, N5, N2, N3, M0, M5, M2, M3,         &
                                K0, K5, K2, K3,  NV, dt, Q, Q1, D_Flux )

#if defined (WITH_SPECIES_EQUATION)
  call PS_Adjust_Mass_Fraction (N0, N5, N2, N3, M0, M5, M2, M3,            &
                                K0, K5, K2, K3,                            &
                                NV, NF, NS, Q ,                            &
                                Adjust_Mass_Fraction, Check_Mass_Fraction, &
                                F_Zero) 
#endif

  call BC   (Q,Time_n+dt,dt/3.0d0)

#if defined (HYBRID) 
  call Central_Smoothing (Q)
  call BC                (Q,Time_n+dt,dt/3.0d0)
#endif

  End Subroutine Runge_Kutta

#if defined (MULTI_RESOLUTION)
#undef _RK_MR_FLAG
#undef _D_FLUX_FLAG
#endif

