#if defined (VARIABLE_GAMMA)
#define _EULER_GAMMA	Specific_Heat
#else
#define _EULER_GAMMA	Gamma
#endif

  Subroutine Time_Step (dt, Euler_dt, NS_dt)

  REALTYPE :: dt
  REALTYPE :: dt_Euler, dt_NS
  INTEGER :: i,j,k,l

  REALTYPE, OPTIONAL :: Euler_dt, NS_dt

  dt_Euler = ZERO 

#if defined (VARIABLE_GAMMA)
  call PS_Specific_Heat_Gamma (Gamma_Formulation, Q, Gamma_Gas, Specific_Heat)
#endif

    call PS_Euler_Time_Step (N0, N2, N3, N5,                                &
                             M0, M2, M3, M5,                                &
                             K0, K2, K3, K5, NV,                            &
                             x, y, z, Q, _EULER_GAMMA, CFL_Euler, dt_Euler  )

  dt = dt_Euler

  if (PRESENT(Euler_dt)) Euler_dt = dt_Euler

#if defined (NAVIER_STOKES) 
  if (Navier_Stokes) then

    call PS_Navier_Stokes_Time_Step (N0, N2, N3, N5,                        &
                                     M0, M2, M3, M5,                        &
                                     K0, K2, K3, K5, NV,                    &
                                     x, y, z, Q, _EULER_GAMMA, CFL_NS, dt_NS)

    dt = MIN(dt_Euler, dt_NS)

    if (dt_NS == HUGE(ONE)) dt_NS = ZERO
  endif

  if (PRESENT(NS_dt)) NS_dt = dt_NS
#else
  if (PRESENT(NS_dt)) NS_dt = ZERO
#endif

!---------------------------------------------
          ! KAUSHIK ADDED THIS

!  if (I_Am == First_Processor) then
!       open(345,file='time_step',form='formatted',access='append')
!          write(345,*) dt 
!       close(345)
!  endif
!---------------------------------------------

  END Subroutine Time_Step
!
!======================================================================
!
  Subroutine Stability_Check

  if (dt >  EPSILON(ONE)) RETURN

  if (I_Am == First_Processor) then
    do lid_loop = 1,lid_loop_Max
      if (lid_loop == 1) lid = lid60
      if (lid_loop == 2) lid = lid6 

      write (lid,*) '================  STOP ! Unstable  ====================='
      write (lid,*)
    enddo
  endif

    TecPlot_Header = .TRUE. ; call TecPlot (Plot_BlowUp_File)

  call PS_STOP ('Stability_Check',' ')

  END Subroutine Stability_Check

#if defined (_EULER_GAMMA)
#undef _EULER_GAMMA
#endif
