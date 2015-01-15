  if ((Adjust_D_Flux >= 0) .AND. (Adjust_D_Flux <=3)) then

    call PS_Adjust_D_F_Epsilon (Adjust_D_Flux, Navier_Stokes, D_F_Zero, &
                                D_F_Epsilon, dt, dx, dy, dz)

    if ((I_Am == First_Processor) .AND. (Step == 1)) &
    write (6,*) 'Roundoff error control : ', D_F_Epsilon

  endif
