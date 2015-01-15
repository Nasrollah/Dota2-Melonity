if (I_Am == First_Processor) then
  do lid_loop = 1,lid_loop_Max
    if (lid_loop == 1) lid = lid60
    if (lid_loop == 2) lid = lid6

    call PS_Convert_Second_To_Clock (0, CPU_Used)

    write (lid,200) INT((Final_Time-Time)/dt)
    write (lid,201)
    write (lid,202) INT((Final_Time-Time)/dt), Step, dt, Time, ZERO, CPU_Used
  enddo
endif

200 format (1x/1x,'Estimate Number of Time Step : ', i20//)
201 format (1x,' Total_Step    Step    dt       Time        CPU Time'/  &
            1x,72('-'))
202 format ('*',2x, 2(1x,i7), 3(1x,G12.4), 5i3 )
