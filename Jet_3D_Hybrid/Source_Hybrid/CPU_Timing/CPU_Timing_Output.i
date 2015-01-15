  if (I_Am == First_Processor) then
    CPU_Usage = CPU_RK_END-CPU_Start

    call PS_Convert_Second_To_Clock (INT(CPU_Usage), CPU_Used)

    do lid_loop = 1,lid_loop_Max
      if (lid_loop == 1) lid = lid60

      if (lid_loop == 2) then
        lid = lid6 
        if ((Max_Output_Step >= 0) .AND. (Step >  Max_Output_Step) .AND.  &
            (Step/CPU_Step_Factor*CPU_Step_Factor /= Step)) then
             Max_Output_Step = 2*Max_Output_Step 
             CPU_Step_Factor = MIN(2*CPU_Step_Factor,100)
             EXIT
        endif
      endif

      write (lid,202) INT((Final_Time-Time)/dt_Original), Step,  &
                      dt_Original*Time_Factor, Time*Time_Factor, &
                      CPU_RK_End-CPU_RK_Begin, CPU_Used
    enddo
  endif
