  if (MR_Activation == 1) then
!   if (I_Am == First_Processor) then
      if (MR_Output_Init == 1) then
        call PS_Make_Filename (MR_Flag_IC_File, '.dat', MR_Flag_Filename, &
                               Step, I_Am)
        OPEN (UNIT=79, FILE=MR_Flag_Filename, STATUS='UNKNOWN', IOSTAT=IOS)

        write (79,300) ; MR_Output_Init = 0
!     endif
    endif
  endif
