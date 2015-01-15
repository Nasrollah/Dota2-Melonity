if (I_Am == First_Processor) then
  lid62=62

  OPEN (UNIT=lid62, FILE=TRIM(Restart_List_File), STATUS='UNKNOWN', &
                    POSITION='APPEND', IOSTAT=IOS)

  write (lid62,'(i20,2(1x,g21.14),1x,i10)') Step, Time, dt, N_Processor

  CLOSE (lid62)
endif
