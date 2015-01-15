if (I_Am == First_Processor) then
  lid61=61

  OPEN (UNIT=lid61, FILE=TRIM(Auto_Restart_File), STATUS='UNKNOWN', &
                    POSITION='REWIND', IOSTAT=IOS)

#if defined (PARTICLE)
  write (lid61,*) Step, Time, dt, N_Processor, P(1)%N
#else
  write (lid61,*) Step, Time, dt, N_Processor
#endif

  CLOSE (lid61)
endif
