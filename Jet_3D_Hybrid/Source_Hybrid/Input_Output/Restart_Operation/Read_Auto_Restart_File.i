if (Auto_Restart == 1) then
  if (I_Am == First_Processor) &
    write (6,*) 'Reading .... ', TRIM(Auto_Restart_File)

  lid61=61

  OPEN (UNIT=lid61, FILE=TRIM(Auto_Restart_File), STATUS='OLD', &
                    POSITION='REWIND', IOSTAT=IOS)

  if (IOS /= 0) then
    Restart = 0 
  else
    read (lid61,*) Restart
  endif

  CLOSE (lid61)

  if (I_Am == First_Processor) &
    write (6,*) '....... Done ', TRIM(Auto_Restart_File)
endif
