if (I_Am == First_Processor) then
  lid60=60

  CLOSE (lid60) 

  call PS_Make_Filename (Output_File, '.dat', Output_Filename, Restart)

  OPEN (UNIT=lid60, FILE=TRIM(Output_Filename), STATUS='UNKNOWN', &
                    POSITION='APPEND', IOSTAT=IOS)
endif
