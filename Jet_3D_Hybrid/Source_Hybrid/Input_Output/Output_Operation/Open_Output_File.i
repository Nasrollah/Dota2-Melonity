  lid60=60

  call PS_Make_Filename (Output_File, '.dat', Output_Filename, Restart)

  OPEN (UNIT=lid60, FILE=TRIM(Output_Filename), STATUS='UNKNOWN', &
                    POSITION='APPEND', IOSTAT=IOS)
