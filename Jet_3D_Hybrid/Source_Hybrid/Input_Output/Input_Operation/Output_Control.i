  OPEN (UNIT=LID, FILE=TRIM(Output_Input_File), STATUS='OLD', IOSTAT=IOS)

  read (LID,'(//)')
  read (LID,*) Graphic_Format
  read (LID,*)    Data_Format

  read (LID,*)
  read (LID,*) Write_IC  
  read (LID,*) Write_Initial_Movie
  read (LID,*) Write_Initial_MR_Flag
  read (LID,*) Write_Final

  read (LID,*)
  read (LID,*) Lid_Loop_Max

  read (LID,*)
  read (LID,*)   Movie_Resolution
  read (LID,*) Restart_Resolution, Restart_Distributed_Read,   &
                                   Restart_Distributed_Write

  CLOSE (LID)
