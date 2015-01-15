  if (I_Am == First_Processor) &
    write (6,*) 'Reading .... ', TRIM(MR_Input_File)

  OPEN (UNIT=LID, FILE=TRIM(MR_Input_File), STATUS='OLD', IOSTAT=IOS)

  read (LID,'(//)')
  read (LID,*) MR_Activation
  read (LID,*) MR_Activate(1:Dimensions)
  read (LID,*) 

  read (LID,*) MR_Order
  read (LID,*) MR_Index
  read (LID,*) MR_Tolerance
  read (LID,*) MR_Jump_Factor
  read (LID,*) MR_i0, MR_i1
  read (LID,*) 
  read (LID,*) MR_Output_Screen

  read (LID,*) 
  read (LID,*) Hybrid_Start_Time 

  if (MR_Activation == 0) MR_Activate(1:Dimensions) = 0

  if (ANY(MR_Activate(1:Dimensions) == 1)) MR_Activation = 1

  CLOSE (LID)

  if (I_Am == First_Processor) &
    write (6,*) '....... Done ', TRIM(MR_Input_File)
