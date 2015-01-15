MODULE File_OPEN

#if defined (DEBUG)
  USE MPI_Abort_Routine
#endif

implicit NONE

PRIVATE  

PUBLIC  :: PS_OPEN_File
#if defined (DEBUG)
PUBLIC  :: PS_IO_Status_Check
#endif

CONTAINS
!
!=======================================================================
!
  Subroutine PS_OPEN_File (Unit, File, Form, Position, Status, IO_Status)

  integer                    :: Unit, IOS
  character(LEN=*)           :: File
  character(LEN=*), OPTIONAL :: Form, Position, Status
  integer         , OPTIONAL :: IO_Status
  logical                    :: Opened, Named
  character(LEN=80)          :: Filename

#if defined (DEBUG)
  if (.NOT. PRESENT(IO_Status)) then
    Inquire (UNIT=Unit, OPENED=Opened, NAMED=Named, NAME=Filename)

    if (Opened) then
      if (Named .AND. (TRIM(Filename) == TRIM(File))) then
        RETURN
      else
        call PS_IO_Status_Check (Unit, File, -1)
      endif
    endif
  endif
#endif

  if (PRESENT(Form)) then
    if (PRESENT(Position)) then
      if (PRESENT(Status)) then
        OPEN (UNIT=Unit, FILE=TRIM(File), FORM=TRIM(Form),  &
              POSITION=TRIM(Position), STATUS=TRIM(Status), IOSTAT=IOS)
      else
        OPEN (UNIT=Unit, FILE=TRIM(File), FORM=TRIM(Form),  &
              POSITION=TRIM(Position),                      IOSTAT=IOS)
      endif
    else
      if (PRESENT(Status)) then
        OPEN (UNIT=Unit, FILE=TRIM(File), FORM=TRIM(Form),  &
                                       STATUS=TRIM(Status), IOSTAT=IOS)
      else
        OPEN (UNIT=Unit, FILE=TRIM(File), FORM=TRIM(Form),  &
                                                            IOSTAT=IOS)
      endif
    endif
  else
    if (PRESENT(Position)) then
      if (PRESENT(Status)) then
        OPEN (UNIT=Unit, FILE=TRIM(File),                   &
              POSITION=TRIM(Position), STATUS=TRIM(Status), IOSTAT=IOS)
      else
        OPEN (UNIT=Unit, FILE=TRIM(File),                   &
              POSITION=TRIM(Position),                      IOSTAT=IOS)
      endif
    else
      if (PRESENT(Status)) then
        OPEN (UNIT=Unit, FILE=TRIM(File),                   &
                                       STATUS=TRIM(Status), IOSTAT=IOS)
      else
        OPEN (UNIT=Unit, FILE=TRIM(File),                   &
                                                            IOSTAT=IOS)
      endif
    endif
  endif

  if (PRESENT(IO_Status)) then
    IO_Status = IOS
#if defined (DEBUG)
  else
    call PS_IO_Status_Check (Unit, File, IOS)
#endif
  endif

  END Subroutine PS_OPEN_File

#if defined (DEBUG)
!
!=======================================================================
!
  Subroutine PS_IO_Status_Check (Unit, Filename, IO_Status)

    USE PS_IO_Unit, ONLY: lid6, lid99
    USE Processor

  implicit NONE

  integer           :: Unit, IO_Status
  character(LEN=*)  :: Filename
  logical           :: Named
  character(LEN=80) :: Existing_Filename

  SELECT CASE (IO_Status)
    CASE (-1)
      Inquire (UNIT=Unit, NAMED=Named, NAME=Existing_Filename)

      write (lid6 ,101) Unit, TRIM(Filename), TRIM(Existing_Filename)
      write (lid99,101) Unit, TRIM(Filename), TRIM(Existing_Filename)

    CASE (0)
      RETURN

    CASE DEFAULT
      write (lid6 ,100) I_Am, TRIM(Filename), Unit, IO_Status
      write (lid99,100) I_Am, TRIM(Filename), Unit, IO_Status

  END SELECT

  call PS_STOP ('Open_File', 'MPI_Abort')

 100 format (1x,72('*')/                                         &
             1x,'<-- FATAL ERROR : Unable to Open File -->'//    &
             1x,'                  Processor   : ',i4/           &
             1x,'                  Filename    : ',(a)/          &
             1x,'                  Unit Number : ',i4/           &
             1x,'                  IO Status   : ',i4)

 101 format (1x,72('*')/                                         &
             1x,'<-- FATAL ERROR : Inconsistent Filename -->'//  &
             1x,'                  Unit Number       : ',i4/     &
             1x,'                  Read in  Filename : ',(a)/    &
             1x,'                  Existing Filename : ',(a))

  END Subroutine PS_IO_Status_Check
#endif

END MODULE File_OPEN
