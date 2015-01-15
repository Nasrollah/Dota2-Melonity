MODULE String_Maker

   USE PS_IO_Unit, ONLY: lid6, lid99

implicit NONE

INTERFACE PS_Make_String
  MODULE PROCEDURE PS_Make_String_0
  MODULE PROCEDURE PS_Make_String_1
END INTERFACE

PRIVATE

PUBLIC  :: PS_Make_String
PUBLIC  :: PS_Convert_Integer_To_String
PUBLIC  :: PS_Convert_String_To_Integer

CONTAINS
!
!=======================================================================
!
  Subroutine PS_Make_String_0 (String_Header, String_Trailer, String, &
                               Number_1, Number_2, Number_3)

  integer , OPTIONAL :: Number_1, Number_2, Number_3

  character(LEN=*)        :: String_Header, String_Trailer, String
  character(LEN=RANGE(1)) :: String_ID

        String = TRIM(String_Header)

   call Form_String (Number_1)
   call Form_String (Number_2)
   call Form_String (Number_3)

        String = TRIM(String)//TRIM(String_Trailer)

  CONTAINS

    Subroutine Form_String (Number)
  
    integer , OPTIONAL :: Number
    integer            :: N_Digits
 
    if (PRESENT(Number)) then
      call PS_Convert_Integer_To_String (Number, N_Digits, String_ID)

      if (LEN_TRIM(String_ID(1:N_Digits)) /= 0) then
        String = TRIM(String)//'_'

        if (Number <  0) String = TRIM(String)//'-'

        String = TRIM(String)//String_ID(1:N_Digits)
      endif
    endif

    END Subroutine Form_String 

  END Subroutine PS_Make_String_0
!
!=======================================================================
!
  Subroutine PS_Make_String_1 (String_Header, String_Trailer, String, &
                               Number)

  integer , dimension(:)           :: Number
  integer                          :: n

  character(LEN=*)        :: String_Header, String_Trailer, String
  character(LEN=RANGE(1)) :: String_ID

        String = TRIM(String_Header)

     do n = 1,SIZE(Number)
       call Form_String (Number(n))
     enddo

        String = TRIM(String)//TRIM(String_Trailer)

  CONTAINS

    Subroutine Form_String (Number)
  
    integer  :: Number
    integer  :: N_Digits
 
      call PS_Convert_Integer_To_String (Number, N_Digits, String_ID)

      if (LEN_TRIM(String_ID(1:N_Digits)) /= 0) then
        String = TRIM(String)//'_'

        if (Number <  0) String = TRIM(String)//'-'

        String = TRIM(String)//String_ID(1:N_Digits)
      endif

    END Subroutine Form_String 

  END Subroutine PS_Make_String_1
!
!=======================================================================
!
  Subroutine PS_Convert_Integer_To_String (Number, N_Digits, String)

  integer          :: N_Digits, Number
  character(LEN=*) :: String 

  integer                       :: k, Digit, N
  integer , dimension(RANGE(1)) :: Digit_i

  N = ABS(Number)

  if (N <  10) then
    N_Digits = 1 ; Digit_i(1) = N
  else
    do k = 1, RANGE(1)
      if (N/10**k == 0) then ; N_Digits = k ; EXIT ; endif
    enddo

      Digit      = 0
    do k = 1, N_Digits
      Digit_i(k) = (N-Digit)    /10**(N_Digits-k)

      Digit      = Digit + Digit_i(k)*10**(N_Digits-k) 
    enddo
  endif

  write (String,'(10i1)') Digit_i

  END Subroutine PS_Convert_Integer_To_String
!
!=======================================================================
!
  Subroutine PS_Convert_String_To_Integer (String, N_Digits, Number)

  integer          :: N_Digits, Number
  character(LEN=*) :: String

  integer                   :: k, L, ABS_Number
  integer                   :: Leading_Number, Leading_String
  character(LEN=RANGE(1)+1) :: AdjustL_String

  do k = 2,LEN(String)-1
    if ( (String(k-1:k-1) /= ' ') .AND.         &
         (String(k  :k  ) == ' ') .AND.         &
         (String(k+1:k+1) /= ' ') ) then
      write (lid6 ,100) TRIM(ADJUSTL(String))
      write (lid99,100) TRIM(ADJUSTL(String))
      STOP
    endif
  enddo

  AdjustL_String = ' '

  L = MIN(RANGE(1)+1,LEN_TRIM(String))

  AdjustL_String(1:L) = ADJUSTL(String(1:L))

  read (AdjustL_String, '(i1)') Leading_String

  Leading_Number = HUGE(1)/10**RANGE(1)

  if ( (LEN_TRIM(AdjustL_String) >  RANGE(1)) .AND. &
       (Leading_String >  Leading_Number) ) then
    write (lid6 ,101) TRIM(String), HUGE(1)
    write (lid99,101) TRIM(String), HUGE(1)
    STOP
  endif

  read (AdjustL_String, '(i10)') Number

  ABS_Number = ABS(Number)

  if (ABS_Number <  10) then
    N_Digits = 1
  else
    do k = 1, RANGE(1)+1
      if (ABS_Number/10**k == 0) then ; N_Digits = k ; EXIT ; endif
    enddo
  endif

 100 format (1x/1x,72('*')/ &
             15x,'STOP! Subroutine Convert_String_To_Integer'// &
             10x,'Character String "',(a),'" must contains no empty space!'/ &
             1x,72('*')/)
 101 format (1x/1x,72('*')/ &
             15x,'STOP! Subroutine Convert_String_To_Integer'//  &
             8x,'Overflow in converting Character String ',(a),' to Integer!'/ &
             1x/10x,'Maximum alllowed range for integer are +- ',i12/ &
             1x,72('*')/)

  END Subroutine PS_Convert_String_To_Integer

END MODULE String_Maker






#if 0
  if (PRESENT(Number_1) .OR. PRESENT(Number_2)) then

    if (PRESENT(Number_1)) then
      call PS_Convert_Integer_To_String (Number_1, N_Digits, String_ID)

      if (LEN_TRIM(String_ID(1:N_Digits)) /= 0) then
        String = TRIM(String)//'_'

        if (Number_1 <  0) String = TRIM(String)//'-'

        String = TRIM(String)//String_ID(1:N_Digits)
      endif
    endif

    if (PRESENT(Number_2)) then
      call PS_Convert_Integer_To_String (Number_2, N_Digits, String_ID)

      if (LEN_TRIM(String_ID(1:N_Digits)) /= 0) then
        String = TRIM(String)//'_'

        if (Number_2 <  0) String = TRIM(String)//'-'

        String = TRIM(String)//String_ID(1:N_Digits)
      endif
    endif

  endif
#endif
