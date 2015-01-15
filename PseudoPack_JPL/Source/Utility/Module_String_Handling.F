! This module contains String Handling procedures.
! Built and tested with Compaq Visual Fortran v6.1a
! Test platform: Windows98 second edition
! Usage: see 'test.f90'
!
! Please use or modify freely.
! Please send me feedback or improved versions of the module. Thank you.
! I'll repost the latest best version to comp.lang.fortran.
!
! Author: fortranuser@aol.com
! release 001: spring 2001

! Cleanup by Wai Sun Don

MODULE String_Handling

  IMPLICIT NONE

PRIVATE

PUBLIC :: OPERATOR(+), OPERATOR(.LTU.), OPERATOR(.UTL.)
PUBLIC :: ASSIGNMENT(=)

INTERFACE OPERATOR(+)
  MODULE PROCEDURE Concatenation
END INTERFACE

INTERFACE OPERATOR(+)
  MODULE PROCEDURE String_Plus_Integer
END INTERFACE

INTERFACE OPERATOR(.LTU.)
  MODULE PROCEDURE Lower_to_Upper_Case
END INTERFACE

INTERFACE OPERATOR(.UTL.)
  MODULE PROCEDURE Upper_to_Lower_Case
END INTERFACE

INTERFACE ASSIGNMENT(=)
! MODULE PROCEDURE Integer_to_String
  MODULE PROCEDURE Integer_to_String_Padded
END INTERFACE

CONTAINS

  FUNCTION Concatenation (a, b) RESULT (x)

  character(LEN=*)            , INTENT(IN) :: a, b
  character(LEN=LEN(a)+LEN(b))             :: x

  if ((LEN(a) == 1) .AND. (a == ' ')) then
    x = a//b
  else
    x = TRIM(a)//b
  endif

  END FUNCTION Concatenation
!
! ===================================================================
!
  FUNCTION String_Plus_Integer (a, b) RESULT (x)

  character(LEN=*)       , INTENT(IN) :: a
  integer ,                INTENT(IN) :: b
  character(LEN=LEN(a)+8)             :: x
  character(LEN=8)                    :: Integer_String

  call Integer_to_String (Integer_String, b)

  x = a+Integer_String

  END FUNCTION String_Plus_Integer
!
! ===================================================================
!
  FUNCTION Lower_to_Upper_Case (String) RESULT (x)

  character(LEN=*)          , INTENT(IN) :: String
  character(LEN=LEN(String))             :: x

  integer  :: agcollat, apcollat, collat, offset, i, zpcollat
  INTRINSIC CHAR, ICHAR, LEN
  
  apcollat = ICHAR ('a')
  zpcollat = ICHAR ('z')
  agcollat = ICHAR ('A')

  offset = apcollat - agcollat

  do i = 1, LEN(String)

    collat = ICHAR(String(i:i))
    if ((collat >= apcollat) .AND. (collat <= zpcollat)) then
      x(i:i) = CHAR(collat - offset)
    else
      x(i:i) = String(i:i)
    endif

  enddo

  END FUNCTION Lower_to_Upper_Case
!
! ===================================================================
!
  FUNCTION Upper_to_Lower_Case (String) RESULT (x)

  character(LEN=*)          , INTENT(IN) :: String
  character(LEN=LEN(String))             :: x

  integer  :: agcollat, apcollat, collat, offset, i, zgcollat
  INTRINSIC CHAR, ICHAR, LEN
  
  apcollat = ICHAR('a')
  agcollat = ICHAR('A')
  zgcollat = ICHAR('Z')

  offset = apcollat - agcollat

  do i = 1, LEN(String)

    collat = ICHAR(String(i:i))
    if ((collat >= agcollat) .AND. (collat <= zgcollat)) then
      x(i:i) = CHAR(collat + offset)
    else
      x(i:i) = String(i:i)
    endif

  enddo

  END FUNCTION Upper_to_Lower_Case
!
! ===================================================================
!
  SUBROUTINE Integer_to_String (String, i)

  integer ,         INTENT(IN)  :: i
  character(LEN=*), INTENT(OUT) :: String
  character(LEN=5)              :: Integer_Format = '(Inn)'

  if (LEN(String) <= 0) RETURN

  write (Integer_Format(3:4),'(I2.2)') LEN(String)
  write (String,Integer_Format       ) ABS(i)

  String = ADJUSTL(String)

  if (i <  0) String = '-'//String

  END SUBROUTINE Integer_to_String
!
! ===================================================================
!
  SUBROUTINE Integer_to_String_Padded (String, i)

  integer ,         INTENT(IN)  :: i
  character(LEN=*), INTENT(OUT) :: String
  character(LEN=8)              :: Integer_Format = '(Inn.nn)'

  if (LEN(String) <= 0) RETURN

  write (Integer_Format(3:4),'(I2.2)') LEN(String)

  Integer_Format(6:7) = Integer_Format(3:4)
  write (String,Integer_Format) ABS(i)

  if (i <  0) String = '-'//String

  END SUBROUTINE Integer_to_String_Padded

END MODULE String_Handling





#if 0
!Program to test the above Module is the following:

PROGRAM String_Test

  USE String_Handling

  integer ::i
  character(LEN=10) :: iString, one, two, three
  character(LEN=256) :: filename

! Integer to String test
!-----------------------
  print '(A)','>>>>> Integer to String test'
  i=-12345
  iString=i
  print *,iString

  i=123456789
  iString=i
  print *,iString
  pause

! Concatenation test
!-------------------
  print '(A)','>>>>> Concatenation test'
  one='abc'
  two='def'
  three=one+two
  print *,three

  three=one + ":" + two
  print *,three
  three=one + (" " + two)
  print *,three
  pause

! Upper/Lower Case test
!----------------------
  print '(A)','>>>>> Upper/Lower Case test'
  one='abc'
  three=.LTU. one
  print *,three

  three=.UTL. three
  print *,three
  pause

! Other tests
  print '(A)',"Other tests"
  print *, .LTU. 'abc123defg'
  print *, .UTL.(.LTU. 'abc123defg')
  print *, .LTU.(.UTL.(.LTU. 'abc123defg'))

  print *,"'abcd  ' + '1234'"
  print *,'abcd  ' + '1234'
  print *,"'abcd  ' + ' 1234'"
  print *,'abcd  ' + ' 1234'
  print *,'abc'+(' ' + (' ' +'def'))
  pause

  do i=1,10
    filename='file'+'_'+i+'.dat'+'_'+(i+1)
    print *,'Processing ',TRIM(filename)
  enddo

  pause

END PROGRAM String_Test
#endif
