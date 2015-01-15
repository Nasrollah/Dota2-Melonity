#include "pseudopack.h"
 

!***************************************************************************
! MODULE      :: Timer
!                 
! PURPOSE     :: CPU Timer and utility subroutine for timing.
!
! VARIABLE    :: Timing     (Logical
!                            Decided on whether timing should be done
!
!                CPU_Timer  (REALTYPE Array of size (1:32)
!
!                CLOCK      (Integer Array of size (5)
!                            Store the unit (Week, Day, Hour, Minute, Second)
!
! SUBROUTINE  :: PS_Convert_Second_To_Clock (Second, Clock)
!                PS_Convert_Clock_To_Second (Clock, Second)
!
!                Note :  These two subroutine convert the time in Second
!                        format to (Week, Day, Hour, Minute, Second) format.
!
! ARGUMENTS   :: Second  (Integer
!                               Time in Second format
!
!                Clock   (Integer Array of size (5)
!                               Time in (Week, Day, Hour, Minute, Second) format
!
! FUNCTION    :: CPU_TIME 
!
! ARGUMENTS   :: Dummy    (Integer, OPTIONAL, INPUT
!
! USAGE       :: Standard
! 
! MODULE USED :: NONE
!
! AUTHOR      :: Wai Sun Don, PhD
! DATE        :: 3-14-1999
! VERSION     :: 2.0
!***************************************************************************


MODULE Timer

#if defined (PARALLEL_MPI)
! USE MPIF, ONLY: MPI_WTIME
  USE Processor, ONLY: MPI_WTIME
#endif

implicit NONE

logical  :: Timing = .TRUE.

integer , parameter :: Number_Of_Timer = 32

REALTYPE, dimension(Number_Of_Timer) :: CPU_Timer

integer , parameter :: N_Clock = 5

integer , dimension(N_Clock), parameter :: &
   Clock_Unit = (/ ONE_WEEK, ONE_DAY, ONE_HOUR, ONE_MINUTE, ONE_SECOND /)

PRIVATE

PUBLIC :: Timing, Number_Of_Timer, CPU_Timer, N_Clock, Clock_Unit

PUBLIC :: PS_Convert_Second_To_Clock
PUBLIC :: PS_Convert_Clock_To_Second

PUBLIC :: CPU_TIME

#if defined (PARALLEL_MPI)
! Intel compiler said that it can't have attribute.
#if 0
PUBLIC :: MPI_WTIME
#endif
#endif

CONTAINS
!
!======================================================================
!
  Subroutine PS_Convert_Second_To_Clock (Second, Clock)
       
  integer                      :: Second, Second_Now, i
  integer , dimension(N_Clock) ::  Clock

        Second_Now = Second
  do i = 1,N_Clock
      Clock(i) = FLOOR(ONE*Second_Now/Clock_Unit(i))
    Second_Now = Second_Now-Clock(i)*Clock_Unit(i)
  enddo

  END Subroutine PS_Convert_Second_To_Clock 
!
!======================================================================
!
  Subroutine PS_Convert_Clock_To_Second (Clock, Second)
       
  integer                      :: Second
  integer , dimension(N_Clock) ::  Clock

  Second = SUM(Clock*Clock_Unit)

  END Subroutine PS_Convert_Clock_To_Second 
!
!======================================================================
!
!-----
#if 0
!-----
  REALTYPE Function CPU_TIME (Dummy)

  integer ,              OPTIONAL :: Dummy
  integer ,          dimension(8) :: Date_Time
  character(LEN=10), dimension(3) :: Big_Ben
  
  call Date_And_Time (Big_Ben(1), Big_Ben(2), Big_Ben(3), Date_Time)

  CPU_Time = Date_Time(5)*ONE_HOUR              +  &  ! Hour
             Date_Time(6)*ONE_MINUTE            +  &  ! Min
             Date_Time(7)                       +  &  ! Second
             Date_Time(8)*ONE/ONE_MILLISECOND         ! MilliSecond

  END Function CPU_TIME 
!-----
#else
!-----
  REALTYPE Function CPU_TIME (Dummy)

  integer ,              OPTIONAL :: Dummy
  integer                         :: CPU_Counter, CPU_Count_Rate
  
  call SYSTEM_CLOCK (CPU_Counter, CPU_Count_Rate)

  CPU_TIME = (ONE*CPU_Counter)/CPU_Count_Rate

  END Function CPU_TIME 
!-----
#endif
!-----

END MODULE Timer
