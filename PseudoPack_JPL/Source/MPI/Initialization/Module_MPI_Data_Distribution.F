#include "pseudopack.h"


!--------------------------------------------------------------------
! Module Name : Data Distribution 
!
! Subroutines : PS_MPI_Block_Distribution
!               PS_MPI_Data_Distribution
!
! Author      : Wai Sun Don
! Description : Distribute M data points to N_Processor 
!               See bottom for details.
!
! Module Used : Processor
!
! M         integer (Input)
! M_Remain  integer (Output)
! M_0       integer (Output)
! M_Star    integer (Output)
! M_Block   integer (Output)
! MM        integer (Output)
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------------
#if 0
!---------------------------------------------------------------------------

MODULE MPI_Data_Distribution
  USE Processor

implicit NONE

PRIVATE 
PUBLIC  :: PS_MPI_Block_Distribution
PUBLIC  :: PS_MPI_Data_Distribution

CONTAINS
!
!=======================================================================
!
  Subroutine PS_MPI_Block_Distribution (M, M_Remain, M_0, M_Star, &
                                               M_Block, MM)

  integer :: M, M_Remain, M_0, M_Star, M_Block, MM

  call PS_MPI_Data_Distribution (M, M_Remain, M_0, M_Star)

  M_Block = CEILING(ONE*M/N_Processor)
  MM      = M_Block*N_Processor

  END Subroutine PS_MPI_Block_Distribution
!
!=======================================================================
!
  Subroutine PS_MPI_Data_Distribution (M, M_Remain, M_0, M_Star)

  integer :: M, M_Remain, M_0, M_Star, L

  L        = M/N_Processor
  M_Remain = M-L*N_Processor

  if (I_Am <= M_Remain-1) then
    M_Star = L+1  ; M_0 = ID_Processor*M_Star
  else
    M_Star = L    ; M_0 = ID_Processor*M_Star+M_Remain
  endif

  END Subroutine PS_MPI_Data_Distribution

END MODULE MPI_Data_Distribution

!---------------------------------------------------------------------------
#else
!---------------------------------------------------------------------------

MODULE MPI_Data_Distribution
  USE Processor
  USE MPI_Processor_Topology

implicit NONE

PRIVATE
PUBLIC  :: PS_MPI_Block_Distribution
PUBLIC  :: PS_MPI_Data_Distribution

CONTAINS
!
!=======================================================================
!
  Subroutine PS_MPI_Block_Distribution (M, M_Remain, M_0, M_Star, &
                                           M_Block, MM, Index, MPI_Communicator)

  integer :: M, M_Remain, M_0, M_Star, M_Block, MM
  integer :: I_Am, Last_Processor, Number_Of_Processor

  integer , OPTIONAL                  :: Index, MPI_Communicator

#if 0
  integer                             :: N_Dims
  integer , dimension(:), ALLOCATABLE :: DIMS, COORDS
  logical , dimension(:), ALLOCATABLE :: PERIODS

  if (PRESENT(MPI_Communicator) .AND. PRESENT(Index)) then
    call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

    ALLOCATE (DIMS(N_Dims), PERIODS(N_Dims), COORDS(N_Dims))

    call MPI_CART_GET    (MPI_Communicator, N_Dims, DIMS, PERIODS, COORDS, &
                          MPI_Error_Status)

    Number_Of_Processor = DIMS(Index)

    DEALLOCATE (DIMS, PERIODS, COORDS)
  else
    Number_Of_Processor = N_Processor
  endif
#else
  call PS_MPI_Processor_Local_Info (I_Am, Last_Processor, &
                                    Index, MPI_Communicator)

    Number_Of_Processor = Last_Processor+1
#endif

  call PS_MPI_Data_Distribution (M, M_Remain, M_0, M_Star,    &
                                    Index, MPI_Communicator)

  M_Block = CEILING(ONE*M/Number_Of_Processor)
  MM      = M_Block*Number_Of_Processor

  call PS_MPI_Processor (I_Am, Last_Processor)

  END Subroutine PS_MPI_Block_Distribution
!
!=======================================================================
!
  Subroutine PS_MPI_Data_Distribution (M, M_Remain, M_0, M_Star, &
                                          Index, MPI_Communicator)

  integer :: M, M_Remain, M_0, M_Star, L
  integer :: I_Am, Last_Processor, Number_Of_Processor

  integer , OPTIONAL                  :: Index, MPI_Communicator

#if 0
  integer                             :: N_Dims
  integer , dimension(:), ALLOCATABLE :: DIMS, COORDS
  logical , dimension(:), ALLOCATABLE :: PERIODS

  if (PRESENT(MPI_Communicator) .AND. PRESENT (Index)) then
    call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

    ALLOCATE (DIMS(N_Dims), PERIODS(N_Dims), COORDS(N_Dims))

    call MPI_CART_GET    (MPI_Communicator, N_Dims, DIMS, PERIODS, COORDS, &
                          MPI_Error_Status)

    Number_Of_Processor = DIMS(Index)
    I_Am                = COORDS(Index)

    DEALLOCATE (DIMS, PERIODS, COORDS)
  else
    Number_Of_Processor = N_Processor
    I_Am                = ID_Processor
  endif
#else
  call PS_MPI_Processor_Local_Info (I_Am, Last_Processor, &
                                    Index, MPI_Communicator)

    Number_Of_Processor = Last_Processor+1
#endif

  L        = M/Number_Of_Processor
  M_Remain = M-L*Number_Of_Processor

  if (I_Am <= M_Remain-1) then
    M_Star = L+1  ; M_0 = I_Am*M_Star
  else
    M_Star = L    ; M_0 = I_Am*M_Star+M_Remain
  endif

  call PS_MPI_Processor (I_Am, Last_Processor)

  END Subroutine PS_MPI_Data_Distribution

END MODULE MPI_Data_Distribution

!---------------------------------------------------------------------------
#endif
!---------------------------------------------------------------------------




! ----------------------------------------------------------------------
! FUNCTION NAME: PS_MPI_Block_Distribution
! Author       : Wai Sun Don
! Descritpion  : Setup paramters for block distribution of the data array
!
!  M          Number of data points which will be distributed across
!               N_Processor in block of size M_Star 
!  M_Remain   Number of data points remained after an uniform
!               distribution of block size INT(M/N_Processor) 
!               = M-INT(M/N_Processor)*L
!  M_0        Starting location of the global data 
!               of each local data block starting from 0
!  M_Star     Size of local data block 
!  M_Block    Maximum block size = CEILING(M/N_Processor)
!  MM         = M_Block*N_Processor
!
!  Consider the Global data set of size M 
!     f(0),....,f(i),.....,f(M-1)
!  and distributed across N_Processor=NP .
!
!  Each local data block has size of M_Star
!
!  For the first M_Remain=M-INT(M/N_Processor)*N_Processor Processors,
!    it will have M_Star=M_Block number of data points.
!  For the rest of the N_Processor-M_Remain  Processors,
!    it will have M_Star=INT(M/N_Processor)=L=M_Block-1 number of data points.
!
!  For example, let Pn denote Processor n  and M_Remain=2, then
!    (P0, P1) will have M_Star=M_Block   number of data points and the rest 
!    (P2, P3) will have M_Star=M_Block-1 number of data points.
!                 
!             Global                                       Local 
! 
!   P0 : f(        0  ),.,f(  M_Block    -1)  g(1),.,g(M_Block)   
!   P1 : f(  M_Block  ),.,f(2*M_Block    -1)  g(1),.,g(M_Block)   
!   P2 : f(2*M_Block  ),.,f(2*M_Block+  L-1)  g(1),.,g(M_Block-1) 
!   P3 : f(2*M_Block+L),.,f(2*M_Block+2*L-1)  g(1),.,g(M_Block-1) 
!
! and 
!   P0 : M_0 = 0             M_Star = M_Block
!   P1 : M_0 =   M_Block     M_Star = M_Block
!   P2 : M_0 = 2*M_Block     M_Star = M_Block-1
!   P3 : M_0 = 2*M_Block+L   M_Star = M_Block-1
!
!  Let say, N_Processor = 4, M = 10, M_Remain=2, M_Block=3, L=2
!
!             Global              Local 
! 
!   P0 :   f(0),f(1),f(2)     g(1),g(2),g(3)   M_0=0  M_Star=3
!   P1 :   f(3),f(4),f(5)     g(1),g(2),g(3)   M_0=3  M_Star=3
!   P2 :   f(6),f(7)          g(1),g(2)        M_0=6  M_Star=2
!   P3 :   f(8),f(9)          g(1),g(2)        M_0=8  M_Star=2
!   
! ----------------------------------------------------------------------
