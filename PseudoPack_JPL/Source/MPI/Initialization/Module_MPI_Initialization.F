#include "pseudopack.h"



! ----------------------------------------------------------------------
! Module Name  : Processor
! Author       : Wai Sun Don
! Descritpion  : Initialize_CPU_And_Comunication and MPI Initialization
!  
! Subroutines  : PS_MPI_Initialize_Processors
! Module Used  : MPIF   (PARALLEL_MPI)
!
!  N_Processor : Total number of Processor used -- Output from the call to
!                  call MPI_comm_size (MPI_Comm_World,  N_Processor, IErr)
!                Initialized as 1
! ID_Processor : Processor ID                   -- output from the call to
!                  call MPI_comm_rank (MPI_Comm_World, ID_Processor, IErr)
!                Initialized as 0
!
!  First_Processor = 0
!   Last_Processor = N_Processor-1
!             I_Am = ID_Processor
! ----------------------------------------------------------------------

MODULE Processor

#if 1
  USE MPIF
#endif

implicit NONE

#if defined (PARALLEL_MPI)
  include "mpif.h"
#endif

integer ::     M_Processor = 1
integer ::     N_Processor = 1,   ID_Processor = 0
integer :: First_Processor = 0, Last_Processor = 0
integer ::    Processor_ID = 0,           I_Am = 0 , Number_Of_Processor = 1

integer  :: MPI_Communicator        = 1
integer  :: MPI_Error_Status        = 0
integer  :: MPI_Error_Status_Global = 0

PUBLIC 

PUBLIC  :: M_Processor
PUBLIC  :: N_Processor, ID_Processor, First_Processor, Last_Processor
PUBLIC  :: Processor_ID, I_Am, Number_Of_Processor

PUBLIC  :: MPI_Communicator, MPI_Error_Status, MPI_Error_Status_Global

PUBLIC  :: PS_MPI_Initialization
PUBLIC  :: PS_MPI_Processor
PUBLIC  :: PS_MPI_Processor_Info
PUBLIC  :: PS_Initialization

CONTAINS
!
!=======================================================================
!
  Subroutine PS_MPI_Initialization (MPI_Comm_Type)

  integer , OPTIONAL :: MPI_Comm_Type

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

#if defined (PARALLEL_MPI)
  call MPI_init      (                                MPI_Error_Status)
  if (MPI_Error_Status /= MPI_SUCCESS) STOP

  call MPI_comm_rank (MPI_Communicator, ID_Processor, MPI_Error_Status)
  call MPI_comm_size (MPI_Communicator,  N_Processor, MPI_Error_Status)
  call MPI_comm_size (MPI_Communicator,  M_Processor, MPI_Error_Status)
#else
  MPI_Communicator = -1
#endif

  END Subroutine PS_MPI_Initialization
!
!=======================================================================
!
  Subroutine PS_MPI_Processor_Info (N_Processor, ID_Processor,           &
                                    First_Processor, Last_Processor, I_Am)

  integer  :: N_Processor, ID_Processor
  integer  :: First_Processor, Last_Processor, I_Am

  First_Processor = 0
   Last_Processor =  N_Processor-1
     Processor_ID = ID_Processor

             I_Am = Processor_ID 

  END Subroutine PS_MPI_Processor_Info
!
!=======================================================================
!
  Subroutine PS_MPI_Processor (I_Am, Last_Processor)

  integer  :: I_Am, Last_Processor

  First_Processor = 0
   Last_Processor =  N_Processor-1
     Processor_ID = ID_Processor

             I_Am = Processor_ID 

  END Subroutine PS_MPI_Processor 
!
!=======================================================================
!
  Subroutine PS_Initialization

  call PS_MPI_Initialization
  call PS_MPI_Processor       (I_Am, Last_Processor)

  END Subroutine PS_Initialization

END MODULE Processor
