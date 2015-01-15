MODULE MPI_Processor_Rank

  USE Processor, ONLY : MPI_Comm_World, MPI_Error_Status

implicit NONE

INTERFACE PS_MPI_Processor_Rank
  MODULE PROCEDURE Processor_Rank
END INTERFACE

PRIVATE
PUBLIC  :: PS_MPI_Processor_Rank

CONTAINS

#if 0
!
!=======================================================================
!
  Subroutine Processor_Rank (I_Am, First_Processor, Last_Processor, &
                                     MPI_Communicator)

  integer            :: First_Processor, Last_Processor, N_Processor
  integer            :: I_Am, ID_Processor
  integer , OPTIONAL :: MPI_Communicator

  integer            :: MPI_Comm_Type

#if defined (PARALLEL_MPI)
                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  call MPI_comm_size (MPI_Comm_Type,  N_Processor, MPI_Error_Status)
  call MPI_comm_rank (MPI_Comm_Type, ID_Processor, MPI_Error_Status)
#else
   N_Processor = 1 ; ID_Processor = 0 
#endif

  First_Processor = 0
   Last_Processor =  N_Processor-1
             I_Am = ID_Processor

  END Subroutine Processor_Rank
#else
!
!=======================================================================
!
  Subroutine Processor_Rank (I_Am, First_Processor, Last_Processor, &
                                       N_Processor, MPI_Communicator)

  integer            :: First_Processor, Last_Processor, N_Processor
  integer            :: I_Am, ID_Processor
  integer , OPTIONAL :: MPI_Communicator

  integer            :: MPI_Comm_Type

#if defined (PARALLEL_MPI)
                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  call MPI_comm_size (MPI_Comm_Type,  N_Processor, MPI_Error_Status)
  call MPI_comm_rank (MPI_Comm_Type, ID_Processor, MPI_Error_Status)
#else
   N_Processor = 1 ; ID_Processor = 0 
#endif

  First_Processor = 0
   Last_Processor =  N_Processor-1
             I_Am = ID_Processor

  END Subroutine Processor_Rank
#endif

END MODULE MPI_Processor_Rank
