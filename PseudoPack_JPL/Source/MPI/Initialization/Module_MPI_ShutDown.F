MODULE MPI_ShutDown

#if defined (PARALLEL_MPI)
  USE Processor, ONLY : MPI_Comm_World, MPI_Error_Status
#endif

implicit NONE

PRIVATE
PUBLIC  :: PS_MPI_ShutDown

CONTAINS
!
!======================================================================
!
  Subroutine PS_MPI_ShutDown

#if defined (PARALLEL_MPI)
  call MPI_Barrier  (MPI_Comm_World, MPI_Error_Status)
  call MPI_Finalize (                MPI_Error_Status)
#endif

  END Subroutine PS_MPI_ShutDown

END MODULE MPI_ShutDown
