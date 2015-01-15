MODULE MPI_Synchronize

  USE Processor, ONLY: MPI_Comm_World, MPI_Error_Status

implicit NONE

PRIVATE 
PUBLIC  :: PS_MPI_Synchronize

CONTAINS
!
!======================================================================
!
  Subroutine PS_MPI_Synchronize (MPI_Comm_Type)

  Integer, OPTIONAL :: MPI_Comm_Type

#if defined (PARALLEL_MPI)
  if (PRESENT(MPI_Comm_Type)) then
    call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)
  else
    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
  endif
#endif

  END Subroutine PS_MPI_Synchronize

END MODULE MPI_Synchronize
