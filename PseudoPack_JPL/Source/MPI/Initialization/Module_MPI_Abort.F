MODULE MPI_Abort_Routine

  USE Processor,  ONLY: MPI_Comm_World, MPI_Error_Status
  USE Processor,  ONLY: I_Am, First_Processor, Last_Processor
  USE PS_IO_Unit, ONLY: lid6, lid99

implicit NONE

PRIVATE 
PUBLIC  :: PS_STOP

CONTAINS
!
!======================================================================
!
  Subroutine PS_STOP (Sub_Name, MPI_Call, Comment)

  character(LEN=256)         :: C_Sub_Name, C_MPI_Call, C_Comment
  character(LEN=*), OPTIONAL ::   Sub_Name,   MPI_Call,   Comment
  
  if (I_Am == First_Processor) then
    C_Sub_Name = ' ' ; if (PRESENT(Sub_Name)) C_Sub_Name = TRIM(Sub_Name)
    C_MPI_Call = ' ' ; if (PRESENT(MPI_Call)) C_MPI_Call = TRIM(MPI_Call)
    C_Comment  = ' ' ; if (PRESENT(Comment )) C_Comment  = TRIM(Comment )

    write (lid6 ,100) I_Am, TRIM(C_Sub_Name), TRIM(C_Comment), &
                            TRIM(C_MPI_Call), MPI_Error_Status
    write (lid99,100) I_Am, TRIM(C_Sub_Name), TRIM(C_Comment),  &
                            TRIM(C_MPI_Call), MPI_Error_Status
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
  call MPI_ABORT   (MPI_Comm_World, MPI_Error_Status)
#else
  STOP
#endif

 100 format (1x,72('*')/                     &
             1x,'  Processor ID    : ',i5/   &
             1x,'  Subroutine Name : ',(a)/  &     
             1x,'    Comment       : ',(a)/  &     
             1x,'    MPI_Call      : ',(a)/  &
             1x,'  Error Code      : ',i5/)  

  END Subroutine PS_STOP

END MODULE MPI_Abort_Routine
