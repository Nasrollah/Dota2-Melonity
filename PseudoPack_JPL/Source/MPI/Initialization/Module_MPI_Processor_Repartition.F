#include "pseudopack.h"

#if defined (PARALLEL_MPI)

MODULE MPI_Processor_Repartition

  USE Processor
  USE MPI_Processor_Topology
  USE MPI_Processor_Rank
  USE MPI_Synchronize

implicit NONE

INTERFACE PS_MPI_Processor_Repartition
  MODULE PROCEDURE Repartition_1D
END INTERFACE

PRIVATE
PUBLIC  :: PS_MPI_Processor_Repartition

CONTAINS

  Subroutine Repartition_1D (M, Topology, MPI_Communicator,   &
                             Processor_Topology_File)

  TYPE (Processor_Topology) :: Topology

  integer                 :: M, n
  integer                 :: MPI_Communicator
  integer                 :: MPI_World_Group, MPI_Group

  character(LEN=*), OPTIONAL :: Processor_Topology_File

  integer , dimension(:), ALLOCATABLE :: Group_Rank

  call PS_MPI_Synchronize    ()

  M_Processor = MIN(N_Processor,M)

  ALLOCATE (Group_Rank(0:M_Processor-1))

  do n = 0, M_Processor-1
    Group_Rank(n) = n
  enddo

  call MPI_Comm_Group  (MPI_Comm_World , MPI_World_Group,                 &
                                                    MPI_Error_Status)
  call MPI_Group_Incl  (MPI_World_Group, M_Processor, Group_Rank,         &
                                         MPI_Group, MPI_Error_Status)
  call MPI_Comm_Create (MPI_Comm_World , MPI_Group, MPI_Communicator,     &
                                                    MPI_Error_Status)

  DEALLOCATE (Group_Rank)

  if (MPI_Communicator /= MPI_Comm_NULL) then
    call MPI_Topology_Setup_1D (Topology, MPI_Communicator,  &
                                Processor_Topology_File)
  endif

  call PS_MPI_Synchronize    ()

  END Subroutine Repartition_1D 
!
! ========================================================================
!                  
  Subroutine MPI_Topology_Setup_1D (Topology, MPI_Communicator, &
                                    Processor_Topology_File)

  TYPE (Processor_Topology) :: Topology

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type

  character(LEN=*), OPTIONAL :: Processor_Topology_File
  character(LEN=256)         :: Processor_Topology_File_1D
                                        
  integer               :: Dimensions
  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid

                                 MPI_Comm_type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_type = MPI_Communicator

  Dimensions     = 1
  Distributed    = .FALSE. ; Periodicity    = .FALSE. ; Grid    = 0
  Distributed(1) = .TRUE.  ; Periodicity(1) = .FALSE. ; Grid(1) = 0

  if (MPI_Comm_Type /= MPI_Comm_Null) then
    call PS_MPI_Processor_Topology_Setup (Topology, Dimensions,           &
                                          Distributed, Periodicity, Grid, &
                                          MPI_Comm_Type)
  endif

#if defined (PARALLEL_MPI)
  if (MPI_Comm_Type /= MPI_Comm_Null) then
    call MPI_CART_GET    (Topology%MPI_Communicator,                 &
                          Topology%Dimension       ,                 &
                          Topology%Grid            ,                 &
                          Topology%Periodicity     ,                 &
                          Topology%Coordination    , MPI_Error_Status)
  endif
#endif

  if (.NOT. PRESENT(Processor_Topology_File)) RETURN

  Processor_Topology_File_1D = TRIM(Processor_Topology_File)//'_1D'

  if (MPI_Comm_Type /= MPI_Comm_Null) then
    call PS_MPI_Write_Processor_Topology (Topology%Dimension     ,     &
                                          Topology%Grid          ,     &
                                          Topology%Periodicity   ,     &
                                          Topology%Coordination  ,     &
                                          Processor_Topology_File_1D,  &
                                          MPI_Comm_Type                )
  endif

  END Subroutine MPI_Topology_Setup_1D 

END MODULE MPI_Processor_Repartition

#else

MODULE MPI_Processor_Repartition

  USE MPI_Processor_Topology

implicit NONE

INTERFACE PS_MPI_Processor_Repartition
  MODULE PROCEDURE Repartition_1D
END INTERFACE

PRIVATE
PUBLIC  :: PS_MPI_Processor_Repartition

CONTAINS

  Subroutine Repartition_1D (M, Topology, MPI_Communicator,   &
                             Processor_Topology_File)

  TYPE (Processor_Topology) :: Topology

  integer , OPTIONAL      :: M
  integer , OPTIONAL      :: MPI_Communicator

  character(LEN=*), OPTIONAL :: Processor_Topology_File

  RETURN

  END Subroutine Repartition_1D 

END MODULE MPI_Processor_Repartition

#endif
