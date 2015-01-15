Module MPI_Communicator_Sub_Group

  USE Processor

implicit NONE

logical  :: Processor_Active_Status = .FALSE.

INTERFACE PS_MPI_Sub_Communicator
  MODULE PROCEDURE MPI_Sub_Communicator_2D 
  MODULE PROCEDURE MPI_Sub_Communicator_3D 
END INTERFACE

INTERFACE PS_MPI_Processor_Active_Status
  MODULE PROCEDURE MPI_Processor_Active_Status
END INTERFACE

INTERFACE PS_MPI_Processor_Reset
  MODULE PROCEDURE MPI_Processor_Reset
END INTERFACE

PRIVATE
PUBLIC  :: Processor_Active_Status
PUBLIC  :: PS_MPI_Sub_Communicator
PUBLIC  :: PS_MPI_Processor_Active_Status
PUBLIC  :: PS_MPI_Processor_Reset

CONTAINS

  Subroutine MPI_Sub_Communicator_2D (MPI_Communicator, MPI_Communicator_1D, &
                                      MPI_Done)

  integer                             :: MPI_Communicator
  integer , dimension(2)              :: MPI_Communicator_1D
  logical               , OPTIONAL    :: MPI_Done

  integer                             :: N_Dims, Status,  i, j
  logical , dimension(2)              :: Sub_Grid_1D

#if defined (PARALLEL_MPI)
  if (PRESENT(MPI_Done)) then
    if (MPI_Done) then
      do i = 1,2
        call MPI_COMM_FREE (MPI_Communicator_1D(i), MPI_Error_Status)
      enddo

      RETURN
    endif
  endif

  call MPI_TOPO_Test (MPI_Communicator, Status, MPI_Error_Status)

  if (Status /= MPI_CART) STOP 'MPI_Sub_Communicator : Non-Cartesian Topology'

  call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

  do i = 1,2
    do j = 1,2 
      Sub_Grid_1D(j) = (i == j) 
    enddo

    call MPI_CART_SUB (MPI_Communicator, Sub_Grid_1D, MPI_Communicator_1D(i), &
                                                      MPI_Error_Status)
  enddo
#endif

  END Subroutine MPI_Sub_Communicator_2D
!
! --------------------------------------------------------------------------
!
  Subroutine MPI_Sub_Communicator_3D (MPI_Communicator, MPI_Communicator_1D, &
                                      MPI_Communicator_2D, MPI_Done)

  integer                             :: MPI_Communicator
  integer , dimension(3)              :: MPI_Communicator_1D
  integer , dimension(3)              :: MPI_Communicator_2D
  logical               , OPTIONAL    :: MPI_Done

  integer                             :: N_Dims, Status,  i, j
  logical , dimension(3)              :: Sub_Grid_1D, Sub_Grid_2D

#if defined (PARALLEL_MPI)
  if (PRESENT(MPI_Done)) then
    if (MPI_Done) then
      do i = 1,3
        call MPI_COMM_FREE (MPI_Communicator_1D(i), MPI_Error_Status)
        call MPI_COMM_FREE (MPI_Communicator_2D(i), MPI_Error_Status)
      enddo

      RETURN
    endif
  endif

  call MPI_TOPO_Test (MPI_Communicator, Status, MPI_Error_Status)

  if (Status /= MPI_CART) STOP 'MPI_Sub_Communicator : Non-Cartesian Topology'

  call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

  do i = 1,3
    do j = 1,3 
      Sub_Grid_1D(j) = (i == j) 
      Sub_Grid_2D(j) = (i /= j) 
    enddo

    write (6,*) 'Sub_Grid', i,j
    call MPI_CART_SUB (MPI_Communicator, Sub_Grid_1D, MPI_Communicator_1D(i), &
                                                      MPI_Error_Status)
    write (6,*) 'Sub_Grid 1D : ', I_Am, Sub_Grid_1D(j), MPI_Error_Status
    call MPI_CART_SUB (MPI_Communicator, Sub_Grid_2D, MPI_Communicator_2D(i), &
                                                      MPI_Error_Status)
    write (6,*) 'Sub_Grid 2D : ', I_Am, Sub_Grid_2D(j), MPI_Error_Status
  enddo
#endif

  END Subroutine MPI_Sub_Communicator_3D
!
! --------------------------------------------------------------------------
!
  Subroutine MPI_Processor_Active_Status (I_Am, Processor, Active_Status) 

  integer , OPTIONAL  :: I_Am, Processor
  logical , OPTIONAL  :: Active_Status
  logical             :: Current_Status

#if defined (PARALLEL_MPI)
                              Current_Status = .TRUE.
  if (PRESENT(Active_Status)) Current_Status = Active_Status

                           Processor_Active_Status = Current_Status
  if (PRESENT(I_Am) .AND. (PRESENT(Processor))) then
    if (I_Am == Processor) Processor_Active_Status = Current_Status
  endif
#else
                           Processor_Active_Status = .TRUE.
#endif

  END Subroutine MPI_Processor_Active_Status 
!
! --------------------------------------------------------------------------
!
  Subroutine MPI_Processor_Reset      (Active_Status)

  logical , OPTIONAL :: Active_Status

#if defined (PARALLEL_MPI)
                              Processor_Active_Status = .TRUE.
  if (PRESENT(Active_Status)) Processor_Active_Status = Active_Status
#else
                              Processor_Active_Status = .TRUE.
#endif

  END Subroutine MPI_Processor_Reset

END Module MPI_Communicator_Sub_Group
