#include "pseudopack.h"


MODULE Multi_Resolution_Variable

  USE Processor
  USE MPI_Abort_Routine
  USE MPI_Processor_Topology
  USE MPI_Ghost_Cell_Exchange
  USE PS_IO_Unit

IMPLICIT NONE

Interface PS_MR_Variable
  Module Procedure MR_Variable_1D
  Module Procedure MR_Variable_2D
  Module Procedure MR_Variable_3D
END Interface

PRIVATE

PUBLIC  :: PS_MR_Variable
!PUBLIC  :: PS_MPI_Processor_Local_Index

CONTAINS

  Subroutine MR_Variable_1D (Index_x, Distributed_x,                     &
                                N0,    N2,    N3,    N5,       NV,    Q, &
                             MR_N0, MR_N2, MR_N3, MR_N5, MR_Index, MR_Q, &
                             MPI_Communicator)

  integer  :: Index_x,    N0,    N2,    N3,    N5,    NV
  integer  ::          MR_N0, MR_N2, MR_N3, MR_N5, MR_NV, MR_Index
  logical  :: Distributed_x

  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(MR_N0:MR_N5,1,1,1) :: MR_Q

  integer  , OPTIONAL :: MPI_Communicator

  integer  :: i
  integer  :: I_Am_i, First_Processor_i, Last_Processor_i

#if defined (MR_DEBUG)
  call Data_Check
#endif

  MR_Q = ZERO ; MR_NV = 1 

  MR_Q(N0:N5,1,1,1) = Q(N0:N5,MR_Index) 

  call PS_MPI_Ghost_Cell_Exchange (Index_x, Distributed_x,                    &
                                     MR_N0, MR_N2, MR_N3, MR_N5,              &
                                   MR_NV, MR_Q(:,:,1,1), MPI_Communicator)

  call PS_MPI_Processor_Local_Index (Index_x, I_Am_i, First_Processor_i,      &
                                        Last_Processor_i, MPI_Communicator)

  if ((I_Am_i == First_Processor_i) .AND. (MR_N0 <  N0)) then 
    do i = MR_N0,N0-1
      MR_Q(i,1,1,1) = Q(N0,MR_Index) 
    enddo
  endif
    
  if ((I_Am_i ==  Last_Processor_i) .AND. (MR_N5 >  N5)) then 
    do i = N5+1,MR_N5
      MR_Q(i,1,1,1) = Q(N5,MR_Index) 
    enddo
  endif
    
#if defined (MR_DEBUG)
  CONTAINS

    Subroutine Data_Check 

    logical  :: Error

    Error = .FALSE.

    if ((MR_N0 > N0) .OR. (MR_N5 < N5)) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_N0 > N0 or MR_N5 < N5'
      write (lid99,*) 'MR_Variable, MR_N0 > N0 or MR_N5 < N5'
    endif

    if  (MR_Index > NV) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_Index > NV'
      write (lid99,*) 'MR_Variable, MR_Index > NV'
    endif

    if  (MR_N5-MR_N3 /= MR_N2-MR_N0) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_N5-MR_N3 /= MR_N2-MR_N0'
      write (lid99,*) 'MR_Variable, MR_N5-MR_N3 /= MR_N2-MR_N0'
    endif

    if (Error) then
      write (lid6 ,*) 'My processor ID is ', I_Am
      write (lid99,*) 'My processor ID is ', I_Am

      call PS_STOP
    endif

    END Subroutine Data_Check 
#endif

  END Subroutine MR_Variable_1D
!
! ==============================================================================
!
  Subroutine MR_Variable_2D (Index_x, Distributed_x,                     &
                                N0,    N2,    N3,    N5,                 &
                             Index_y, Distributed_y,                     &
                                M0,    M2,    M3,    M5,       NV,    Q, &
                             MR_N0, MR_N2, MR_N3, MR_N5,                 &
                             MR_M0, MR_M2, MR_M3, MR_M5, MR_Index, MR_Q, &
                             MPI_Communicator)

  integer  :: Index_x,    N0,    N2,    N3,    N5
  integer  :: Index_y,    M0,    M2,    M3,    M5,    NV
  integer  ::          MR_N0, MR_N2, MR_N3, MR_N5
  integer  ::          MR_M0, MR_M2, MR_M3, MR_M5, MR_NV, MR_Index
  logical  :: Distributed_x
  logical  :: Distributed_y

  integer , OPTIONAL :: MPI_Communicator

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(MR_N0:MR_N5,MR_M0:MR_M5,1,1) :: MR_Q

  integer  :: i, j
  integer  :: I_Am_i, First_Processor_i, Last_Processor_i

#if defined (MR_DEBUG)
  call Data_Check
#endif

  MR_Q = ZERO ; MR_NV = 1 

  MR_Q(N0:N5,M0:M5,1,1) = Q(:,:,MR_Index) 

  call PS_MPI_Ghost_Cell_Exchange (Index_x, Distributed_x,                    &
                                     MR_N0, MR_N2, MR_N3, MR_N5,              &
                                   Index_y, Distributed_y,                    &
                                     MR_M0, MR_M2, MR_M3, MR_M5,              &
                                   MR_NV, MR_Q(:,:,:,1), MPI_Communicator)

  call PS_MPI_Processor_Local_Index (Index_x, I_Am_i, First_Processor_i,      &
                                        Last_Processor_i, MPI_Communicator)

  if ((I_Am_i == First_Processor_i) .AND. (MR_N0 <  N0)) then 
    do i = MR_N0,N0-1
      MR_Q(i,M0:M5,1,1) = Q(N0,:,MR_Index) 
    enddo
  endif
    
  if ((I_Am_i ==  Last_Processor_i) .AND. (MR_N5 >  N5)) then 
    do i = N5+1,MR_N5
      MR_Q(i,M0:M5,1,1) = Q(N5,:,MR_Index) 
    enddo
  endif
    
  call PS_MPI_Processor_Local_Index (Index_y, I_Am_i, First_Processor_i,      &
                                        Last_Processor_i, MPI_Communicator)

  if ((I_Am_i == First_Processor_i) .AND. (MR_M0 <  M0)) then 
    do j = MR_M0,M0-1
      MR_Q(N0:N5,j,1,1) = Q(:,M0,MR_Index) 
    enddo
  endif

  if ((I_Am_i ==  Last_Processor_i) .AND. (MR_M5 >  M5)) then 
    do j = M5+1,MR_M5
      MR_Q(N0:N5,j,1,1) = Q(:,M5,MR_Index) 
    enddo
  endif

#if defined (MR_DEBUG)
  CONTAINS

    Subroutine Data_Check 

    logical  :: Error

    Error = .FALSE.

    if ((MR_N0 > N0) .OR. (MR_N5 < N5)) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_N0 > N0 or MR_N5 < N5'
      write (lid99,*) 'MR_Variable, MR_N0 > N0 or MR_N5 < N5'
    endif

    if ((MR_M0 > M0) .OR. (MR_M5 < M5)) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_M0 > M0 or MR_M5 < M5'
      write (lid99,*) 'MR_Variable, MR_M0 > M0 or MR_M5 < M5'
    endif

    if  (MR_Index > NV) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_Index > NV'
      write (lid99,*) 'MR_Variable, MR_Index > NV'
    endif

    if  (MR_N5-MR_N3 /= MR_N2-MR_N0) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_N5-MR_N3 /= MR_N2-MR_N0'
      write (lid99,*) 'MR_Variable, MR_N5-MR_N3 /= MR_N2-MR_N0'
    endif

    if  (MR_M5-MR_M3 /= MR_M2-MR_M0) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_M5-MR_M3 /= MR_M2-MR_M0'
      write (lid99,*) 'MR_Variable, MR_M5-MR_M3 /= MR_M2-MR_M0'
    endif

    if (Error) then
      write (lid6 ,*) 'My processor ID is ', I_Am
      write (lid99,*) 'My processor ID is ', I_Am

      call PS_STOP
    endif

    END Subroutine Data_Check 
#endif

  END Subroutine MR_Variable_2D 
!
! ==============================================================================
!
  Subroutine MR_Variable_3D (Index_x, Distributed_x,                     &
                                N0,    N2,    N3,    N5,                 &
                             Index_y, Distributed_y,                     &
                                M0,    M2,    M3,    M5,                 &
                             Index_z, Distributed_z,                     &
                                K0,    K2,    K3,    K5,       NV,    Q, &
                             MR_N0, MR_N2, MR_N3, MR_N5,                 &
                             MR_M0, MR_M2, MR_M3, MR_M5,                 &
                             MR_K0, MR_K2, MR_K3, MR_K5, MR_Index, MR_Q, &
                             MPI_Communicator)

  integer  :: Index_x,    N0,    N2,    N3,    N5
  integer  :: Index_y,    M0,    M2,    M3,    M5
  integer  :: Index_z,    K0,    K2,    K3,    K5,    NV
  integer  ::          MR_N0, MR_N2, MR_N3, MR_N5
  integer  ::          MR_M0, MR_M2, MR_M3, MR_M5
  integer  ::          MR_K0, MR_K2, MR_K3, MR_K5, MR_NV, MR_Index
  logical  :: Distributed_x
  logical  :: Distributed_y
  logical  :: Distributed_z

  integer , OPTIONAL :: MPI_Communicator

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,1) :: MR_Q

  integer  :: i, j, k
  integer  :: I_Am_i, First_Processor_i, Last_Processor_i

#if defined (MR_DEBUG)
  call Data_Check
#endif

  MR_Q = ZERO ; MR_NV = 1 

  MR_Q(N0:N5,M0:M5,K0:K5,1) = Q(N0:N5,M0:M5,K0:K5,MR_Index) 

  call PS_MPI_Ghost_Cell_Exchange (Index_x, Distributed_x,                    &
                                     MR_N0, MR_N2, MR_N3, MR_N5,              &
                                   Index_y, Distributed_y,                    &
                                     MR_M0, MR_M2, MR_M3, MR_M5,              &
                                   Index_z, Distributed_z,                    &
                                     MR_K0, MR_K2, MR_K3, MR_K5,              &
                                   MR_NV, MR_Q(:,:,:,:), MPI_Communicator)

  call PS_MPI_Processor_Local_Index (Index_x, I_Am_i, First_Processor_i,      &
                                        Last_Processor_i, MPI_Communicator)

  if ((I_Am_i == First_Processor_i) .AND. (MR_N0 <  N0)) then 
    do i = MR_N0,N0-1
      MR_Q(i,M0:M5,K0:K5,1) = Q(N0,:,:,MR_Index) 
    enddo
  endif
    
  if ((I_Am_i ==  Last_Processor_i) .AND. (MR_N5 >  N5)) then 
    do i = N5+1,MR_N5
      MR_Q(i,M0:M5,K0:K5,1) = Q(N5,:,:,MR_Index) 
    enddo
  endif
    
  call PS_MPI_Processor_Local_Index (Index_y, I_Am_i, First_Processor_i,      &
                                        Last_Processor_i, MPI_Communicator)

  if ((I_Am_i == First_Processor_i) .AND. (MR_M0 <  M0)) then 
    do j = MR_M0,M0-1
      MR_Q(N0:N5,j,K0:K5,1) = Q(:,M0,:,MR_Index) 
    enddo
  endif

  if ((I_Am_i ==  Last_Processor_i) .AND. (MR_M5 >  M5)) then 
    do j = M5+1,MR_M5
      MR_Q(N0:N5,j,K0:K5,1) = Q(:,M5,:,MR_Index) 
    enddo
  endif
    
  call PS_MPI_Processor_Local_Index (Index_z, I_Am_i, First_Processor_i,      &
                                        Last_Processor_i, MPI_Communicator)

  if ((I_Am_i == First_Processor_i) .AND. (MR_K0 <  K0)) then
    do k = MR_K0,K0-1
      MR_Q(N0:N5,M0:M5,k,1) = Q(:,:,K0,MR_Index)
    enddo
  endif

  if ((I_Am_i ==  Last_Processor_i) .AND. (MR_K5 >  K5)) then
    do k = K5+1,MR_K5
      MR_Q(N0:N5,M0:M5,k,1) = Q(:,:,K5,MR_Index)
    enddo
  endif

#if defined (MR_DEBUG)
  CONTAINS

    Subroutine Data_Check 

    logical  :: Error

    Error = .FALSE.

    if ((MR_N0 > N0) .OR. (MR_N5 < N5)) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_N0 > N0 or MR_N5 < N5'
      write (lid99,*) 'MR_Variable, MR_N0 > N0 or MR_N5 < N5'
    endif

    if ((MR_M0 > M0) .OR. (MR_M5 < M5)) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_M0 > M0 or MR_M5 < M5'
      write (lid99,*) 'MR_Variable, MR_M0 > M0 or MR_M5 < M5'
    endif

    if ((MR_K0 > K0) .OR. (MR_K5 < K5)) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_K0 > K0 or MR_K5 < K5'
      write (lid99,*) 'MR_Variable, MR_K0 > K0 or MR_K5 < K5'
    endif

    if  (MR_Index > NV) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_Index > NV'
      write (lid99,*) 'MR_Variable, MR_Index > NV'
    endif

    if  (MR_N5-MR_N3 /= MR_N2-MR_N0) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_N5-MR_N3 /= MR_N2-MR_N0'
      write (lid99,*) 'MR_Variable, MR_N5-MR_N3 /= MR_N2-MR_N0'
    endif

    if  (MR_M5-MR_M3 /= MR_M2-MR_M0) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_M5-MR_M3 /= MR_M2-MR_M0'
      write (lid99,*) 'MR_Variable, MR_M5-MR_M3 /= MR_M2-MR_M0'
    endif

    if  (MR_K5-MR_K3 /= MR_K2-MR_K0) then
      Error = .TRUE.
      write (lid6 ,*) 'MR_Variable, MR_K5-MR_K3 /= MR_K2-MR_K0'
      write (lid99,*) 'MR_Variable, MR_K5-MR_K3 /= MR_K2-MR_K0'
    endif

    if (Error) then
      write (lid6 ,*) 'My processor ID is ', I_Am
      write (lid99,*) 'My processor ID is ', I_Am

      call PS_STOP
    endif

    END Subroutine Data_Check 
#endif

  END Subroutine MR_Variable_3D
#if 0
!
! -------------------------------------------------------------------------
!
#if defined (PARALLEL_MPI)
    Subroutine PS_MPI_Processor_Local_Index (Index, I_Am_i, First_Processor_i, &
                                             Last_Processor_i, MPI_Communicator)

    USE MPI_Processor_Topology

    integer            :: Index
    integer            :: I_Am, First_Processor, Last_Processor
    integer , OPTIONAL :: MPI_Communicator

    integer                             :: MPI_Communicator_1D
    integer                             :: N_Dims
    logical , dimension(:), ALLOCATABLE :: Sub_Grid

               I_Am_i =     I_Am
    First_Processor_i = First_Processor
     Last_Processor_i =  Last_Processor

    if (.NOT. PRESENT(MPI_Communicator)) RETURN

    call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

    ALLOCATE (Sub_Grid(N_Dims))

    Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

    call MPI_CART_SUB  (MPI_Communicator, Sub_Grid, MPI_Communicator_1D,  &
                                                    MPI_Error_Status)

    call PS_MPI_Processor_Local_Info (I_Am_i, Last_Processor_i, Index,    &
                                      MPI_Communicator_1D)

    DEALLOCATE (Sub_Grid)

    call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)

    END Subroutine PS_MPI_Processor_Local_Index 
#else
    Subroutine PS_MPI_Processor_Local_Index (Index, I_Am_i, First_Processor_i, &
                                             Last_Processor_i, MPI_Communicator)

    integer            :: Index
    integer            :: I_Am, First_Processor, Last_Processor
    integer , OPTIONAL :: MPI_Communicator

               I_Am_i =     I_Am
    First_Processor_i = First_Processor
     Last_Processor_i =  Last_Processor

    END Subroutine PS_MPI_Processor_Local_Index 
#endif
#endif

END MODULE Multi_Resolution_Variable
