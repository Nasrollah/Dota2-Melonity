#include "pseudopack.h"


! ----------------------------------------------------------------------
! Module Name  : Processor
! Author       : Wai Sun Don
! Descritpion  : 
!  
! Subroutines  : 
! Module Used  : MPIF   (PARALLEL_MPI)
!
! ----------------------------------------------------------------------

MODULE MPI_Processor_Topology

  USE Processor
  USE PS_IO_Unit
  USE MPI_ShutDown

implicit NONE

TYPE Processor_Topology
  integer  :: MPI_Communicator

  integer  :: Dimension
  logical  :: Reorder

  integer , dimension(3) :: Grid, Coordination
  logical , dimension(3) :: Periodicity
END TYPE Processor_Topology

INTERFACE PS_MPI_Processor_Topology_Setup
  MODULE PROCEDURE MPI_Processor_Topology_Setup
END INTERFACE

INTERFACE PS_MPI_Topology_Setup
  MODULE PROCEDURE Setup_MPI_Topology_1D 
  MODULE PROCEDURE Setup_MPI_Topology_2D 
  MODULE PROCEDURE Setup_MPI_Topology_3D 
END INTERFACE

INTERFACE PS_MPI_Processor_Local_Info
  MODULE PROCEDURE MPI_Processor_Local_Info_0
  MODULE PROCEDURE MPI_Processor_Local_Info_1
  MODULE PROCEDURE MPI_Processor_Local_Info_2
END INTERFACE

INTERFACE PS_MPI_Read_Processor_Topology
  MODULE PROCEDURE MPI_Read_Processor_Topology_0
  MODULE PROCEDURE MPI_Read_Processor_Topology_1
END INTERFACE

INTERFACE PS_MPI_Write_Processor_Topology
  MODULE PROCEDURE MPI_Write_Processor_Topology_0
  MODULE PROCEDURE MPI_Write_Processor_Topology_1
END INTERFACE

PRIVATE

PUBLIC  :: Processor_Topology
PUBLIC  :: PS_MPI_Processor_Topology_Setup
PUBLIC  :: PS_MPI_Processor_Local_Info
PUBLIC  :: PS_MPI_Read_Processor_Topology
PUBLIC  :: PS_MPI_Write_Processor_Topology
PUBLIC  :: PS_MPI_Processor_Local_Index

CONTAINS

  Subroutine MPI_Processor_Topology_Setup (Topology, Dimension,             &
                                           Distributed, Periodicity, Grid,  &
                                           MPI_Communicator)

  TYPE (Processor_Topology) :: Topology

  integer                          :: Dimension, N_Dimension, n
  logical , dimension(:)           :: Distributed
  logical , dimension(:), OPTIONAL :: Periodicity
  integer , dimension(:), OPTIONAL :: Grid

  integer                          :: M_Processor
  integer ,               OPTIONAL :: MPI_Communicator
  integer                          :: MPI_Comm_Type

                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  N_Dimension = MIN(Dimension, 3)

  Topology%Dimension           = N_Dimension
  Topology%Reorder             = .FALSE.
  Topology%Periodicity         = .FALSE.

  do n = 1,N_Dimension
                                Topology%Periodicity(n) = .FALSE.
      if (PRESENT(Periodicity)) Topology%Periodicity(n) = Periodicity(n)

    if (Distributed(n)) then
                                Topology%Grid(n) = 0
      if (PRESENT(Grid))        Topology%Grid(n) = Grid(n)
    else
                                Topology%Grid(n) = 1
    endif
  enddo
    
#if defined (PARALLEL_MPI)
  Topology%MPI_Communicator    = MPI_Comm_Type

  call MPI_Comm_Size   (MPI_Comm_Type, M_Processor, MPI_Error_Status)

  call MPI_DIMS_CREATE (M_Processor, Topology%Dimension,                  &
                                     Topology%Grid, MPI_Error_Status)

  call MPI_CART_CREATE (MPI_Comm_Type, Topology%Dimension,                &
                        Topology%Grid,                                    &
                        Topology%Periodicity, Topology%Reorder,           &
                        Topology%MPI_Communicator, MPI_Error_Status)

#else

  Topology%MPI_Communicator    = 0
  Topology%Grid                = 1
  Topology%Coordination        = 1

#endif

  do n = 1,N_Dimension
    if (Topology%Grid(n) == 1) Distributed(n) = .FALSE.
  enddo

  END Subroutine MPI_Processor_Topology_Setup

!
!=======================================================================
!
  Subroutine MPI_Processor_Local_Info_0 (I_Am, Last_Processor,        &
                                         Index, MPI_Communicator)

  integer                             :: I_Am, Last_Processor
  integer , OPTIONAL                  :: Index, MPI_Communicator

  integer                             :: N_Dims, i
  integer , dimension(:), ALLOCATABLE :: Grid, COORDS
  logical , dimension(:), ALLOCATABLE :: PERIODS

#if defined (PARALLEL_MPI)
  if (PRESENT(MPI_Communicator)) then
    call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

    ALLOCATE (Grid(N_Dims), PERIODS(N_Dims), COORDS(N_Dims))

    call MPI_CART_GET    (MPI_Communicator, N_Dims, Grid, PERIODS, COORDS, &
                          MPI_Error_Status)

    i = 1 ; if (PRESENT(Index)) i = Index

    I_Am           =  COORDS(i)
    Last_Processor =    Grid(i)-1

    DEALLOCATE (Grid, PERIODS, COORDS)
  else
              I_Am = Processor_ID 
    Last_Processor = N_Processor-1
  endif
#else
              I_Am = Processor_ID 
    Last_Processor = N_Processor-1
#endif

  END Subroutine MPI_Processor_Local_Info_0
!
!=======================================================================
!
  Subroutine MPI_Processor_Local_Info_1 (Periodicity,                      &
                                         Index, MPI_Communicator)

  logical                             :: Periodicity
  integer , OPTIONAL                  :: Index, MPI_Communicator

  integer                             :: N_Dims, i
  integer , dimension(:), ALLOCATABLE :: Grid, COORDS
  logical , dimension(:), ALLOCATABLE :: PERIODS

#if defined (PARALLEL_MPI)
  if (PRESENT(MPI_Communicator)) then
    call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

    ALLOCATE (Grid(N_Dims), PERIODS(N_Dims), COORDS(N_Dims))

    call MPI_CART_GET    (MPI_Communicator, N_Dims, Grid, PERIODS, COORDS, &
                          MPI_Error_Status)

    i = 1 ; if (PRESENT(Index)) i = Index

    Periodicity    = PERIODS(i)

    DEALLOCATE (Grid, PERIODS, COORDS)
  else
    Periodicity         = .FALSE.
  endif
#else
    Periodicity         = .FALSE.
#endif

  END Subroutine MPI_Processor_Local_Info_1
!
!=======================================================================
!
  Subroutine MPI_Processor_Local_Info_2 (I_Am, Last_Processor, Periodicity, &
                                         Index, MPI_Communicator)

  integer                             :: I_Am, Last_Processor
  logical                             :: Periodicity
  integer , OPTIONAL                  :: Index, MPI_Communicator

  integer                             :: N_Dims, i
  integer , dimension(:), ALLOCATABLE :: Grid, COORDS
  logical , dimension(:), ALLOCATABLE :: PERIODS

#if defined (PARALLEL_MPI)
  if (PRESENT(MPI_Communicator)) then
    call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

    ALLOCATE (Grid(N_Dims), PERIODS(N_Dims), COORDS(N_Dims))

    call MPI_CART_GET    (MPI_Communicator, N_Dims, Grid, PERIODS, COORDS, &
                          MPI_Error_Status)

    i = 1 ; if (PRESENT(Index)) i = Index

    I_Am           =  COORDS(i)
    Last_Processor =    Grid(i)-1
    Periodicity    = PERIODS(i)

    DEALLOCATE (Grid, PERIODS, COORDS)
  else
              I_Am = Processor_ID 
    Last_Processor = N_Processor-1
    Periodicity    = .FALSE.
  endif
#else
              I_Am = Processor_ID 
    Last_Processor = N_Processor-1
    Periodicity    = .FALSE.
#endif

  END Subroutine MPI_Processor_Local_Info_2
!
!=======================================================================
!
  Subroutine PS_MPI_Processor_Local_Index (Index, I_Am_i, First_Processor_i, &
                                           Last_Processor_i, MPI_Communicator)

  integer            :: Index
  integer            :: I_Am_i, First_Processor_i, Last_Processor_i
  integer , OPTIONAL :: MPI_Communicator

#if defined (PARALLEL_MPI)
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
#else
             I_Am_i =     I_Am
  First_Processor_i = First_Processor
   Last_Processor_i =  Last_Processor
#endif

  END Subroutine PS_MPI_Processor_Local_Index 

!#if defined (PARALLEL_MPI)
!
!=======================================================================
!
  Subroutine MPI_Read_Processor_Topology_0 (Dimension, Grid,              &
                                            Periodicity, lid, Error)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  integer                          :: lid
  logical , OPTIONAL               :: Error
  logical                          :: Error_p, Error_l
  integer                          :: M_Processor, Dims

  Error_p = .FALSE.

  read (lid,100) M_Processor
  if (M_Processor /= N_Processor) then
    write (lid6 ,200) N_Processor, M_Processor
    write (lid99,200) N_Processor, M_Processor
    Error_p = .TRUE.
  endif

  read (lid,100) Dims
  if (Dims /= Dimension) then
    write (lid6 ,201) Dimension, Dims
    write (lid99,201) Dimension, Dims
    Error_p = .TRUE.
  endif

#if defined (PARALLEL_MPI)
  call MPI_ALLREDUCE (Error_p, Error_l, 1, MPI_LOGICAL, MPI_LOR,       &
                               MPI_Comm_World, MPI_Error_Status)
#else
  Error_l = Error_p
#endif

  if (Error_l) then
    if (PRESENT(Error)) then
      CLOSE (lid) ; Error = Error_l ; RETURN 
    else
      call PS_MPI_ShutDown
    endif
  endif

  read (lid,101)        Grid(1:Dimension)
  read (lid,102) Periodicity(1:Dimension)

  CLOSE (lid)

 100 format (23x,i10)
 101 format (23x,3(i10,:))
 102 format (23x,3(l10,:))

 200 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Number of Processor Specified : ',i10        / &
             1x,'  Number of Processor Read in   : ',i10        / &
             1x,72(1x,'*')/)
 201 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Topology Dimension Specified  : ',i10        / &
             1x,'  Topology Dimension Read in    : ',i10        / &
             1x,72(1x,'*')/)

  END Subroutine MPI_Read_Processor_Topology_0
!
!=======================================================================
!
  Subroutine MPI_Write_Processor_Topology_0 (Dimension, Grid,                &
                                             Periodicity, Coordination, lid, &
                                             MPI_Communicator                )

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  integer , dimension(1:Dimension) :: Coordination
  integer                          :: lid

  integer , dimension(0:Dimension) :: A 

  integer                          :: M_Processor

#if defined (PARALLEL_MPI)
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)
#endif

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type

                                 MPI_Comm_type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_type = MPI_Communicator

#if defined (PARALLEL_MPI)
  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)
#else
  M_Processor = N_Processor
#endif

  A(0) = I_Am ; A(1:Dimension) = Coordination(1:Dimension) 

  if (I_Am  == First_Processor) then
    write (lid,100) M_Processor, Dimension
    write (lid,101) Grid 
    write (lid,102) Periodicity
  endif

  if (I_Am == First_Processor) write (lid,103) A 

#if defined (PARALLEL_MPI)
  do K_Processor = First_Processor+1,M_Processor-1

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

    if (I_Am ==     K_Processor)   &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,  &
                                 MPI_Comm_Type ,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,  &
                                 MPI_Comm_Type , Status, MPI_Error_Status)

      write (lid,103) A
    endif
  enddo

  call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)
#endif

100 format (1x,'Number of Processor : ', i10      / &
            1x,'Cartesian Dimension : ', i10      )
101 format (1x,'Cartesian Grid Size : ', 3(i10,:) )
102 format (1x,'Periodicity         : ', 3(l10,:) )
103 format (1x,'Rank/Coorindate     : ', i10, 3(i10,:))

  END Subroutine MPI_Write_Processor_Topology_0

!---------------------------------
#if 0
!---------------------------------
!
!=======================================================================
!
  Subroutine MPI_Read_Processor_Topology_0 (Dimension, Grid,              &
                                            Periodicity, lid)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  integer                          :: lid

  END Subroutine MPI_Read_Processor_Topology_0
!
!=======================================================================
!
  Subroutine MPI_Write_Processor_Topology_0 (Dimension, Grid,              &
                                             Periodicity, Coordination, lid)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  integer , dimension(1:Dimension) :: Coordination
  integer                          :: lid

  END Subroutine MPI_Write_Processor_Topology_0

!---------------------------------
#endif
!---------------------------------

!#if defined (PARALLEL_MPI)
!
!=======================================================================
!
  Subroutine MPI_Read_Processor_Topology_1 (Dimension, Grid,              &
                                            Periodicity, Filename, Error)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  character(LEN=*)                 :: Filename
  integer                          :: lid, IOS, IOS_p
  integer                          :: M_Processor, Dims
  logical                          :: Error_p, Error_l
  logical , OPTIONAL               :: Error

  IOS_p = 0 ; Error_p = .FALSE.

    lid = 77

    OPEN (UNIT=lid, FILE=TRIM(Filename), STATUS='OLD'    , POSITION='REWIND',  &
                    IOSTAT=IOS_p)

    if (IOS_p /= 0) then
      write (lid6 ,202) I_Am, lid, IOS_p, TRIM(Filename)
      write (lid99,202) I_Am, lid, IOS_p, TRIM(Filename)

      Error_p = .TRUE.
    endif

#if defined (PARALLEL_MPI)
  call MPI_ALLREDUCE (Error_p, Error_l, 1, MPI_LOGICAL, MPI_LOR,       &
                               MPI_Comm_World, MPI_Error_Status)
#else
  Error_l = Error_p
#endif

  if (Error_l) call PS_MPI_ShutDown

  Error_p = .FALSE.

  read (lid,100) M_Processor
  if (M_Processor /= N_Processor) then
    write (lid6 ,200) N_Processor, M_Processor
    write (lid99,200) N_Processor, M_Processor
    Error_p = .TRUE.
  endif

  read (lid,100) Dims
  if (Dims /= Dimension) then
    write (lid6 ,201) Dimension, Dims
    write (lid99,201) Dimension, Dims
    Error_p = .TRUE.
  endif

#if defined (PARALLEL_MPI)
  call MPI_ALLREDUCE (Error_p, Error_l, 1, MPI_LOGICAL, MPI_LOR,       &
                               MPI_Comm_World, MPI_Error_Status)
#else
  Error_l = Error_p
#endif

  if (Error_l) then
    if (PRESENT(Error)) then
      CLOSE (lid) ; Error = Error_l ; RETURN 
    else
      call PS_MPI_ShutDown
    endif
  endif

  read (lid,101)        Grid(1:Dimension)
  read (lid,102) Periodicity(1:Dimension)

  CLOSE (lid)

 100 format (23x,i10)
 101 format (23x,3(i10,:))
 102 format (23x,3(l10,:))

 200 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Number of Processor Specified : ',i10        / &
             1x,'  Number of Processor Read in   : ',i10        / &
             1x,72(1x,'*')/)
 201 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Topology Dimension Specified  : ',i10        / &
             1x,'  Topology Dimension Read in    : ',i10        / &
             1x,72(1x,'*')/)
 202 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Processor ID                  : ',i10        / &
             1x,'  File Unit, Status, Filename   : ',i5,i5,a32  )

  END Subroutine MPI_Read_Processor_Topology_1
!
!=======================================================================
!
  Subroutine MPI_Write_Processor_Topology_1 (Dimension, Grid,              &
                                             Periodicity, Coordination,    &
                                             Filename, MPI_Communicator)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  integer , dimension(1:Dimension) :: Coordination

  character(LEN=*)                 :: Filename
  integer                          :: lid, IOS, IOS_p
  logical                          :: Error_p, Error_l

  integer , dimension(0:Dimension) :: A 

  integer                          :: M_Processor

#if defined (PARALLEL_MPI)
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)
#endif

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type

                                 MPI_Comm_type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_type = MPI_Communicator

#if defined (PARALLEL_MPI)
  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)
#else
  M_Processor = N_Processor
#endif

  A(0) = I_Am ; A(1:Dimension) = Coordination(1:Dimension) 

  IOS_p = 0 ; Error_p = .FALSE.
  if (I_Am  == First_Processor) then
    lid = 78

    OPEN (UNIT=lid, FILE=TRIM(Filename), STATUS='UNKNOWN', POSITION='REWIND',  &
                    IOSTAT=IOS_p)

    if (IOS_p /= 0) then
      write (lid6 ,202) I_Am, lid, IOS_p, TRIM(Filename)
      write (lid99,202) I_Am, lid, IOS_p, TRIM(Filename)

      Error_p = .TRUE.
    endif
  endif

#if defined (PARALLEL_MPI)
  call MPI_ALLREDUCE (Error_p, Error_l, 1, MPI_LOGICAL, MPI_LOR,       &
                               MPI_Comm_Type , MPI_Error_Status)
#else
  Error_l = Error_p
#endif

  if (Error_l) call PS_MPI_ShutDown

  if (I_Am  == First_Processor) then
    write (lid,100) M_Processor, Dimension
    write (lid,101) Grid 
    write (lid,102) Periodicity
  endif

  if (I_Am == First_Processor) write (lid,103) A 

#if defined (PARALLEL_MPI)
  do K_Processor = First_Processor+1,M_Processor-1

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

    if (I_Am ==     K_Processor)   &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,  &
                                 MPI_Comm_Type ,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,  &
                                 MPI_Comm_Type , Status, MPI_Error_Status)


      write (lid,103) A
    endif
  enddo

  call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)
#endif

  if (I_Am == First_Processor) CLOSE (lid)

 100 format (1x,'Number of Processor : ', i10      / &
             1x,'Cartesian Dimension : ', i10      )
 101 format (1x,'Cartesian Grid Size : ', 3(i10,:) )
 102 format (1x,'Periodicity         : ', 3(l10,:) )
 103 format (1x,'Rank/Coorindate     : ', i10, 3(i10,:))

 202 format (1x,'FATAL ERROR : PS_MPI_Write_Processor_Topology '/ &
             1x,'  Processor ID                  : ',i10        / &
             1x,'  File Unit, Status, Filename   : ',i5,i5,a32  )

  END Subroutine MPI_Write_Processor_Topology_1

!----------------
#if 0
!----------------
!
!=======================================================================
!
  Subroutine MPI_Read_Processor_Topology_1 (Dimension, Grid,              &
                                            Periodicity, Filename)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  character(LEN=*)                 :: Filename

  END Subroutine MPI_Read_Processor_Topology_1
!
!=======================================================================
!
  Subroutine MPI_Write_Processor_Topology_1 (Dimension, Grid,              &
                                             Periodicity, Coordination,    &
                                             Filename)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  integer , dimension(1:Dimension) :: Coordination
  character(LEN=*)                 :: Filename

  END Subroutine MPI_Write_Processor_Topology_1
!----------------
#endif
!----------------




!
! ========================================================================
!  
  Subroutine Setup_MPI_Topology_1D (                        &
    Index_x, Distributed_x, Periodicity_x, Topology_Grid_x, &
    Topology, Processor_Topology_File, Read_Topology_From_File)
 
  logical :: Distributed_x, Periodicity_x
  integer :: Index_x, Topology_Grid_x

  integer :: Dimension

  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid
  
  TYPE (Processor_Topology) :: Topology

  character(LEN=*), OPTIONAL :: Processor_Topology_File
  logical         , OPTIONAL :: Read_Topology_From_File
 
  Dimension = 1

                   Distributed(Index_x) =   Distributed_x
                   Periodicity(Index_x) =   Periodicity_x
                          Grid(Index_x) = Topology_Grid_x
 
                   Distributed(2      ) = .FALSE.
                   Periodicity(2      ) = .FALSE.
                          Grid(2      ) = 0

                   Distributed(3      ) = .FALSE.
                   Periodicity(3      ) = .FALSE.
                          Grid(3      ) = 0

  if (PRESENT(Processor_Topology_File) .AND.   &
      PRESENT(Read_Topology_From_File)) then
    if (Read_Topology_From_File) then
      call PS_MPI_Read_Processor_Topology (Dimension, Grid, Periodicity,     &
                                           Processor_Topology_File           )
    endif
  endif

  call PS_MPI_Processor_Topology_Setup (Topology, Dimension, &
                                        Distributed, Periodicity, Grid)

  Distributed_x = Distributed(Index_x)

#if defined (PARALLEL_MPI)
  call MPI_CART_GET    (Topology%MPI_Communicator, Topology%Dimension,     &
                        Topology%Grid            , Topology%Periodicity,   &
                        Topology%Coordination    , MPI_Error_Status)
#endif

  if (.NOT. PRESENT(Processor_Topology_File)) RETURN
                                        
  call PS_MPI_Write_Processor_Topology (Topology%Dimension, Topology%Grid, &
                                        Topology%Periodicity,              &
                                        Topology%Coordination,             &
                                        Processor_Topology_File,           &
                                        Topology%MPI_Communicator          )

  END Subroutine Setup_MPI_Topology_1D
!
! ========================================================================
!  
  Subroutine Setup_MPI_Topology_2D (                        &
    Index_x, Distributed_x, Periodicity_x, Topology_Grid_x, &
    Index_y, Distributed_y, Periodicity_y, Topology_Grid_y, &
    Topology, Processor_Topology_File, Read_Topology_From_File)
 
  logical :: Distributed_x, Periodicity_x
  logical :: Distributed_y, Periodicity_y

  integer :: Index_x, Topology_Grid_x
  integer :: Index_y, Topology_Grid_y

  integer :: Dimension

  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid
  
  TYPE (Processor_Topology) :: Topology

  character(LEN=*), OPTIONAL :: Processor_Topology_File
  logical         , OPTIONAL :: Read_Topology_From_File
 
  Dimension = 2

                   Distributed(Index_x) =   Distributed_x
                   Periodicity(Index_x) =   Periodicity_x
                          Grid(Index_x) = Topology_Grid_x
 
                   Distributed(Index_y) =   Distributed_y
                   Periodicity(Index_y) =   Periodicity_y
                          Grid(Index_y) = Topology_Grid_y

                   Distributed(3      ) = .FALSE.
                   Periodicity(3      ) = .FALSE.
                          Grid(3      ) = 0

  if (PRESENT(Processor_Topology_File) .AND.   &
      PRESENT(Read_Topology_From_File)) then
    if (Read_Topology_From_File) then
      call PS_MPI_Read_Processor_Topology (Dimension, Grid, Periodicity,     &
                                           Processor_Topology_File           )
    endif
  endif

  call PS_MPI_Processor_Topology_Setup (Topology, Dimension, &
                                        Distributed, Periodicity, Grid)

  Distributed_x = Distributed(Index_x)
  Distributed_y = Distributed(Index_y)

#if defined (PARALLEL_MPI)
  call MPI_CART_GET    (Topology%MPI_Communicator, Topology%Dimension,     &
                        Topology%Grid            , Topology%Periodicity,   &
                        Topology%Coordination    , MPI_Error_Status)
#endif

  if (.NOT. PRESENT(Processor_Topology_File)) RETURN
                                        
  call PS_MPI_Write_Processor_Topology (Topology%Dimension, Topology%Grid, &
                                        Topology%Periodicity,              &
                                        Topology%Coordination,             &
                                        Processor_Topology_File,           &
                                        Topology%MPI_Communicator          )

  END Subroutine Setup_MPI_Topology_2D
!
! ========================================================================
!  
  Subroutine Setup_MPI_Topology_3D (                        &
    Index_x, Distributed_x, Periodicity_x, Topology_Grid_x, &
    Index_y, Distributed_y, Periodicity_y, Topology_Grid_y, &
    Index_z, Distributed_z, Periodicity_z, Topology_Grid_z, &
    Topology, Processor_Topology_File, Read_Topology_From_File)
 
  logical :: Distributed_x, Periodicity_x
  logical :: Distributed_y, Periodicity_y
  logical :: Distributed_z, Periodicity_z

  integer :: Index_x, Topology_Grid_x
  integer :: Index_y, Topology_Grid_y
  integer :: Index_z, Topology_Grid_z

  integer :: Dimension

  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid
  
  TYPE (Processor_Topology) :: Topology

  character(LEN=*), OPTIONAL :: Processor_Topology_File
  logical         , OPTIONAL :: Read_Topology_From_File
 
  Dimension = 3

                   Distributed(Index_x) =   Distributed_x
                   Periodicity(Index_x) =   Periodicity_x
                          Grid(Index_x) = Topology_Grid_x
 
                   Distributed(Index_y) =   Distributed_y
                   Periodicity(Index_y) =   Periodicity_y
                          Grid(Index_y) = Topology_Grid_y

                   Distributed(Index_z) =   Distributed_z
                   Periodicity(Index_z) =   Periodicity_z
                          Grid(Index_z) = Topology_Grid_z

  if (PRESENT(Processor_Topology_File) .AND.   &
      PRESENT(Read_Topology_From_File)) then
    if (Read_Topology_From_File) then
      call PS_MPI_Read_Processor_Topology (Dimension, Grid, Periodicity,     &
                                           Processor_Topology_File           )
    endif
  endif

  call PS_MPI_Processor_Topology_Setup (Topology, Dimension, &
                                        Distributed, Periodicity, Grid)

  Distributed_x = Distributed(Index_x)
  Distributed_y = Distributed(Index_y)
  Distributed_z = Distributed(Index_z)

#if defined (PARALLEL_MPI)
  call MPI_CART_GET    (Topology%MPI_Communicator, Topology%Dimension,     &
                        Topology%Grid            , Topology%Periodicity,   &
                        Topology%Coordination    , MPI_Error_Status)
#endif

  if (.NOT. PRESENT(Processor_Topology_File)) RETURN
                                        
  call PS_MPI_Write_Processor_Topology (Topology%Dimension, Topology%Grid, &
                                        Topology%Periodicity,              &
                                        Topology%Coordination,             &
                                        Processor_Topology_File,           &
                                        Topology%MPI_Communicator          )

  END Subroutine Setup_MPI_Topology_3D

END MODULE MPI_Processor_Topology






!-------------------------------
#if 0
!-------------------------------
!
!=======================================================================
!
  Subroutine MPI_Processor_Topology_Setup_0D (Topology, Periodicity)

  TYPE (Processor_Topology) :: Topology

  logical  :: Periodicity

  Topology%MPI_Communicator    = MPI_Comm_World
  Topology%Dimension           = 1
  Topology%Reorder             = .FALSE.
  Topology%Grid                = N_Processor
  Topology%Periodicity         = Periodicity

  call MPI_Grid_CREATE (N_Processor, Topology%Dimension,                  &
                                     Topology%Grid, MPI_Error_Status)

  call MPI_CART_CREATE (MPI_Comm_World, Topology%Dimension,               &
                        Topology%Grid,                                    &
                        Topology%Periodicity, Topology%Reorder,           &
                        Topology%MPI_Communicator, MPI_Error_Status)

  END Subroutine MPI_Processor_Topology_Setup_0D 
!
!=======================================================================
!
  Subroutine MPI_Processor_Topology_Setup_1D (Topology,                        &
                                        Index_x, Distributed_x, Periodicity_x  )

  TYPE (Processor_Topology) :: Topology

  integer  :: Index_x, Distributed_x, Periodicity_x

  Topology%MPI_Communicator    = MPI_Comm_World
  Topology%Dimension           = 1
  Topology%Reorder             = .FALSE.
  Topology%Periodicity         = .FALSE.

  if (Distributed_x) then
                            Topology%Periodicity(Index_x) = .FALSE.
    if (Periodicity_x == 1) Topology%Periodicity(Index_x) = .TRUE.

    Topology%Grid(Index_x) = 0
  else
    Topology%Grid(Index_x) = 1
  endif
    
  call MPI_Grid_CREATE (N_Processor, Topology%Dimension,                   &
                                     Topology%Grid, MPI_Error_Status)

  call MPI_CART_CREATE (MPI_Comm_World, Topology%Dimension,               &
                        Topology%Grid,                             &
                        Topology%Periodicity, Topology%Reorder,           &
                        Topology%MPI_Communicator, MPI_Error_Status)

  END Subroutine MPI_Processor_Topology_Setup_1D 
!
!=======================================================================
!
  Subroutine MPI_Processor_Topology_Setup_2D (Topology,                        &
                                        Index_x, Distributed_x, Periodicity_x, &
                                        Index_y, Distributed_y, Periodicity_y )

  TYPE (Processor_Topology) :: Topology

  integer  :: Index_x, Distributed_x, Periodicity_x
  integer  :: Index_y, Distributed_y, Periodicity_y

  Topology%MPI_Communicator    = MPI_Comm_World
  Topology%Dimension           = 2
  Topology%Reorder             = .FALSE.
  Topology%Periodicity         = .FALSE.

  if (Distributed_x) then
                            Topology%Periodicity(Index_x) = .FALSE.
    if (Periodicity_x == 1) Topology%Periodicity(Index_x) = .TRUE.

    Topology%Grid(Index_x) = 0
  else
    Topology%Grid(Index_x) = 1
  endif
    
  if (Distributed_y) then
                            Topology%Periodicity(Index_y) = .FALSE.
    if (Periodicity_y == 1) Topology%Periodicity(Index_y) = .TRUE.

    Topology%Grid(Index_y) = 0
  else
    Topology%Grid(Index_y) = 1
  endif
    
  call MPI_Grid_CREATE (N_Processor, Topology%Dimension,                   &
                                     Topology%Grid, MPI_Error_Status)

  call MPI_CART_CREATE (MPI_Comm_World, Topology%Dimension,               &
                        Topology%Grid,                             &
                        Topology%Periodicity, Topology%Reorder,           &
                        Topology%MPI_Communicator, MPI_Error_Status)

  END Subroutine MPI_Processor_Topology_Setup_2D 
!
!=======================================================================
!
  Subroutine MPI_Processor_Topology_Setup_3D (Topology,                        &
                                        Index_x, Distributed_x, Periodicity_x, &
                                        Index_y, Distributed_y, Periodicity_y, &
                                        Index_z, Distributed_z, Periodicity_z )

  TYPE (Processor_Topology) :: Topology

  integer  :: Index_x, Distributed_x, Periodicity_x
  integer  :: Index_y, Distributed_y, Periodicity_y
  integer  :: Index_z, Distributed_z, Periodicity_z 

  Topology%MPI_Communicator    = MPI_Comm_World
  Topology%Dimension           = 3
  Topology%Reorder             = .FALSE.
  Topology%Periodicity         = .FALSE.

  if (Distributed_x) then
                            Topology%Periodicity(Index_x) = .FALSE.
    if (Periodicity_x == 1) Topology%Periodicity(Index_x) = .TRUE.

    Topology%Grid(Index_x) = 0
  else
    Topology%Grid(Index_x) = 1
  endif
    
  if (Distributed_y) then
                            Topology%Periodicity(Index_y) = .FALSE.
    if (Periodicity_y == 1) Topology%Periodicity(Index_y) = .TRUE.

    Topology%Grid(Index_y) = 0
  else
    Topology%Grid(Index_y) = 1
  endif
    
  if (Distributed_z) then
                            Topology%Periodicity(Index_z) = .FALSE.
    if (Periodicity_z == 1) Topology%Periodicity(Index_z) = .TRUE.

    Topology%Grid(Index_z) = 0
  else
    Topology%Grid(Index_z) = 1
  endif
    
  call MPI_Grid_CREATE (N_Processor, Topology%Dimension,                   &
                                     Topology%Grid, MPI_Error_Status)

  call MPI_CART_CREATE (MPI_Comm_World, Topology%Dimension,               &
                        Topology%Grid,                             &
                        Topology%Periodicity, Topology%Reorder,           &
                        Topology%MPI_Communicator, MPI_Error_Status)

  END Subroutine MPI_Processor_Topology_Setup_3D 
!-------------------------------
#endif
!-------------------------------


!-------------------------------
#if 0
#include "pseudopack.h"


! ----------------------------------------------------------------------
! Module Name  : Processor
! Author       : Wai Sun Don
! Descritpion  : 
!  
! Subroutines  : 
! Module Used  : MPIF   (PARALLEL_MPI)
!
! ----------------------------------------------------------------------

MODULE MPI_Processor_Topology

  USE Processor
  USE PS_IO_Unit
  USE MPI_ShutDown

implicit NONE

TYPE Processor_Topology
  integer  :: MPI_Communicator

  integer  :: Dimension
  logical  :: Reorder

  integer , dimension(3) :: Grid, Coordination
  logical , dimension(3) :: Periodicity
END TYPE Processor_Topology

INTERFACE PS_MPI_Processor_Topology_Setup
  MODULE PROCEDURE MPI_Processor_Topology_Setup
END INTERFACE

INTERFACE PS_MPI_Topology_Setup
  MODULE PROCEDURE Setup_MPI_Topology_1D 
  MODULE PROCEDURE Setup_MPI_Topology_2D 
  MODULE PROCEDURE Setup_MPI_Topology_3D 
END INTERFACE

INTERFACE PS_MPI_Processor_Local_Info
  MODULE PROCEDURE MPI_Processor_Local_Info_0
  MODULE PROCEDURE MPI_Processor_Local_Info_1
  MODULE PROCEDURE MPI_Processor_Local_Info_2
END INTERFACE

INTERFACE PS_MPI_Read_Processor_Topology
  MODULE PROCEDURE MPI_Read_Processor_Topology_0
  MODULE PROCEDURE MPI_Read_Processor_Topology_1
END INTERFACE

INTERFACE PS_MPI_Write_Processor_Topology
  MODULE PROCEDURE MPI_Write_Processor_Topology_0
  MODULE PROCEDURE MPI_Write_Processor_Topology_1
END INTERFACE

PRIVATE

PUBLIC  :: Processor_Topology
PUBLIC  :: PS_MPI_Processor_Topology_Setup
PUBLIC  :: PS_MPI_Processor_Local_Info
PUBLIC  :: PS_MPI_Read_Processor_Topology
PUBLIC  :: PS_MPI_Write_Processor_Topology
PUBLIC  :: PS_MPI_Processor_Local_Index

CONTAINS

  Subroutine MPI_Processor_Topology_Setup (Topology, Dimension,             &
                                           Distributed, Periodicity, Grid,  &
                                           MPI_Communicator)

  TYPE (Processor_Topology) :: Topology

  integer                          :: Dimension, N_Dimension, n
  logical , dimension(:)           :: Distributed
  logical , dimension(:), OPTIONAL :: Periodicity
  integer , dimension(:), OPTIONAL :: Grid

  integer                          :: M_Processor

  integer ,               OPTIONAL :: MPI_Communicator
  integer                          :: MPI_Comm_Type

                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  N_Dimension = MIN(Dimension, 3)

  Topology%Dimension           = N_Dimension
  Topology%Reorder             = .FALSE.
  Topology%Periodicity         = .FALSE.

  do n = 1,N_Dimension
                                Topology%Periodicity(n) = .FALSE.
      if (PRESENT(Periodicity)) Topology%Periodicity(n) = Periodicity(n)

    if (Distributed(n)) then
                                Topology%Grid(n) = 0
      if (PRESENT(Grid))        Topology%Grid(n) = Grid(n)
    else
                                Topology%Grid(n) = 1
    endif
  enddo
    
#if defined (PARALLEL_MPI)
  Topology%MPI_Communicator    = MPI_Comm_Type

  call MPI_Comm_Size   (MPI_Comm_Type, M_Processor, MPI_Error_Status)

  call MPI_DIMS_CREATE (M_Processor, Topology%Dimension,                  &
                                     Topology%Grid, MPI_Error_Status)

  call MPI_CART_CREATE (MPI_Comm_Type, Topology%Dimension,                &
                        Topology%Grid,                                    &
                        Topology%Periodicity, Topology%Reorder,           &
                        Topology%MPI_Communicator, MPI_Error_Status)

#else

  Topology%MPI_Communicator    = 0
  Topology%Grid                = 1
  Topology%Coordination        = 1

#endif

  do n = 1,N_Dimension
    if (Topology%Grid(n) == 1) Distributed(n) = .FALSE.
  enddo

  END Subroutine MPI_Processor_Topology_Setup

!
!=======================================================================
!
  Subroutine MPI_Processor_Local_Info_0 (I_Am, Last_Processor,        &
                                         Index, MPI_Communicator)

  integer                             :: I_Am, Last_Processor
  integer , OPTIONAL                  :: Index, MPI_Communicator

  integer                             :: N_Dims, i
  integer , dimension(:), ALLOCATABLE :: Grid, COORDS
  logical , dimension(:), ALLOCATABLE :: PERIODS

#if defined (PARALLEL_MPI)
  if (PRESENT(MPI_Communicator)) then
    call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

    ALLOCATE (Grid(N_Dims), PERIODS(N_Dims), COORDS(N_Dims))

    call MPI_CART_GET    (MPI_Communicator, N_Dims, Grid, PERIODS, COORDS, &
                          MPI_Error_Status)

    i = 1 ; if (PRESENT(Index)) i = Index

    I_Am           =  COORDS(i)
    Last_Processor =    Grid(i)-1

    DEALLOCATE (Grid, PERIODS, COORDS)
  else
              I_Am = Processor_ID 
    Last_Processor = N_Processor-1
  endif
#else
              I_Am = Processor_ID 
    Last_Processor = N_Processor-1
#endif

  END Subroutine MPI_Processor_Local_Info_0
!
!=======================================================================
!
  Subroutine MPI_Processor_Local_Info_1 (Periodicity,                      &
                                         Index, MPI_Communicator)

  logical                             :: Periodicity
  integer , OPTIONAL                  :: Index, MPI_Communicator

  integer                             :: N_Dims, i
  integer , dimension(:), ALLOCATABLE :: Grid, COORDS
  logical , dimension(:), ALLOCATABLE :: PERIODS

#if defined (PARALLEL_MPI)
  if (PRESENT(MPI_Communicator)) then
    call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

    ALLOCATE (Grid(N_Dims), PERIODS(N_Dims), COORDS(N_Dims))

    call MPI_CART_GET    (MPI_Communicator, N_Dims, Grid, PERIODS, COORDS, &
                          MPI_Error_Status)

    i = 1 ; if (PRESENT(Index)) i = Index

    Periodicity    = PERIODS(i)

    DEALLOCATE (Grid, PERIODS, COORDS)
  else
    Periodicity         = .FALSE.
  endif
#else
    Periodicity         = .FALSE.
#endif

  END Subroutine MPI_Processor_Local_Info_1
!
!=======================================================================
!
  Subroutine MPI_Processor_Local_Info_2 (I_Am, Last_Processor, Periodicity, &
                                         Index, MPI_Communicator)

  integer                             :: I_Am, Last_Processor
  logical                             :: Periodicity
  integer , OPTIONAL                  :: Index, MPI_Communicator

  integer                             :: N_Dims, i
  integer , dimension(:), ALLOCATABLE :: Grid, COORDS
  logical , dimension(:), ALLOCATABLE :: PERIODS

#if defined (PARALLEL_MPI)
  if (PRESENT(MPI_Communicator)) then
    call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

    ALLOCATE (Grid(N_Dims), PERIODS(N_Dims), COORDS(N_Dims))

    call MPI_CART_GET    (MPI_Communicator, N_Dims, Grid, PERIODS, COORDS, &
                          MPI_Error_Status)

    i = 1 ; if (PRESENT(Index)) i = Index

    I_Am           =  COORDS(i)
    Last_Processor =    Grid(i)-1
    Periodicity    = PERIODS(i)

    DEALLOCATE (Grid, PERIODS, COORDS)
  else
              I_Am = Processor_ID 
    Last_Processor = N_Processor-1
    Periodicity    = .FALSE.
  endif
#else
              I_Am = Processor_ID 
    Last_Processor = N_Processor-1
    Periodicity    = .FALSE.
#endif

  END Subroutine MPI_Processor_Local_Info_2
!
!=======================================================================
!
  Subroutine PS_MPI_Processor_Local_Index (Index, I_Am_i, First_Processor_i, &
                                           Last_Processor_i, MPI_Comm_Type)

  integer            :: Index
  integer            :: I_Am_i, First_Processor_i, Last_Processor_i
  integer , OPTIONAL :: MPI_Comm_Type

#if defined (PARALLEL_MPI)
  integer                             :: MPI_Communicator_1D
  integer                             :: N_Dims
  logical , dimension(:), ALLOCATABLE :: Sub_Grid

             I_Am_i =     I_Am
  First_Processor_i = First_Processor
   Last_Processor_i =  Last_Processor

  if (.NOT. PRESENT(MPI_Comm_Type)) RETURN

  call MPI_CARTDIM_GET (MPI_Comm_Type, N_Dims, MPI_Error_Status)

  ALLOCATE (Sub_Grid(N_Dims))

  Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

  call MPI_CART_SUB  (MPI_Comm_Type, Sub_Grid, MPI_Communicator_1D,  &
                                                  MPI_Error_Status)

  call PS_MPI_Processor_Local_Info (I_Am_i, Last_Processor_i, Index,    &
                                    MPI_Communicator_1D)

  DEALLOCATE (Sub_Grid)

  call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
#else
             I_Am_i =     I_Am
  First_Processor_i = First_Processor
   Last_Processor_i =  Last_Processor
#endif

  END Subroutine PS_MPI_Processor_Local_Index 

!#if defined (PARALLEL_MPI)
!
!=======================================================================
!
  Subroutine MPI_Read_Processor_Topology_0 (Dimension, Grid,              &
                                            Periodicity, lid, Error)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  integer                          :: lid
  logical , OPTIONAL               :: Error
  logical                          :: Error_p, Error_l
  integer                          :: M_Processor, Dims

  Error_p = .FALSE.

  read (lid,100) M_Processor
  if (M_Processor /= N_Processor) then
    write (lid6 ,200) N_Processor, M_Processor
    write (lid99,200) N_Processor, M_Processor
    Error_p = .TRUE.
  endif

  read (lid,100) Dims
  if (Dims /= Dimension) then
    write (lid6 ,201) Dimension, Dims
    write (lid99,201) Dimension, Dims
    Error_p = .TRUE.
  endif

#if defined (PARALLEL_MPI)
  call MPI_ALLREDUCE (Error_p, Error_l, 1, MPI_LOGICAL, MPI_LOR,       &
                               MPI_Comm_World, MPI_Error_Status)
#else
  Error_l = Error_p
#endif

  if (Error_l) then
    if (PRESENT(Error)) then
      CLOSE (lid) ; Error = Error_l ; RETURN 
    else
      call PS_MPI_ShutDown
    endif
  endif

  read (lid,101)        Grid(1:Dimension)
  read (lid,102) Periodicity(1:Dimension)

  CLOSE (lid)

 100 format (23x,i10)
 101 format (23x,3(i10,:))
 102 format (23x,3(l10,:))

 200 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Number of Processor Specified : ',i10        / &
             1x,'  Number of Processor Read in   : ',i10        / &
             1x,72(1x,'*')/)
 201 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Topology Dimension Specified  : ',i10        / &
             1x,'  Topology Dimension Read in    : ',i10        / &
             1x,72(1x,'*')/)

  END Subroutine MPI_Read_Processor_Topology_0
!
!=======================================================================
!
  Subroutine MPI_Write_Processor_Topology_0 (Dimension, Grid,              &
                                             Periodicity, Coordination, lid)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  integer , dimension(1:Dimension) :: Coordination
  integer                          :: lid

  integer , dimension(0:Dimension) :: A 

#if defined (PARALLEL_MPI)
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)
#endif

  A(0) = I_Am ; A(1:Dimension) = Coordination(1:Dimension) 

  if (I_Am  == First_Processor) then
    write (lid,100) N_Processor, Dimension
    write (lid,101) Grid 
    write (lid,102) Periodicity
  endif

  if (I_Am == First_Processor) write (lid,103) A 

#if defined (PARALLEL_MPI)
  do K_Processor = First_Processor+1,M_Processor-1

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

    if (I_Am ==     K_Processor)   &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,  &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,  &
                                 MPI_Comm_World, Status, MPI_Error_Status)

      write (lid,103) A
    endif
  enddo

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

100 format (1x,'Number of Processor : ', i10      / &
            1x,'Cartesian Dimension : ', i10      )
101 format (1x,'Cartesian Grid Size : ', 3(i10,:) )
102 format (1x,'Periodicity         : ', 3(l10,:) )
103 format (1x,'Rank/Coorindate     : ', i10, 3(i10,:))

  END Subroutine MPI_Write_Processor_Topology_0
!
!=======================================================================
!
  Subroutine MPI_Read_Processor_Topology_1 (Dimension, Grid,              &
                                            Periodicity, Filename, Error)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  character(LEN=*)                 :: Filename
  integer                          :: lid, IOS, IOS_p
  integer                          :: M_Processor, Dims
  logical                          :: Error_p, Error_l
  logical , OPTIONAL               :: Error

  IOS_p = 0 ; Error_p = .FALSE.

    lid = 77

    OPEN (UNIT=lid, FILE=TRIM(Filename), STATUS='OLD'    , POSITION='REWIND',  &
                    IOSTAT=IOS_p)

    if (IOS_p /= 0) then
      write (lid6 ,202) I_Am, lid, IOS_p, TRIM(Filename)
      write (lid99,202) I_Am, lid, IOS_p, TRIM(Filename)

      Error_p = .TRUE.
    endif

#if defined (PARALLEL_MPI)
  call MPI_ALLREDUCE (Error_p, Error_l, 1, MPI_LOGICAL, MPI_LOR,       &
                               MPI_Comm_World, MPI_Error_Status)
#else
  Error_l = Error_p
#endif

  if (Error_l) call PS_MPI_ShutDown

  Error_p = .FALSE.

  read (lid,100) M_Processor
  if (M_Processor /= N_Processor) then
    write (lid6 ,200) N_Processor, M_Processor
    write (lid99,200) N_Processor, M_Processor
    Error_p = .TRUE.
  endif

  read (lid,100) Dims
  if (Dims /= Dimension) then
    write (lid6 ,201) Dimension, Dims
    write (lid99,201) Dimension, Dims
    Error_p = .TRUE.
  endif

#if defined (PARALLEL_MPI)
  call MPI_ALLREDUCE (Error_p, Error_l, 1, MPI_LOGICAL, MPI_LOR,       &
                               MPI_Comm_World, MPI_Error_Status)
#else
  Error_l = Error_p
#endif

  if (Error_l) then
    if (PRESENT(Error)) then
      CLOSE (lid) ; Error = Error_l ; RETURN 
    else
      call PS_MPI_ShutDown
    endif
  endif

  read (lid,101)        Grid(1:Dimension)
  read (lid,102) Periodicity(1:Dimension)

  CLOSE (lid)

 100 format (23x,i10)
 101 format (23x,3(i10,:))
 102 format (23x,3(l10,:))

 200 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Number of Processor Specified : ',i10        / &
             1x,'  Number of Processor Read in   : ',i10        / &
             1x,72(1x,'*')/)
 201 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Topology Dimension Specified  : ',i10        / &
             1x,'  Topology Dimension Read in    : ',i10        / &
             1x,72(1x,'*')/)
 202 format (1x,'FATAL ERROR : PS_MPI_Read_Processor_Topology ' / &
             1x,'  Processor ID                  : ',i10        / &
             1x,'  File Unit, Status, Filename   : ',i5,i5,a32  )

  END Subroutine MPI_Read_Processor_Topology_1
!
!=======================================================================
!
  Subroutine MPI_Write_Processor_Topology_1 (Dimension, Grid,              &
                                             Periodicity, Coordination,    &
                                             Filename)

  integer                          :: Dimension
  integer , dimension(1:Dimension) :: Grid
  logical , dimension(1:Dimension) :: Periodicity
  integer , dimension(1:Dimension) :: Coordination

  character(LEN=*)                 :: Filename
  integer                          :: lid, IOS, IOS_p
  logical                          :: Error_p, Error_l

  integer , dimension(0:Dimension) :: A 

#if defined (PARALLEL_MPI)
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)
#endif

  A(0) = I_Am ; A(1:Dimension) = Coordination(1:Dimension) 

  IOS_p = 0 ; Error_p = .FALSE.
  if (I_Am  == First_Processor) then
    lid = 78

    OPEN (UNIT=lid, FILE=TRIM(Filename), STATUS='UNKNOWN', POSITION='REWIND',  &
                    IOSTAT=IOS_p)

    if (IOS_p /= 0) then
      write (lid6 ,202) I_Am, lid, IOS_p, TRIM(Filename)
      write (lid99,202) I_Am, lid, IOS_p, TRIM(Filename)

      Error_p = .TRUE.
    endif
  endif

#if defined (PARALLEL_MPI)
  call MPI_ALLREDUCE (Error_p, Error_l, 1, MPI_LOGICAL, MPI_LOR,       &
                               MPI_Comm_World, MPI_Error_Status)
#else
  Error_l = Error_p
#endif

  if (Error_l) call PS_MPI_ShutDown

  if (I_Am  == First_Processor) then
    write (lid,100) N_Processor, Dimension
    write (lid,101) Grid 
    write (lid,102) Periodicity
  endif

  if (I_Am == First_Processor) write (lid,103) A 

#if defined (PARALLEL_MPI)
  do K_Processor = First_Processor+1,M_Processor-1

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

    if (I_Am ==     K_Processor)   &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,  &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,  &
                                 MPI_Comm_World, Status, MPI_Error_Status)


      write (lid,103) A
    endif
  enddo

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  if (I_Am == First_Processor) CLOSE (lid)

 100 format (1x,'Number of Processor : ', i10      / &
             1x,'Cartesian Dimension : ', i10      )
 101 format (1x,'Cartesian Grid Size : ', 3(i10,:) )
 102 format (1x,'Periodicity         : ', 3(l10,:) )
 103 format (1x,'Rank/Coorindate     : ', i10, 3(i10,:))

 202 format (1x,'FATAL ERROR : PS_MPI_Write_Processor_Topology '/ &
             1x,'  Processor ID                  : ',i10        / &
             1x,'  File Unit, Status, Filename   : ',i5,i5,a32  )

  END Subroutine MPI_Write_Processor_Topology_1

!
! ========================================================================
!  
  Subroutine Setup_MPI_Topology_1D (                        &
    Index_x, Distributed_x, Periodicity_x, Topology_Grid_x, &
    Topology, Processor_Topology_File, Read_Topology_From_File)
 
  logical :: Distributed_x, Periodicity_x
  integer :: Index_x, Topology_Grid_x

  integer :: Dimension

  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid
  
  TYPE (Processor_Topology) :: Topology

  character(LEN=*), OPTIONAL :: Processor_Topology_File
  logical         , OPTIONAL :: Read_Topology_From_File
 
  Dimension = 1

                   Distributed(Index_x) =   Distributed_x
                   Periodicity(Index_x) =   Periodicity_x
                          Grid(Index_x) = Topology_Grid_x
 
                   Distributed(2      ) = .FALSE.
                   Periodicity(2      ) = .FALSE.
                          Grid(2      ) = 0

                   Distributed(3      ) = .FALSE.
                   Periodicity(3      ) = .FALSE.
                          Grid(3      ) = 0

  if (PRESENT(Processor_Topology_File) .AND.   &
      PRESENT(Read_Topology_From_File)) then
    if (Read_Topology_From_File) then
      call PS_MPI_Read_Processor_Topology (Dimension, Grid, Periodicity,     &
                                           Processor_Topology_File           )
    endif
  endif

  call PS_MPI_Processor_Topology_Setup (Topology, Dimension, &
                                        Distributed, Periodicity, Grid)

  Distributed_x = Distributed(Index_x)

#if defined (PARALLEL_MPI)
  call MPI_CART_GET    (Topology%MPI_Communicator, Topology%Dimension,     &
                        Topology%Grid            , Topology%Periodicity,   &
                        Topology%Coordination    , MPI_Error_Status)
#endif

  if (.NOT. PRESENT(Processor_Topology_File)) RETURN
                                        
  call PS_MPI_Write_Processor_Topology (Topology%Dimension, Topology%Grid, &
                                        Topology%Periodicity,              &
                                        Topology%Coordination,             &
                                        Processor_Topology_File,           &
                                        Topology%MPI_Communicator          )

  END Subroutine Setup_MPI_Topology_1D
!
! ========================================================================
!  
  Subroutine Setup_MPI_Topology_2D (                        &
    Index_x, Distributed_x, Periodicity_x, Topology_Grid_x, &
    Index_y, Distributed_y, Periodicity_y, Topology_Grid_y, &
    Topology, Processor_Topology_File, Read_Topology_From_File)
 
  logical :: Distributed_x, Periodicity_x
  logical :: Distributed_y, Periodicity_y

  integer :: Index_x, Topology_Grid_x
  integer :: Index_y, Topology_Grid_y

  integer :: Dimension

  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid
  
  TYPE (Processor_Topology) :: Topology

  character(LEN=*), OPTIONAL :: Processor_Topology_File
  logical         , OPTIONAL :: Read_Topology_From_File
 
  Dimension = 2

                   Distributed(Index_x) =   Distributed_x
                   Periodicity(Index_x) =   Periodicity_x
                          Grid(Index_x) = Topology_Grid_x
 
                   Distributed(Index_y) =   Distributed_y
                   Periodicity(Index_y) =   Periodicity_y
                          Grid(Index_y) = Topology_Grid_y

                   Distributed(3      ) = .FALSE.
                   Periodicity(3      ) = .FALSE.
                          Grid(3      ) = 0

  if (PRESENT(Processor_Topology_File) .AND.   &
      PRESENT(Read_Topology_From_File)) then
    if (Read_Topology_From_File) then
      call PS_MPI_Read_Processor_Topology (Dimension, Grid, Periodicity,     &
                                           Processor_Topology_File           )
    endif
  endif

  call PS_MPI_Processor_Topology_Setup (Topology, Dimension, &
                                        Distributed, Periodicity, Grid)

  Distributed_x = Distributed(Index_x)
  Distributed_y = Distributed(Index_y)

#if defined (PARALLEL_MPI)
  call MPI_CART_GET    (Topology%MPI_Communicator, Topology%Dimension,     &
                        Topology%Grid            , Topology%Periodicity,   &
                        Topology%Coordination    , MPI_Error_Status)
#endif

  if (.NOT. PRESENT(Processor_Topology_File)) RETURN
                                        
  call PS_MPI_Write_Processor_Topology (Topology%Dimension, Topology%Grid, &
                                        Topology%Periodicity,              &
                                        Topology%Coordination,             &
                                        Processor_Topology_File,           &
                                        Topology%MPI_Communicator          )

  END Subroutine Setup_MPI_Topology_2D
!
! ========================================================================
!  
  Subroutine Setup_MPI_Topology_3D (                        &
    Index_x, Distributed_x, Periodicity_x, Topology_Grid_x, &
    Index_y, Distributed_y, Periodicity_y, Topology_Grid_y, &
    Index_z, Distributed_z, Periodicity_z, Topology_Grid_z, &
    Topology, Processor_Topology_File, Read_Topology_From_File)
 
  logical :: Distributed_x, Periodicity_x
  logical :: Distributed_y, Periodicity_y
  logical :: Distributed_z, Periodicity_z

  integer :: Index_x, Topology_Grid_x
  integer :: Index_y, Topology_Grid_y
  integer :: Index_z, Topology_Grid_z

  integer :: Dimension

  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid
  
  TYPE (Processor_Topology) :: Topology

  character(LEN=*), OPTIONAL :: Processor_Topology_File
  logical         , OPTIONAL :: Read_Topology_From_File
 
  Dimension = 3

                   Distributed(Index_x) =   Distributed_x
                   Periodicity(Index_x) =   Periodicity_x
                          Grid(Index_x) = Topology_Grid_x
 
                   Distributed(Index_y) =   Distributed_y
                   Periodicity(Index_y) =   Periodicity_y
                          Grid(Index_y) = Topology_Grid_y

                   Distributed(Index_z) =   Distributed_z
                   Periodicity(Index_z) =   Periodicity_z
                          Grid(Index_z) = Topology_Grid_z

  if (PRESENT(Processor_Topology_File) .AND.   &
      PRESENT(Read_Topology_From_File)) then
    if (Read_Topology_From_File) then
      call PS_MPI_Read_Processor_Topology (Dimension, Grid, Periodicity,     &
                                           Processor_Topology_File           )
    endif
  endif

  call PS_MPI_Processor_Topology_Setup (Topology, Dimension, &
                                        Distributed, Periodicity, Grid)

  Distributed_x = Distributed(Index_x)
  Distributed_y = Distributed(Index_y)
  Distributed_z = Distributed(Index_z)

#if defined (PARALLEL_MPI)
  call MPI_CART_GET    (Topology%MPI_Communicator, Topology%Dimension,     &
                        Topology%Grid            , Topology%Periodicity,   &
                        Topology%Coordination    , MPI_Error_Status)
#endif

  if (.NOT. PRESENT(Processor_Topology_File)) RETURN
                                        
  call PS_MPI_Write_Processor_Topology (Topology%Dimension, Topology%Grid, &
                                        Topology%Periodicity,              &
                                        Topology%Coordination,             &
                                        Processor_Topology_File,           &
                                        Topology%MPI_Communicator          )

  END Subroutine Setup_MPI_Topology_3D

END MODULE MPI_Processor_Topology

#endif
!-------------------------------
