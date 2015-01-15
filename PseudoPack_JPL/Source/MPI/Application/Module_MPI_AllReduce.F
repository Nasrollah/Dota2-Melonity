#include "pseudopack.h"


MODULE MPI_AllReduce_Operation

  USE Processor
  USE PS_IO_Unit, ONLY: lid6, lid99

implicit NONE

INTERFACE PS_MPI_AllReduce
  MODULE PROCEDURE AllReduce_0D
  MODULE PROCEDURE AllReduce_1D
  MODULE PROCEDURE AllReduce_2D
  MODULE PROCEDURE AllReduce_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_MPI_AllReduce

CONTAINS
!
!=======================================================================
!
  Subroutine AllReduce_0D (Operation, Index, N, f_p, f, MPI_Comm_Type)

  character(LEN=32)                   ::     Operation
  integer                             :: MPI_Operation

  integer                             :: Index, N
  REALTYPE, dimension(*)              :: f_p, f

  integer ,                  OPTIONAL :: MPI_Comm_Type

#if defined (PARALLEL_MPI)
  integer                             :: MPI_Communicator
  integer                             :: MPI_Communicator_0D
  integer                             :: N_Dims
  logical , dimension(:), ALLOCATABLE :: Sub_Grid

  SELECT CASE (TRIM(Operation))
    CASE ('MAX')
      MPI_Operation = MPI_MAX

    CASE ('MIN')
      MPI_Operation = MPI_MIN

    CASE ('SUM')
      MPI_Operation = MPI_SUM

    CASE ('PROD')
      MPI_Operation = MPI_PROD

    CASE ('MAXLOC')
      MPI_Operation = MPI_MAXLOC

    CASE ('MINLOC')
      MPI_Operation = MPI_MINLOC

    CASE DEFAULT
      write (lid6 ,*) 'STOP! MPI_AllReduce : Unknown Operation'
      write (lid99,*) 'STOP! MPI_AllReduce : Unknown Operation'
      STOP

  END SELECT

  if (N_Processor == 1) then
    f(1:N) = f_p(1:N)
  else
    if (PRESENT(MPI_Comm_Type)) then
      MPI_Communicator = MPI_Comm_Type

      call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

      ALLOCATE (Sub_Grid(N_Dims))

      Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

      call MPI_CART_SUB  (MPI_Communicator, Sub_Grid, MPI_Communicator_0D,  &
                                                      MPI_Error_Status)

      call MPI_ALLREDUCE (f_p, f, N,                                        &
                          MPI_REALTYPE, MPI_Operation, MPI_Communicator_0D, &
                          MPI_Error_Status)

      DEALLOCATE (Sub_Grid)

      call MPI_COMM_FREE (MPI_Communicator_0D, MPI_Error_Status)
    else
      MPI_Communicator = MPI_Comm_World

      call MPI_ALLREDUCE (f_p, f, N,                                        &
                          MPI_REALTYPE, MPI_Operation, MPI_Communicator,    &
                          MPI_Error_Status)
    endif
  endif
#else
    f(1:N) = f_p(1:N)
#endif

  END Subroutine AllReduce_0D
!
!=======================================================================
!
#if 0
  Subroutine AllReduce_1D (Operation, Index, N, f_p, f, MPI_Comm_Type)
#else
  Subroutine AllReduce_1D (Operation, Index, f_p, f, MPI_Comm_Type)
#endif

  character(LEN=32)                   ::     Operation
  integer                             :: MPI_Operation

  integer                             :: Index, N
  REALTYPE, dimension(:)              :: f_p, f

  integer ,                  OPTIONAL :: MPI_Comm_Type

#if defined (PARALLEL_MPI)
  integer                             :: MPI_Communicator
  integer                             :: MPI_Communicator_1D
  integer                             :: N_Dims
  logical , dimension(:), ALLOCATABLE :: Sub_Grid

  SELECT CASE (TRIM(Operation))
    CASE ('MAX')
      MPI_Operation = MPI_MAX

    CASE ('MIN')
      MPI_Operation = MPI_MIN

    CASE ('SUM')
      MPI_Operation = MPI_SUM

    CASE ('PROD')
      MPI_Operation = MPI_PROD

    CASE ('MAXLOC')
      MPI_Operation = MPI_MAXLOC

    CASE ('MINLOC')
      MPI_Operation = MPI_MINLOC

    CASE DEFAULT
      write (lid6 ,*) 'STOP! MPI_AllReduce : Unknown Operation'
      write (lid99,*) 'STOP! MPI_AllReduce : Unknown Operation'
      STOP

  END SELECT

  if (N_Processor == 1) then
    f = f_p
  else
    if (PRESENT(MPI_Comm_Type)) then
      MPI_Communicator = MPI_Comm_Type

      call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

      ALLOCATE (Sub_Grid(N_Dims))

      Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

      call MPI_CART_SUB  (MPI_Communicator, Sub_Grid, MPI_Communicator_1D,  &
                                                      MPI_Error_Status)

      call MPI_ALLREDUCE (f_p, f, SIZE(f),                                  &
                          MPI_REALTYPE, MPI_Operation, MPI_Communicator_1D, &
                          MPI_Error_Status)

      DEALLOCATE (Sub_Grid)

      call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
    else
      MPI_Communicator = MPI_Comm_World

      call MPI_ALLREDUCE (f_p, f, SIZE(f),                                  &
                          MPI_REALTYPE, MPI_Operation, MPI_Communicator,    &
                          MPI_Error_Status)
    endif
  endif
#else
    f = f_p
#endif

  END Subroutine AllReduce_1D
!
!=======================================================================
!
#if 0
  Subroutine AllReduce_2D (Operation, Index, N, f_p, f, MPI_Comm_Type)
#else
  Subroutine AllReduce_2D (Operation, Index, f_p, f, MPI_Comm_Type)
#endif

  character(LEN=32)                   ::     Operation
  integer                             :: MPI_Operation

  integer                             :: Index, N
  REALTYPE, dimension(:,:)            :: f_p, f

  integer ,                  OPTIONAL :: MPI_Comm_Type

#if defined (PARALLEL_MPI)
  integer                             :: MPI_Communicator
  integer                             :: MPI_Communicator_1D
  integer                             :: N_Dims
  logical , dimension(:), ALLOCATABLE :: Sub_Grid

  SELECT CASE (TRIM(Operation))
    CASE ('MAX')
      MPI_Operation = MPI_MAX

    CASE ('MIN')
      MPI_Operation = MPI_MIN

    CASE ('SUM')
      MPI_Operation = MPI_SUM

    CASE ('PROD')
      MPI_Operation = MPI_PROD

    CASE ('MAXLOC')
      MPI_Operation = MPI_MAXLOC

    CASE ('MINLOC')
      MPI_Operation = MPI_MINLOC

    CASE DEFAULT
      write (lid6 ,*) 'STOP! MPI_AllReduce : Unknown Operation'
      write (lid99,*) 'STOP! MPI_AllReduce : Unknown Operation'
      STOP

  END SELECT

  if (N_Processor == 1) then
    f = f_p
  else
    if (PRESENT(MPI_Comm_Type)) then
      MPI_Communicator = MPI_Comm_Type

      call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

      ALLOCATE (Sub_Grid(N_Dims))

      Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

      call MPI_CART_SUB  (MPI_Communicator, Sub_Grid, MPI_Communicator_1D,  &
                                                      MPI_Error_Status)

      call MPI_ALLREDUCE (f_p, f, SIZE(f),                                  &
                          MPI_REALTYPE, MPI_Operation, MPI_Communicator_1D, &
                          MPI_Error_Status)

      DEALLOCATE (Sub_Grid)

      call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
    else
      MPI_Communicator = MPI_Comm_World

      call MPI_ALLREDUCE (f_p, f, SIZE(f),                                  &
                          MPI_REALTYPE, MPI_Operation, MPI_Communicator,    &
                          MPI_Error_Status)
    endif
  endif
#else
    f = f_p
#endif

  END Subroutine AllReduce_2D
!
!=======================================================================
!
#if 0
  Subroutine AllReduce_3D (Operation, Index, N, f_p, f, MPI_Comm_Type)
#else
  Subroutine AllReduce_3D (Operation, Index, f_p, f, MPI_Comm_Type)
#endif

  character(LEN=32)                   ::     Operation
  integer                             :: MPI_Operation

  integer                             :: Index, N
  REALTYPE, dimension(:,:,:)          :: f_p, f

  integer ,                  OPTIONAL :: MPI_Comm_Type

#if defined (PARALLEL_MPI)
  integer                             :: MPI_Communicator
  integer                             :: MPI_Communicator_1D
  integer                             :: N_Dims
  logical , dimension(:), ALLOCATABLE :: Sub_Grid

  SELECT CASE (TRIM(Operation))
    CASE ('MAX')
      MPI_Operation = MPI_MAX

    CASE ('MIN')
      MPI_Operation = MPI_MIN

    CASE ('SUM')
      MPI_Operation = MPI_SUM

    CASE ('PROD')
      MPI_Operation = MPI_PROD

    CASE ('MAXLOC')
      MPI_Operation = MPI_MAXLOC

    CASE ('MINLOC')
      MPI_Operation = MPI_MINLOC

    CASE DEFAULT
      write (lid6 ,*) 'STOP! MPI_AllReduce : Unknown Operation'
      write (lid99,*) 'STOP! MPI_AllReduce : Unknown Operation'
      STOP

  END SELECT

  if (N_Processor == 1) then
    f = f_p
  else
    if (PRESENT(MPI_Comm_Type)) then
      MPI_Communicator = MPI_Comm_Type

      call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

      ALLOCATE (Sub_Grid(N_Dims))

      Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

      call MPI_CART_SUB  (MPI_Communicator, Sub_Grid, MPI_Communicator_1D,  &
                                                      MPI_Error_Status)

      call MPI_ALLREDUCE (f_p, f, SIZE(f),                                  &
                          MPI_REALTYPE, MPI_Operation, MPI_Communicator_1D, &
                          MPI_Error_Status)

      DEALLOCATE (Sub_Grid)

      call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
    else
      MPI_Communicator = MPI_Comm_World

      call MPI_ALLREDUCE (f_p, f, SIZE(f),                                  &
                          MPI_REALTYPE, MPI_Operation, MPI_Communicator,    &
                          MPI_Error_Status)
    endif
  endif
#else
    f = f_p
#endif

  END Subroutine AllReduce_3D

END MODULE MPI_AllReduce_Operation
