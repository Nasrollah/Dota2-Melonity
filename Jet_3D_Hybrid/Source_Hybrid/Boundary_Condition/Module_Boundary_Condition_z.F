#include "pseudopack.h"


MODULE Boundary_Condition_z

  USE PseudoPack

IMPLICIT NONE

INTERFACE PS_Boundary_Condition_z
  MODULE PROCEDURE BC_z
END INTERFACE 

PRIVATE
PUBLIC  :: PS_Boundary_Condition_z

CONTAINS
!
! =========================================================================
!
SUBROUTINE BC_z          (Direction,                                  &
                          N0, N5, N2, N3,                             &
                          M0, M5, M2, M3,                             &
                          K0, K5, K2, K3,  NV, Q,                     &
                          MPI_Comm_Type)

integer  :: Direction, Index
integer  :: N0, N5, N2, N3
integer  :: M0, M5, M2, M3
integer  :: K0, K5, K2, K3,  NV

REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q

integer ,                  OPTIONAL :: MPI_Comm_Type
#if defined (PARALLEL_MPI)
integer                             :: MPI_Communicator
integer                             :: MPI_Communicator_1D
integer                             :: N_Dims
logical , dimension(:), ALLOCATABLE :: Sub_Grid

                            MPI_Communicator = MPI_Comm_World
if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

ALLOCATE (Sub_Grid(N_Dims))
#endif

Index = 3

#if defined (PARALLEL_MPI)
  call PS_MPI_Processor_Local_Info (I_Am, Last_Processor, &
                                          Index, MPI_Communicator)

  Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

  call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                 MPI_Error_Status)
#endif


      if (Index == 3) call BC

#if defined (PARALLEL_MPI)
  call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)

DEALLOCATE (Sub_Grid)

call PS_MPI_Processor (I_Am, Last_Processor)
#endif

CONTAINS

  SUBROUTINE BC



  SELECT CASE (Direction)
    CASE (-1)
#if defined (PARALLEL_MPI)
      if (I_Am == First_Processor) then
#endif

!--------------------------------------------------
           ! DANGER
        ! KAUSHIK COMMENTED THE FOLLOWING LINE
!#include "z_Boundary_Bottom.i"
!--------------------------------------------------

#if defined (PARALLEL_MPI)
       endif
#endif

   CASE DEFAULT

   CASE (1)
#if defined (PARALLEL_MPI)
      if (I_Am ==  Last_Processor) then
#endif

!--------------------------------------------------
           ! DANGER
        ! KAUSHIK COMMENTED THE FOLLOWING LINE
!#include "z_Boundary_Top.i"
!--------------------------------------------------

#if defined (PARALLEL_MPI)
       endif
#endif

  END SELECT

  END SUBROUTINE BC

END SUBROUTINE BC_z

END MODULE Boundary_Condition_z
