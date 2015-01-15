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

MODULE MPI_Processor_Topology_Setup

#if 0
  USE MPIF
#else
  USE Processor
#endif
  USE MPI_Processor_Topology
  USE PseudoPack_Definition

implicit NONE

PRIVATE

PUBLIC  :: PS_MPI_Topology_Setup_1D
PUBLIC  :: PS_MPI_Topology_Setup_2D
PUBLIC  :: PS_MPI_Topology_Setup_3D

CONTAINS

!
! ========================================================================
!  
  Subroutine PS_MPI_Topology_Setup_1D (Topolgy, Processor_Topology_File)
 
  TYPE (Processor_Toplogy) :: Toplogy

  character(LEN=*), OPTIONAL :: Processor_Topology_File
 
  integer :: Dimension

  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid
  
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
                                        Processor_Topology_File            )

  END Subroutine PS_MPI_Topology_Setup_1D
!
! ========================================================================
!  
  Subroutine PS_MPI_Topology_Setup_2D (Topolgy, Processor_Topology_File)
 
  TYPE (Processor_Toplogy) :: Toplogy
 
  character(LEN=*), OPTIONAL :: Processor_Topology_File
 
  integer :: Dimension

  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid
  
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
                                        Processor_Topology_File            )

  END Subroutine PS_MPI_Topology_Setup_2D
!
! ========================================================================
!  
  Subroutine PS_MPI_Topology_Setup_2D (Topolgy, Processor_Topology_File)
 
  TYPE (Processor_Toplogy) :: Toplogy
 
  character(LEN=*), OPTIONAL :: Processor_Topology_File
 
  integer :: Dimension

  logical, dimension(3) :: Distributed, Periodicity
  integer, dimension(3) :: Grid
  
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
                                        Processor_Topology_File            )

  END Subroutine PS_MPI_Topology_Setup_3D

END MODULE MPI_Processor_Topology_Setup
