#include "pseudopack.h"
#include "grid.h"

MODULE MPI_Spectral_Grid_Indexing

  USE PseudoPack_Utility
  USE Processor
  USE MPI_Processor_Topology
  USE MPI_Data_Distribution

implicit NONE

INTERFACE PS_MPI_Spectral_Grid_Index
  MODULE PROCEDURE MPI_Spectral_Grid_Index_1D
  MODULE PROCEDURE MPI_Spectral_Grid_Index_2D
  MODULE PROCEDURE MPI_Spectral_Grid_Index_3D
END INTERFACE

INTERFACE PS_Spectral_Grid_Index
  MODULE PROCEDURE     Spectral_Grid_Index_1D
  MODULE PROCEDURE     Spectral_Grid_Index_2D
  MODULE PROCEDURE     Spectral_Grid_Index_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_MPI_Spectral_Grid_Index
PUBLIC  :: PS_Spectral_Grid_Index

CONTAINS
!
! ============================================================================
!
  Subroutine MPI_Spectral_Grid_Index (_PS_SPECTRAL_GRID_X_1, &
                                      _PS_SPECTRAL_GRID_X_2, &
                                      _PS_SPECTRAL_GRID_X_3, &
                                      _PS_SPECTRAL_GRID_X_4, &
                                      MPI_Communicator       )

  logical  :: Distributed_x
  integer  :: _PS_SPECTRAL_GRID_X_0
  integer  :: _PS_SPECTRAL_GRID_X_2
  integer  :: _PS_SPECTRAL_GRID_X_3
  integer  :: _PS_SPECTRAL_GRID_X_4

  integer , OPTIONAL :: MPI_Communicator

!Set Global Index
  call Spectral_Grid_Index (_PS_SPECTRAL_GRID_X_0, &
                            _PS_SPECTRAL_GRID_X_2  )
  
  if (Distributed_x) then
    call PS_MPI_Data_Distribution (N3_g-N2_g+1, _PS_SPECTRAL_GRID_X_4, &
                                   Index_x, MPI_Communicator)
  
    NN = n_Star
    N0 = 0 ; N1 = 0 ; N2 = 0 ; N3 = n_Star ; N4 = N3 ; N5 = N3
  else
    n_Remain = 0 ; n_0 = 0 ; n_Star = N3_g-N2_g+1

    N  = N_g  
    N0 = N0_g ; N1 = N1_g ; N2 = N2_g ; N3 = N3_g ; N4 = N4_g ; N5 = N5_g
  endif

  END Subroutine MPI_Spectral_Grid_Index
!
! ============================================================================
!
  Subroutine MPI_SPECTRAL_Grid_Index_1D (N_x, _PS_SPECTRAL_GRID_X_1,  &
                                              _PS_SPECTRAL_GRID_X_2,  &
                                              _PS_SPECTRAL_GRID_X_3,  &
                                              _PS_SPECTRAL_GRID_X_4,  &
                                         MPI_Communicator             )

  logical  :: Distributed_x
  integer  :: N_x
  integer  :: _PS_SPECTRAL_GRID_X_0
  integer  :: _PS_SPECTRAL_GRID_X_2
  integer  :: _PS_SPECTRAL_GRID_X_3
  integer  :: _PS_SPECTRAL_GRID_X_4

  integer , OPTIONAL :: MPI_Communicator

  NN_g = N_x

  call MPI_Spectral_Grid_Index (_PS_SPECTRAL_GRID_X_1,  &
                                _PS_SPECTRAL_GRID_X_2,  &
                                _PS_SPECTRAL_GRID_X_3,  &
                                _PS_SPECTRAL_GRID_X_4,  &
                                MPI_Communicator        )

  END Subroutine MPI_Spectral_Grid_Index_1D
!
! ============================================================================
!
  Subroutine MPI_SPECTRAL_Grid_Index_2D (N_x, _PS_SPECTRAL_GRID_X_1,  &
                                              _PS_SPECTRAL_GRID_X_2,  &
                                              _PS_SPECTRAL_GRID_X_3,  &
                                              _PS_SPECTRAL_GRID_X_4,  &
                                         N_y, _PS_SPECTRAL_GRID_Y_1,  &
                                              _PS_SPECTRAL_GRID_Y_2,  &
                                              _PS_SPECTRAL_GRID_Y_3,  &
                                              _PS_SPECTRAL_GRID_Y_4,  &
                                         MPI_Communicator             )

  logical  :: Distributed_x
  integer  :: N_x
  integer  :: _PS_SPECTRAL_GRID_X_0
  integer  :: _PS_SPECTRAL_GRID_X_2
  integer  :: _PS_SPECTRAL_GRID_X_3
  integer  :: _PS_SPECTRAL_GRID_X_4

  logical  :: Distributed_y
  integer  :: N_y
  integer  :: _PS_SPECTRAL_GRID_Y_0
  integer  :: _PS_SPECTRAL_GRID_Y_2
  integer  :: _PS_SPECTRAL_GRID_Y_3
  integer  :: _PS_SPECTRAL_GRID_Y_4

  integer , OPTIONAL :: MPI_Communicator

  NN_g = N_x ; MM_g = N_y

  call MPI_Spectral_Grid_Index (_PS_SPECTRAL_GRID_X_1,  &
                                _PS_SPECTRAL_GRID_X_2,  &
                                _PS_SPECTRAL_GRID_X_3,  &
                                _PS_SPECTRAL_GRID_X_4,  &
                                MPI_Communicator        )

  call MPI_Spectral_Grid_Index (_PS_SPECTRAL_GRID_Y_1,  &
                                _PS_SPECTRAL_GRID_Y_2,  &
                                _PS_SPECTRAL_GRID_Y_3,  &
                                _PS_SPECTRAL_GRID_Y_4,  &
                                MPI_Communicator        )

  END Subroutine MPI_Spectral_Grid_Index_2D
!
! ============================================================================
!
  Subroutine MPI_SPECTRAL_Grid_Index_3D (N_x, _PS_SPECTRAL_GRID_X_1,  &
                                              _PS_SPECTRAL_GRID_X_2,  &
                                              _PS_SPECTRAL_GRID_X_3,  &
                                              _PS_SPECTRAL_GRID_X_4,  &
                                         N_y, _PS_SPECTRAL_GRID_Y_1,  &
                                              _PS_SPECTRAL_GRID_Y_2,  &
                                              _PS_SPECTRAL_GRID_Y_3,  &
                                              _PS_SPECTRAL_GRID_Y_4,  &
                                         N_z, _PS_SPECTRAL_GRID_Z_1,  &
                                              _PS_SPECTRAL_GRID_Z_2,  &
                                              _PS_SPECTRAL_GRID_Z_3,  &
                                              _PS_SPECTRAL_GRID_Z_4,  &
                                         MPI_Communicator             )

  logical  :: Distributed_x
  integer  :: N_x
  integer  :: _PS_SPECTRAL_GRID_X_0
  integer  :: _PS_SPECTRAL_GRID_X_2
  integer  :: _PS_SPECTRAL_GRID_X_3
  integer  :: _PS_SPECTRAL_GRID_X_4

  logical  :: Distributed_y
  integer  :: N_y
  integer  :: _PS_SPECTRAL_GRID_Y_0
  integer  :: _PS_SPECTRAL_GRID_Y_2
  integer  :: _PS_SPECTRAL_GRID_Y_3
  integer  :: _PS_SPECTRAL_GRID_Y_4

  logical  :: Distributed_z
  integer  :: N_z
  integer  :: _PS_SPECTRAL_GRID_Z_0
  integer  :: _PS_SPECTRAL_GRID_Z_2
  integer  :: _PS_SPECTRAL_GRID_Z_3
  integer  :: _PS_SPECTRAL_GRID_Z_4

  integer , OPTIONAL :: MPI_Communicator

  NN_g = N_x ; MM_g = N_y ; KK_g = N_z

  call MPI_Spectral_Grid_Index (_PS_SPECTRAL_GRID_X_1,  &
                                _PS_SPECTRAL_GRID_X_2,  &
                                _PS_SPECTRAL_GRID_X_3,  &
                                _PS_SPECTRAL_GRID_X_4,  &
                                MPI_Communicator        )

  call MPI_Spectral_Grid_Index (_PS_SPECTRAL_GRID_Y_1,  &
                                _PS_SPECTRAL_GRID_Y_2,  &
                                _PS_SPECTRAL_GRID_Y_3,  &
                                _PS_SPECTRAL_GRID_Y_4,  &
                                MPI_Communicator        )

  call MPI_Spectral_Grid_Index (_PS_SPECTRAL_GRID_Z_1,  &
                                _PS_SPECTRAL_GRID_Z_2,  &
                                _PS_SPECTRAL_GRID_Z_3,  &
                                _PS_SPECTRAL_GRID_Z_4,  &
                                MPI_Communicator        )

  END Subroutine MPI_Spectral_Grid_Index_3D
!
! ============================================================================
!
  Subroutine Spectral_Grid_Index_1D (N_x, _PS_SPECTRAL_GRID_X_0,  &
                                             _PS_SPECTRAL_GRID_X_3   )
  
  integer  :: N_x
  integer  :: _PS_SPECTRAL_GRID_X_0
  integer  :: _PS_SPECTRAL_GRID_X_3

  NN = N_x

  call Spectral_Grid_Index (_PS_SPECTRAL_GRID_X_0,  &
                            _PS_SPECTRAL_GRID_X_3   )
  
  END Subroutine Spectral_Grid_Index_1D
!
! ============================================================================
!
  Subroutine Spectral_Grid_Index_2D (N_x, _PS_SPECTRAL_GRID_X_0,  &
                                          _PS_SPECTRAL_GRID_X_3,  &
                                     N_y, _PS_SPECTRAL_GRID_Y_0,  &
                                          _PS_SPECTRAL_GRID_Y_3   )

  integer  :: N_x
  integer  :: _PS_SPECTRAL_GRID_X_0
  integer  :: _PS_SPECTRAL_GRID_X_3

  integer  :: N_y
  integer  :: _PS_SPECTRAL_GRID_Y_0
  integer  :: _PS_SPECTRAL_GRID_Y_3

  NN = N_x  ; MM = N_y

  call Spectral_Grid_Index (_PS_SPECTRAL_GRID_X_0,  &
                            _PS_SPECTRAL_GRID_X_3   )

  call Spectral_Grid_Index (_PS_SPECTRAL_GRID_Y_0,  &
                            _PS_SPECTRAL_GRID_Y_3   )

  END Subroutine Spectral_Grid_Index_2D
!
! ============================================================================
!
  Subroutine Spectral_Grid_Index_3D (N_x, _PS_SPECTRAL_GRID_X_0,  &
                                          _PS_SPECTRAL_GRID_X_3,  &
                                     N_y, _PS_SPECTRAL_GRID_Y_0,  &
                                          _PS_SPECTRAL_GRID_Y_3,  &
                                     N_z, _PS_SPECTRAL_GRID_Z_0,  &
                                          _PS_SPECTRAL_GRID_Z_3   )
  integer  :: N_x
  integer  :: _PS_SPECTRAL_GRID_X_0
  integer  :: _PS_SPECTRAL_GRID_X_3

  integer  :: N_y
  integer  :: _PS_SPECTRAL_GRID_Y_0
  integer  :: _PS_SPECTRAL_GRID_Y_3

  integer  :: N_z
  integer  :: _PS_SPECTRAL_GRID_Z_0
  integer  :: _PS_SPECTRAL_GRID_Z_3

  NN = N_x ; MM = N_y ; KK = N_z

  call Spectral_Grid_Index (_PS_SPECTRAL_GRID_X_0,  &
                            _PS_SPECTRAL_GRID_X_3   )

  call Spectral_Grid_Index (_PS_SPECTRAL_GRID_Y_0,  &
                            _PS_SPECTRAL_GRID_Y_3   )

  call Spectral_Grid_Index (_PS_SPECTRAL_GRID_Z_0,  &
                            _PS_SPECTRAL_GRID_Z_3   )

  END Subroutine Spectral_Grid_Index_3D
!
!======================================================================
!
  Subroutine Spectral_Grid_Index (Index, Method, Symmetry, Point_Type, &
                                     N_Pnts, N, N0, N1, N2, N3, N4, N5)

  integer  :: Index, Method, Symmetry, Point_Type, N_Pnts
  integer  :: N, N0, N1, N2, N3, N4, N5

  call PS_Get_Grid_Size (N_Pnts, Method, Point_Type, Symmetry, N)

  N0 = 0 ; N1 = 0 ; N2 = 0 ; N3 = N-1 ; N4 = N-1 ; N5 = N-1

  END Subroutine Spectral_Grid_Index

END MODULE MPI_Spectral_Grid_Indexing
