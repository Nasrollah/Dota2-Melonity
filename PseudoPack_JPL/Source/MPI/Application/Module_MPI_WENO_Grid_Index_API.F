#include "pseudopack.h"
#include "grid.h"

MODULE MPI_WENO_Grid_Indexing

  USE Processor
  USE MPI_Processor_Topology
  USE MPI_Data_Distribution

implicit NONE

INTERFACE PS_MPI_WENO_Grid_Index
  MODULE PROCEDURE MPI_WENO_Grid_Index_1D
  MODULE PROCEDURE MPI_WENO_Grid_Index_2D
  MODULE PROCEDURE MPI_WENO_Grid_Index_3D
END INTERFACE

INTERFACE PS_WENO_Grid_Index
  MODULE PROCEDURE     WENO_Grid_Index_1D
  MODULE PROCEDURE     WENO_Grid_Index_2D
  MODULE PROCEDURE     WENO_Grid_Index_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_MPI_WENO_Grid_Index
PUBLIC  :: PS_WENO_Grid_Index

CONTAINS
!
! ============================================================================
!
  Subroutine MPI_WENO_Grid_Index (_PS_WENO_GRID_X_1,  &
                                  _PS_WENO_GRID_X_2,  &
                                  _PS_WENO_GRID_X_3,  &
                                  _PS_WENO_GRID_X_4,  &
                                  MPI_Communicator    )

  logical  :: Distributed_x
  integer  :: _PS_WENO_GRID_X_0
  integer  :: _PS_WENO_GRID_X_2
  integer  :: _PS_WENO_GRID_X_3
  integer  :: _PS_WENO_GRID_X_4

  integer , OPTIONAL :: MPI_Communicator

!Set Global Index
  call WENO_Grid_Index (_PS_WENO_GRID_X_0,  &
                        _PS_WENO_GRID_X_2)
  
  if (Distributed_x) then
    call PS_MPI_Data_Distribution (N3_g-N2_g+1, _PS_WENO_GRID_X_4, &
                                   Index_x, MPI_Communicator)
  
! Adjust n_0 and n_Star for the WENO Cell Center scheme
    n_0 = n_0 + N2_g

    call PS_MPI_Processor_Local_Info (I_Am, Last_Processor, &
                                      Index_x, MPI_Communicator)

    NN = n_Star+1
  
!Set Local Index
    call WENO_Grid_Index (Index_x, 1, 0, WENO_Order_x, &
                          _PS_WENO_GRID_X_3)

    call PS_MPI_Processor (I_Am, Last_Processor)
  else
    n_Remain = 0 ; n_0 = 0 ; n_Star = N3_g-N2_g+1

    N  = N_g  
    N0 = N0_g ; N1 = N1_g ; N2 = N2_g ; N3 = N3_g ; N4 = N4_g ; N5 = N5_g
  endif

! write (6,*) Index_x, NN_g, N_g, N, N0, N1, N2, N3, N4, N5

  END Subroutine MPI_WENO_Grid_Index
!
! ============================================================================
!
  Subroutine MPI_WENO_Grid_Index_1D (N_x, _PS_WENO_GRID_X_1,  &
                                          _PS_WENO_GRID_X_2,  &
                                          _PS_WENO_GRID_X_3,  &
                                          _PS_WENO_GRID_X_4,  &
                                     MPI_Communicator         )

  logical  :: Distributed_x
  integer  :: N_x
  integer  :: _PS_WENO_GRID_X_0
  integer  :: _PS_WENO_GRID_X_2
  integer  :: _PS_WENO_GRID_X_3
  integer  :: _PS_WENO_GRID_X_4

  integer , OPTIONAL :: MPI_Communicator
 
  NN_g = N_x

  call MPI_WENO_Grid_Index (_PS_WENO_GRID_X_1,  &
                            _PS_WENO_GRID_X_2,  &
                            _PS_WENO_GRID_X_3,  &
                            _PS_WENO_GRID_X_4,  &
                            MPI_Communicator    )

  END Subroutine MPI_WENO_Grid_Index_1D 
!
! ============================================================================
!
  Subroutine MPI_WENO_Grid_Index_2D (N_x, _PS_WENO_GRID_X_1,  &
                                          _PS_WENO_GRID_X_2,  &
                                          _PS_WENO_GRID_X_3,  &
                                          _PS_WENO_GRID_X_4,  &
                                     N_y, _PS_WENO_GRID_Y_1,  &
                                          _PS_WENO_GRID_Y_2,  &
                                          _PS_WENO_GRID_Y_3,  &
                                          _PS_WENO_GRID_Y_4,  &
                                     MPI_Communicator         )

  logical  :: Distributed_x
  integer  :: N_x
  integer  :: _PS_WENO_GRID_X_0
  integer  :: _PS_WENO_GRID_X_2
  integer  :: _PS_WENO_GRID_X_3
  integer  :: _PS_WENO_GRID_X_4

  logical  :: Distributed_y
  integer  :: N_y
  integer  :: _PS_WENO_GRID_Y_0
  integer  :: _PS_WENO_GRID_Y_2
  integer  :: _PS_WENO_GRID_Y_3
  integer  :: _PS_WENO_GRID_Y_4

  integer , OPTIONAL :: MPI_Communicator
 
  NN_g = N_x ; MM_g = N_y

  call MPI_WENO_Grid_Index (_PS_WENO_GRID_X_1,  &
                            _PS_WENO_GRID_X_2,  &
                            _PS_WENO_GRID_X_3,  &
                            _PS_WENO_GRID_X_4,  &
                            MPI_Communicator    )

  call MPI_WENO_Grid_Index (_PS_WENO_GRID_Y_1,  &
                            _PS_WENO_GRID_Y_2,  &
                            _PS_WENO_GRID_Y_3,  &
                            _PS_WENO_GRID_Y_4,  &
                            MPI_Communicator    )

  END Subroutine MPI_WENO_Grid_Index_2D 
!
! ============================================================================
!
  Subroutine MPI_WENO_Grid_Index_3D (N_x, _PS_WENO_GRID_X_1,  &
                                          _PS_WENO_GRID_X_2,  &
                                          _PS_WENO_GRID_X_3,  &
                                          _PS_WENO_GRID_X_4,  &
                                     N_y, _PS_WENO_GRID_Y_1,  &
                                          _PS_WENO_GRID_Y_2,  &
                                          _PS_WENO_GRID_Y_3,  &
                                          _PS_WENO_GRID_Y_4,  &
                                     N_z, _PS_WENO_GRID_Z_1,  &
                                          _PS_WENO_GRID_Z_2,  &
                                          _PS_WENO_GRID_Z_3,  &
                                          _PS_WENO_GRID_Z_4,  &
                                     MPI_Communicator         )

  logical  :: Distributed_x
  integer  :: N_x
  integer  :: _PS_WENO_GRID_X_0
  integer  :: _PS_WENO_GRID_X_2
  integer  :: _PS_WENO_GRID_X_3
  integer  :: _PS_WENO_GRID_X_4

  logical  :: Distributed_y
  integer  :: N_y
  integer  :: _PS_WENO_GRID_Y_0
  integer  :: _PS_WENO_GRID_Y_2
  integer  :: _PS_WENO_GRID_Y_3
  integer  :: _PS_WENO_GRID_Y_4

  logical  :: Distributed_z
  integer  :: N_z
  integer  :: _PS_WENO_GRID_Z_0
  integer  :: _PS_WENO_GRID_Z_2
  integer  :: _PS_WENO_GRID_Z_3
  integer  :: _PS_WENO_GRID_Z_4

  integer , OPTIONAL :: MPI_Communicator

  NN_g = N_x ; MM_g = N_y ; KK_g = N_z

  call MPI_WENO_Grid_Index (_PS_WENO_GRID_X_1,  &
                            _PS_WENO_GRID_X_2,  &
                            _PS_WENO_GRID_X_3,  &
                            _PS_WENO_GRID_X_4,  &
                            MPI_Communicator    )

  call MPI_WENO_Grid_Index (_PS_WENO_GRID_Y_1,  &
                            _PS_WENO_GRID_Y_2,  &
                            _PS_WENO_GRID_Y_3,  &
                            _PS_WENO_GRID_Y_4,  &
                            MPI_Communicator    )

  call MPI_WENO_Grid_Index (_PS_WENO_GRID_Z_1,  &
                            _PS_WENO_GRID_Z_2,  &
                            _PS_WENO_GRID_Z_3,  &
                            _PS_WENO_GRID_Z_4,  &
                            MPI_Communicator    )

  END Subroutine MPI_WENO_Grid_Index_3D 
!
! ============================================================================
!
  Subroutine WENO_Grid_Index_1D (N_x, _PS_WENO_GRID_X_0,  &
                                      _PS_WENO_GRID_X_3   )
  
  integer  :: N_x
  integer  :: _PS_WENO_GRID_X_0
  integer  :: _PS_WENO_GRID_X_3

  NN = N_x 

  call WENO_Grid_Index (_PS_WENO_GRID_X_0,  &
                        _PS_WENO_GRID_X_3   )

  END Subroutine WENO_Grid_Index_1D 
!
! ============================================================================
!
  Subroutine WENO_Grid_Index_2D (N_x, _PS_WENO_GRID_X_0,  &
                                      _PS_WENO_GRID_X_3,  &
                                 N_y, _PS_WENO_GRID_Y_0,  &
                                      _PS_WENO_GRID_Y_3   )

  integer  :: N_x
  integer  :: _PS_WENO_GRID_X_0
  integer  :: _PS_WENO_GRID_X_3
  
  integer  :: N_y
  integer  :: _PS_WENO_GRID_Y_0
  integer  :: _PS_WENO_GRID_Y_3

  NN = N_x  ; MM = N_y

  call WENO_Grid_Index (_PS_WENO_GRID_X_0,  &
                        _PS_WENO_GRID_X_3   )
  
  call WENO_Grid_Index (_PS_WENO_GRID_Y_0,  &
                        _PS_WENO_GRID_Y_3   )

  END Subroutine WENO_Grid_Index_2D 
!
! ============================================================================
!
  Subroutine WENO_Grid_Index_3D (N_x, _PS_WENO_GRID_X_0,  &
                                      _PS_WENO_GRID_X_3,  &
                                 N_y, _PS_WENO_GRID_Y_0,  &
                                      _PS_WENO_GRID_Y_3,  &
                                 N_z, _PS_WENO_GRID_Z_0,  &
                                      _PS_WENO_GRID_Z_3   )

  integer  :: N_x
  integer  :: _PS_WENO_GRID_X_0
  integer  :: _PS_WENO_GRID_X_3
  
  integer  :: N_y
  integer  :: _PS_WENO_GRID_Y_0
  integer  :: _PS_WENO_GRID_Y_3

  integer  :: N_z
  integer  :: _PS_WENO_GRID_Z_0
  integer  :: _PS_WENO_GRID_Z_3

  NN = N_x ; MM = N_y ; KK = N_z

  call WENO_Grid_Index (_PS_WENO_GRID_X_0,  &
                        _PS_WENO_GRID_X_3   )
  
  call WENO_Grid_Index (_PS_WENO_GRID_Y_0,  &
                        _PS_WENO_GRID_Y_3   )

  call WENO_Grid_Index (_PS_WENO_GRID_Z_0,  &
                        _PS_WENO_GRID_Z_3   )

  END Subroutine WENO_Grid_Index_3D 
!
!======================================================================
!
  Subroutine WENO_Grid_Index (Index, Method, Symmetry, Order, &
                              N_Pnts, N, N0, N1, N2, N3, N4, N5, Tight_Bound)

  integer  :: Index, Method, Symmetry, Order, N_Pnts
  integer  :: N, N0, N1, N2, N3, N4, N5, M_s

  integer , OPTIONAL :: Tight_Bound

  M_s = (Order+1)/2

  N   = N_Pnts-2 ; if (Symmetry == 1) N = (N-1)/2

  N0 = -M_s ; N1 = -1 ; N2 = 0 ; N3 = N ; N4 = N3+1 ; N5 = N3+M_s

! write (6,*) Index, Method, Symmetry, Order, N_Pnts
! write (6,*) N, N0, N1, N2, N3, N4, N5, M_s

  if (.NOT. PRESENT(Tight_Bound)) RETURN

  if (Tight_Bound == 1) N5 = N5-1

  END Subroutine WENO_Grid_Index

END MODULE MPI_WENO_Grid_Indexing
