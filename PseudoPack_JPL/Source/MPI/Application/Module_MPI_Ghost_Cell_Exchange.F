#include "pseudopack.h"


MODULE MPI_Ghost_Cell_Exchange

  USE Processor
  USE MPI_Data_Exchange

implicit NONE

#if 0
INTERFACE PS_MPI_Ghost_Cell_Exchange
  MODULE PROCEDURE Ghost_Cell_Exchange_0D_A
  MODULE PROCEDURE Ghost_Cell_Exchange_1D_A
  MODULE PROCEDURE Ghost_Cell_Exchange_2D_A
  MODULE PROCEDURE Ghost_Cell_Exchange_3D_A

  MODULE PROCEDURE Ghost_Cell_Exchange_0D_B
  MODULE PROCEDURE Ghost_Cell_Exchange_1D_B
  MODULE PROCEDURE Ghost_Cell_Exchange_2D_B
  MODULE PROCEDURE Ghost_Cell_Exchange_3D_B
END INTERFACE
#endif

INTERFACE PS_MPI_Ghost_Cell_Exchange
  MODULE PROCEDURE Ghost_Cell_Exchange_0D_C
  MODULE PROCEDURE Ghost_Cell_Exchange_1D_C
  MODULE PROCEDURE Ghost_Cell_Exchange_2D_C
  MODULE PROCEDURE Ghost_Cell_Exchange_3D_C
END INTERFACE

PRIVATE
PUBLIC  :: PS_MPI_Ghost_Cell_Exchange

CONTAINS

#if 0
  Subroutine Ghost_Cell_Exchange_0D_A (Index_x, Distributed_x, M_s_x,    &
                                       Q, MPI_Communicator)  

  integer  :: Index_x, M_s_x
  logical  :: Distributed_x

  REALTYPE, dimension(:) :: Q

  integer , OPTIONAL :: MPI_Communicator

  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1, &
                                      MPI_Communicator)

  END Subroutine Ghost_Cell_Exchange_0D_A 
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_1D_A (Index_x, Distributed_x, M_s_x,    &
                                       Q, MPI_Communicator)  

  integer  :: Index_x, M_s_x
  logical  :: Distributed_x

  integer , OPTIONAL :: MPI_Communicator

  REALTYPE, dimension(:,:) :: Q

  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1, &
                                      MPI_Communicator)

  END Subroutine Ghost_Cell_Exchange_1D_A
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_2D_A (Index_x, Distributed_x, M_s_x,    &
                                       Index_y, Distributed_y, M_s_y,    &
                                       Q, MPI_Communicator)  

  integer  :: Index_x, M_s_x
  integer  :: Index_y, M_s_y
  logical  :: Distributed_x
  logical  :: Distributed_y

  REALTYPE, dimension(:,:,:) :: Q

  integer , OPTIONAL :: MPI_Communicator

  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1, &
                                      MPI_Communicator)

  call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q, -1, &
                                      MPI_Communicator)

  END Subroutine Ghost_Cell_Exchange_2D_A 
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_3D_A (Index_x, Distributed_x, M_s_x,    &
                                       Index_y, Distributed_y, M_s_y,    &
                                       Index_z, Distributed_z, M_s_z,    &
                                       Q, MPI_Communicator)  

  integer  :: Index_x, M_s_x
  integer  :: Index_y, M_s_y
  integer  :: Index_z, M_s_z
  logical  :: Distributed_x
  logical  :: Distributed_y
  logical  :: Distributed_z

  REALTYPE, dimension(:,:,:,:) :: Q

  integer , OPTIONAL :: MPI_Communicator

  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1, &
                                      MPI_Communicator)

  call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q, -1, &
                                      MPI_Communicator)

  call PS_MPI_Exchange_Boundary_Data (Index_z, Distributed_z, M_s_z, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_z, Distributed_z, M_s_z, Q, -1, &
                                      MPI_Communicator)

  END Subroutine Ghost_Cell_Exchange_3D_A 
!                                      
! -------------------------------------------------------------------------
!                                      
  Subroutine Ghost_Cell_Exchange_0D_B (Index_x, Distributed_x, N0, N2, N3, N5, &
                                           Q, MPI_Communicator)
  
  integer  :: Index_x, N0, N2, N3, N5, M_s_x
  logical  :: Distributed_x
  
  REALTYPE, dimension(N0:N5) :: Q

  integer , OPTIONAL :: MPI_Communicator

  M_s_x = N2-N0

  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1, &
                                      MPI_Communicator)
  
  END Subroutine Ghost_Cell_Exchange_0D_B 
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_1D_B (Index_x, Distributed_x, N0, N2, N3, N5, &
                                       NV, Q, MPI_Communicator)  

  integer  :: Index_x, N0, N2, N3, N5, M_s_x
  logical  :: Distributed_x
  integer  :: NV

  REALTYPE, dimension(N0:N5,NV) :: Q

  integer , OPTIONAL :: MPI_Communicator

  M_s_x = N2-N0

  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1, &
                                      MPI_Communicator)

  END Subroutine Ghost_Cell_Exchange_1D_B 
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_2D_B (Index_x, Distributed_x, N0, N2, N3, N5, &
                                       Index_y, Distributed_y, M0, M2, M3, M5, &
                                       NV, Q, MPI_Communicator)  

  integer  :: Index_x, N0, N2, N3, N5, M_s_x
  integer  :: Index_y, M0, M2, M3, M5, M_s_y
  logical  :: Distributed_x
  logical  :: Distributed_y
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q

  integer , OPTIONAL :: MPI_Communicator

  M_s_x = N2-N0
  M_s_y = M2-M0

  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1, &
                                      MPI_Communicator)

  call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q, -1, &
                                      MPI_Communicator)

  END Subroutine Ghost_Cell_Exchange_2D_B
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_3D_B (Index_x, Distributed_x, N0, N2, N3, N5, &
                                       Index_y, Distributed_y, M0, M2, M3, M5, &
                                       Index_z, Distributed_z, K0, K2, K3, K5, &
                                       NV, Q, MPI_Communicator)  

  integer  :: Index_x, N0, N2, N3, N5, M_s_x
  integer  :: Index_y, M0, M2, M3, M5, M_s_y
  integer  :: Index_z, K0, K2, K3, K5, M_s_z
  logical  :: Distributed_x
  logical  :: Distributed_y
  logical  :: Distributed_z
  integer  :: NV

  integer , OPTIONAL :: MPI_Communicator

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q

  M_s_x = N2-N0
  M_s_y = M2-M0
  M_s_z = K2-K0

  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1, &
                                      MPI_Communicator)

  call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q, -1, &
                                      MPI_Communicator)

  call PS_MPI_Exchange_Boundary_Data (Index_z, Distributed_z, M_s_z, Q,  1, &
                                      MPI_Communicator)
  call PS_MPI_Exchange_Boundary_Data (Index_z, Distributed_z, M_s_z, Q, -1, &
                                      MPI_Communicator)

  END Subroutine Ghost_Cell_Exchange_3D_B
#endif
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_0D_C (Index_x, Distributed_x,    &
                                         N0, N2, N3, N5,          &
                                           Q, MPI_Communicator)  

  integer  :: Index_x, N0, N2, N3, N5, M_s_x
  logical  :: Distributed_x

  REALTYPE, dimension(N0:N5)    :: Q

  integer            :: n, M_s_0, M_s

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type
  
                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  if (Distributed_x) then
    M_s_x = N2-N0 ; M_s = MIN(N3-N2+1,M_s_x)

#if defined (PARALLEL_MPI)
    call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                             MPI_Comm_Type, MPI_Error_Status  )
#else
    M_s_0 = M_s
#endif
  
    do n = M_s_0, M_s_x                  
      call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1,&
                                          MPI_Comm_Type)
      call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1,&
                                          MPI_Comm_Type)
    enddo
  endif

  END Subroutine Ghost_Cell_Exchange_0D_C
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_1D_C (Index_x, Distributed_x,    &
                                         N0, N2, N3, N5,          &
                                       NV, Q, MPI_Communicator)  

  integer  :: Index_x, N0, N2, N3, N5, M_s_x
  logical  :: Distributed_x
  integer  :: NV

  REALTYPE, dimension(N0:N5,NV) :: Q

  integer            :: n, M_s_0, M_s

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type
  
                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  if (Distributed_x) then
    M_s_x = N2-N0 ; M_s = MIN(N3-N2+1,M_s_x)

#if defined (PARALLEL_MPI)
    call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                             MPI_Comm_Type, MPI_Error_Status  )
#else
    M_s_0 = M_s
#endif
  
    do n = M_s_0, M_s_x                  
      call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1,&
                                          MPI_Comm_Type)
      call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1,&
                                          MPI_Comm_Type)
    enddo
  endif

  END Subroutine Ghost_Cell_Exchange_1D_C
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_2D_C (Index_x, Distributed_x,    &
                                         N0, N2, N3, N5,          &
                                       Index_y, Distributed_y,    &
                                         M0, M2, M3, M5,          &
                                       NV, Q, MPI_Communicator)  

  integer  :: Index_x, N0, N2, N3, N5, M_s_x
  integer  :: Index_y, M0, M2, M3, M5, M_s_y
  logical  :: Distributed_x
  logical  :: Distributed_y
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q

  integer            :: n, M_s_0, M_s

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type
  
                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  if (Distributed_x) then
    M_s_x = N2-N0 ; M_s = MIN(N3-N2+1,M_s_x)

#if defined (PARALLEL_MPI)
    call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                             MPI_Comm_Type, MPI_Error_Status  )
#else
    M_s_0 = M_s
#endif

    do n = M_s_0, M_s_x
      call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1,&
                                          MPI_Comm_Type)
      call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1,&
                                          MPI_Comm_Type)
    enddo
  endif

  if (Distributed_y) then
    M_s_y = M2-M0 ; M_s = MIN(M3-M2+1,M_s_y)

#if defined (PARALLEL_MPI)
    call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                             MPI_Comm_Type, MPI_Error_Status  )
#else
    M_s_0 = M_s
#endif

    do n = M_s_0, M_s_y
      call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q,  1,&
                                          MPI_Comm_Type)
      call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q, -1,&
                                          MPI_Comm_Type)
    enddo
  endif

  END Subroutine Ghost_Cell_Exchange_2D_C
!
! -------------------------------------------------------------------------
!
  Subroutine Ghost_Cell_Exchange_3D_C (Index_x, Distributed_x,    &
                                         N0, N2, N3, N5,          &
                                       Index_y, Distributed_y,    &
                                         M0, M2, M3, M5,          &
                                       Index_z, Distributed_z,    &
                                         K0, K2, K3, K5,          &
                                       NV, Q, MPI_Communicator)  

  integer  :: Index_x, N0, N2, N3, N5, M_s_x
  integer  :: Index_y, M0, M2, M3, M5, M_s_y
  integer  :: Index_z, K0, K2, K3, K5, M_s_z
  logical  :: Distributed_x
  logical  :: Distributed_y
  logical  :: Distributed_z
  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q

  integer            :: n, M_s_0, M_s

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type
  
                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  if (Distributed_x) then
    M_s_x = N2-N0 ; M_s = MIN(N3-N2+1,M_s_x)

#if defined (PARALLEL_MPI)
    call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                             MPI_Comm_Type, MPI_Error_Status  )
#else
    M_s_0 = M_s
#endif

    do n = M_s_0, M_s_x
      call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1,&
                                          MPI_Comm_Type)
      call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1,&
                                          MPI_Comm_Type)
    enddo
  endif

  if (Distributed_y) then
    M_s_y = M2-M0 ; M_s = MIN(M3-M2+1,M_s_y)

#if defined (PARALLEL_MPI)
    call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                             MPI_Comm_Type, MPI_Error_Status  )
#else
    M_s_0 = M_s
#endif

    do n = M_s_0, M_s_y
      call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q,  1,&
                                          MPI_Comm_Type)
      call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q, -1,&
                                          MPI_Comm_Type)
    enddo
  endif

  if (Distributed_z) then
    M_s_z = K2-K0 ; M_s = MIN(K3-K2+1,M_s_z)

#if defined (PARALLEL_MPI)
    call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                             MPI_Comm_Type, MPI_Error_Status  )
#else
    M_s_0 = M_s
#endif

    do n = M_s_0, M_s_z
      call PS_MPI_Exchange_Boundary_Data (Index_z, Distributed_z, M_s_z, Q,  1,&
                                          MPI_Comm_Type)
      call PS_MPI_Exchange_Boundary_Data (Index_z, Distributed_z, M_s_z, Q, -1,&
                                          MPI_Comm_Type)
    enddo
  endif

  END Subroutine Ghost_Cell_Exchange_3D_C

END MODULE MPI_Ghost_Cell_Exchange
