#include "pseudopack.h"

#define SENDRECV_TOP	SIZE(P_Top ), MPI_REALTYPE, Down, Send_Tag
#define SENDRECV_DOWN	SIZE(P_Down), MPI_REALTYPE, Top , Recv_Tag
#define SEND_STATUS	MPI_Communicator,             MPI_Error_Status
#define RECV_STATUS	MPI_Communicator, MPI_Status, MPI_Error_Status
#define TOP_DOWN_STATUS	Top_Or_Down, Top, Down, MPI_Error_Status

MODULE MPI_Data_Exchange

  USE Processor

implicit NONE

INTERFACE PS_MPI_Exchange_Boundary_Data
  MODULE PROCEDURE Exchange_Boundary_Data_0D
  MODULE PROCEDURE Exchange_Boundary_Data_1D
  MODULE PROCEDURE Exchange_Boundary_Data_2D
  MODULE PROCEDURE Exchange_Boundary_Data_3D

  MODULE PROCEDURE Exchange_Boundary_Data_0DG
  MODULE PROCEDURE Exchange_Boundary_Data_1DG
  MODULE PROCEDURE Exchange_Boundary_Data_2DG
  MODULE PROCEDURE Exchange_Boundary_Data_3DG
END INTERFACE

  integer  :: Status

PRIVATE
PUBLIC  :: PS_MPI_Exchange_Boundary_Data

CONTAINS
#if defined (PARALLEL_MPI)
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_0D  (Index, Distributed, M, Q, Direction,  &
                                                MPI_Comm_Type)

  integer  :: Index, M
  logical  :: Distributed

  REALTYPE, dimension(:), TARGET  :: Q
  REALTYPE, dimension(:), POINTER :: P_Top, P_Down

  integer , OPTIONAL :: Direction

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer  :: N0, N5, M0, M5
  integer  :: L0, L1, L4, L5, P0, P1, P4, P5
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag
  integer  :: M_Processor, Processor_Last
  logical  :: Topology

  integer , dimension(MPI_Status_Size) :: MPI_Status

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)

  Nullify (P_Top, P_Down)

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1)

      P0 = L4 ; P1 = L5 ; P4 = L0 ; P5 = L1

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5)

      P0 = L0 ; P1 = L1 ; P4 = L4 ; P5 = L5

  END SELECT 

  Send_Tag = 5555 ;  Recv_Tag = 5555

  if (PRESENT(MPI_Comm_Type))  &
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)

  Topology = PRESENT(MPI_Comm_Type) .AND. (Status == MPI_CART)

#if defined (USE_POINTER)
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (P_Top   , SENDRECV_TOP ,            &
                           P_Down  , SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (P_Top   , SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (P_Down  , SENDRECV_DOWN, RECV_STATUS)
  endif
#else
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (Q(P0:P1), SENDRECV_TOP ,            &
                           Q(P4:P5), SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (Q(P0:P1), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(P4:P5), SENDRECV_DOWN, RECV_STATUS)
  endif
#endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_0D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_1D  (Index, Distributed, M, Q, Direction,  &
                                                MPI_Comm_Type)

  integer  :: Index, M
  logical  :: Distributed

  REALTYPE, dimension(:,:), TARGET  :: Q
  REALTYPE, dimension(:,:), POINTER :: P_Top, P_Down

  integer , OPTIONAL :: Direction

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer  :: N0, N5, M0, M5
  integer  :: L0, L1, L4, L5, P0, P1, P4, P5
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag
  integer  :: M_Processor, Processor_Last
  logical  :: Topology

  integer , dimension(MPI_Status_Size) :: MPI_Status

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)

  Nullify (P_Top, P_Down)

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5,:)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1,:)

      P0 = L4 ; P1 = L5 ; P4 = L0 ; P5 = L1

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1,:)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5,:)

      P0 = L0 ; P1 = L1 ; P4 = L4 ; P5 = L5

  END SELECT 

  Send_Tag = 5555 ;  Recv_Tag = 5555

  if (PRESENT(MPI_Comm_Type))  &
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)

  Topology = PRESENT(MPI_Comm_Type) .AND. (Status == MPI_CART)

#if defined (USE_POINTER)
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (P_Top     , SENDRECV_TOP ,            &
                           P_Down    , SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (P_Top     , SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (P_Down    , SENDRECV_DOWN, RECV_STATUS)
  endif
#else
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (Q(P0:P1,:), SENDRECV_TOP ,            &
                           Q(P4:P5,:), SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (Q(P0:P1,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(P4:P5,:), SENDRECV_DOWN, RECV_STATUS)
  endif
#endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_1D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_2D  (Index, Distributed, M, Q, Direction,  &
                                                MPI_Comm_Type)

  integer  :: Index, M
  logical  :: Distributed

  REALTYPE, dimension(:,:,:), TARGET  :: Q
  REALTYPE, dimension(:,:,:), POINTER :: P_Top, P_Down

  integer , OPTIONAL :: Direction

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer  :: N0, N5, M0, M5
  integer  :: L0, L1, L4, L5, P0, P1, P4, P5
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag
  integer  :: M_Processor, Processor_Last
  logical  :: Topology

  integer , dimension(MPI_Status_Size) :: MPI_Status

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)
  M0 = LBOUND(Q,DIM=2) ; M5 = UBOUND(Q,DIM=2)

  Nullify (P_Top, P_Down)

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

      SELECT CASE (Index)
        CASE (1)
          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5,:,:)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1,:,:)

        CASE (2)
          L4 = M5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(:,L4:L5,:)
          L0 = M0       ; L1 = L0+M-1 ; P_Down => Q(:,L0:L1,:)

      END SELECT

      P0 = L4 ; P1 = L5 ; P4 = L0 ; P5 = L1

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

      SELECT CASE (Index)
        CASE (1)
          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1,:,:)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5,:,:)

        CASE (2)
          L0 = M0+M     ; L1 = L0+M-1 ; P_Top  => Q(:,L0:L1,:)
          L4 = M5-M+1   ; L5 = M5     ; P_Down => Q(:,L4:L5,:)

      END SELECT

      P0 = L0 ; P1 = L1 ; P4 = L4 ; P5 = L5

  END SELECT 

  Send_Tag = 5555 ;  Recv_Tag = 5555

  if (PRESENT(MPI_Comm_Type))  &
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)

  Topology = PRESENT(MPI_Comm_Type) .AND. (Status == MPI_CART)

#if defined (USE_POINTER)
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (P_Top        , SENDRECV_TOP ,            &
                           P_Down       , SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (P_Top        , SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (P_Down       , SENDRECV_DOWN, RECV_STATUS)
  endif
#else
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

    SELECT CASE (Index)
      CASE (1)
        call MPI_SendRecv (Q(P0:P1,:,:), SENDRECV_TOP ,            &
                           Q(P4:P5,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (2)
        call MPI_SendRecv (Q(:,P0:P1,:), SENDRECV_TOP ,            &
                           Q(:,P4:P5,:), SENDRECV_DOWN, RECV_STATUS)

    END SELECT
  else
    SELECT CASE (Index)
      CASE (1)
        call MPI_Send     (Q(P0:P1,:,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(P4:P5,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (2)
        call MPI_Send     (Q(:,P0:P1,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(:,P4:P5,:), SENDRECV_DOWN, RECV_STATUS)

    END SELECT
  endif
#endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_2D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_3D  (Index, Distributed, M, Q, Direction,  &
                                                MPI_Comm_Type)

  integer  :: Index, M
  logical  :: Distributed

  REALTYPE, dimension(:,:,:,:), TARGET  :: Q
  REALTYPE, dimension(:,:,:,:), POINTER :: P_Top, P_Down

  integer , OPTIONAL :: Direction

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer  :: N0, N5, M0, M5, K0, K5, NV
  integer  :: L0, L1, L4, L5, P0, P1, P4, P5
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag
  integer  :: M_Processor, Processor_Last
  logical  :: Topology

  integer , dimension(MPI_Status_Size) :: MPI_Status

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)
  M0 = LBOUND(Q,DIM=2) ; M5 = UBOUND(Q,DIM=2)
  K0 = LBOUND(Q,DIM=3) ; K5 = UBOUND(Q,DIM=3)

  Nullify (P_Top, P_Down)

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

      SELECT CASE (Index)
        CASE (1)
          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5,:,:,:)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1,:,:,:)

        CASE (2)
          L4 = M5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(:,L4:L5,:,:)
          L0 = M0       ; L1 = L0+M-1 ; P_Down => Q(:,L0:L1,:,:)

        CASE (3)
          L4 = K5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(:,:,L4:L5,:)
          L0 = K0       ; L1 = L0+M-1 ; P_Down => Q(:,:,L0:L1,:)

      END SELECT

      P0 = L4 ; P1 = L5 ; P4 = L0 ; P5 = L1

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

      SELECT CASE (Index)
        CASE (1)
          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1,:,:,:)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5,:,:,:)

        CASE (2)
          L0 = M0+M     ; L1 = L0+M-1 ; P_Top  => Q(:,L0:L1,:,:)
          L4 = M5-M+1   ; L5 = M5     ; P_Down => Q(:,L4:L5,:,:)

        CASE (3)
          L0 = K0+M     ; L1 = L0+M-1 ; P_Top  => Q(:,:,L0:L1,:)
          L4 = K5-M+1   ; L5 = K5     ; P_Down => Q(:,:,L4:L5,:)

      END SELECT

      P0 = L0 ; P1 = L1 ; P4 = L4 ; P5 = L5

  END SELECT 

  Send_Tag = 5555 ;  Recv_Tag = 5555

  if (PRESENT(MPI_Comm_Type))  &
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)

  Topology = PRESENT(MPI_Comm_Type) .AND. (Status == MPI_CART)

#if defined (USE_POINTER)
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (P_Top         , SENDRECV_TOP ,            &
                           P_Down        , SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (P_Top         , SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (P_Down        , SENDRECV_DOWN, RECV_STATUS)
  endif
#else
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

    SELECT CASE (Index)
      CASE (1)
        call MPI_SendRecv (Q(P0:P1,:,:,:), SENDRECV_TOP ,            &
                           Q(P4:P5,:,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (2)
        call MPI_SendRecv (Q(:,P0:P1,:,:), SENDRECV_TOP ,            &
                           Q(:,P4:P5,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (3)
        call MPI_SendRecv (Q(:,:,P0:P1,:), SENDRECV_TOP ,            &
                           Q(:,:,P4:P5,:), SENDRECV_DOWN, RECV_STATUS)

    END SELECT
  else
    SELECT CASE (Index)
      CASE (1)
        call MPI_Send     (Q(P0:P1,:,:,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(P4:P5,:,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (2)
        call MPI_Send     (Q(:,P0:P1,:,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(:,P4:P5,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (3)
        call MPI_Send     (Q(:,:,P0:P1,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(:,:,P4:P5,:), SENDRECV_DOWN, RECV_STATUS)
    END SELECT
  endif
#endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_3D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_0DG (Index, Distributed,            &
                                         N0,N5,                   M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0, N5, M
  logical  :: Distributed

  REALTYPE, dimension(N0:N5), TARGET  :: Q
  REALTYPE, dimension(:)    , POINTER :: P_Top, P_Down

  integer , OPTIONAL :: Direction

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer  :: L0, L1, L4, L5, P0, P1, P4, P5
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag
  integer  :: M_Processor, Processor_Last
  logical  :: Topology

  integer , dimension(MPI_Status_Size) :: MPI_Status

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)

  Nullify (P_Top, P_Down)

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1)

      P0 = L4 ; P1 = L5 ; P4 = L0 ; P5 = L1

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5)

      P0 = L0 ; P1 = L1 ; P4 = L4 ; P5 = L5

  END SELECT 

  Send_Tag = 5555 ;  Recv_Tag = 5555

  if (PRESENT(MPI_Comm_Type))  &
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)

  Topology = PRESENT(MPI_Comm_Type) .AND. (Status == MPI_CART)

#if defined (USE_POINTER)
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (P_Top   , SENDRECV_TOP ,            &
                           P_Down  , SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (P_Top   , SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (P_Down  , SENDRECV_DOWN, RECV_STATUS)
  endif
#else
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (Q(P0:P1), SENDRECV_TOP ,            &
                           Q(P4:P5), SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (Q(P0:P1), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(P4:P5), SENDRECV_DOWN, RECV_STATUS)
  endif
#endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_0DG
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_1DG (Index, Distributed,            &
                                         N0,N5,               NV, M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0, N5, NV, M
  logical  :: Distributed

  REALTYPE, dimension(N0:N5,NV), TARGET  :: Q
  REALTYPE, dimension(:,:)     , POINTER :: P_Top, P_Down

  integer , OPTIONAL :: Direction

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer  :: L0, L1, L4, L5, P0, P1, P4, P5
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag
  integer  :: M_Processor, Processor_Last
  logical  :: Topology

  integer , dimension(MPI_Status_Size) :: MPI_Status

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)

  Nullify (P_Top, P_Down)

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5,:)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1,:)

      P0 = L4 ; P1 = L5 ; P4 = L0 ; P5 = L1

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1,:)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5,:)

      P0 = L0 ; P1 = L1 ; P4 = L4 ; P5 = L5

  END SELECT 

  Send_Tag = 5555 ;  Recv_Tag = 5555

  if (PRESENT(MPI_Comm_Type))  &
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)

  Topology = PRESENT(MPI_Comm_Type) .AND. (Status == MPI_CART)

#if defined (USE_POINTER)
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (P_Top   , SENDRECV_TOP ,            &
                           P_Down  , SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (P_Top   , SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (P_Down  , SENDRECV_DOWN, RECV_STATUS)
  endif
#else
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (Q(P0:P1,:), SENDRECV_TOP ,            &
                           Q(P4:P5,:), SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (Q(P0:P1,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(P4:P5,:), SENDRECV_DOWN, RECV_STATUS)
  endif
#endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_1DG
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_2DG (Index, Distributed,            &
                                         N0,N5, M0,M5,        NV, M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0, N5, M0, M5, NV, M
  logical  :: Distributed

  REALTYPE, dimension(N0:N5,M0:M5,NV), TARGET  :: Q
  REALTYPE, dimension(:,:,:)         , POINTER :: P_Top, P_Down

  integer , OPTIONAL :: Direction

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer  :: L0, L1, L4, L5, P0, P1, P4, P5
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag
  integer  :: M_Processor, Processor_Last
  logical  :: Topology

  integer , dimension(MPI_Status_Size) :: MPI_Status

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)
  M0 = LBOUND(Q,DIM=2) ; M5 = UBOUND(Q,DIM=2)

  Nullify (P_Top, P_Down)

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

      SELECT CASE (Index)
        CASE (1)
          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5,:,:)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1,:,:)

        CASE (2)
          L4 = M5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(:,L4:L5,:)
          L0 = M0       ; L1 = L0+M-1 ; P_Down => Q(:,L0:L1,:)

      END SELECT

      P0 = L4 ; P1 = L5 ; P4 = L0 ; P5 = L1

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

      SELECT CASE (Index)
        CASE (1)
          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1,:,:)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5,:,:)

        CASE (2)
          L0 = M0+M     ; L1 = L0+M-1 ; P_Top  => Q(:,L0:L1,:)
          L4 = M5-M+1   ; L5 = M5     ; P_Down => Q(:,L4:L5,:)

      END SELECT

      P0 = L0 ; P1 = L1 ; P4 = L4 ; P5 = L5

  END SELECT 

  Send_Tag = 5555 ;  Recv_Tag = 5555

  if (PRESENT(MPI_Comm_Type))  &
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)

  Topology = PRESENT(MPI_Comm_Type) .AND. (Status == MPI_CART)

#if defined (USE_POINTER)
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (P_Top   , SENDRECV_TOP ,            &
                           P_Down  , SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (P_Top   , SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (P_Down  , SENDRECV_DOWN, RECV_STATUS)
  endif
#else
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

    SELECT CASE (Index)
      CASE (1)
        call MPI_SendRecv (Q(P0:P1,:,:), SENDRECV_TOP ,            &
                           Q(P4:P5,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (2)
        call MPI_SendRecv (Q(:,P0:P1,:), SENDRECV_TOP ,            &
                           Q(:,P4:P5,:), SENDRECV_DOWN, RECV_STATUS)

    END SELECT
  else
    SELECT CASE (Index)
      CASE (1)
        call MPI_Send     (Q(P0:P1,:,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(P4:P5,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (2)
        call MPI_Send     (Q(:,P0:P1,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(:,P4:P5,:), SENDRECV_DOWN, RECV_STATUS)

    END SELECT
  endif
#endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_2DG
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_3DG (Index, Distributed,            &
                                         N0,N5, M0,M5, K0,K5, NV, M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0, N5, M0, M5, K0, K5, NV, M
  logical  :: Distributed

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET  :: Q
  REALTYPE, dimension(:,:,:,:)             , POINTER :: P_Top, P_Down

  integer , OPTIONAL :: Direction

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer  :: L0, L1, L4, L5, P0, P1, P4, P5
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag
  integer  :: M_Processor, Processor_Last
  logical  :: Topology

  integer , dimension(MPI_Status_Size) :: MPI_Status

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)
  M0 = LBOUND(Q,DIM=2) ; M5 = UBOUND(Q,DIM=2)
  K0 = LBOUND(Q,DIM=3) ; K5 = UBOUND(Q,DIM=3)

  Nullify (P_Top, P_Down)

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

      SELECT CASE (Index)
        CASE (1)
          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5,:,:,:)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1,:,:,:)

        CASE (2)
          L4 = M5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(:,L4:L5,:,:)
          L0 = M0       ; L1 = L0+M-1 ; P_Down => Q(:,L0:L1,:,:)

        CASE (3)
          L4 = K5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(:,:,L4:L5,:)
          L0 = K0       ; L1 = L0+M-1 ; P_Down => Q(:,:,L0:L1,:)

      END SELECT

      P0 = L4 ; P1 = L5 ; P4 = L0 ; P5 = L1

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

      SELECT CASE (Index)
        CASE (1)
          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1,:,:,:)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5,:,:,:)

        CASE (2)
          L0 = M0+M     ; L1 = L0+M-1 ; P_Top  => Q(:,L0:L1,:,:)
          L4 = M5-M+1   ; L5 = M5     ; P_Down => Q(:,L4:L5,:,:)

        CASE (3)
          L0 = K0+M     ; L1 = L0+M-1 ; P_Top  => Q(:,:,L0:L1,:)
          L4 = K5-M+1   ; L5 = K5     ; P_Down => Q(:,:,L4:L5,:)

      END SELECT

      P0 = L0 ; P1 = L1 ; P4 = L4 ; P5 = L5

  END SELECT 

  Send_Tag = 5555 ;  Recv_Tag = 5555

  if (PRESENT(MPI_Comm_Type))  &
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)

  Topology = PRESENT(MPI_Comm_Type) .AND. (Status == MPI_CART)

#if defined (USE_POINTER)
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

        call MPI_SendRecv (P_Top   , SENDRECV_TOP ,            &
                           P_Down  , SENDRECV_DOWN, RECV_STATUS)
  else
        call MPI_Send     (P_Top   , SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (P_Down  , SENDRECV_DOWN, RECV_STATUS)
  endif
#else
  if (Topology) then
    call MPI_Cart_Shift (MPI_Communicator, Index-1, TOP_DOWN_STATUS)

    SELECT CASE (Index)
      CASE (1)
        call MPI_SendRecv (Q(P0:P1,:,:,:), SENDRECV_TOP ,            &
                           Q(P4:P5,:,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (2)
        call MPI_SendRecv (Q(:,P0:P1,:,:), SENDRECV_TOP ,            &
                           Q(:,P4:P5,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (3)
        call MPI_SendRecv (Q(:,:,P0:P1,:), SENDRECV_TOP ,            &
                           Q(:,:,P4:P5,:), SENDRECV_DOWN, RECV_STATUS)

    END SELECT
  else
    SELECT CASE (Index)
      CASE (1)
        call MPI_Send     (Q(P0:P1,:,:,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(P4:P5,:,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (2)
        call MPI_Send     (Q(:,P0:P1,:,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(:,P4:P5,:,:), SENDRECV_DOWN, RECV_STATUS)

      CASE (3)
        call MPI_Send     (Q(:,:,P0:P1,:), SENDRECV_TOP , SEND_STATUS)
        call MPI_Recv     (Q(:,:,P4:P5,:), SENDRECV_DOWN, RECV_STATUS)
    END SELECT
  endif
#endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_3DG
#else
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_0D  (Index, Distributed, M, Q, Direction,  &
                                                MPI_Comm_Type)

  integer  :: Index, M

  logical            :: Distributed
  integer , OPTIONAL :: Direction
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:)          :: Q

  RETURN

  END Subroutine Exchange_Boundary_Data_0D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_1D  (Index, Distributed, M, Q, Direction,  &
                                                MPI_Comm_Type)

  integer  :: Index, M

  logical            :: Distributed
  integer , OPTIONAL :: Direction
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:)          :: Q

  RETURN

  END Subroutine Exchange_Boundary_Data_1D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_2D  (Index, Distributed, M, Q, Direction, &
                                                MPI_Comm_Type)

  integer  :: Index, M

  logical            :: Distributed
  integer , OPTIONAL :: Direction
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:,:), TARGET  :: Q

  RETURN

  END Subroutine Exchange_Boundary_Data_2D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_3D  (Index, Distributed, M, Q, Direction, &
                                                MPI_Comm_Type)

  integer  :: Index, M

  logical            :: Distributed
  integer , OPTIONAL :: Direction
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(:,:,:,:)          :: Q

  RETURN

  END Subroutine Exchange_Boundary_Data_3D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_0DG (Index, Distributed,            &
                                         N0,N5,                   M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0,N5, M

  logical            :: Distributed
  integer , OPTIONAL :: Direction
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5)          :: Q

  RETURN

  END Subroutine Exchange_Boundary_Data_0DG
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_1DG (Index, Distributed,            &
                                         N0,N5,               NV, M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0,N5, NV, M

  logical            :: Distributed
  integer , OPTIONAL :: Direction
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,NV)          :: Q

  RETURN

  END Subroutine Exchange_Boundary_Data_1DG
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_2DG (Index, Distributed,            &
                                         N0,N5, M0,M5,        NV, M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0,N5, M0,M5, NV, M

  logical            :: Distributed
  integer , OPTIONAL :: Direction
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,NV)          :: Q

  RETURN

  END Subroutine Exchange_Boundary_Data_2DG
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_3DG (Index, Distributed,            &
                                         N0,N5, M0,M5, K0,K5, NV, M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0,N5, M0,M5, K0,K5, NV, M

  logical            :: Distributed
  integer , OPTIONAL :: Direction
  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)          :: Q

  RETURN

  END Subroutine Exchange_Boundary_Data_3DG
#endif

END MODULE MPI_Data_Exchange

#undef SENDRECV_TOP
#undef SENDRECV_DOWN
#undef SEND_STATUS
#undef RECV_STATUS
#undef TOP_DOWN_STATUS
