#include "pseudopack.h"


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

  integer  :: Index, N0,N5, L0,L1,L4,L5, M
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag

  integer , dimension(MPI_Status_Size) :: MPI_Status

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer , OPTIONAL :: Direction

  logical            :: Distributed

  REALTYPE, dimension(:), TARGET  :: Q
  REALTYPE, dimension(:), POINTER :: P_Top, P_Down

  integer  :: M_Processor, Processor_Last

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)

  Nullify (P_Top, P_Down)

  Send_Tag = 5555 ;  Recv_Tag = 5555

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1)

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5)

  END SELECT

  if (PRESENT(MPI_Comm_Type)) then
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)
  endif

  if ((PRESENT(MPI_Comm_Type)) .AND. (Status == MPI_CART)) then
    call MPI_CART_SHIFT (MPI_Communicator, Index-1, Top_Or_Down,           &
                                           Top, Down, MPI_Error_Status)

    call MPI_SendRecv (P_Top , SIZE(P_Top ), MPI_REALTYPE, Down, Send_Tag, &
                       P_Down, SIZE(P_Down), MPI_REALTYPE, Top , Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  else
    call MPI_Send     (P_Top , SIZE(P_Top ), MPI_REALTYPE, Top , Send_Tag, &
                       MPI_Communicator,             MPI_Error_Status)
    call MPI_Recv     (P_Down, SIZE(P_Down), MPI_REALTYPE, Down, Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_0D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_1D  (Index, Distributed, M, Q, Direction,  &
                                                MPI_Comm_Type)

  integer  :: Index, N0,N5, NV, L0,L1,L4,L5, M
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag

  integer , dimension(MPI_Status_Size) :: MPI_Status

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer , OPTIONAL :: Direction

  logical            :: Distributed

  REALTYPE, dimension(:,:), TARGET  :: Q
  REALTYPE, dimension(:,:), POINTER :: P_Top, P_Down

  integer  :: M_Processor, Processor_Last

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)

  Nullify (P_Top, P_Down)

  Send_Tag = 5555 ;  Recv_Tag = 5555

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5,:)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1,:)

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1,:)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5,:)

  END SELECT

  if (PRESENT(MPI_Comm_Type)) then
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)
  endif

  if ((PRESENT(MPI_Comm_Type)) .AND. (Status == MPI_CART)) then
    call MPI_CART_SHIFT (MPI_Communicator, Index-1, Top_Or_Down,           &
                                           Top, Down, MPI_Error_Status)

    call MPI_SendRecv (P_Top , SIZE(P_Top ), MPI_REALTYPE, Down, Send_Tag, &
                       P_Down, SIZE(P_Down), MPI_REALTYPE, Top , Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  else
    call MPI_Send     (P_Top , SIZE(P_Top ), MPI_REALTYPE, Top , Send_Tag, &
                       MPI_Communicator,             MPI_Error_Status)
    call MPI_Recv     (P_Down, SIZE(P_Down), MPI_REALTYPE, Down, Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_1D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_2D  (Index, Distributed, M, Q, Direction,  &
                                                MPI_Comm_Type)

  integer  :: Index, N0,N5, M0,M5, NV, L0,L1,L4,L5, M
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag

  integer , dimension(MPI_Status_Size) :: MPI_Status

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer , OPTIONAL :: Direction

  logical            :: Distributed

  REALTYPE, dimension(:,:,:), TARGET  :: Q
  REALTYPE, dimension(:,:,:), POINTER :: P_Top, P_Down

  integer  :: M_Processor, Processor_Last

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  N0 = LBOUND(Q,DIM=1) ; N5 = UBOUND(Q,DIM=1)
  M0 = LBOUND(Q,DIM=2) ; M5 = UBOUND(Q,DIM=2)

  Nullify (P_Top, P_Down)

  Send_Tag = 5555 ;  Recv_Tag = 5555

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

  END SELECT 

  if (PRESENT(MPI_Comm_Type)) then
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)
  endif

  if ((PRESENT(MPI_Comm_Type)) .AND. (Status == MPI_CART)) then
    call MPI_CART_SHIFT (MPI_Communicator, Index-1, Top_Or_Down,           &
                                           Top, Down, MPI_Error_Status)

    call MPI_SendRecv (P_Top , SIZE(P_Top ), MPI_REALTYPE, Down, Send_Tag, &
                       P_Down, SIZE(P_Down), MPI_REALTYPE, Top , Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  else
    call MPI_Send     (P_Top , SIZE(P_Top ), MPI_REALTYPE, Top , Send_Tag, &
                       MPI_Communicator,             MPI_Error_Status)
    call MPI_Recv     (P_Down, SIZE(P_Down), MPI_REALTYPE, Down, Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_2D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_3D  (Index, Distributed, M, Q, Direction,  &
                                                MPI_Comm_Type)

  integer  :: Index, N0,N5, M0,M5, K0,K5, NV, L0,L1,L4,L5, M
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag

  integer , dimension(MPI_Status_Size) :: MPI_Status

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer , OPTIONAL :: Direction

  logical            :: Distributed

  REALTYPE, dimension(:,:,:,:), TARGET  :: Q
  REALTYPE, dimension(:,:,:,:), POINTER :: P_Top, P_Down

  integer  :: M_Processor, Processor_Last

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

  Send_Tag = 5555 ;  Recv_Tag = 5555

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

  END SELECT 

  if (PRESENT(MPI_Comm_Type)) then
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)
  endif

  if ((PRESENT(MPI_Comm_Type)) .AND. (Status == MPI_CART)) then
    call MPI_CART_SHIFT (MPI_Communicator, Index-1, Top_Or_Down,           &
                                           Top, Down, MPI_Error_Status)

    call MPI_SendRecv (P_Top , SIZE(P_Top ), MPI_REALTYPE, Down, Send_Tag, &
                       P_Down, SIZE(P_Down), MPI_REALTYPE, Top , Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  else
    call MPI_Send     (P_Top , SIZE(P_Top ), MPI_REALTYPE, Top , Send_Tag, &
                       MPI_Communicator,             MPI_Error_Status)
    call MPI_Recv     (P_Down, SIZE(P_Down), MPI_REALTYPE, Down, Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_3D
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_0DG (Index, Distributed,            &
                                         N0,N5,                   M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0,N5, L0,L1,L4,L5, M
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag

  integer , dimension(MPI_Status_Size) :: MPI_Status

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer , OPTIONAL :: Direction

  logical            :: Distributed

  REALTYPE, dimension(N0:N5), TARGET  :: Q
  REALTYPE, dimension(:)    , POINTER :: P_Top, P_Down

  integer  :: M_Processor, Processor_Last

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  Nullify (P_Top, P_Down)

  Send_Tag = 5555 ;  Recv_Tag = 5555

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1)

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5)

  END SELECT

  if (PRESENT(MPI_Comm_Type)) then
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)
  endif

  if ((PRESENT(MPI_Comm_Type)) .AND. (Status == MPI_CART)) then
    call MPI_CART_SHIFT (MPI_Communicator, Index-1, Top_Or_Down,           &
                                           Top, Down, MPI_Error_Status)

    call MPI_SendRecv (P_Top , SIZE(P_Top ), MPI_REALTYPE, Down, Send_Tag, &
                       P_Down, SIZE(P_Down), MPI_REALTYPE, Top , Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  else
    call MPI_Send     (P_Top , SIZE(P_Top ), MPI_REALTYPE, Top , Send_Tag, &
                       MPI_Communicator,             MPI_Error_Status)
    call MPI_Recv     (P_Down, SIZE(P_Down), MPI_REALTYPE, Down, Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_0DG
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_1DG (Index, Distributed,            &
                                         N0,N5,               NV, M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0,N5, NV, L0,L1,L4,L5, M
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag

  integer , dimension(MPI_Status_Size) :: MPI_Status

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer , OPTIONAL :: Direction

  logical            :: Distributed

  REALTYPE, dimension(N0:N5,NV), TARGET  :: Q
  REALTYPE, dimension(:,:)     , POINTER :: P_Top, P_Down

  integer  :: M_Processor, Processor_Last

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  Nullify (P_Top, P_Down)

  Send_Tag = 5555 ;  Recv_Tag = 5555

  Top_Or_Down = 1 ; if (PRESENT(Direction)) Top_Or_Down = Direction

  SELECT CASE (Top_Or_Down)
    CASE DEFAULT
      Top  = I_Am + 1 ; if (I_Am ==  Processor_Last) Top  = MPI_PROC_NULL
      Down = I_Am - 1 ; if (I_Am == First_Processor) Down = MPI_PROC_NULL

          L4 = N5-2*M+1 ; L5 = L4+M-1 ; P_Top  => Q(L4:L5,:)
          L0 = N0       ; L1 = L0+M-1 ; P_Down => Q(L0:L1,:)

    CASE (-1)
      Top  = I_Am - 1 ; if (I_Am == First_Processor) Top  = MPI_PROC_NULL
      Down = I_Am + 1 ; if (I_Am ==  Processor_Last) Down = MPI_PROC_NULL

          L0 = N0+M     ; L1 = L0+M-1 ; P_Top  => Q(L0:L1,:)
          L4 = N5-M+1   ; L5 = N5     ; P_Down => Q(L4:L5,:)

  END SELECT

  if (PRESENT(MPI_Comm_Type)) then
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)
  endif

  if ((PRESENT(MPI_Comm_Type)) .AND. (Status == MPI_CART)) then
    call MPI_CART_SHIFT (MPI_Communicator, Index-1, Top_Or_Down,           &
                                           Top, Down, MPI_Error_Status)

    call MPI_SendRecv (P_Top , SIZE(P_Top ), MPI_REALTYPE, Down, Send_Tag, &
                       P_Down, SIZE(P_Down), MPI_REALTYPE, Top , Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  else
    call MPI_Send     (P_Top , SIZE(P_Top ), MPI_REALTYPE, Top , Send_Tag, &
                       MPI_Communicator,             MPI_Error_Status)
    call MPI_Recv     (P_Down, SIZE(P_Down), MPI_REALTYPE, Down, Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_1DG
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_2DG (Index, Distributed,            &
                                         N0,N5, M0,M5,        NV, M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0,N5, M0,M5, NV, L0,L1,L4,L5, M
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag

  integer , dimension(MPI_Status_Size) :: MPI_Status

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer , OPTIONAL :: Direction

  logical            :: Distributed

  REALTYPE, dimension(N0:N5,M0:M5,NV), TARGET  :: Q
  REALTYPE, dimension(:,:,:)         , POINTER :: P_Top, P_Down

  integer  :: M_Processor, Processor_Last

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  Nullify (P_Top, P_Down)

  Send_Tag = 5555 ;  Recv_Tag = 5555

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

  END SELECT 

  if (PRESENT(MPI_Comm_Type)) then
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)
  endif

  if ((PRESENT(MPI_Comm_Type)) .AND. (Status == MPI_CART)) then
    call MPI_CART_SHIFT (MPI_Communicator, Index-1, Top_Or_Down,           &
                                           Top, Down, MPI_Error_Status)

    call MPI_SendRecv (P_Top , SIZE(P_Top ), MPI_REALTYPE, Down, Send_Tag, &
                       P_Down, SIZE(P_Down), MPI_REALTYPE, Top , Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  else
    call MPI_Send     (P_Top , SIZE(P_Top ), MPI_REALTYPE, Top , Send_Tag, &
                       MPI_Communicator,             MPI_Error_Status)
    call MPI_Recv     (P_Down, SIZE(P_Down), MPI_REALTYPE, Down, Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  endif

  Nullify (P_Top, P_Down)

  END Subroutine Exchange_Boundary_Data_2DG
!
! ============================================================================
!
  Subroutine Exchange_Boundary_Data_3DG (Index, Distributed,            &
                                         N0,N5, M0,M5, K0,K5, NV, M, Q, &
                                         Direction, MPI_Comm_Type)

  integer  :: Index, N0,N5, M0,M5, K0,K5, NV, L0,L1,L4,L5, M
  integer  :: Top_Or_Down, Top, Down, Send_Tag, Recv_Tag

  integer , dimension(MPI_Status_Size) :: MPI_Status

  integer , OPTIONAL :: MPI_Comm_Type
  integer            :: MPI_Communicator

  integer , OPTIONAL :: Direction

  logical            :: Distributed

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET  :: Q
  REALTYPE, dimension(:,:,:,:)             , POINTER :: P_Top, P_Down

  integer  :: M_Processor, Processor_Last

  if ((M == 0) .OR. (N_Processor == 1) .OR. (.NOT. Distributed)) RETURN

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)

                                       Processor_Last = Last_Processor 
  if (M_Processor-1 <  Last_Processor) Processor_Last = M_Processor-1

  Nullify (P_Top, P_Down)

  Send_Tag = 5555 ;  Recv_Tag = 5555

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

  END SELECT 

  if (PRESENT(MPI_Comm_Type)) then
    call MPI_TOPO_TEST (MPI_Comm_Type, Status, MPI_Error_Status)
  endif

  if ((PRESENT(MPI_Comm_Type)) .AND. (Status == MPI_CART)) then
    call MPI_CART_SHIFT (MPI_Communicator, Index-1, Top_Or_Down,           &
                                           Top, Down, MPI_Error_Status)

    call MPI_SendRecv (P_Top , SIZE(P_Top ), MPI_REALTYPE, Down, Send_Tag, &
                       P_Down, SIZE(P_Down), MPI_REALTYPE, Top , Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  else
    call MPI_Send     (P_Top , SIZE(P_Top ), MPI_REALTYPE, Top , Send_Tag, &
                       MPI_Communicator,             MPI_Error_Status)
    call MPI_Recv     (P_Down, SIZE(P_Down), MPI_REALTYPE, Down, Recv_Tag, &
                       MPI_Communicator, MPI_Status, MPI_Error_Status)
  endif

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
