#include "pseudopack.h"


MODULE Navier_Stokes_Flux_Boundary_Condition

  USE PseudoPack

implicit NONE

INTERFACE PS_Flux_Boundary_Condition
  MODULE PROCEDURE Flux_BC_0D
  MODULE PROCEDURE Flux_BC_1D
  MODULE PROCEDURE Flux_BC_2D
  MODULE PROCEDURE Flux_BC_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Flux_Boundary_Condition

CONTAINS
!
!======================================================================
!
  Subroutine Flux_BC_0D (N0, N5, N2, N3, Q,                                    &
                         BC_On_Off, BC_Distributed, BC_Order,                  &
                         BC_Method, BC_Symmetry   , BC_Type , BC_Even_or_Odd,  &
                         MPI_Comm_Type)

  integer  :: i, Index
  integer  :: N0, N1, N2, N3, N4, N5

  integer  :: II, N_Vector

  REALTYPE, dimension(N0:N5)               , TARGET  :: Q
  REALTYPE, dimension(:)                   , POINTER :: P0, P1

  integer                             :: BC_On_Off
  integer                             :: BC_Order
  logical                             :: BC_Distributed
  integer                             :: BC_Method
  integer                  , OPTIONAL :: BC_Symmetry
  integer , dimension(2)   , OPTIONAL :: BC_Type
  REALTYPE                 , OPTIONAL :: BC_Even_or_Odd

  integer                             :: Symmetry
  integer , dimension(2)              :: BC_IO_Type
  REALTYPE                            :: IEO

  integer ,                  OPTIONAL :: MPI_Comm_Type
#if defined (PARALLEL_MPI)
  integer                             :: MPI_Communicator
  integer                             :: MPI_Communicator_1D
  integer                             :: N_Dims
  logical , dimension(:), ALLOCATABLE :: Sub_Grid

                              MPI_Communicator = MPI_Comm_World
  if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type
#endif

                            Symmetry = 0 
  if (PRESENT(BC_Symmetry)) Symmetry = BC_Symmetry

                        BC_IO_Type = 1
  if (PRESENT(BC_Type)) BC_IO_Type = BC_Type

                                                          IEO = ONE
  if (PRESENT(BC_Symmetry) .AND. PRESENT(BC_Even_or_Odd)) IEO = BC_Even_or_Odd

  if (BC_On_Off /= 1) RETURN

#if defined (PARALLEL_MPI)
    call MPI_COMM_DUP (MPI_Communicator, MPI_Communicator_1D, MPI_Error_Status)
#endif

    N1 = N2-1 ; N4 = N3+1

    SELECT CASE (BC_Method)
      CASE DEFAULT                  ! Periodical BC
#if defined (PGI704)
        call Flux_BC_0D_Periodic (Q, P0, P1)
#else
        call Flux_BC_0D_Periodic 
#endif
 
      CASE (1)                      ! Freestream BC
        call Flux_BC_0D_Freestream

    END SELECT 

#if defined (PARALLEL_MPI)
    call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
#endif

  CONTAINS

#if defined (PGI704)
    Subroutine Flux_BC_0D_Periodic (Q, P0, P1)

    REALTYPE, dimension(N0:N5)               , TARGET  :: Q
    REALTYPE, dimension(:)                   , POINTER :: P0, P1
#else
    Subroutine Flux_BC_0D_Periodic
#endif

#if defined (PARALLEL_MPI)
    integer                              :: Send_Tag, Recv_Tag
    integer , dimension(MPI_Status_Size) :: MPI_Status

    if (BC_Distributed) then

      if (Symmetry == 0) then
        Nullify (P0, P1) ; Send_Tag = 5555 ; Recv_Tag = 5555

        if (I_Am ==  Last_Processor) then
          N_Vector = ABS(N0-N1)+1 ; II = N3-N_Vector+1
          P0 => Q(II:N3)

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE, First_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am == First_Processor) then
          P1 => Q(N0:N1)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE,  Last_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1) ; Send_Tag = 6666 ; Recv_Tag = 6666

        if (I_Am == First_Processor) then
          N_Vector = ABS(N5-N4)+1 ; II = N1+N_Vector
          P0 => Q(N2:II)

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE,  Last_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am ==  Last_Processor) then
          P1 => Q(N4:N5)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE, First_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1)
      else
        if (I_Am == First_Processor) then
            do i = 0,ABS(N1-N0)
              Q(N1-i) = IEO*Q(N2+i)
            enddo
        endif

        if (I_Am ==  Last_Processor) then
            do i = 0,ABS(N5-N4)
              Q(N4+i) = IEO*Q(N3-i)
            enddo
        endif
      endif

    else
#endif

      Nullify (P0, P1)

      if (Symmetry == 0) then
          N_Vector = ABS(N0-N1)+1 ; II = N3-N_Vector+1
          P1 => Q(N0:N1) ; P0 => Q(II:N3) ; P1 = P0

          N_Vector = ABS(N5-N4)+1 ; II = N2+N_Vector-1
          P1 => Q(N4:N5) ; P0 => Q(N2:II) ; P1 = P0
      else
            do i = 0,ABS(N1-N0)
              Q(N1-i) = IEO*Q(N2+i)
            enddo
            do i = 0,ABS(N5-N4)
              Q(N4+i) = IEO*Q(N3-i)
            enddo
      endif

#if defined (PARALLEL_MPI)
    endif
#endif

      Nullify (P0, P1)

    END Subroutine Flux_BC_0D_Periodic
!
!======================================================================
!
    Subroutine Flux_BC_0D_Freestream

#if defined (PARALLEL_MPI)
      if (I_Am == First_Processor) then
#endif
        SELECT CASE (BC_IO_Type(1))
          CASE (-1)     ! Outflow
            do i = N0,N1
              Q(i) = ZERO
            enddo
 
          CASE ( 0)     ! No BC

          CASE DEFAULT  ! Freestream
            do i = N0,N1
              Q(i) = ZERO
            enddo

        END SELECT

#if defined (PARALLEL_MPI)
      endif

      if (I_Am ==  Last_Processor) then
#endif
        if (Symmetry == 1) then

            N_Vector = ABS(N5-N4)+1 ; II = N3-N_Vector+1
              Q(N4:N5) = IEO*Q(N3:II:-1)

        else

          SELECT CASE (BC_IO_Type(2))
            CASE (-1)     !  Outflow
              do i = N4,N5
                Q(i) = ZERO
              enddo

            CASE ( 0)     ! No BC

            CASE DEFAULT  ! Freestream
              do i = N4,N5
                Q(i) = ZERO
              enddo

          END SELECT

      endif

#if defined (PARALLEL_MPI)
    endif
#endif

    END Subroutine Flux_BC_0D_Freestream

  END Subroutine Flux_BC_0D 
!
!======================================================================
!
  Subroutine Flux_BC_1D (N0, N5, N2, N3,  NV, Q,                               &
                         BC_On_Off, BC_Distributed, BC_Order,                  &
                         BC_Method, BC_Symmetry   , BC_Type , BC_Even_or_Odd,  &
                         MPI_Comm_Type)

  integer  :: i, n, Index
  integer  :: N0, N1, N2, N3, N4, N5, NV

  integer  :: II, N_Vector

  REALTYPE, dimension(N0:N5,NV)            , TARGET  :: Q
  REALTYPE, dimension(:,:)                 , POINTER :: P0, P1

  integer , dimension(1)              :: BC_On_Off
  integer , dimension(1)              :: BC_Order
  logical , dimension(1)              :: BC_Distributed
  integer , dimension(1)              :: BC_Method
  integer , dimension(1)   , OPTIONAL :: BC_Symmetry
  integer , dimension(3,2) , OPTIONAL :: BC_Type
  REALTYPE, dimension(3,NV), OPTIONAL :: BC_Even_or_Odd

  integer , dimension(1)              :: Symmetry
  integer , dimension(3,2)            :: BC_IO_Type
  REALTYPE, dimension(3,NV)           :: IEO

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

                            Symmetry = 0 
  if (PRESENT(BC_Symmetry)) Symmetry = BC_Symmetry

                        BC_IO_Type = 1
  if (PRESENT(BC_Type)) BC_IO_Type = BC_Type

                                                          IEO = ONE
  if (PRESENT(BC_Symmetry) .AND. PRESENT(BC_Even_or_Odd)) IEO = BC_Even_or_Odd

  do Index = 1,SIZE(BC_On_Off)
    if (BC_On_Off(Index) /= 1) CYCLE

#if defined (PARALLEL_MPI)
    call PS_MPI_Processor_Local_Info (I_Am, Last_Processor, &
                                            Index, MPI_Communicator)

    Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

    call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                   MPI_Error_Status)
#endif

    N1 = N2-1 ; N4 = N3+1

    SELECT CASE (BC_Method(Index))
      CASE DEFAULT                  ! Periodical BC
#if defined (PGI704)
        call Flux_BC_1D_x_Periodic (Q, P0, P1)
#else
        call Flux_BC_1D_x_Periodic 
#endif
 
      CASE (1)                      ! Freestream BC
        call Flux_BC_1D_x_Freestream

    END SELECT 

#if defined (PARALLEL_MPI)
    call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
  enddo

  DEALLOCATE (Sub_Grid)

  call PS_MPI_Processor (I_Am, Last_Processor)
#else
  enddo
#endif

  CONTAINS

#if defined (PGI704)
    Subroutine Flux_BC_1D_x_Periodic (Q, P0, P1)

    REALTYPE, dimension(N0:N5,NV)            , TARGET  :: Q
    REALTYPE, dimension(:,:)                 , POINTER :: P0, P1
#else
    Subroutine Flux_BC_1D_x_Periodic
#endif

#if defined (PARALLEL_MPI)
    integer                              :: Send_Tag, Recv_Tag
    integer , dimension(MPI_Status_Size) :: MPI_Status

    if (BC_Distributed(Index)) then

      if (Symmetry(Index) == 0) then
        Nullify (P0, P1) ; Send_Tag = 5555 ; Recv_Tag = 5555

        if (I_Am ==  Last_Processor) then
          N_Vector = ABS(N0-N1)+1 ; II = N3-N_Vector+1
          P0 => Q(II:N3,:)

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE, First_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am == First_Processor) then
          P1 => Q(N0:N1,:)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE,  Last_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1) ; Send_Tag = 6666 ; Recv_Tag = 6666

        if (I_Am == First_Processor) then
          N_Vector = ABS(N5-N4)+1 ; II = N2+N_Vector-1
          P0 => Q(N2:II,:)

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE,  Last_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am ==  Last_Processor) then
          P1 => Q(N4:N5,:)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE, First_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1)
      else
        if (I_Am == First_Processor) then
          do n = 1,NV
            do i = 0,ABS(N1-N0)
              Q(N1-i,n) = IEO(Index,n)*Q(N2+i,n)
            enddo
          enddo
        endif

        if (I_Am ==  Last_Processor) then
          do n = 1,NV
            do i = 0,ABS(N5-N4)
              Q(N4+i,n) = IEO(Index,n)*Q(N3-i,n)
            enddo
          enddo
        endif
      endif

    else
#endif

      Nullify (P0, P1)

      if (Symmetry(Index) == 0) then
          N_Vector = ABS(N0-N1)+1 ; II = N3-N_Vector+1
          P1 => Q(N0:N1,:) ; P0 => Q(II:N3,:) ; P1 = P0

          N_Vector = ABS(N5-N4)+1 ; II = N2+N_Vector-1
          P1 => Q(N4:N5,:) ; P0 => Q(N2:II,:) ; P1 = P0
      else
          do n = 1,NV
            do i = 0,ABS(N1-N0)
              Q(N1-i,n) = IEO(Index,n)*Q(N2+i,n)
            enddo
            do i = 0,ABS(N5-N4)
              Q(N4+i,n) = IEO(Index,n)*Q(N3-i,n)
            enddo
          enddo
      endif

#if defined (PARALLEL_MPI)
    endif
#endif

      Nullify (P0, P1)

    END Subroutine Flux_BC_1D_x_Periodic
!
!======================================================================
!
    Subroutine Flux_BC_1D_x_Freestream

#if defined (PARALLEL_MPI)
      if (I_Am == First_Processor) then
#endif
        SELECT CASE (BC_IO_Type(Index,1))
          CASE (-1)     ! Outflow
            do i = N0,N1
              Q(i,:) = ZERO
            enddo
 
          CASE ( 0)     ! No BC

          CASE DEFAULT  ! Freestream
            do i = N0,N1
              Q(i,:) = ZERO
            enddo

        END SELECT

#if defined (PARALLEL_MPI)
      endif

      if (I_Am ==  Last_Processor) then
#endif
        if (Symmetry(Index) == 1) then

            N_Vector = ABS(N5-N4)+1 ; II = N3-N_Vector+1
            do n = 1,NV
              Q(N4:N5,n) = IEO(Index,n)*Q(N3:II:-1,n)
            enddo

        else

          SELECT CASE (BC_IO_Type(Index,2))
            CASE (-1)     ! Outflow
              do i = N4,N5
                Q(i,:) = ZERO
              enddo

            CASE ( 0)     ! No BC

            CASE DEFAULT  ! Freestream
              do i = N4,N5
                Q(i,:) = ZERO
              enddo

          END SELECT

        endif

#if defined (PARALLEL_MPI)
      endif
#endif

    END Subroutine Flux_BC_1D_x_Freestream

  END Subroutine Flux_BC_1D 
!
!======================================================================
!
  Subroutine Flux_BC_2D (N0, N5, N2, N3,                                       &
                         M0, M5, M2, M3,  NV, Q,                               &
                         BC_On_Off, BC_Distributed, BC_Order,                  &
                         BC_Method, BC_Symmetry   , BC_Type , BC_Even_or_Odd,  &
                         MPI_Comm_Type)

  integer  :: i, j, n, Index

  integer  :: N0, N1, N2, N3, N4, N5
  integer  :: M0, M1, M2, M3, M4, M5,  NV

  integer  :: II, N_Vector

  REALTYPE, dimension(N0:N5,M0:M5,NV)      , TARGET  :: Q
  REALTYPE, dimension(:,:,:)               , POINTER :: P0, P1

  integer , dimension(2)              :: BC_On_Off
  integer , dimension(2)              :: BC_Order
  logical , dimension(2)              :: BC_Distributed
  integer , dimension(2)              :: BC_Method
  integer , dimension(2)   , OPTIONAL :: BC_Symmetry
  integer , dimension(3,2) , OPTIONAL :: BC_Type
  REALTYPE, dimension(3,NV), OPTIONAL :: BC_Even_or_Odd

  integer , dimension(2)              :: Symmetry
  integer , dimension(3,2)            :: BC_IO_Type
  REALTYPE, dimension(3,NV)           :: IEO

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

                            Symmetry = 0 
  if (PRESENT(BC_Symmetry)) Symmetry = BC_Symmetry

                        BC_IO_Type = 1
  if (PRESENT(BC_Type)) BC_IO_Type = BC_Type

                                                          IEO = ONE
  if (PRESENT(BC_Symmetry) .AND. PRESENT(BC_Even_or_Odd)) IEO = BC_Even_or_Odd

  do Index = 1,SIZE(BC_On_Off)
    if (BC_On_Off(Index) /= 1) CYCLE

#if defined (PARALLEL_MPI)
    call PS_MPI_Processor_Local_Info (I_Am, Last_Processor, &
                                            Index, MPI_Communicator)

    Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

    call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                   MPI_Error_Status)
#endif

    N1 = N2-1 ; N4 = N3+1
    M1 = M2-1 ; M4 = M3+1

    SELECT CASE (BC_Method(Index))
      CASE DEFAULT                  ! Periodical BC
#if defined (PGI704)
        if (Index == 1) call Flux_BC_2D_x_Periodic (Q, P0, P1)
        if (Index == 2) call Flux_BC_2D_y_Periodic (Q, P0, P1)
#else
        if (Index == 1) call Flux_BC_2D_x_Periodic 
        if (Index == 2) call Flux_BC_2D_y_Periodic 
#endif
 
      CASE (1)                      ! Freestream BC
        if (Index == 1) call Flux_BC_2D_x_Freestream
        if (Index == 2) call Flux_BC_2D_y_Freestream

    END SELECT 

#if defined (PARALLEL_MPI)
    call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
  enddo

  DEALLOCATE (Sub_Grid)

  call PS_MPI_Processor (I_Am, Last_Processor)
#else
  enddo
#endif

  CONTAINS

#if defined (PGI704)
    Subroutine Flux_BC_2D_x_Periodic (Q, P0, P1)

    REALTYPE, dimension(N0:N5,M0:M5,NV)      , TARGET  :: Q
    REALTYPE, dimension(:,:,:)               , POINTER :: P0, P1
#else
    Subroutine Flux_BC_2D_x_Periodic
#endif

#if defined (PARALLEL_MPI)
    integer                              :: Send_Tag, Recv_Tag
    integer , dimension(MPI_Status_Size) :: MPI_Status

    if (BC_Distributed(Index)) then

      if (Symmetry(Index) == 0) then

        Nullify (P0, P1) ; Send_Tag = 5555 ; Recv_Tag = 5555

        if (I_Am ==  Last_Processor) then
          N_Vector = ABS(N0-N1)+1 ; II = N3-N_Vector+1
          P0 => Q(II:N3,:,:)

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE, First_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am == First_Processor) then
          P1 => Q(N0:N1,:,:)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE,  Last_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1) ; Send_Tag = 6666 ; Recv_Tag = 6666

        if (I_Am == First_Processor) then
          N_Vector = ABS(N5-N4)+1 ; II = N2+N_Vector-1
          P0 => Q(N2:II,:,:)

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE,  Last_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am ==  Last_Processor) then
          P1 => Q(N4:N5,:,:)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE, First_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1)
      else
        if (I_Am == First_Processor) then
          do n = 1,NV
            do i = 0,ABS(N1-N0)
              Q(N1-i,:,n) = IEO(Index,n)*Q(N2+i,:,n)
            enddo
          enddo
        endif

        if (I_Am ==  Last_Processor) then
          do n = 1,NV
            do i = 0,ABS(N5-N4)
              Q(N4+i,:,n) = IEO(Index,n)*Q(N3-i,:,n)
            enddo
          enddo
        endif
      endif

    else
#endif

      Nullify (P0, P1)

      if (Symmetry(Index) == 0) then
          N_Vector = ABS(N0-N1)+1 ; II = N3-N_Vector+1
          P1 => Q(N0:N1,:,:) ; P0 => Q(II:N3,:,:) ; P1 = P0

          N_Vector = ABS(N5-N4)+1 ; II = N2+N_Vector-1
          P1 => Q(N4:N5,:,:) ; P0 => Q(N2:II,:,:) ; P1 = P0
      else
          do n = 1,NV
            do i = 0,ABS(N1-N0)
              Q(N1-i,:,n) = IEO(Index,n)*Q(N2+i,:,n)
            enddo
            do i = 0,ABS(N5-N4)
              Q(N4+i,:,n) = IEO(Index,n)*Q(N3-i,:,n)
            enddo
          enddo
      endif

#if defined (PARALLEL_MPI)
    endif
#endif

      Nullify (P0, P1)

    END Subroutine Flux_BC_2D_x_Periodic
!
!======================================================================
!
    Subroutine Flux_BC_2D_x_Freestream

#if defined (PARALLEL_MPI)
      if (I_Am == First_Processor) then
#endif
        SELECT CASE (BC_IO_Type(Index,1))
          CASE (-1)     ! Outflow
            do i = N0,N1
              Q(i,:,:) = ZERO
            enddo
 
          CASE ( 0)     ! No BC

          CASE DEFAULT  ! Freestream
            do i = N0,N1
              Q(i,:,:) = ZERO
            enddo

        END SELECT

#if defined (PARALLEL_MPI)
      endif

      if (I_Am ==  Last_Processor) then
#endif
        if (Symmetry(Index) == 1) then

            N_Vector = ABS(N5-N4)+1 ; II = N3-N_Vector+1
            do n = 1,NV
              Q(N4:N5,:,n) = IEO(Index,n)*Q(N3:II:-1,:,n)
            enddo

        else

          SELECT CASE (BC_IO_Type(Index,2))
            CASE (-1)     ! Outflow
              do i = N4,N5
                Q(i,:,:) = ZERO
              enddo

            CASE ( 0)     ! No BC

            CASE DEFAULT  ! Freestream
              do i = N4,N5
                Q(i,:,:) = ZERO
              enddo

          END SELECT

        endif

#if defined (PARALLEL_MPI)
      endif
#endif

    END Subroutine Flux_BC_2D_x_Freestream
!
!======================================================================
!
#if defined (PGI704)
    Subroutine Flux_BC_2D_y_Periodic (Q, P0, P1)

    REALTYPE, dimension(N0:N5,M0:M5,NV)      , TARGET  :: Q
    REALTYPE, dimension(:,:,:)               , POINTER :: P0, P1
#else
    Subroutine Flux_BC_2D_y_Periodic 
#endif

#if defined (PARALLEL_MPI)
    integer                              :: Send_Tag, Recv_Tag
    integer , dimension(MPI_Status_Size) :: MPI_Status

    if (BC_Distributed(Index)) then

      if (Symmetry(Index) == 0) then
        Nullify (P0, P1) ; Send_Tag = 5555 ; Recv_Tag = 5555

        if (I_Am ==  Last_Processor) then
          N_Vector = ABS(M0-M1)+1 ; II = M3-N_Vector+1
          P0 => Q(:,II:M3,:)

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE, First_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am == First_Processor) then
          P1 => Q(:,M0:M1,:)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE,  Last_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1) ; Send_Tag = 6666 ; Recv_Tag = 6666

        if (I_Am == First_Processor) then
          N_Vector = ABS(M5-M4)+1 ; II = M2+N_Vector-1
          P0 => Q(:,M2:II,:)

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE,  Last_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am ==  Last_Processor) then
          P1 => Q(:,M4:M5,:)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE, First_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1)
      else
        if (I_Am == First_Processor) then
          do n = 1,NV
            do j = 0,ABS(M1-M0)
              Q(:,M1-j,n) = IEO(Index,n)*Q(:,M2+j,n)
            enddo
          enddo
        endif

        if (I_Am ==  Last_Processor) then
          do n = 1,NV
            do j = 0,ABS(M5-M4)
              Q(:,M4+j,n) = IEO(Index,n)*Q(:,M3-j,n)
            enddo
          enddo
        endif
      endif

    else
#endif

      Nullify (P0, P1)

      if (Symmetry(Index) == 0) then
          N_Vector = ABS(M0-M1)+1 ; II = M3-N_Vector+1
          P1 => Q(:,M0:M1,:) ; P0 => Q(:,II:M3,:) ; P1 = P0

          N_Vector = ABS(M5-M4)+1 ; II = M2+N_Vector-1
          P1 => Q(:,M4:M5,:) ; P0 => Q(:,M2:II,:) ; P1 = P0
      else
          do n = 1,NV
            do j = 0,ABS(M1-M0)
              Q(:,M1-j,n) = IEO(Index,n)*Q(:,M2+j,n)
            enddo
            do j = 0,ABS(M5-M4)
              Q(:,M4+j,n) = IEO(Index,n)*Q(:,M3-j,n)
            enddo
          enddo
      endif

#if defined (PARALLEL_MPI)
    endif
#endif

      Nullify (P0, P1)

    END Subroutine Flux_BC_2D_y_Periodic
!
!======================================================================
!
    Subroutine Flux_BC_2D_y_Freestream

#if defined (PARALLEL_MPI)
      if (I_Am == First_Processor) then
#endif
        SELECT CASE (BC_IO_Type(Index,1))
          CASE (-1)      ! Outflow
            do j = M0,M1
              Q(:,j,:) = ZERO
            enddo
 
          CASE ( 0)     ! No BC

          CASE DEFAULT  ! Freestream
            do j = M0,M1
              Q(:,j,:) = ZERO
            enddo

        END SELECT

#if defined (PARALLEL_MPI)
      endif

      if (I_Am ==  Last_Processor) then
#endif
        if (Symmetry(Index) == 1) then

            N_Vector = ABS(M5-M4)+1 ; II = M3-N_Vector+1
            do n = 1,NV
              Q(:,M4:M5,n) = IEO(Index,n)*Q(:,M3:II:-1,n)
            enddo

        else

          SELECT CASE (BC_IO_Type(Index,2))
            CASE (-1)     ! Outflow
              do j = M4,M5
                Q(:,j,:) = ZERO
              enddo

            CASE ( 0)     ! No BC

            CASE DEFAULT  ! Freestream
              do j = M4,M5
                Q(:,j,:) = ZERO
              enddo

          END SELECT

        endif

#if defined (PARALLEL_MPI)
      endif
#endif

    END Subroutine Flux_BC_2D_y_Freestream

  END Subroutine Flux_BC_2D 
!
!======================================================================
!
  Subroutine Flux_BC_3D (N0, N5, N2, N3,                                       &
                         M0, M5, M2, M3,                                       &
                         K0, K5, K2, K3,  NV, Q,                               &
                         BC_On_Off, BC_Distributed, BC_Order,                  &
                         BC_Method, BC_Symmetry   , BC_Type , BC_Even_or_Odd,  &
                         MPI_Comm_Type)

  integer  :: i, j, k, n, Index
  integer  :: N0, N1, N2, N3, N4, N5
  integer  :: M0, M1, M2, M3, M4, M5
  integer  :: K0, K1, K2, K3, K4, K5,  NV

  integer  :: II, N_Vector

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET  :: Q
  REALTYPE, dimension(:,:,:,:)             , POINTER :: P0, P1

  integer , dimension(3)              :: BC_On_Off
  integer , dimension(3)              :: BC_Order
  logical , dimension(3)              :: BC_Distributed
  integer , dimension(3)              :: BC_Method
  integer , dimension(3)   , OPTIONAL :: BC_Symmetry
  integer , dimension(3,2) , OPTIONAL :: BC_Type
  REALTYPE, dimension(3,NV), OPTIONAL :: BC_Even_or_Odd

  integer , dimension(3)              :: Symmetry
  integer , dimension(3,2)            :: BC_IO_Type
  REALTYPE, dimension(3,NV)           :: IEO

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

                            Symmetry = 0 
  if (PRESENT(BC_Symmetry)) Symmetry = BC_Symmetry

                        BC_IO_Type = 1
  if (PRESENT(BC_Type)) BC_IO_Type = BC_Type

                                                          IEO = ONE
  if (PRESENT(BC_Symmetry) .AND. PRESENT(BC_Even_or_Odd)) IEO = BC_Even_or_Odd

  do Index = 1,SIZE(BC_On_Off)
    if (BC_On_Off(Index) /= 1) CYCLE

#if defined (PARALLEL_MPI)
    call PS_MPI_Processor_Local_Info (I_Am, Last_Processor, &
                                            Index, MPI_Communicator)

    Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

    call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                   MPI_Error_Status)
#endif

    N1 = N2-1 ; N4 = N3+1
    M1 = M2-1 ; M4 = M3+1
    K1 = K2-1 ; K4 = K3+1

    SELECT CASE (BC_Method(Index))
      CASE DEFAULT                  ! Periodical BC
#if defined (PGI704)
        if (Index == 1) call Flux_BC_3D_x_Periodic (Q)
        if (Index == 2) call Flux_BC_3D_y_Periodic (Q)
        if (Index == 3) call Flux_BC_3D_z_Periodic (Q)
#else
        if (Index == 1) call Flux_BC_3D_x_Periodic 
        if (Index == 2) call Flux_BC_3D_y_Periodic 
        if (Index == 3) call Flux_BC_3D_z_Periodic 
#endif
 
      CASE (1)                      ! Freestream BC
        if (Index == 1) call Flux_BC_3D_x_Freestream
        if (Index == 2) call Flux_BC_3D_y_Freestream
        if (Index == 3) call Flux_BC_3D_z_Freestream

    END SELECT 

#if defined (PARALLEL_MPI)
    call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
  enddo

  DEALLOCATE (Sub_Grid)

  call PS_MPI_Processor (I_Am, Last_Processor)
#else
  enddo
#endif

  CONTAINS

#if defined (PGI704)
    Subroutine Flux_BC_3D_x_Periodic (Q)

    REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET  :: Q
#else
    Subroutine Flux_BC_3D_x_Periodic
#endif

#if defined (PARALLEL_MPI)
    integer                              :: Send_Tag, Recv_Tag
    integer , dimension(MPI_Status_Size) :: MPI_Status

    if (BC_Distributed(Index)) then

      if (Symmetry(Index) == 0) then

        Nullify (P0, P1) ; Send_Tag = 5555 ; Recv_Tag = 5555

        if (I_Am ==  Last_Processor) then
          N_Vector = ABS(N0-N1)+1 ; II = N3-N_Vector+1
          P0 => Q(II:N3,:,:,:) 

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE, First_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am == First_Processor) then
          P1 => Q(N0:N1,:,:,:)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE,  Last_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1) ; Send_Tag = 6666 ; Recv_Tag = 6666

        if (I_Am == First_Processor) then
          N_Vector = ABS(N5-N4)+1 ; II = N1+N_Vector
          P0 => Q(N2:II,:,:,:) 

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE,  Last_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am ==  Last_Processor) then
          P1 => Q(N4:N5,:,:,:) 

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE, First_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1)
      else
        if (I_Am == First_Processor) then
          do n = 1,NV
            do i = 0,ABS(N1-N0)
              Q(N1-i,:,:,n) = IEO(Index,n)*Q(N2+i,:,:,n)
            enddo
          enddo
        endif

        if (I_Am ==  Last_Processor) then
          do n = 1,NV
            do i = 0,ABS(N5-N4)
              Q(N4+i,:,:,n) = IEO(Index,n)*Q(N3-i,:,:,n)
            enddo
          enddo
        endif
      endif

    else
#endif

      Nullify (P0, P1)

      if (Symmetry(Index) == 0) then
          N_Vector = ABS(N0-N1)+1 ; II = N3-N_Vector+1
          P1 => Q(N0:N1,:,:,:) ; P0 => Q(II:N3,:,:,:) ; P1 = P0
            
          N_Vector = ABS(N5-N4)+1 ; II = N1+N_Vector
          P1 => Q(N4:N5,:,:,:) ; P0 => Q(N2:II,:,:,:) ; P1 = P0
      else
          do n = 1,NV
            do i = 0,ABS(N1-N0)
              Q(N1-i,:,:,n) = IEO(Index,n)*Q(N2+i,:,:,n)
            enddo
            do i = 0,ABS(N5-N4)
              Q(N4+i,:,:,n) = IEO(Index,n)*Q(N3-i,:,:,n)
            enddo
          enddo
      endif

#if defined (PARALLEL_MPI)
    endif
#endif

      Nullify (P0, P1)

    END Subroutine Flux_BC_3D_x_Periodic
!
!======================================================================
!
    Subroutine Flux_BC_3D_x_Freestream

#if defined (PARALLEL_MPI)
      if (I_Am == First_Processor) then
#endif
        SELECT CASE (BC_IO_Type(Index,1))
          CASE (-1)     ! Outflow
            do i = N0,N1
              Q(i,:,:,:) = ZERO
            enddo
 
          CASE ( 0)     ! No BC

          CASE DEFAULT  ! Freestream
            do i = N0,N1
              Q(i,:,:,:) = ZERO
            enddo

        END SELECT

#if defined (PARALLEL_MPI)
      endif

      if (I_Am ==  Last_Processor) then
#endif
        if (Symmetry(Index) == 1) then

            N_Vector = ABS(N5-N4)+1 ; II = N3-N_Vector+1
            do n = 1,NV
              Q(N4:N5,:,:,n) = IEO(Index,n)*Q(N3:II:-1,:,:,n)
            enddo

        else

          SELECT CASE (BC_IO_Type(Index,2))
            CASE (-1)     ! Outflow
              do i = N4,N5
                Q(i,:,:,:) = ZERO
              enddo

            CASE ( 0)     ! No BC

            CASE DEFAULT  ! Freestream
              do i = N4,N5
                Q(i,:,:,:) = ZERO
              enddo

          END SELECT

        endif

#if defined (PARALLEL_MPI)
      endif
#endif

    END Subroutine Flux_BC_3D_x_Freestream
!
!======================================================================
!
#if defined (PGI704)
    Subroutine Flux_BC_3D_y_Periodic (Q)

    REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET  :: Q
#else
    Subroutine Flux_BC_3D_y_Periodic
#endif

#if defined (PARALLEL_MPI)
    integer                              :: Send_Tag, Recv_Tag
    integer , dimension(MPI_Status_Size) :: MPI_Status

    if ((BC_Distributed(Index)) .AND. (First_Processor /= Last_Processor)) then

      if (Symmetry(Index) == 0) then

        Nullify (P0, P1) ; Send_Tag = 5555 ; Recv_Tag = 5555

        if (I_Am ==  Last_Processor) then
          N_Vector = ABS(M0-M1)+1 ; II = M3-N_Vector+1 
          P0 => Q(:,II:M3,:,:) 

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE, First_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am == First_Processor) then
          P1 => Q(:,M0:M1,:,:)

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE,  Last_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1) ; Send_Tag = 6666 ; Recv_Tag = 6666

        if (I_Am == First_Processor) then
          N_Vector = ABS(M5-M4)+1 ; II = M2+N_Vector-1
          P0 => Q(:,M2:II,:,:) 

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE,  Last_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am ==  Last_Processor) then
          P1 => Q(:,M4:M5,:,:) 

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE, First_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1)
      else
        if (I_Am == First_Processor) then
          do n = 1,NV
            do j = 0,ABS(M1-M0)
              Q(:,M1-j,:,n) = IEO(Index,n)*Q(:,M2+j,:,n)
            enddo
          enddo
        endif

        if (I_Am ==  Last_Processor) then
          do n = 1,NV
            do j = 0,ABS(M5-M4)
              Q(:,M4+j,:,n) = IEO(Index,n)*Q(:,M3-j,:,n)
            enddo
          enddo
        endif
      endif

    else
#endif

      Nullify (P0, P1)

      if (Symmetry(Index) == 0) then
          N_Vector = ABS(M0-M1)+1 ; II = M3-N_Vector+1
          P1 => Q(:,M0:M1,:,:) ; P0 => Q(:,II:M3,:,:) ; P1 = P0
            
          N_Vector = ABS(M5-M4)+1 ; II = M2+N_Vector-1
          P1 => Q(:,M4:M5,:,:) ; P0 => Q(:,M2:II,:,:) ; P1 = P0
      else
          do n = 1,NV
            do j = 0,ABS(M1-M0)
              Q(:,M1-j,:,n) = IEO(Index,n)*Q(:,M2+j,:,n)
            enddo
            do j = 0,ABS(M5-M4)
              Q(:,M4+j,:,n) = IEO(Index,n)*Q(:,M3-j,:,n)
            enddo
          enddo
      endif

#if defined (PARALLEL_MPI)
    endif
#endif

      Nullify (P0, P1)

    END Subroutine Flux_BC_3D_y_Periodic
!
!======================================================================
!
    Subroutine Flux_BC_3D_y_Freestream

#if defined (PARALLEL_MPI)
      if (I_Am == First_Processor) then
#endif
        SELECT CASE (BC_IO_Type(Index,1))
          CASE (-1)     ! Outflow
            do j = M0,M1
              Q(:,j,:,:) = ZERO
            enddo
 
          CASE ( 0)     ! No BC

          CASE DEFAULT  ! Freestream
            do j = M0,M1
              Q(:,j,:,:) = ZERO
            enddo

        END SELECT

#if defined (PARALLEL_MPI)
      endif

      if (I_Am ==  Last_Processor) then
#endif
        if (Symmetry(Index) == 1) then

            N_Vector = ABS(M5-M4)+1 ; II = M3-N_Vector+1
            do n = 1,NV
              Q(:,M4:M5,:,n) = IEO(Index,n)*Q(:,M3:II:-1,:,n)
            enddo

        else

          SELECT CASE (BC_IO_Type(Index,2))
            CASE (-1)     ! Outflow
              do j = M4,M5
                Q(:,j,:,:) = ZERO
              enddo

            CASE ( 0)     ! No BC

            CASE DEFAULT  ! Freestream
              do j = M4,M5
                Q(:,j,:,:) = ZERO
              enddo

          END SELECT

        endif

#if defined (PARALLEL_MPI)
      endif
#endif

    END Subroutine Flux_BC_3D_y_Freestream
!
!======================================================================
!
#if defined (PGI704)
    Subroutine Flux_BC_3D_z_Periodic (Q)

    REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV), TARGET  :: Q
#else
    Subroutine Flux_BC_3D_z_Periodic
#endif

#if defined (PARALLEL_MPI)
    integer                              :: Send_Tag, Recv_Tag
    integer , dimension(MPI_Status_Size) :: MPI_Status

    if (BC_Distributed(Index)) then

      if (Symmetry(Index) == 0) then

        Nullify (P0, P1) ; Send_Tag = 5555 ; Recv_Tag = 5555

        if (I_Am ==  Last_Processor) then
          N_Vector = ABS(K0-K1)+1 ; II = K3-N_Vector+1
          P0 => Q(:,:,II:K3,:) 

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE, First_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am == First_Processor) then
          P1 => Q(:,:,K0:K1,:) 

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE,  Last_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1) ; Send_Tag = 6666 ; Recv_Tag = 6666

        if (I_Am == First_Processor) then
          N_Vector = ABS(K5-K4)+1 ; II = K2+N_Vector-1
          P0 => Q(:,:,K2:II,:)

          call MPI_SEND (P0, SIZE(P0), MPI_REALTYPE,  Last_Processor, Send_Tag,&
                             MPI_Communicator_1D,             MPI_Error_Status)
        endif
        if (I_Am ==  Last_Processor) then
          P1 => Q(:,:,K4:K5,:) 

          call MPI_RECV (P1, SIZE(P1), MPI_REALTYPE, First_Processor, Recv_Tag,&
                             MPI_Communicator_1D, MPI_Status, MPI_Error_Status)
        endif

        Nullify (P0, P1)
      else
        if (I_Am == First_Processor) then
          do n = 1,NV
            do k = 0,ABS(K1-K0)
              Q(:,:,K1-k,n) = IEO(Index,n)*Q(:,:,K2+k,n)
            enddo
          enddo
        endif

        if (I_Am ==  Last_Processor) then
          do n = 1,NV
            do k = 0,ABS(K5-K4)
              Q(:,:,K4+k,n) = IEO(Index,n)*Q(:,:,K3-k,n)
            enddo
          enddo
        endif
      endif

    else
#endif

      Nullify (P0, P1)

      if (Symmetry(Index) == 0) then
          N_Vector = ABS(K0-K1)+1 ; II = K3-N_Vector+1
          P1 => Q(:,:,K0:K1,:) ; P0 => Q(:,:,II:K3,:) ; P1 = P0
            
          N_Vector = ABS(K5-K4)+1 ; II = K2+N_Vector-1
          P1 => Q(:,:,K4:K5,:) ; P0 => Q(:,:,K2:II,:) ; P1 = P0
      else
          do n = 1,NV
            do k = 0,ABS(K1-K0)
              Q(:,:,K1-k,n) = IEO(Index,n)*Q(:,:,K2+k,n)
            enddo
            do k = 0,ABS(K5-K4)
              Q(:,:,K4+k,n) = IEO(Index,n)*Q(:,:,K3-k,n)
            enddo
          enddo
      endif

#if defined (PARALLEL_MPI)
    endif
#endif

      Nullify (P0, P1)

    END Subroutine Flux_BC_3D_z_Periodic
!
!======================================================================
!
    Subroutine Flux_BC_3D_z_Freestream

#if defined (PARALLEL_MPI)
      if (I_Am == First_Processor) then
#endif
        SELECT CASE (BC_IO_Type(Index,1))
          CASE (-1)     ! Outflow
            do k = K0,K1
              Q(:,:,k,:) = ZERO
            enddo
 
          CASE ( 0)     ! No BC

          CASE DEFAULT  ! Freestream
            do k = K0,K1
              Q(:,:,k,:) = ZERO
            enddo

        END SELECT

#if defined (PARALLEL_MPI)
      endif

      if (I_Am ==  Last_Processor) then
#endif
        if (Symmetry(Index) == 1) then

            N_Vector = ABS(K5-K4)+1 ; II = K3-N_Vector+1
            do n = 1,NV
              Q(:,:,K4:K5,n) = IEO(Index,n)*Q(:,:,K3:II:-1,n)
            enddo

        else

          SELECT CASE (BC_IO_Type(Index,2))
            CASE (-1)     ! Outflow
              do k = K4,K5
                Q(:,:,k,:) = ZERO
              enddo

            CASE ( 0)     ! No BC

            CASE DEFAULT  ! Freestream
              do k = K4,K5
                Q(:,:,k,:) = ZERO
              enddo

          END SELECT

        endif

#if defined (PARALLEL_MPI)
      endif
#endif

    END Subroutine Flux_BC_3D_z_Freestream

  END Subroutine Flux_BC_3D 

END MODULE Navier_Stokes_Flux_Boundary_Condition
