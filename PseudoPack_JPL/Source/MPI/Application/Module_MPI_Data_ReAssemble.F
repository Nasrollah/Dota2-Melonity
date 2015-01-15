!------------------
#if 1
!------------------

#include "pseudopack.h"


MODULE MPI_Data_ReAssemble

#if defined (PARALLEL_MPI)
  USE Processor 
#endif

implicit NONE

INTERFACE PS_MPI_Data_ReAssemble
  MODULE PROCEDURE MPI_Data_ReAssemble_1D
  MODULE PROCEDURE MPI_Data_ReAssemble_2D
  MODULE PROCEDURE MPI_Data_ReAssemble_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_MPI_Data_ReAssemble

CONTAINS
#if defined (PARALLEL_MPI)
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_1D (f  , N0  , N5  ,   &
                                     f_g, N0_g, N5_g,   &
                                     n_0, MPI_Communicator)

  integer  :: N0_g, N5_g
  integer  :: N0  , N5      
  integer  :: n_0

  REALTYPE, dimension(N0  :N5  )          :: f
  REALTYPE, dimension(N0_g:N5_g), TARGET  :: f_g
  REALTYPE, dimension(:)        , POINTER :: p

  integer , dimension(-1:1) :: A 
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)

  integer            :: M_Processor

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type

                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)
 
  if (M_Processor == 1) then ; f_g = f ; RETURN ; endif

  call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

                A(0) = 0
  A(-1) = n_0 ; A(1) = SIZE(f,DIM=1)-1

  if (I_Am == First_Processor) &
    f_g(A(-1):A(-1)+A(1)) = f

  do K_Processor = First_Processor+1,M_Processor-1

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,   &
                                 MPI_Comm_Type ,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,   &
                                 MPI_Comm_Type , Status, MPI_Error_Status)

      Nullify (p) ; p => f_g(A(-1):A(-1)+A(1))
    endif

    call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (f, SIZE(f), MPI_REALTYPE, First_Processor, Tag_2,   &
                                 MPI_Comm_Type ,         MPI_Error_Status)

    if (I_Am == First_Processor)                                         &
      call MPI_Recv (p, SIZE(p), MPI_REALTYPE,     K_Processor, Tag_2,   &
                                 MPI_Comm_Type , Status, MPI_Error_Status)
  enddo

  call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

  END Subroutine MPI_Data_ReAssemble_1D
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_2D (f  , N0  , N5  , M0  , M5  ,  &
                                     f_g, N0_g, N5_g, M0_g, M5_g,  &
                                     n_0, m_0, MPI_Communicator)

  integer  :: N0_g, N5_g, M0_g, M5_g
  integer  :: N0  , N5  , M0  , M5    
  integer  :: n_0, m_0

  REALTYPE, dimension(N0  :N5  ,M0  :M5  )          :: f
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g), TARGET  :: f_g
  REALTYPE, dimension(:,:)                , POINTER :: p

  integer , dimension(-2:2) :: A 
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)

  integer            :: M_Processor

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type

                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)
 
  if (M_Processor == 1) then ; f_g = f ; RETURN ; endif

  call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

                A(0) = 0
  A(-1) = n_0 ; A(1) = SIZE(f,DIM=1)-1
  A(-2) = m_0 ; A(2) = SIZE(f,DIM=2)-1

  if (I_Am == First_Processor)                                           &
    f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2)) = f

  do K_Processor = First_Processor+1,M_Processor-1

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,   &
                                 MPI_Comm_Type ,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,   &
                                 MPI_Comm_Type , Status, MPI_Error_Status)

      Nullify (p) ; p => f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2))
    endif

    call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (f, SIZE(f), MPI_REALTYPE, First_Processor, Tag_2,   &
                                 MPI_Comm_Type ,         MPI_Error_Status)

    if (I_Am == First_Processor)                                         &
      call MPI_Recv (p, SIZE(p), MPI_REALTYPE,     K_Processor, Tag_2,   &
                                 MPI_Comm_Type , Status, MPI_Error_Status)

  enddo

  call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

  END Subroutine MPI_Data_ReAssemble_2D
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_3D (f  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                     f_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                     n_0, m_0, k_0, MPI_Communicator)

  integer  :: N0_g, N5_g, M0_g, M5_g, K0_g, K5_g
  integer  :: N0  , N5  , M0  , M5  , K0  , K5  
  integer  :: n_0, m_0, k_0

  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )          :: f
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g), TARGET  :: f_g
  REALTYPE, dimension(:,:,:)                        , POINTER :: p

  integer , dimension(-3:3) :: A 
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)

  integer            :: M_Processor

  integer , OPTIONAL :: MPI_Communicator
  integer            :: MPI_Comm_Type

                                 MPI_Comm_Type = MPI_Comm_World
  if (PRESENT(MPI_Communicator)) MPI_Comm_Type = MPI_Communicator

  call MPI_Comm_Size (MPI_Comm_Type, M_Processor, MPI_Error_Status)
 
  if (M_Processor == 1) then ; f_g = f ; RETURN ; endif

  call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

                A(0) = 0
  A(-1) = n_0 ; A(1) = SIZE(f,DIM=1)-1
  A(-2) = m_0 ; A(2) = SIZE(f,DIM=2)-1
  A(-3) = k_0 ; A(3) = SIZE(f,DIM=3)-1

  if (I_Am == First_Processor)                                           &
    f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2),A(-3):A(-3)+A(3)) = f

  do K_Processor = First_Processor+1,M_Processor-1

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,   &
                                 MPI_Comm_Type ,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,   &
                                 MPI_Comm_Type , Status, MPI_Error_Status)

      Nullify (p) ; p => f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2),A(-3):A(-3)+A(3))
    endif

    call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (f, SIZE(f), MPI_REALTYPE, First_Processor, Tag_2,   &
                                 MPI_Comm_Type ,         MPI_Error_Status)

    if (I_Am == First_Processor)                                         &
      call MPI_Recv (p, SIZE(p), MPI_REALTYPE,     K_Processor, Tag_2,   &
                                 MPI_Comm_Type , Status, MPI_Error_Status)
  enddo

  call MPI_Barrier (MPI_Comm_Type , MPI_Error_Status)

  END Subroutine MPI_Data_ReAssemble_3D
#else
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_1D (f  , N0  , N5  ,   &
                                     f_g, N0_g, N5_g,   &
                                     n_0)

  integer                   :: N0_g, N5_g
  integer                   :: N0  , N5      
  integer , OPTIONAL        :: n_0

  REALTYPE, dimension(N0  :N5  )          :: f
  REALTYPE, dimension(N0_g:N5_g)          :: f_g

  f_g = f

  END Subroutine MPI_Data_ReAssemble_1D 
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_2D (f  , N0  , N5  , M0  , M5  ,  &
                                     f_g, N0_g, N5_g, M0_g, M5_g,  &
                                     n_0, m_0)

  integer                   :: N0_g, N5_g, M0_g, M5_g
  integer                   :: N0  , N5  , M0  , M5    
  integer , OPTIONAL        :: n_0, m_0

  REALTYPE, dimension(N0  :N5  ,M0  :M5  )          :: f
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g)          :: f_g

  f_g = f

  END Subroutine MPI_Data_ReAssemble_2D 
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_3D (f  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                     f_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                     n_0, m_0, k_0)

  integer                   :: N0_g, N5_g, M0_g, M5_g, K0_g, K5_g
  integer                   :: N0  , N5  , M0  , M5  , K0  , K5  
  integer , OPTIONAL        :: n_0, m_0, k_0

  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )          :: f
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g)          :: f_g

  f_g = f

  END Subroutine MPI_Data_ReAssemble_3D 
#endif

END MODULE MPI_Data_ReAssemble

!------------------
#else
!------------------

#include "pseudopack.h"


MODULE MPI_Data_ReAssemble

#if defined (PARALLEL_MPI)
  USE Processor 
#endif

implicit NONE

INTERFACE PS_MPI_Data_ReAssemble
  MODULE PROCEDURE MPI_Data_ReAssemble_1D
  MODULE PROCEDURE MPI_Data_ReAssemble_2D
  MODULE PROCEDURE MPI_Data_ReAssemble_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_MPI_Data_ReAssemble

CONTAINS
#if defined (PARALLEL_MPI)
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_1D (f  , N0  , N5  ,   &
                                     f_g, N0_g, N5_g,   &
                                     n_0)

  integer  :: N0_g, N5_g
  integer  :: N0  , N5      
  integer  :: n_0

  REALTYPE, dimension(N0  :N5  )          :: f
  REALTYPE, dimension(N0_g:N5_g), TARGET  :: f_g
  REALTYPE, dimension(:)        , POINTER :: p

  integer , dimension(-1:1) :: A 
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)

  if (N_Processor == 1) then ; f_g = f ; RETURN ; endif

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

                A(0) = 0
  A(-1) = n_0 ; A(1) = SIZE(f,DIM=1)-1

  if (I_Am == First_Processor) &
    f_g(A(-1):A(-1)+A(1)) = f

  do K_Processor = First_Processor+1,Last_Processor

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,   &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,   &
                                 MPI_Comm_World, Status, MPI_Error_Status)

      Nullify (p) ; p => f_g(A(-1):A(-1)+A(1))
    endif

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (f, SIZE(f), MPI_REALTYPE, First_Processor, Tag_2,   &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor)                                         &
      call MPI_Recv (p, SIZE(p), MPI_REALTYPE,     K_Processor, Tag_2,   &
                                 MPI_Comm_World, Status, MPI_Error_Status)
  enddo

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

  END Subroutine MPI_Data_ReAssemble_1D
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_2D (f  , N0  , N5  , M0  , M5  ,  &
                                     f_g, N0_g, N5_g, M0_g, M5_g,  &
                                     n_0, m_0)

  integer  :: N0_g, N5_g, M0_g, M5_g
  integer  :: N0  , N5  , M0  , M5    
  integer  :: n_0, m_0

  REALTYPE, dimension(N0  :N5  ,M0  :M5  )          :: f
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g), TARGET  :: f_g
  REALTYPE, dimension(:,:)                , POINTER :: p

  integer , dimension(-2:2) :: A 
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)

  if (N_Processor == 1) then ; f_g = f ; RETURN ; endif

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

                A(0) = 0
  A(-1) = n_0 ; A(1) = SIZE(f,DIM=1)-1
  A(-2) = m_0 ; A(2) = SIZE(f,DIM=2)-1

  if (I_Am == First_Processor)                                           &
    f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2)) = f

  do K_Processor = First_Processor+1,Last_Processor

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,   &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,   &
                                 MPI_Comm_World, Status, MPI_Error_Status)

      Nullify (p) ; p => f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2))
    endif

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (f, SIZE(f), MPI_REALTYPE, First_Processor, Tag_2,   &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor)                                         &
      call MPI_Recv (p, SIZE(p), MPI_REALTYPE,     K_Processor, Tag_2,   &
                                 MPI_Comm_World, Status, MPI_Error_Status)

  enddo

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

  END Subroutine MPI_Data_ReAssemble_2D
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_3D (f  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                     f_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                     n_0, m_0, k_0)

  integer  :: N0_g, N5_g, M0_g, M5_g, K0_g, K5_g
  integer  :: N0  , N5  , M0  , M5  , K0  , K5  
  integer  :: n_0, m_0, k_0

  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )          :: f
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g), TARGET  :: f_g
  REALTYPE, dimension(:,:,:)                        , POINTER :: p

  integer , dimension(-3:3) :: A 
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)

  if (N_Processor == 1) then ; f_g = f ; RETURN ; endif

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

                A(0) = 0
  A(-1) = n_0 ; A(1) = SIZE(f,DIM=1)-1
  A(-2) = m_0 ; A(2) = SIZE(f,DIM=2)-1
  A(-3) = k_0 ; A(3) = SIZE(f,DIM=3)-1

  if (I_Am == First_Processor)                                           &
    f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2),A(-3):A(-3)+A(3)) = f

  do K_Processor = First_Processor+1,Last_Processor

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,   &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,   &
                                 MPI_Comm_World, Status, MPI_Error_Status)

      Nullify (p) ; p => f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2),A(-3):A(-3)+A(3))
    endif

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (f, SIZE(f), MPI_REALTYPE, First_Processor, Tag_2,   &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor)                                         &
      call MPI_Recv (p, SIZE(p), MPI_REALTYPE,     K_Processor, Tag_2,   &
                                 MPI_Comm_World, Status, MPI_Error_Status)
  enddo

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

  END Subroutine MPI_Data_ReAssemble_3D
#else
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_1D (f  , N0  , N5  ,   &
                                     f_g, N0_g, N5_g,   &
                                     n_0)

  integer                   :: N0_g, N5_g
  integer                   :: N0  , N5      
  integer , OPTIONAL        :: n_0

  REALTYPE, dimension(N0  :N5  )          :: f
  REALTYPE, dimension(N0_g:N5_g)          :: f_g

  f_g = f

  END Subroutine MPI_Data_ReAssemble_1D 
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_2D (f  , N0  , N5  , M0  , M5  ,  &
                                     f_g, N0_g, N5_g, M0_g, M5_g,  &
                                     n_0, m_0)

  integer                   :: N0_g, N5_g, M0_g, M5_g
  integer                   :: N0  , N5  , M0  , M5    
  integer , OPTIONAL        :: n_0, m_0

  REALTYPE, dimension(N0  :N5  ,M0  :M5  )          :: f
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g)          :: f_g

  f_g = f

  END Subroutine MPI_Data_ReAssemble_2D 
!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_3D (f  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                     f_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                     n_0, m_0, k_0)

  integer                   :: N0_g, N5_g, M0_g, M5_g, K0_g, K5_g
  integer                   :: N0  , N5  , M0  , M5  , K0  , K5  
  integer , OPTIONAL        :: n_0, m_0, k_0

  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )          :: f
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g)          :: f_g

  f_g = f

  END Subroutine MPI_Data_ReAssemble_3D 
#endif

END MODULE MPI_Data_ReAssemble
!------------------
#endif
!------------------

