!
!=======================================================================
!
  Subroutine MPI_Data_ReAssemble_2D (f  , N0  , N5  , M0  , M5  ,  &
                                     f_g, N0_g, N5_g, M0_g, M5_g,  &
                                     n_0, m_0)

  integer  :: N0_g, N5_g, M0_g, M5_g
  integer  :: N0  , N5  , M0  , M5    
  integer  :: n_0, m_0

  DOUBLE PRECISION, dimension(N0  :N5  ,M0  :M5  )          :: f
  DOUBLE PRECISION, dimension(N0_g:N5_g,M0_g:M5_g), TARGET  :: f_g
  DOUBLE PRECISION, dimension(:,:)                , POINTER :: p

  DOUBLE PRECISION, dimension(:,:),                 ALLOCATABLE :: g

  integer , dimension(-2:2) :: A 
  integer                   :: K_Processor, Tag_1, Tag_2
  integer                   :: Status(MPI_Status_Size)

#if 1
  integer :: i,j, NN
#endif
#if 0
  if (N_Processor == 1) then ; f_g = f ; RETURN ; endif
#endif

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

                A(0) = 0
  A(-1) = n_0 ; A(1) = SIZE(f,DIM=1)-1
  A(-2) = m_0 ; A(2) = SIZE(f,DIM=2)-1

#if 1
write (6,*) N0_g, N5_g, M0_g, M5_g
write (6,*) N0  , N5  , M0  , M5
write (6,*) I_Am, A

  do j = M0,M5
    do i = N0,N5
      f(i,j) = (10*I_Am+1)*(i*j)
    enddo
  enddo
#endif

#if 0
  if (I_Am == First_Processor)                                           &
    f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2)) = f

  do K_Processor = First_Processor+1,Last_Processor
#else
  do K_Processor = First_Processor  ,Last_Processor
#endif

    Tag_1 = K_Processor+5555 ; Tag_2 = K_Processor+6666

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

    if (I_Am ==     K_Processor)                                         &
      call MPI_Send (A, SIZE(A), MPI_Integer , First_Processor, Tag_1,   &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor) then
      call MPI_Recv (A, SIZE(A), MPI_Integer ,     K_Processor, Tag_1,   &
                                 MPI_Comm_World, Status, MPI_Error_Status)

#if 0
      Nullify (p) ; p => f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2))
write (6,*) 'p :', SIZE(p), A
#else
!     Allocate (g(A(1)+1,A(2)+1))
!     Allocate (g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2)))
write (6,*) 'g :', SIZE(g), A
#endif
    endif
      Allocate (g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2)))

#if 1
    if (I_Am ==     K_Processor) then
      do j = M0,M5
        write (6,100) 'f',j, f(:,j)
      enddo
    endif
#endif

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

#if 1
      NN = SIZE(g)
      write (6,*) 'NN ', SIZE(g), SIZE(f)
#endif

    if (I_Am ==     K_Processor)                                                 &
      call MPI_Send (f, NN, MPI_DOUBLE_PRECISION, First_Processor, Tag_2,   &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor) then
#if 0
      call MPI_Recv (p, SIZE(p), MPI_DOUBLE_PRECISION,     K_Processor, Tag_2,   &
                                 MPI_Comm_World, Status, MPI_Error_Status)

      do j = LBOUND(p,DIM=2),UBOUND(p,DIM=2)
        write (6,100) 'p',j, p(:,j)
      enddo
100 format (a,1x,i5,16(1x,f5.0,:))
#else
      call MPI_Recv (g, NN, MPI_DOUBLE_PRECISION,     K_Processor, Tag_2,   &
                                 MPI_Comm_World, Status, MPI_Error_Status)

      do j = A(-2),A(-2)+A(2)
        write (6,100) 'g',j, g(:,j)
      enddo
100 format (a,1x,i5,16(1x,f5.0,:))

      f_g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2)) = g(A(-1):A(-1)+A(1),A(-2):A(-2)+A(2))
#endif
    endif

    call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
  enddo

  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)

#if 1
  STOP
#endif

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

  DOUBLE PRECISION, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )          :: f
  DOUBLE PRECISION, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g), TARGET  :: f_g
  DOUBLE PRECISION, dimension(:,:,:)                        , POINTER :: p

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
      call MPI_Send (f, SIZE(f), MPI_DOUBLE_PRECISION, First_Processor, Tag_2,   &
                                 MPI_Comm_World,         MPI_Error_Status)

    if (I_Am == First_Processor)                                         &
      call MPI_Recv (p, SIZE(p), MPI_DOUBLE_PRECISION,     K_Processor, Tag_2,   &
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
