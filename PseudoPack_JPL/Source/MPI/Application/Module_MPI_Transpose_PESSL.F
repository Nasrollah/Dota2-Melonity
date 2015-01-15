#include "pseudopack.h"



MODULE MPI_Transpose_PESSL
  USE Processor, ONLY: N_Processor
  USE MPI_Data_Distribution
  USE Memory_Allocation_Status

implicit NONE

PRIVATE
PUBLIC  :: P_Transpose

CONTAINS

#if defined (PARALLEL_MPI)
# if defined (PESSL) || defined (SCALAPACK)
!
!=======================================================================
!
  Subroutine P_Transpose   (N, M, f, g)

  integer                :: N, M
  REALTYPE, dimension(*) :: f, g

  if (N_Processor == 1) then
    call Transpose_S  (N, M, f, g)
  else  
    if (MOD(M,N_Processor)+MOD(N,N_Processor) == 0) then
      call Transpose_NM (N, M, f, g)
    else
      call Transpose_G  (N, M, f, g)
    endif
  endif

  END Subroutine P_Transpose

!
!=======================================================================
!
  Subroutine Transpose_S  (N, M, f, g)

  integer                :: i, j
  integer                :: N, M
  REALTYPE               :: f(N,M), u(M,N), g(*)

  do i = 1,N
    do j = 1,M
      u(j,i) = f(i,j)
    enddo
  enddo

  call Copy_u_to_f (M, N, u, g)

  END Subroutine Transpose_S
!
!=======================================================================
!
  Subroutine Transpose_NM (N, M, f, g)

  integer                :: N, M

  integer  :: n_Block, m_Block
  integer  :: NP_Row, NP_Col, My_Row, My_Col, IContxt
  integer  :: NumRoC

  integer , dimension(9)                :: Desc_f, Desc_u
  REALTYPE, dimension(*)                :: f, g
  REALTYPE, dimension(:,:), ALLOCATABLE :: u

  n_Block = N/N_Processor; m_Block = M/N_Processor

  if (.NOT. ALLOCATED(u)) then
    ALLOCATE (u(M,n_Block), STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = M*n_Block
    Subroutine_Name  = TRIM(Subroutine_Call)//'Transpose_NM : u'

    call Memory_Error_Check (0)
#endif
  endif

  NP_Row = 1 ; NP_Col = N_Processor

  call BLACS_Get      (0, 0, IContxt)
  call BLACS_GridInit (IContxt, 'R', NP_Row,NP_Col)
  call BLACS_GridInfo (IContxt,      NP_Row,NP_Col, My_Row,My_Col)

  Desc_f(1) = 1
  Desc_f(2) = IContxt
  Desc_f(3) = N
  Desc_f(4) = M
  Desc_f(5) = n_Block
  Desc_f(6) = m_Block
  Desc_f(7) = 0
  Desc_f(8) = 0
  Desc_f(9) = MAX(1, NumRoC(Desc_f(3),Desc_f(5),My_Row,Desc_f(7),NP_Row))

  Desc_u(1) = 1
  Desc_u(2) = IContxt
  Desc_u(3) = M
  Desc_u(4) = N
  Desc_u(5) = m_Block
  Desc_u(6) = n_Block
  Desc_u(7) = 0
  Desc_u(8) = 0
  Desc_u(9) = MAX(1, NumRoC(Desc_u(3),Desc_u(5),My_Row,Desc_u(7),NP_Row))

  call PDTRAN (M, N, ONE, f, 1,1, Desc_f, ZERO, u, 1,1, Desc_u)

  call Copy_u_to_f (M, n_Block, u, g)

  if (ALLOCATED(u)) then
    DEALLOCATE (u, STAT=M_Error) 

#if defined (DEBUG)
    Memory_Requested = M*n_Block
    Subroutine_Name  = TRIM(Subroutine_Call)//'Transpose_NM : u'

    call Memory_Error_Check (1)
#endif 
  endif

  call BLACS_GridEXIT (IContxt)

  END Subroutine Transpose_NM 
!
!=======================================================================
!
  Subroutine Transpose_G  (N, M, f, g)

  integer                :: N, M

  integer  :: n_Remain, n_0, n_Star, n_Block, NN
  integer  :: m_Remain, m_0, m_Star, m_Block, MM
  integer  :: NP_Row, NP_Col, My_Row, My_Col, IContxt
  integer  :: NumRoC, NM

  integer , dimension(9)                :: Desc_f, Desc_u
  REALTYPE, dimension(*)                :: f, g
  REALTYPE, dimension(:,:), ALLOCATABLE :: uu, ff

  call PS_MPI_Block_Distribution (N, n_Remain, n_0, n_Star, n_Block, NN)
  call PS_MPI_Block_Distribution (M, m_Remain, m_0, m_Star, m_Block, MM)

  if (MM >= N) then; NM = MM; else; NM = NN; endif

  if (.NOT. ALLOCATED(ff)) then
    ALLOCATE (ff(NM,m_Block), STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = NM*n_Block
    Subroutine_Name  = TRIM(Subroutine_Call)//'Transpose_G : ff'

    call Memory_Error_Check (0)
#endif
  endif

  if (.NOT. ALLOCATED(uu)) then
    ALLOCATE (uu(MM,n_Block), STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = MM*n_Block
    Subroutine_Name  = TRIM(Subroutine_Call)//'Transpose_G : uu'

    call Memory_Error_Check (0)
#endif
  endif

  NP_Row = 1 ; NP_Col = N_Processor

  call BLACS_Get      (0, 0, IContxt)
  call BLACS_GridInit (IContxt, 'R', NP_Row,NP_Col)
  call BLACS_GridInfo (IContxt,      NP_Row,NP_Col, My_Row,My_Col)

  Desc_f(1) = 1
  Desc_f(2) = IContxt
  Desc_f(3) = NM
  Desc_f(4) = MM
  Desc_f(5) = n_Block
  Desc_f(6) = m_Block
  Desc_f(7) = 0
  Desc_f(8) = 0
  Desc_f(9) = MAX(1, NumRoC(Desc_f(3),Desc_f(5),My_Row,Desc_f(7),NP_Row))

  Desc_u(1) = 1
  Desc_u(2) = IContxt
  Desc_u(3) = MM
  Desc_u(4) = NM
  Desc_u(5) = m_Block
  Desc_u(6) = n_Block
  Desc_u(7) = 0
  Desc_u(8) = 0
  Desc_u(9) = MAX(1, NumRoC(Desc_u(3),Desc_u(5),My_Row,Desc_u(7),NP_Row))

  call Expand  (N ,NM, m_Star, n_Block, n_Remain, f,ff)

  call PDTRAN  (MM,NN, ONE, ff, 1,1, Desc_f, ZERO, uu, 1,1, Desc_u)

  call Restore (MM,M , n_Star, m_Block, m_Remain, uu,g)

  if (ALLOCATED(uu)) then
    DEALLOCATE (uu, STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = MM*n_Block
    Subroutine_Name  = TRIM(Subroutine_Call)//'Transpose_G : uu'

    call Memory_Error_Check (1)
#endif 
  endif

  if (ALLOCATED(ff)) then
    DEALLOCATE (ff, STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = NM*n_Block
    Subroutine_Name  = TRIM(Subroutine_Call)//'Transpose_G : ff'

    call Memory_Error_Check (1)
#endif 
  endif

  call BLACS_GridEXIT (IContxt)

  CONTAINS
!
!=======================================================================
!
    Subroutine Expand  (N ,NN, m_Star, n_Block, n_Remain, f,ff)

    integer  :: N, NN, m_Star, n_Block, n_Remain, NP
    integer  :: i, k, i_Block, ii

    REALTYPE :: f(N,m_Star), ff(NN,m_Star)

    if (n_Remain == 0) then
      NP = N_Processor ; else ; NP = n_Remain
    endif

    i = NP*n_Block

    ff(1:i,:) = f(1:i,:)

    if (i >= N) RETURN

    do k = NP, N_Processor-1

      do i_Block = 1,n_Block-1
          ii = i_Block+k*n_Block ; i = i+1 
       
        ff(ii,:) = f(i,:)
      enddo

      if (i >= N) EXIT
    enddo

    END Subroutine Expand 
!
!=======================================================================
!
    Subroutine Restore (MM,M , n_Star, m_Block, m_Remain, uu,u)

    integer  :: M, MM, n_Star, m_Block, m_Remain, NP
    integer  :: j, k, j_Block, jj

    REALTYPE :: uu(MM,n_Star), u(M,n_Star)

    if (m_Remain == 0) then
      NP = N_Processor ; else ; NP = m_Remain
    endif

    j = NP*m_Block

    u(1:j,:) = uu(1:j,:)

    if (j >= M) RETURN

    do k = NP, N_Processor-1

      do j_Block = 1,m_Block-1
        j = j+1 ; jj = j_Block+k*m_Block

        u(j,:) = uu(jj,:)

        if (j >= M) RETURN
      enddo

    enddo

    END Subroutine Restore 

  END Subroutine Transpose_G 
# else
!
!=======================================================================
!
  Subroutine P_Transpose   (N, M, f, g)

  integer                :: i, j
  integer                :: N, M
  REALTYPE               :: f(N,M), u(M,N), g(M,N)

  do i = 1,N
    do j = 1,M
      u(j,i) = f(i,j)
    enddo
  enddo

  call Copy_u_to_f (M, N, u, g)
  
  END Subroutine P_Transpose
# endif
#else
!
!=======================================================================
!
  Subroutine P_Transpose   (N, M, f, g)

  integer                :: i, j
  integer                :: N, M
  REALTYPE               :: f(N,M), u(M,N), g(M,N)

  do i = 1,N
    do j = 1,M
      u(j,i) = f(i,j)
    enddo
  enddo

  call Copy_u_to_f (M, N, u, g)
  
  END Subroutine P_Transpose
#endif
!
!=======================================================================
!
  Subroutine Copy_u_to_f (N, M, u, f)

  integer  :: N, M

  REALTYPE, dimension(N*M) :: u, f

  f = u

  END Subroutine Copy_u_to_f 

END MODULE MPI_Transpose_PESSL
