#include "pseudopack.h"


MODULE Matrix_Multiply_MPI
  USE Processor
  USE MPI_Data_Distribution
  USE Memory_Allocation_Status

#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE
 
PRIVATE
PUBLIC  :: MPI_GEMM

CONTAINS

#if defined (PARALLEL_MPI) 
! ----------------------------------------------------------------------
! FUNCTION NAME: MPI_GEMM
! Author       : Wai Sun Don
! Description  : Perform the matrix multiply of F with D :
!
!                if (Index = 1) DF = alpha*D*F+beta*DF
!                  D(N,N)*F(N,M) -> DF(N,M)
!
!                if (Index = 2) DF = alpha*TRANSPOSE(D*TRANSPOSE(F))+beta*DF
!                                  = alpha*F*TRANSPOSE(D)+beta*DF
!                  D(M,M)*F(N,M) -> DF(N,M)
! ----------------------
!
! if (Index = 1) then
!
!   Data D         is distributed along the Second index (IMPORTANT)
!   Data F and DF are distributed along the FIRST  index
!
!   Global : D(K_D,N     ), F(K_F,M), DF(K_DF,M)
!   Local  : D(L_D,n_Star), F(L_F,M), DF(L_DF,M)
!
!   INIT is not used in current implementation
!
!   N      <= K_D 
!   N      <= K_F
!   N      <= K_DF
!   N      <= L_D 
!   n_Star <= L_F
!   n_Star <= L_DF
!
! ----------------------
!
! if (Index = 2) then
!
!   Data D, F and DF are distributed along the SECOND index
!
!   Global : D(K_D,M)     , F(K_F,M)     , DF(K_DF,M)
!   Local  : D(L_D,m_Star), F(L_F,m_Star), DF(L_DF,m_Star)
!
!   M      <= K_D
!   N      <= K_F
!   N      <= K_DF
!   M      <= L_D
!   N      <= L_F
!   N      <= L_DF
!
! ----------------------
!
!   K_D  is the Global leading dimension of Square array D
!   K_F  is the Global leading dimension of the data array F 
!   K_DF is the Global leading dimension of the data array DF 
!   L_D  is the Local  leading dimension of Square array D
!   L_F  is the Local  leading dimension of the data array F
!   L_DF is the Local  leading dimension of the data array DF
!
!   N      is the number of Rows for arrays D, F and DF (Global)
!   M      is the number of columns(vectors) (Global)
!
!   n_Star is the number of rows distributed across N_Processors
!             used for calculation if (Index = 1) (Local)
!   m_Star is the number of Cols distributed across N_Processors
!             used for calculation if (Index = 2) (Local)
!
! ----------------------
!
!  Note :  Should be run under SP2 User Space Communication Subsystem library
!          with environment setting
!               setenv MP_EUILIB       us
!          for timing purpose.
! ----------------------------------------------------------------------
!
! =====================================================================
 
  Subroutine MPI_GEMM (Index,   &
                       INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF) 
 
  integer            :: Index, INIT, N, M, L_D, L_F, L_DF
  REALTYPE           :: alpha, beta

  REALTYPE, dimension(*) :: D, F, DF

  if ((alpha == ZERO) .AND. (beta == ZERO)) RETURN

  SELECT CASE (Index)
    CASE (1)
      call MPI_GEMM_1 (INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF) 

    CASE (2)
      call MPI_GEMM_2 (INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF) 

  END SELECT

  END Subroutine MPI_GEMM 

!
!=======================================================================
!
  Subroutine MPI_GEMM_1 (INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF) 
 
  integer            :: INIT, N, M, L_D, L_F, L_DF
  REALTYPE           :: alpha, beta

  integer  :: n_Remain, n_0, n_Star
  integer  :: i, j, k
  REALTYPE :: S

  Integer , dimension(N_Processor)  :: RecvCounts
  REALTYPE, dimension(L_D ,*)       :: D
  REALTYPE, dimension(L_F ,*)       :: F
  REALTYPE, dimension(L_DF,*)       :: DF
  REALTYPE, dimension(N,M)          :: Local_Sum
#if 0
  REALTYPE, dimension(:), ALLOCATABLE :: DG
#else
  REALTYPE, dimension(INT(CEILING(ONE*N/N_Processor))) :: DG
#endif

  call PS_MPI_Data_Distribution (N, n_Remain, n_0, n_Star)

  if (alpha == ZERO) then
    if ((beta /= ZERO) .OR. (beta /= ONE)) then
      DF(1:n_Star,1:M) = beta*DF(1:n_Star,1:M)
    endif

    RETURN
  endif

  call MPI_ALLGATHER (n_Star, 1, MPI_Integer, RecvCounts, 1, MPI_Integer, &
                                 MPI_Comm_World, MPI_Error_Status)

#if 0
  Local_Sum = MATMUL(D(1:N,1:n_Star), F(1:n_Star,1:M))
#endif
#if 0
    do i = 1,N
      do j = 1,M

        S = ZERO
        do k = 1,n_Star
          S = S + D(i,k)*F(k,j)
        enddo
        Local_Sum(i,j) = S 

      enddo
    enddo
  endif
#endif
      
  call GEMM ('N','N', N, M, n_Star, ONE, D, L_D, F, L_F, ZERO, Local_Sum, N)

  if (beta == ZERO) then
    do j = 1,M
      call MPI_Reduce_Scatter (Local_Sum(1,j), DF(1,j), RecvCounts,    &
                                               MPI_REALTYPE, MPI_SUM,  &
                                               MPI_Comm_World, MPI_Error_Status)
    enddo

    if (alpha /= ONE) then
      DF(1:n_Star,1:M) = alpha*DF(1:n_Star,1:M)
    endif
  else
#if 0
    if (.NOT. ALLOCATED(DG)) then
      ALLOCATE (DG(n_Star), STAT=M_Error)

#if defined (DEBUG)
      Memory_Requested = n_Star
      Subroutine_Name  = TRIM(Subroutine_Call)//'MPI_GEMM_1 : DG'

      call Memory_Error_Check (0)
#endif
    endif
#endif

    do j = 1,M
      call MPI_Reduce_Scatter (Local_Sum(1,j), DG, RecvCounts,         &
                                               MPI_REALTYPE, MPI_SUM,  &
                                               MPI_Comm_World, MPI_Error_Status)

      if (beta == ONE) then
        if (alpha /= ONE) then
          DF(1:n_Star,j) =      DF(1:n_Star,j)+alpha*DG(1:n_Star)
        else
          DF(1:n_Star,j) =      DF(1:n_Star,j)+      DG(1:n_Star)
        endif
      else
        if (alpha /= ONE) then
          DF(1:n_Star,j) = beta*DF(1:n_Star,j)+alpha*DG(1:n_Star)
        else
          DF(1:n_Star,j) = beta*DF(1:n_Star,j)+      DG(1:n_Star)
        endif
      endif
    enddo
#if 0
    if (ALLOCATED(DG)) then
      DEALLOCATE (DG, STAT=M_Error)

#if defined (DEBUG)
      Memory_Requested = n_Star
      Subroutine_Name  = TRIM(Subroutine_Call)//'MPI_GEMM_1 : DG'

      call Memory_Error_Check (1)
#endif
    endif
#endif
  endif

  END Subroutine MPI_GEMM_1
#if 1
!
!=======================================================================
!
  Subroutine MPI_GEMM_2 (INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF) 
 
  integer            :: INIT, N, M, L_D, L_F, L_DF
  REALTYPE           :: alpha, beta

  integer  :: m_Remain, m_0, m_Star
  integer  :: i, j, k
  REALTYPE :: S

  Integer , dimension(N_Processor)  :: RecvCounts
  REALTYPE, dimension(L_D ,*)       :: D
  REALTYPE, dimension(L_F ,*)       :: F
  REALTYPE, dimension(L_DF,*)       :: DF
  REALTYPE, dimension(M,N)          :: Local_Sum
#if 0
  REALTYPE, dimension(:), ALLOCATABLE :: DG
#else
  REALTYPE, dimension(INT(CEILING(ONE*M/N_Processor))) :: DG
#endif
 
  call PS_MPI_Data_Distribution (M, m_Remain, m_0, m_Star)

  if (alpha == ZERO) then
    if ((beta /= ZERO) .OR. (beta /= ONE)) then
      DF(1:N,1:m_Star) = beta*DF(1:N,1:m_Star)
    endif

    RETURN
  endif

  call MPI_ALLGATHER (m_Star, 1, MPI_Integer, RecvCounts, 1, MPI_Integer, &
                                 MPI_Comm_World, MPI_Error_Status)

#if 0
  Local_Sum = TRANSPOSE(MATMUL(F(1:N,1:m_Star),D(1:m_Star,1:M)))
#endif
#if 0
  do i = 1,N
    do j = 1,M

      S = ZERO
      do k = 1,m_Star
        S = S + F(i,k)*D(k,j)
      enddo
      Local_Sum(j,i) = S 

    enddo
  enddo
#endif
      
  call GEMM ('T','T', M, N, m_Star, ONE, D, L_D, F, L_F, ZERO, Local_Sum, M)

#if 0
  if (.NOT. ALLOCATED(DG)) then
    ALLOCATE (DG(m_Star), STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = m_Star
    Subroutine_Name  = TRIM(Subroutine_Call)//'MPI_GEMM_2 : DG'

    call Memory_Error_Check (0)
#endif
  endif
#endif

  do i = 1,N
    call MPI_Reduce_Scatter (Local_Sum(1,i), DG, RecvCounts,         &
                                             MPI_REALTYPE, MPI_SUM,  &
                                             MPI_Comm_World, MPI_Error_Status)

    if (beta == ZERO) then
      if (alpha /= ONE) then
        DF(i,1:m_Star) =                       alpha*DG(1:m_Star)
      else
        DF(i,1:m_Star) =                             DG(1:m_Star)
      endif
    else
      if (beta == ONE) then
        if (alpha /= ONE) then
          DF(i,1:m_Star) =      DF(i,1:m_Star)+alpha*DG(1:m_Star)
        else
          DF(i,1:m_Star) =      DF(i,1:m_Star)+      DG(1:m_Star)
        endif
      else
        if (alpha /= ONE) then
          DF(i,1:m_Star) = beta*DF(i,1:m_Star)+alpha*DG(1:m_Star)
        else
          DF(i,1:m_Star) = beta*DF(i,1:m_Star)+      DG(1:m_Star)
        endif
      endif
    endif

  enddo

#if 0
  if (ALLOCATED(DG)) then
    DEALLOCATE (DG, STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = m_Star
    Subroutine_Name  = TRIM(Subroutine_Call)//'MPI_GEMM_2 : DG'

    call Memory_Error_Check (1)
#endif
  endif
#endif

  END Subroutine MPI_GEMM_2 
#endif
#if 0
!
!=======================================================================
!
  Subroutine MPI_GEMM_2 (INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF) 
 
  integer            :: INIT, N, M, L_D, L_F, L_DF
  REALTYPE           :: alpha, beta

  integer  :: m_Remain, m_0, m_Star
  integer  :: i, k
  REALTYPE :: S

  Integer , dimension(N_Processor)  :: RecvCounts
  REALTYPE, dimension(L_D ,*)       :: D
  REALTYPE, dimension(L_F ,*)       :: F
  REALTYPE, dimension(L_DF,*)       :: DF
  REALTYPE, dimension(M)            :: Local_Sum
#if 0
  REALTYPE, dimension(:), ALLOCATABLE :: DG
#else
  REALTYPE, dimension(INT(CEILING(ONE*M/N_Processor))) :: DG
#endif

  call PS_MPI_Data_Distribution (M, m_Remain, m_0, m_Star)

  if (alpha == ZERO) then
    if ((beta /= ZERO) .OR. (beta /= ONE)) then
      DF(1:N,1:m_Star) = beta*DF(1:N,1:m_Star)
    endif

    RETURN
  endif

  call MPI_ALLGATHER (m_Star, 1, MPI_Integer, RecvCounts, 1, MPI_Integer, &
                                 MPI_Comm_World, MPI_Error_Status)

#if 0
  if (.NOT. ALLOCATED(DG)) then
    ALLOCATE (DG(m_Star), STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = m_Star
    Subroutine_Name  = TRIM(Subroutine_Call)//'MPI_GEMM_2 : DG'

    call Memory_Error_Check (0)
#endif
  endif
#endif

  do i = 1,N

#if 0
    Local_Sum = MATMUL(F(i,1:m_Star),D(1:m_Star,1:M))
#endif
#if 0
    do j = 1,M
      S = ZERO
      do k = 1,m_Star
        S = S + F(i,k)*D(k,j)
      enddo
      Local_Sum(j) = S 
    enddo
#endif
      
    call GEMV ('T', m_Star, M, ONE, D, L_D, F(i,1), L_F, ZERO, Local_Sum, 1)

    call MPI_Reduce_Scatter (Local_Sum, DG, RecvCounts,         &
                                        MPI_REALTYPE, MPI_SUM,  &
                                        MPI_Comm_World, MPI_Error_Status)

    if (beta == ZERO) then
      if (alpha /= ONE) then
        DF(i,1:m_Star) =                       alpha*DG(1:m_Star)
      else
        DF(i,1:m_Star) =                             DG(1:m_Star)
      endif
    else
      if (beta == ONE) then
        if (alpha /= ONE) then
          DF(i,1:m_Star) =      DF(i,1:m_Star)+alpha*DG(1:m_Star)
        else
          DF(i,1:m_Star) =      DF(i,1:m_Star)+      DG(1:m_Star)
        endif
      else
        if (alpha /= ONE) then
          DF(i,1:m_Star) = beta*DF(i,1:m_Star)+alpha*DG(1:m_Star)
        else
          DF(i,1:m_Star) = beta*DF(i,1:m_Star)+      DG(1:m_Star)
        endif
      endif
    endif

  enddo

#if 0
  if (ALLOCATED(DG)) then
    DEALLOCATE (DG, STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = m_Star
    Subroutine_Name  = TRIM(Subroutine_Call)//'MPI_GEMM_2 : DG'

    call Memory_Error_Check (1)
#endif
  endif

  END Subroutine MPI_GEMM_2 
#endif

# endif
#else
!
!=======================================================================
!
  Subroutine MPI_GEMM (Index,   &
                       INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF)

  integer            :: Index, INIT, N, M, L_D, L_F, L_DF
  REALTYPE           :: alpha, beta

  REALTYPE, dimension(*) :: D, F, DF

  if ((alpha == ZERO) .AND. (beta == ZERO)) RETURN

  SELECT CASE (Index)
    CASE (1)
      call MPI_GEMM_1 (N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF)

    CASE (2)
      call MPI_GEMM_2 (N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF)

  END SELECT

  END Subroutine MPI_GEMM 
! ----------------------------------------------------------------------
! FUNCTION NAME: MPI_GEMM_1
! ----------------------------------------------------------------------
  Subroutine MPI_GEMM_1 (N, M, alpha, A, L_A, B, L_B, beta, C, L_C)

  integer                    :: N, M
  REALTYPE                   :: alpha, beta

  integer                    :: Row_A, Col_A, Col_C
  integer                    ::   L_A,   L_B,   L_C
  REALTYPE, dimension(*)     ::     A,     B,     C

  Row_A = N ; Col_A = N ; Col_C = M

  call   GEMM ('N', 'N', Row_A, Col_C, Col_A, alpha,      &
                             A, L_A, B, L_B, beta, C, L_C)

  END Subroutine MPI_GEMM_1
! ----------------------------------------------------------------------
! FUNCTION NAME: MPI_GEMM_2
! ----------------------------------------------------------------------
  Subroutine MPI_GEMM_2 (N, M, alpha, B, L_B, A, L_A, beta, C, L_C)

  integer                    :: N, M
  REALTYPE                   :: alpha, beta

  integer                    :: Row_A, Col_A, Col_C
  integer                    ::   L_A,   L_B,   L_C
  REALTYPE, dimension(*)     ::     A,     B,     C

  Row_A = N ; Col_A = M ; Col_C = M

  call   GEMM ('N', 'T', Row_A, Col_C, Col_A, alpha,      &
                             A, L_A, B, L_B, beta, C, L_C)

  END Subroutine MPI_GEMM_2
#endif

END MODULE Matrix_Multiply_MPI
