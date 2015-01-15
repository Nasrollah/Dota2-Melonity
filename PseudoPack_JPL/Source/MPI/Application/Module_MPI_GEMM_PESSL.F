#include "pseudopack.h"


#if defined (PESSL) || defined (SCALAPACK)

MODULE Matrix_Multiply_PESSL
  USE Processor, ONLY: N_Processor
  USE MPI_Data_Distribution
  USE Memory_Allocation_Status

implicit NONE

PRIVATE 
PUBLIC  :: P_GEMM_INIT
PUBLIC  :: P_GEMM

CONTAINS

! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_INIT
! Description  : Copy and reorganize the square matrix A into B
!
!  IF (INIT = 1) then
!
!      B =   A 
!    L_B = L_A  Upon RETURN to the calling routine.
!
!  ELSE
!
!  Call this routine to avoid call for INIT=1 inside the main routine P_GEMM
!
!  This is useful because the differentiation matrix usually
!  remain unchanged through out the calculation.
!
! if (Index = 1) then
!
! Data A and B are distributed along the FIRST  index
!
!   A of size (L_A,N)
!   B of size (n_Block,NN)
!   where n_Block=CEILING(N/N_Processor)
!         n_Star =Local number of rows in A < L_A
!         NN     =n_Block*N_Processor
!
!   L_B=n_Block is RETURN as the leading of dimenision of array B (Output)
!   L_A <= L_B
!
!   ONLY needed to be called once if N, N_Processor are not
!     changed during computation.
!
! ----------------
! if (Index = 2) then
!
! Data A and B are distributed along the SECOND index
!
!  A of size (L_A,m_Star)
!  B of size (L_B,m_Block)
!
!  where 
!
!     m_Block = CEILIMG(M/N_Processor)
!     m_Star  = Local number of rows in A
!     MM      =  m_Block*N_Processor
!
!  M  <= L_A
!
!  L_B=MM is RETURN as the leading of dimenision of array B (Output)
!
!  ONLY needed to be called once if M, N_Processor are not
!    changed during calculation.
!
!  Under the following conditions that INIT will be set and return as 1 
!    regardless its initial setting...
!
!    1.)  N_Processor = 1
!    2.)  Index = 1, MOD(N, N_Processor) = 0
!    3.)  Index = 2, MOD(M, N_Processor) = 0
!
!  ERROR = 0    No Error
!        = 1    L_B < L_A
!        = 2    L_B*N_B < L_A*N_A
!        = 3    L_B < L_A  and  L_B*N_B < L_A*N_A
!        = 4    L_A < M
!        
! ----------------------------------------------------------------------
  Subroutine P_GEMM_INIT (Index, INIT, N, M, A, L_A, N_A, B, L_B, N_B,   &
                                                                ERROR)

  integer                 :: Index, INIT, N, M, L_A, N_A, L_B, N_B 
  integer ,      OPTIONAL :: ERROR 
  REALTYPE, dimension(*)  :: A, B

#if defined (DEBUG)
  if (PRESENT(ERROR)) then
    ERROR = 0

    if (INIT == 1) then
      if (L_B     <  L_A    ) ERROR = 1
      if (L_B*N_B <  L_A*N_A) ERROR = ERROR+2
    endif

    if ((Index == 2) .AND. (L_A <  M)) ERROR = 4
    
    if (ERROR /= 0) RETURN
  endif
#endif

  if (N_Processor == 1)                             INIT = 1
  if ((Index == 1) .AND. (MOD(N,N_Processor) == 0)) INIT = 1
  if ((Index == 2) .AND. (MOD(M,N_Processor) == 0)) INIT = 1

  SELECT CASE (INIT)
    CASE (0)
      SELECT CASE (Index)
        CASE (1)
          call P_GEMM_INIT_1 (N, A, L_A, B, L_B)

        CASE (2)
          call P_GEMM_INIT_2 (M, A, L_A, B, L_B)

      END SELECT

    CASE (1)
      B(1:L_A*N_A) = A(1:L_A*N_A) ; L_B = L_A 

  END SELECT

  END Subroutine P_GEMM_INIT 

! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_INIT_1
! Description  : Copy and reorganize the square matrix A into AA
!
!  A  of size (L_A,N) 
!  AA of size (n_Block,NN)
!
!  where 
!    n_Block = CEILING(N/N_Processor)
!    n_Star  = Local number of rows in A < L_A
!    NN      = n_Block*N_Processor
!
!  L_AA=n_Block is the leading of dimenision of array AA
!  
!  ONLY needed to be called once if N, N_Processor are not
!    changed during computation.
!         
! ----------------------------------------------------------------------

  Subroutine P_GEMM_INIT_1 (N, A, L_A, AA, L_AA)

  integer  :: N, L_A, L_AA
  integer  :: n_Remain, n_0, n_Star, n_Block, NN

  REALTYPE, dimension(*) :: A, AA

  call PS_MPI_Block_Distribution (N, n_Remain, n_0, n_Star, n_Block, NN)

  if ((N_Processor == 1) .OR. MOD(N, N_Processor) == 0) then
    call Copy_A_2_AA  (L_A, n_Star, N, n_Block, NN, A, AA)
  else
    call P_Expand_A_1 (L_A, n_Star, N, n_Block, NN, n_Remain, A, AA)
  endif

  L_AA = n_Block

  END Subroutine P_GEMM_INIT_1 

! ----------------------------------------------------------------------
! FUNCTIOM MAME: P_GEMM_INIT_2
! Description  : Copy and reorganize the square matrix B into BB
!
!  B  of size (L_B ,m_Star)
!  BB of size (L_BB,m_Block)
!
!  where 
!
!    m_Block = CEILIMG(M/N_Processor)
!    m_Star  = Local number of rows in B 
!    MM      = m_Block*N_Processor
!
!  M  < L_B
!  L_BB=MM is the leading of dimenision of array BB 
!  
!  ONLY needed to be called once if M, N_Processor are not
!    changed during calculation.
!         
! ----------------------------------------------------------------------

  Subroutine P_GEMM_INIT_2 (M, B, L_B, BB, L_BB)

  integer  :: M, L_B, L_BB
  integer  :: m_Remain, m_0, m_Star, m_Block, MM

  REALTYPE, dimension(*) :: B, BB

  call PS_MPI_Block_Distribution (M, m_Remain, m_0, m_Star, m_Block, MM)

  if ((N_Processor == 1) .OR. MOD(M, N_Processor) == 0) then
    call Copy_B_2_BB  (L_B, M, m_Star, MM, m_Block, B, BB)
  else
    call P_Expand_B_2 (L_B, M, m_Star, MM, m_Block, m_Remain, B, BB)
  endif

  L_BB = MM

  END Subroutine P_GEMM_INIT_2 



! **********************************************************************




! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM
! Author       : Wai Sun Don
! Description  : Perform the matrix multiply of F with D :
!
!                if (Index = 1) DF = alpha*D*F+beta*DF
!                  D(N,N)*F(N,M) -> DF(N,M)
!
!                if (Index = 2) DF = alpha*Transpose(D*Transpose(F))+beta*DF
!                                  = alpha*F*Transpose(D)+beta*DF
!                  D(M,M)*F(N,M) -> DF(N,M)
! ----------------------
!
! if (Index = 1) then
!
!   Data D, F and DF are distributed along the FIRST index
!
!   Global : D(K_D,N), F(K_F,M), DF(K_DF,M)
!   Local  : D(L_D,N), F(L_F,M), DF(L_DF,M)
!
!   INIT : 1  Reorganize array D into DD with proper size
!               if N_Processor>1 and mod(N,N_Processor)>0
!          0  if P_GEMM_INIT_1 was called, or
!               if N_Processor=1 or mod(N,N_Processor)=0
!
!   N      <= K_D 
!   N      <= K_F
!   N      <= K_DF
!   n_Star <= L_D 
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
!   INIT : 1 Yes  Reorganize array D into DD with proper size
!                 if N_Processor>1 and mod(M,N_Processor)>0
!          0 No   if P_GEMM_INIT_2 was called, or
!                 if N_Processor=1 or mod(M,N_Processor)=0
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
! ----------------------------------------------------------------------

  Subroutine P_GEMM (Index,   &
                     INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF)

  integer            :: Index, INIT, N, M, L_D, L_F, L_DF
  REALTYPE           :: alpha, beta

  REALTYPE, dimension(*) :: D, F, DF

  SELECT CASE (Index)
    CASE (1)
      call P_GEMM_1   (INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF)

    CASE (2)
      call P_GEMM_2   (INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF)

  END SELECT

  END Subroutine P_GEMM 

! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_1
! Author       : Wai Sun Don
! Description  : C = alpha*A*B+beta*C
!
!  Global : A(K_A,N), B(K_B,M), C(K_C,M)
!  Local  : A(L_A,N), B(L_B,M), C(L_C,M)
!
!  Data B and C are distributed along the FIRST index
!
!  INIT : 1  Reorganize array A into AA with proper size
!              if N_Processor>1 and mod(N,N_Processor)>0
! 
!         0  if P_GEMM_INIT_1 was called, or
!              if N_Processor=1 or mod(N,N_Processor)=0
!
!  K_A  is the Global leading dimension of Square array A
!  K_B  is the Global leading dimension of the data array B
!  K_C  is the Global leading dimension of the data array C
!  L_A  is the Local  leading dimension of Square array A
!  L_B  is the Local  leading dimension of the data array B
!  L_C  is the Local  leading dimension of the data array C
!  
!  N      is the number of Rows for arrays A, B and C (Global)
!  M      is the number of columns(vectors) (Global)
!  n_Star is the number of rows distributed across N_Processors
!            used for calculation (Local)
!
!  N      <= K_A
!  N      <= K_B
!  N      <= K_C
!  n_Star <= L_A 
!  n_Star <= L_B
!  n_Star <= L_C
!
!  Data A, B, C are distributed in BLOCK of rows
!    see subroutine PS_MPI_Block_Distribution for details.
!              
! ----------------------------------------------------------------------

  Subroutine P_GEMM_1   (INIT, N, M, alpha, A, L_A, B, L_B, beta, C, L_C) 

  integer            :: INIT, N, M, L_A, L_B, L_C
  REALTYPE           :: alpha, beta

  REALTYPE, dimension(*)     :: A, B, C

  if (N_Processor == 1) then
    call S_GEMM_1 (N, M, alpha, A, L_A, B, L_B, beta, C, L_C)
  else
    if (MOD(N, N_Processor) == 0) then
      call P_GEMM_1_N (INIT, N, M, alpha, A, L_A, B, L_B, beta, C, L_C) 
    else
      call P_GEMM_1_G (INIT, N, M, alpha, A, L_A, B, L_B, beta, C, L_C) 
    endif
  endif

  END Subroutine P_GEMM_1 

! ----------------------------------------------------------------------
! FUNCTION NAME: S_GEMM_1
! ----------------------------------------------------------------------
  Subroutine S_GEMM_1 (N, M, alpha, A, L_A, B, L_B, beta, C, L_C)

  integer  :: N, M
  REALTYPE :: alpha, beta

  integer                    :: Row_A, Col_A, Col_C
  integer                    ::   L_A,   L_B,   L_C
  REALTYPE, dimension(*)     ::     A,     B,     C

  Row_A = N ; Col_A = N ; Col_C = M

  call   GEMM ('N', 'N', Row_A, Col_C, Col_A, alpha,      &
                             A, L_A, B, L_B, beta, C, L_C)

  END Subroutine S_GEMM_1 
! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_1_N
! Note : N is divisible by N_Processor
! ----------------------------------------------------------------------
  Subroutine P_GEMM_1_N (INIT, N, M, alpha, A, L_A, B, L_B, beta, C, L_C) 

  integer            :: INIT, N, M, L_A, L_B, L_C 
  REALTYPE           :: alpha, beta

  integer  :: n_Block
  integer  :: NP_Row, NP_Col, My_Row, My_Col, IContxt

  character                  :: Order
  integer                    ::  Row_A,  Col_A,  Col_C
  integer , dimension(8)     :: Desc_A, Desc_B, Desc_C
  REALTYPE, dimension(*)     ::      A,      B,      C

  n_Block = N/N_Processor

  Order = 'C' ; NP_Row = N_Processor ; NP_Col = 1

  call BLACS_Get      (0, 0, IContxt)
  call BLACS_GridInit (IContxt, Order, NP_Row,NP_Col)
  call BLACS_GridInfo (IContxt,        NP_Row,NP_Col, My_Row,My_Col)

  Desc_A(1) = N        ; Desc_B(1) = N        ; Desc_C(1) = N
  Desc_A(2) = N        ; Desc_B(2) = M        ; Desc_C(2) = M
  Desc_A(3) = n_Block  ; Desc_B(3) = n_Block  ; Desc_C(3) = n_Block
  Desc_A(4) = n_Block  ; Desc_B(4) = M        ; Desc_C(4) = M
  Desc_A(5) = 0        ; Desc_B(5) = 0        ; Desc_C(5) = 0
  Desc_A(6) = 0        ; Desc_B(6) = 0        ; Desc_C(6) = 0
  Desc_A(7) = IContxt  ; Desc_B(7) = IContxt  ; Desc_C(7) = IContxt
  Desc_A(8) = L_A      ; Desc_B(8) = L_B      ; Desc_C(8) = L_C

  Row_A = N ; Col_A = N ; Col_C = M

  call PDGEMM ('N', 'N', Row_A, Col_C, Col_A, alpha,          &
                         A, 1,1,Desc_A, B, 1,1,Desc_B, beta,  &
                         C, 1,1,Desc_C)

  call BLACS_GridExit (IContxt)

  END Subroutine P_GEMM_1_N 
! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_1_G
! Note : General Case for N that might or might not be divisible by N_Processor
! ----------------------------------------------------------------------
  Subroutine P_GEMM_1_G (INIT, N, M, alpha, A, L_A, B, L_B, beta, C, L_C) 

  integer            :: INIT, N, M, L_A, L_B, L_C
  REALTYPE           :: alpha, beta

  integer  :: n_Remain, n_0, n_Star, n_Block, NN, MM
  integer  :: NP_Row, NP_Col, My_Row, My_Col, IContxt

  character                  :: Order
  integer                    ::  Row_A,  Col_A,  Col_C
  integer , dimension(8)     :: Desc_A, Desc_B, Desc_C
  REALTYPE, dimension(*)     ::      A,      B,      C

  REALTYPE, dimension(:,:), ALLOCATABLE :: AA, BB, CC

  call PS_MPI_Block_Distribution (N, n_Remain, n_0, n_Star, n_Block, NN)

  MM = M

  Order = 'C' ; NP_Row = N_Processor ; NP_Col = 1

  call BLACS_Get      (0, 0, IContxt)
  call BLACS_GridInit (IContxt, Order, NP_Row,NP_Col)
  call BLACS_GridInfo (IContxt,        NP_Row,NP_Col, My_Row,My_Col)

  Desc_A(1) = NN       ; Desc_B(1) = NN       ; Desc_C(1) = NN
  Desc_A(2) = NN       ; Desc_B(2) = MM       ; Desc_C(2) = MM
  Desc_A(3) = n_Block  ; Desc_B(3) = n_Block  ; Desc_C(3) = n_Block
  Desc_A(4) = n_Block  ; Desc_B(4) = MM       ; Desc_C(4) = MM
  Desc_A(5) = 0        ; Desc_B(5) = 0        ; Desc_C(5) = 0
  Desc_A(6) = 0        ; Desc_B(6) = 0        ; Desc_C(6) = 0
  Desc_A(7) = IContxt  ; Desc_B(7) = IContxt  ; Desc_C(7) = IContxt
  Desc_A(8) = n_Block  ; Desc_B(8) = n_Block  ; Desc_C(8) = n_Block

  Row_A = NN ; Col_A = NN ; Col_C = MM

  ALLOCATE (BB(n_Block,MM), CC(n_Block,MM), STAT=M_Error)

#if defined (DEBUG)
  Memory_Requested = 2*n_Block*MM
  Subroutine_Name  = TRIM(Subroutine_Call)//'P_GEMM_1_G : BB, CC'

  call Memory_Error_Check (0)
#endif

    call P_Expand_B_1 (L_B, n_Star, n_block, MM, B, BB)
  if (beta /= ZERO) then
    call P_Expand_C_1 (L_C, n_Star, n_block, MM, C, CC)
  else
    CC = ZERO
  endif

  if (INIT == 1) then
    ALLOCATE (AA(n_Block,NN), STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = n_Block*NN
    Subroutine_Name  = TRIM(Subroutine_Call)//'P_GEMM_1_G : AA'

    call Memory_Error_Check (0)
#endif

    call P_Expand_A_1 (L_A, n_Star, N, n_Block, NN, n_Remain, A, AA)

    call PDGEMM ('N', 'N', Row_A, Col_C, Col_A, alpha,             &
                           AA, 1,1,Desc_A,  BB, 1,1,Desc_B, beta,  &
                           CC, 1,1,Desc_C)

    DEALLOCATE (AA, STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = n_Block*NN
    Subroutine_Name  = TRIM(Subroutine_Call)//'P_GEMM_1_G : AA'

    call Memory_Error_Check (1)
#endif
  else
    call PDGEMM ('N', 'N', Row_A, Col_C, Col_A, alpha,             &
                            A, 1,1,Desc_A,  BB, 1,1,Desc_B, beta,  &
                           CC, 1,1,Desc_C)
  endif

  call P_Restore_C_1 (L_C, n_Block, n_Star, MM, CC, C)

  DEALLOCATE (CC, BB, STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = 2*n_Block*MM
    Subroutine_Name  = TRIM(Subroutine_Call)//'P_GEMM_1_G : BB, CC'

    call Memory_Error_Check (1)
#endif

  call BLACS_GridExit (IContxt)

  END Subroutine P_GEMM_1_G
!
!=======================================================================
!
  Subroutine P_Expand_A_1 (L_A, n_Star, N, n_Block, NN, n_Remain, A, AA)

  integer  :: L_A, n_Star, N, n_Block, NN, n_Remain
  integer  :: i, k, ii, i_Block, NP

  REALTYPE :: A(L_A,N), AA(n_Block,NN)

  if (N < N_Processor) AA(1:n_Block,:) = ZERO

  if (n_Remain == 0) then
    NP = N_Processor ; else ; NP = n_Remain
  endif

  i = NP*n_Block

  AA(1:n_Star,1:i) = A(:,1:i)

  if (i < N) then
    do k = NP, N_Processor-1
 
      do i_Block = 1,n_Block-1
        i = i+1 ; ii = i_Block+k*n_Block

        AA(1:n_Star,ii) = A(:,i)
      enddo

      AA(:,(k+1)*n_Block) = ZERO
    enddo
  endif

  if (n_Star < n_Block) AA(n_Star+1:n_Block,:) = ZERO

  END Subroutine P_Expand_A_1 
!
!=======================================================================
!
  Subroutine P_Expand_B_1 (L_B, n_Star, n_block, MM, B, BB)

  Integer  ::   L_B,MM,     n_Block,   n_Star
  REALTYPE :: B(L_B,MM), BB(n_Block,MM)

  BB(1:n_Star,:) = B(1:n_Star,:)

  if (n_Star < n_Block) BB(n_Star+1:n_Block,:) = ZERO

  END Subroutine P_Expand_B_1 
!
!=======================================================================
!
  Subroutine P_Expand_C_1 (L_C, n_Star, n_block, MM, C, CC)

  Integer  ::   L_C,MM,     n_Block,   n_Star
  REALTYPE :: C(L_C,MM), CC(n_Block,MM)

  CC(1:n_Star,:) = C(1:n_Star,:)

  if (n_Star < n_Block) CC(n_Star+1:n_Block,:) = ZERO

  END Subroutine P_Expand_C_1 
!
!=======================================================================
!
  Subroutine P_Restore_C_1 (L_C, n_Block, n_Star, MM, CC, C)

  Integer  ::    n_Block,MM,    L_C,   n_Star
  REALTYPE :: CC(n_Block,MM), C(L_C,MM)

  C(1:n_Star,:) = CC(1:n_Star,:)

  END Subroutine P_Restore_C_1 

! ----------------------------------------------------------------------
! FUNCTION NAME: Copy_A_2_AA  
! ----------------------------------------------------------------------

  Subroutine Copy_A_2_AA (L_A, n_Star, N, n_Block, NN, A, AA)

  Integer  ::   L_A,N,     n_Block,NN, n_Star
  REALTYPE :: A(L_A,N), AA(n_Block,NN)

  AA(1:n_Star,1:N) = A(1:n_Star,:)  

  if (N < NN) AA(1:n_Star,N+1:NN) = ZERO

  if (n_Star < n_Block) AA(n_Star+1:n_Block,:) = ZERO

  END Subroutine Copy_A_2_AA 







! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_2
! Author       : Wai Sun Don
! Description  : C = alpha*Transpose(B*Transpose(A))+beta*C
!                  = alpha*A*Transpose(B)+beta*C
!
!  Global : A(K_A,M)     , B(K_B,M)     , C(K_C,M)
!  Local  : A(L_A,m_Star), B(L_B,m_Star), C(L_C,m_Star)
!
!  Data A and C are distributed along the SECOND index
!
!  INIT : 1  Reorganize array B into BB with proper size
!              if N_Processor>1 and mod(N,N_Processor)>0
!         0  if P_GEMM_INIT_2 was called, or
!              if N_Processor=1 or mod(N,N_Processor)=0
!
!  K_B  is the Global leading dimension of Square array B
!  K_A  is the Global leading dimension of the data array A
!  K_C  is the Global leading dimension of the data array C
!  L_B  is the Local  leading dimension of Square array B
!  L_A  is the Local  leading dimension of the data array A
!  L_C  is the Local  leading dimension of the data array C
!
!  N      is the number of Rows for arrays A and C (Global)
!  M      is the number of columns (vectors) (Global)
!  m_Star is the number of Cols distributed across N_Processors
!            used for calculation (Local)
!
!  M      <= K_B
!  N      <= K_A
!  N      <= K_C
!  M      <= L_B
!  N      <= L_A
!  N      <= L_C
!
!  Data A, B, C are distributed in BLOCK of Columns
!    see subroutine PS_MPI_Block_Distribution for details.
!
! ----------------------------------------------------------------------

  Subroutine P_GEMM_2   (INIT, N, M, alpha, B, L_B, A, L_A, beta, C, L_C) 

  integer            :: INIT, N, M, L_A, L_B, L_C
  REALTYPE           :: alpha, beta

  REALTYPE, dimension(*)     :: A, B, C

  if (N_Processor == 1) then
    call S_GEMM_2 (N, M, alpha, B, L_B, A, L_A, beta, C, L_C)
  else
    if (MOD(M, N_Processor) == 0) then
      call P_GEMM_2_M (INIT, N, M, alpha, B, L_B, A, L_A, beta, C, L_C) 
    else
      call P_GEMM_2_G (INIT, N, M, alpha, B, L_B, A, L_A, beta, C, L_C) 
    endif
  endif

  END Subroutine P_GEMM_2 
! ----------------------------------------------------------------------
! FUNCTION NAME: S_GEMM_2
! ----------------------------------------------------------------------
  Subroutine S_GEMM_2 (N, M, alpha, B, L_B, A, L_A, beta, C, L_C)

  integer  :: N, M
  REALTYPE :: alpha, beta

  integer                    :: Row_A, Col_A, Col_C
  integer                    ::   L_A,   L_B,   L_C
  REALTYPE, dimension(*)     ::     A,     B,     C

  Row_A = N ; Col_A = M ; Col_C = M

  call   GEMM ('N', 'T', Row_A, Col_C, Col_A, alpha,      &
                             A, L_A, B, L_B, beta, C, L_C)

  END Subroutine S_GEMM_2 
! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_2_M
! Note : M is divisible by N_Processor
! ----------------------------------------------------------------------
  Subroutine P_GEMM_2_M (INIT, N, M, alpha, B, L_B, A, L_A, beta, C, L_C) 

  integer            :: INIT, N, M, L_A, L_B, L_C
  REALTYPE           :: alpha, beta

  integer  :: m_Block
  integer  :: NP_Row, NP_Col, My_Row, My_Col, IContxt

  character                  :: Order
  integer                    ::  Row_A,  Col_A,  Col_C
  integer , dimension(8)     :: Desc_A, Desc_B, Desc_C
  REALTYPE, dimension(*)     ::      A,      B,      C

  m_Block = M/N_Processor

  Order = 'R' ; NP_Row = 1 ; NP_Col = N_Processor
  
  call BLACS_Get      (0, 0, IContxt)
  call BLACS_GridInit (IContxt, Order, NP_Row,NP_Col)
  call BLACS_GridInfo (IContxt,        NP_Row,NP_Col, My_Row,My_Col)

  Desc_A(1) = L_A      ; Desc_B(1) = M        ; Desc_C(1) = L_C
  Desc_A(2) = M        ; Desc_B(2) = M        ; Desc_C(2) = M
  Desc_A(3) = N        ; Desc_B(3) = m_Block  ; Desc_C(3) = N
  Desc_A(4) = m_Block  ; Desc_B(4) = m_Block  ; Desc_C(4) = m_Block
  Desc_A(5) = 0        ; Desc_B(5) = 0        ; Desc_C(5) = 0
  Desc_A(6) = 0        ; Desc_B(6) = 0        ; Desc_C(6) = 0
  Desc_A(7) = IContxt  ; Desc_B(7) = IContxt  ; Desc_C(7) = IContxt
  Desc_A(8) = L_A      ; Desc_B(8) = L_B      ; Desc_C(8) = L_C

  Row_A = N ; Col_A = M ; Col_C = M

  call PDGEMM ('N', 'T', Row_A, Col_C, Col_A, alpha,          &
                         A, 1,1,Desc_A, B, 1,1,Desc_B, beta,  &
                         C, 1,1,Desc_C)

  call BLACS_GridExit (IContxt)

  END Subroutine P_GEMM_2_M 
! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_2_G
! Note : General M that might or might not be divisible by N_Processor
! ----------------------------------------------------------------------
  Subroutine P_GEMM_2_G (INIT, N, M, alpha, B, L_B, A, L_A, beta, C, L_C) 

  integer            :: INIT, N, M, L_A, L_B, L_C
  REALTYPE           :: alpha, beta

  integer  :: m_Remain, m_0, m_Star, m_Block, MM, NN
  integer  :: NP_Row, NP_Col, My_Row, My_Col, IContxt

  character                  :: Order
  integer                    ::  Row_A,  Col_A,  Col_C
  integer , dimension(8)     :: Desc_A, Desc_B, Desc_C
  REALTYPE, dimension(*)     ::      A,      B,      C

  REALTYPE, dimension(:,:), ALLOCATABLE :: AA, BB, CC

  call PS_MPI_Block_Distribution (M, m_Remain, m_0, m_Star, m_Block, MM)

  NN = N 

  Order = 'R' ; NP_Row = 1 ; NP_Col = N_Processor
  
  call BLACS_Get      (0, 0, IContxt)
  call BLACS_GridInit (IContxt, Order, NP_Row,NP_Col)
  call BLACS_GridInfo (IContxt,        NP_Row,NP_Col, My_Row,My_Col)

  Desc_A(1) = NN       ; Desc_B(1) = MM       ; Desc_C(1) = NN
  Desc_A(2) = MM       ; Desc_B(2) = MM       ; Desc_C(2) = MM
  Desc_A(3) = NN       ; Desc_B(3) = m_Block  ; Desc_C(3) = NN
  Desc_A(4) = m_Block  ; Desc_B(4) = m_Block  ; Desc_C(4) = m_Block
  Desc_A(5) = 0        ; Desc_B(5) = 0        ; Desc_C(5) = 0
  Desc_A(6) = 0        ; Desc_B(6) = 0        ; Desc_C(6) = 0
  Desc_A(7) = IContxt  ; Desc_B(7) = IContxt  ; Desc_C(7) = IContxt
  Desc_A(8) = NN       ; Desc_B(8) = MM       ; Desc_C(8) = NN

  Row_A = NN ; Col_A = MM ; Col_C = MM

  ALLOCATE (AA(NN,m_Block), CC(NN,m_Block), STAT=M_Error)

#if defined (DEBUG)
  Memory_Requested = 2*NN*m_Block
  Subroutine_Name  = TRIM(Subroutine_Call)//'P_GEMM_2_G : AA, CC'

  call Memory_Error_Check (0)
#endif

    call P_Expand_A_2 (L_A, m_Star, NN, m_Block, A, AA)
  if (beta /= ZERO) then
    call P_Expand_C_2 (L_C, m_Star, NN, m_Block, C, CC)
  else
    CC = ZERO
  endif

  if (INIT == 1) then
    ALLOCATE (BB(MM,m_Block), STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = MM*m_Block
    Subroutine_Name  = TRIM(Subroutine_Call)//'P_GEMM_2_G : BB'

    call Memory_Error_Check (0)
#endif

    call P_Expand_B_2 (L_B, M, m_Star, MM, m_Block, m_Remain, B, BB)

    call PDGEMM ('N', 'T', Row_A, Col_C, Col_A, alpha,             &
                           AA, 1,1,Desc_A,  BB, 1,1,Desc_B, beta,  &
                           CC, 1,1,Desc_C)

    DEALLOCATE (BB, STAT=M_Error)

#if defined (DEBUG)
    Memory_Requested = MM*m_Block
    Subroutine_Name  = TRIM(Subroutine_Call)//'P_GEMM_2_G : BB'

    call Memory_Error_Check (1)
#endif
  else
    call PDGEMM ('N', 'T', Row_A, Col_C, Col_A, alpha,             &
                           AA, 1,1,Desc_A,  B , 1,1,Desc_B, beta,  &
                           CC, 1,1,Desc_C)
  endif

  call P_Restore_C_2 (NN , m_Block, L_C, m_Star, CC, C)

  DEALLOCATE (CC, AA, STAT=M_Error)

#if defined (DEBUG)
  Memory_Requested = 2*NN*m_Block
  Subroutine_Name  = TRIM(Subroutine_Call)//'P_GEMM_2_G : AA, CC'

  call Memory_Error_Check (1)
#endif

  call BLACS_GridExit (IContxt)

  END Subroutine P_GEMM_2_G
!
!=======================================================================
!
  Subroutine P_Expand_A_2 (L_A, m_Star, N, m_Block, A, AA)

  Integer  ::   L_A,m_Star,     N,m_Block
  REALTYPE :: A(L_A,m_Star), AA(N,m_Block)

  AA(:,1:m_Star) = A(1:N,:)

  if (m_Star < m_Block) AA(:,m_Star+1:m_Block) = ZERO

  END Subroutine P_Expand_A_2 
!
!=======================================================================
!
  Subroutine P_Expand_B_2 (L_B, M, m_Star, MM, m_Block, m_Remain, B, BB)

  integer  :: L_B, M, m_Star, MM, m_Block, m_Remain
  integer  :: j, k, jj, j_Block, NP

  REALTYPE :: B(L_B,m_Star), BB(MM,m_Block)

  if (M < N_Processor) BB(M+1:MM,1:m_Star) = ZERO

  if (m_Remain == 0) then
    NP = N_Processor ; else ; NP = m_Remain
  endif

  j = NP*m_Block

  BB(1:j,1:m_Star) = B(1:j,:)

  if (j < M) then
    do k = NP, N_Processor-1
 
      do j_Block = 1,m_Block-1
        j = j+1 ; jj = j_Block+k*m_Block

        BB(jj,1:m_Star) = B(j,:)
      enddo

      BB((k+1)*m_Block,:) = ZERO
    enddo
  endif

  if (m_Star < m_Block) BB(:,m_Star+1:m_Block) = ZERO

  END Subroutine P_Expand_B_2 
!
!=======================================================================
!
  Subroutine P_Expand_C_2 (L_C, m_Star, N, m_Block, C, CC)

  Integer  ::   L_C,m_Star,     N,m_Block
  REALTYPE :: C(L_C,m_Star), CC(N,m_Block)

  CC(:,1:m_Star) = C(1:N,:)

  if (m_Star < m_Block) CC(:,m_Star+1:m_Block) = ZERO

  END Subroutine P_Expand_C_2 
!
!=======================================================================
!
  Subroutine P_Restore_C_2 (N, m_Block, L_C, m_Star, CC,C)

  Integer  ::    N,m_Block,    L_C,m_Star
  REALTYPE :: CC(N,m_Block), C(L_C,m_Star)

  C(1:N,:) = CC(:,1:m_Star)

  END Subroutine P_Restore_C_2 
! ----------------------------------------------------------------------
! FUNCTION NAME: Copy_B_2_BB  
! ----------------------------------------------------------------------

  Subroutine Copy_B_2_BB (L_B, M, m_Star, MM, m_Block, B, BB)

  Integer  ::   L_B,m_Star,     MM,m_Block, M
  REALTYPE :: B(L_B,m_Star), BB(MM,m_Block)

  BB(1:M,1:m_Star) = B(1:M,:)

  if (M < MM) BB(M+1:MM,1:m_Star) = ZERO

  if (m_Star < m_Block) BB(:,m_Star+1:m_Block) = ZERO

  END Subroutine Copy_B_2_BB 

  END MODULE Matrix_Multiply_PESSL

!-----
#else
!-----

  MODULE Matrix_Multiply_PESSL

#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

  implicit NONE

  PRIVATE 
  PUBLIC  :: P_GEMM_INIT
  PUBLIC  :: P_GEMM

  CONTAINS
!
!=======================================================================
!
  Subroutine P_GEMM_INIT (Index, INIT, N, M, A, L_A, N_A, B, L_B, N_B,   &
                                                                ERROR)

  integer                 :: Index, INIT, N, M, L_A, N_A, L_B, N_B 
  integer ,      OPTIONAL :: ERROR 
  REALTYPE, dimension(*)  :: A, B

  RETURN

  END Subroutine P_GEMM_INIT 
!
!=======================================================================
!
  Subroutine P_GEMM (Index,   &
                     INIT, N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF)

  integer            :: Index, INIT, N, M, L_D, L_F, L_DF
  REALTYPE           :: alpha, beta

  REALTYPE, dimension(*) :: D, F, DF

  SELECT CASE (Index)
    CASE (1)
      call P_GEMM_1 (N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF)

    CASE (2)
      call P_GEMM_2 (N, M, alpha, D, L_D, F, L_F, beta, DF, L_DF)

  END SELECT

  END Subroutine P_GEMM 
! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_1
! ----------------------------------------------------------------------
  Subroutine P_GEMM_1 (N, M, alpha, A, L_A, B, L_B, beta, C, L_C)

  integer  :: N, M
  REALTYPE :: alpha, beta

  integer                    :: Row_A, Col_A, Col_C
  integer                    ::   L_A,   L_B,   L_C
  REALTYPE, dimension(*)     ::     A,     B,     C

  Row_A = N ; Col_A = N ; Col_C = M

  call   GEMM ('N', 'N', Row_A, Col_C, Col_A, alpha,       &
                             A, L_A, B, L_B, beta, C, L_C)

  END Subroutine P_GEMM_1
! ----------------------------------------------------------------------
! FUNCTION NAME: P_GEMM_2
! ----------------------------------------------------------------------
  Subroutine P_GEMM_2 (N, M, alpha, B, L_B, A, L_A, beta, C, L_C)

  integer  :: N, M
  REALTYPE :: alpha, beta

  integer                    :: Row_A, Col_A, Col_C
  integer                    ::   L_A,   L_B,   L_C
  REALTYPE, dimension(*)     ::     A,     B,     C

  Row_A = N ; Col_A = M ; Col_C = M

  call   GEMM ('N', 'T', Row_A, Col_C, Col_A, alpha,      &
                             A, L_A, B, L_B, beta, C, L_C)

  END Subroutine P_GEMM_2

END MODULE Matrix_Multiply_PESSL
!-----
#endif
!-----
