#include "pseudopack.h"


MODULE Smoothing_EOD_1
#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE

PRIVATE

PUBLIC  :: Smooth_EOD_1

CONTAINS

  Subroutine Smooth_EOD_1      (N, M, y, ldy, Aux_Smth, IEO)

  integer  :: N, M, ldy

  integer , OPTIONAL :: IEO
  logical            :: L_IEO

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call Smooth_EOD_1_BOTH (N, M, y, ldy, Aux_Smth, IEO)
    else
      call Smooth_EOD_1_NONE (N, M, y, ldy, Aux_Smth)
    endif
  else
      call Smooth_EOD_1_NONE (N, M, y, ldy, Aux_Smth)
  endif

  END Subroutine Smooth_EOD_1

! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_EOD_1_NONE
! ----------------------------------------------------------------------
  Subroutine Smooth_EOD_1_NONE (N, M, y, ldy, Aux_Smth)

  integer  :: N, M, ldy
  integer  :: i

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: Method, BlockSize
  integer  :: N0, N1, N2
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_e_s, i_o_s

#if defined (PARALLEL_OPENMP)
  m_Block   = OPERATOR_NUMBER_THREADS(Aux_Smth)
  BlockSize = M/m_Block
#else
  BlockSize = OPERATOR_BLOCKSIZE(Aux_Smth)

  if (BlockSize == 0) BlockSize = M

  BlockSize = MIN(BlockSize, M)
  m_Block   = M/BlockSize
#endif

  m_Left    = M - m_Block*BlockSize

     Method = OPERATOR_METHOD(Aux_Smth)

  SELECT CASE (Method)
    CASE (0)
      N0 = N   ; N1 = N+1 

    CASE (1,2)
      N0 = N-1 ; N1 = N   

  END SELECT

  i_e_s  = EOD_MATRIX_START_EVEN(N1)
  i_o_s  = EOD_MATRIX_START_EVEN(N1) + EOD_HALF_SIZE(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize 
    m_Start  = i*BlockSize+1

    call EOD_Smooth_1 (N0, m_Vector, Aux_Smth(i_e_s), Aux_Smth(i_o_s),  &
                                     y(1,m_Start))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call EOD_Smooth_1 (N0, m_Vector, Aux_Smth(i_e_s), Aux_Smth(i_o_s),  &
                                     y(1,m_Start))
  endif

  CONTAINS

    Subroutine EOD_Smooth_1 (N, M, S_Even, S_Odd, y)

    integer  :: N, M
    integer  :: i, j
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(0:ldy-1,M)   ::  y
    REALTYPE, dimension(0:N/2,0:N/2) ::  S_Even,  S_Odd
    REALTYPE, dimension(0:N/2,M)     ::  y_Even,  y_Odd
    REALTYPE, dimension(0:N/2,M)     :: Sy_Even, Sy_Odd

      N2 = N/2

    LD_A = SIZE( S_Even, DIM=1) ; Row_A = N/2+1 ; Col_A = N/2+1
    LD_B = SIZE( y_Even, DIM=1)
    LD_C = SIZE(Sy_Even, DIM=1) ;                 Col_C = M

    SELECT CASE (Method)
      CASE (0)
        y_Even(0,:) = y(0,:) + y(0,:)

      CASE (1,2)
        y_Even(0,:) = y(0,:) + y(N,:)
        y_Odd (0,:) = y(0,:) - y(N,:)

    END SELECT

    do j = 1, M
      do i = 1,N2
        y_Even(i,j) = y(i,j) + y(N-i,j)
        y_Odd (i,j) = y(i,j) - y(N-i,j)
      enddo
    enddo

    call GEMM ('n', 'n', Row_A, Col_C, Col_A, ONE, S_Even, LD_A,  &
                         y_Even, LD_B, ZERO, Sy_Even, LD_C)
    call GEMM ('n', 'n', Row_A, Col_C, Col_A, ONE, S_Odd , LD_A,  &
                         y_Odd , LD_B, ZERO, Sy_Odd , LD_C)

    SELECT CASE (Method)
      CASE (0)
        y(0,:) = Sy_Even(0,:)

      CASE (1,2)
        y(0,:) =  Sy_Even(0,:) + Sy_Odd(0,:)
        y(N,:) =  Sy_Even(0,:) - Sy_Odd(0,:)

    END SELECT

    if (MOD(N,2) == 0) then
      do j = 1, M
        do i = 1, N2-1
          y(  i,j) =  Sy_Even(i,j) + Sy_Odd(i,j)
          y(N-i,j) =  Sy_Even(i,j) - Sy_Odd(i,j)
        enddo
      enddo

      y(N2,:) = Sy_Even(N2,:)
    else
      do j = 1, M
        do i = 1, N2
          y(  i,j) =  Sy_Even(i,j) + Sy_Odd(i,j)
          y(N-i,j) =  Sy_Even(i,j) - Sy_Odd(i,j)
        enddo
      enddo
    endif

    END Subroutine EOD_Smooth_1 

  END Subroutine Smooth_EOD_1_NONE 
!
! ----------------------------------------------------------------------
!
  Subroutine Smooth_EOD_1_BOTH (N, M, y, ldy, Aux_Smth, IEO)

  integer  :: N, M, ldy, IEO
  integer  :: i

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: Method, BlockSize
  integer  :: N0, N1
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_s

#if defined (PARALLEL_OPENMP)
  m_Block   = OPERATOR_NUMBER_THREADS(Aux_Smth)
  BlockSize = M/m_Block
#else
  BlockSize = OPERATOR_BLOCKSIZE(Aux_Smth)

  if (BlockSize == 0) BlockSize = M

  BlockSize = MIN(BlockSize, M)
  m_Block   = M/BlockSize
#endif

  m_Left    = M - m_Block*BlockSize

     Method = OPERATOR_METHOD(Aux_Smth)

  SELECT CASE (Method)
    CASE (0)
      N0 = N   ; N1 = N+1

    CASE (1,2)
      N0 = N-1 ; N1 = N
        
  END SELECT

  if (IEO ==  1) i_s = EOD_MATRIX_START_EVEN(N1)
  if (IEO == -1) i_s = EOD_MATRIX_START_EVEN(N1) + EOD_HALF_SIZE(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call EOD_Smooth_1 (N0, m_Vector, Aux_Smth(i_s), y(1,m_Start))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call EOD_Smooth_1 (N0, m_Vector, Aux_Smth(i_s), y(1,m_Start))
  endif

  if (IEO == -1) then
    if (Method    == 0) y(1    ,1:M) = ZERO
    if (MOD(N0,2) == 0) y(N/2+1,1:M) = ZERO
  endif

  CONTAINS

    Subroutine EOD_Smooth_1 (N, M, S, y)

    integer  :: N, M, N2
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(0:ldy-1,M)   :: y
    REALTYPE, dimension(0:N/2,0:N/2) :: S
    REALTYPE, dimension(0:N/2,M)     :: Sy

       N2 = N/2

    LD_A = SIZE(S , DIM=1) ; Row_A = N/2+1 ; Col_A = N/2+1
    LD_B = SIZE(Sy, DIM=1)
    LD_C = SIZE(y , DIM=1) ;                 Col_C = M

    if ((IEO == -1) .AND. (MOD(N,2) == 0)) Row_A = Row_A-1

    Sy = y(0:N2,:)

    call GEMM ('n', 'n', Row_A, Col_C, Col_A,  &
                         ONE, S, LD_A, Sy, LD_B, ZERO, y, LD_C)

    END Subroutine EOD_Smooth_1

  END Subroutine Smooth_EOD_1_BOTH 

END MODULE Smoothing_EOD_1

! ----------------------------------------------------------------------

MODULE Smoothing_EOD_2
#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE

PRIVATE

PUBLIC  :: Smooth_EOD_2

CONTAINS

  Subroutine Smooth_EOD_2      (N, M, y, ldy, Aux_Smth, IEO)

  integer  :: N, M, ldy

  integer , OPTIONAL :: IEO
  logical            :: L_IEO

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call Smooth_EOD_2_BOTH (N, M, y, ldy, Aux_Smth, IEO)
    else
      call Smooth_EOD_2_NONE (N, M, y, ldy, Aux_Smth)
    endif
  else
      call Smooth_EOD_2_NONE (N, M, y, ldy, Aux_Smth)
  endif

  END Subroutine Smooth_EOD_2 
! ----------------------------------------------------------------------
  Subroutine Smooth_EOD_2_NONE (N, M, y, ldy, Aux_Smth)

  integer  :: N, M, ldy
  integer  :: i

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: Method, BlockSize
  integer  :: N0, N1, N2
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_e_s, i_o_s

#if defined (PARALLEL_OPENMP)
  m_Block   = OPERATOR_NUMBER_THREADS(Aux_Smth)
  BlockSize = M/m_Block
#else
  BlockSize = OPERATOR_BLOCKSIZE(Aux_Smth)

  if (BlockSize == 0) BlockSize = M

  BlockSize = MIN(BlockSize, M)
  m_Block   = M/BlockSize
#endif

  m_Left    = M - m_Block*BlockSize

     Method = OPERATOR_METHOD(Aux_Smth)

  SELECT CASE (Method)
    CASE (0)
      N0 = N   ; N1 = N+1 

    CASE (1,2)
      N0 = N-1 ; N1 = N   
    
  END SELECT

  i_e_s  = EOD_MATRIX_START_EVEN(N1)
  i_o_s  = EOD_MATRIX_START_EVEN(N1) + EOD_HALF_SIZE(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call EOD_Smooth_2 (N0, m_Vector, Aux_Smth(i_e_s), Aux_Smth(i_o_s),   &
                                     y(m_Start,1))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call EOD_Smooth_2 (N0, m_Vector, Aux_Smth(i_e_s), Aux_Smth(i_o_s),  &
                                     y(m_Start,1))
  endif

  CONTAINS

    Subroutine EOD_Smooth_2 (N, M, S_Even, S_Odd, y)

    integer  :: N, M
    integer  :: i, j
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(ldy,0:N)     ::  y
    REALTYPE, dimension(0:N/2,0:N/2) ::  S_Even,  S_Odd
    REALTYPE, dimension(0:N/2,M)     ::  y_Even,  y_Odd
    REALTYPE, dimension(0:N/2,M)     :: Sy_Even, Sy_Odd

       N2 = N/2

    LD_A = SIZE( S_Even, DIM=1) ; Row_A = N/2+1 ; Col_A = N/2+1
    LD_B = SIZE( y_Even, DIM=1)
    LD_C = SIZE(Sy_Even, DIM=1) ;                 Col_C = M

    SELECT CASE (Method)
      CASE (0)
        y_Even(0,:) = y(1:M,0) + y(1:M,0)
        y_Odd (0,:) = ZERO

      CASE (1,2)
        y_Even(0,:) = y(1:M,0) + y(1:M,N)
        y_Odd (0,:) = y(1:M,0) - y(1:M,N)

    END SELECT

    do i = 1,N2
      y_Even(i,:) = y(1:M,i) + y(1:M,N-i)
      y_Odd (i,:) = y(1:M,i) - y(1:M,N-i)
    enddo

    call GEMM ('n', 'n', Row_A, Col_C, Col_A, ONE, S_Even, LD_A,  &
                         y_Even, LD_B, ZERO, Sy_Even, LD_C)
    call GEMM ('n', 'n', Row_A, Col_C, Col_A, ONE, S_Odd , LD_A,  &
                         y_Odd , LD_B, ZERO, Sy_Odd , LD_C)

    SELECT CASE (Method)
      CASE (0)
        y(1:M,0) = Sy_Even(0,:)

      CASE (1,2)
        y(1:M,0) = Sy_Even(0,:) + Sy_Odd(0,:)
        y(1:M,N) = Sy_Even(0,:) - Sy_Odd(0,:)

    END SELECT

    if (MOD(N,2) == 0) then
      do j = 1,M
        do i = 1,N2-1
          y(j,  i) = Sy_Even(i,j) + Sy_Odd(i,j)
          y(j,N-i) = Sy_Even(i,j) - Sy_Odd(i,j)
        enddo
      enddo

      y(1:M,N2) = Sy_Even(N2,1:M)
    else
      do j = 1,M
        do i = 1,N2
          y(j,  i) = Sy_Even(i,j) + Sy_Odd(i,j)
          y(j,N-i) = Sy_Even(i,j) - Sy_Odd(i,j)
        enddo
      enddo
    endif

    END Subroutine EOD_Smooth_2 

  END Subroutine Smooth_EOD_2_NONE 

! ----------------------------------------------------------------------

  Subroutine Smooth_EOD_2_BOTH (N, M, y, ldy, Aux_Smth, IEO)

  integer  :: N, M, ldy, IEO
  integer  :: i

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: Method, BlockSize
  integer  :: N0, N1, N2
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_s

#if defined (PARALLEL_OPENMP)
  m_Block   = OPERATOR_NUMBER_THREADS(Aux_Smth)
  BlockSize = M/m_Block
#else
  BlockSize = OPERATOR_BLOCKSIZE(Aux_Smth)

  if (BlockSize == 0) BlockSize = M

  BlockSize = MIN(BlockSize, M)
  m_Block   = M/BlockSize
#endif

  m_Left    = M - m_Block*BlockSize

     Method = OPERATOR_METHOD(Aux_Smth)

  SELECT CASE (Method)
    CASE (0)
      N0 = N   ; N1 = N+1

    CASE (1,2)
      N0 = N-1 ; N1 = N

  END SELECT

  if (IEO ==  1) i_s = EOD_MATRIX_START_EVEN(N1)
  if (IEO == -1) i_s = EOD_MATRIX_START_EVEN(N1) + EOD_HALF_SIZE(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call EOD_Smooth_2 (N0, m_Vector, Aux_Smth(i_s), y(m_Start,1))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call EOD_Smooth_2 (N0, m_Vector, Aux_Smth(i_s), y(m_Start,1))
  endif

  if (IEO == -1) then
    if (Method    == 0) y(1:M,1    ) = ZERO
    if (MOD(N0,2) == 0) y(1:M,N/2+1) = ZERO
  endif

  CONTAINS

    Subroutine EOD_Smooth_2 (N, M, S, y)

    integer  :: N, M
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(M,0:N/2)     :: Sy
    REALTYPE, dimension(0:N/2,0:N/2) :: S
    REALTYPE, dimension(ldy,0:N/2)   :: y

       N2 = N/2

    LD_A = SIZE(Sy, DIM=1) ; Row_A = M     ; Col_A = N/2+1
    LD_B = SIZE(S , DIM=1)
    LD_C = SIZE(y , DIM=1) ;                 Col_C = N/2+1

    if ((IEO == -1) .AND. (MOD(N,2) == 0)) Col_C = Col_C-1

    Sy = y(1:M,:)

    call GEMM ('n', 't', Row_A, Col_C, Col_A,  &
                         ONE, Sy, LD_A, S, LD_B, ZERO, y, LD_C)

    END Subroutine EOD_Smooth_2

  END Subroutine Smooth_EOD_2_BOTH 

END MODULE Smoothing_EOD_2
