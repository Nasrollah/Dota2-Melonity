#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_MXM_1
! DESCRIPTION  : Differentiates the vectors in the Matrix f
!                using the plain Matrix-vector algorithm on the first index
!
!  N      = number of grid points
!  M      = number of vectors
!  ldy    = leading dimension of y and of dy
! ----------------------------------------------------------------------

MODULE Differentiation_MXM_1
#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE

PRIVATE

PUBLIC :: Diff_MXM_1

CONTAINS

  Subroutine Diff_MXM_1 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer           :: N, M, O_Of_D, ldy
  integer, OPTIONAL :: IEO

  REALTYPE, dimension(ldy,M) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: BlockSize, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_Offset, i_d_s

#if defined (PARALLEL_OPENMP)
  m_Block   = OPERATOR_NUMBER_THREADS(Aux_Diff)
  BlockSize = M/m_Block
#else
  BlockSize = OPERATOR_BLOCKSIZE(Aux_Diff)

  if (BlockSize == 0) BlockSize = M

  BlockSize = MIN(BlockSize, M)
  m_Block   = M/BlockSize
#endif

  m_Left    = M - m_Block*BlockSize

  i_Offset = (O_Of_D-1)*MXM_FULL_SIZE(N)
  i_d_s    = MXM_MATRIX_START(N) + i_Offset

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call MXM_Diff_1 (N, m_Vector, Aux_Diff(i_d_s), y(1,m_Start), dy(1,m_Start))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call MXM_Diff_1 (N, m_Vector, Aux_Diff(i_d_s), y(1,m_Start), dy(1,m_Start))
  endif

  CONTAINS

    Subroutine MXM_Diff_1 (N, M, D, y, dy)

    integer  :: N, M
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(N,N)   :: D
    REALTYPE, dimension(ldy,M) :: y, dy

    LD_A = SIZE(D , DIM=1) ; Row_A = N ; Col_A = N
    LD_B = SIZE(y , DIM=1)
    LD_C = SIZE(dy, DIM=1) ;             Col_C = M

    call GEMM ('n', 'n', Row_A, Col_C, Col_A,  &
                         ONE, D, LD_A, y, LD_B, ZERO, dy, LD_C)

    END Subroutine MXM_Diff_1 

  END Subroutine Diff_MXM_1 

END MODULE Differentiation_MXM_1

! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_MXM_2
! DESCRIPTION  : Differentiates the vectors in the Matrix f
!                using the plain Matrix-vector algorithm on the second index
!
!  N   = number of grid points
!  M   = number of vectors
!  ldy = leading dimension of y and of dy
! ----------------------------------------------------------------------

MODULE Differentiation_MXM_2
#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE

PRIVATE

PUBLIC :: Diff_MXM_2

CONTAINS

  Subroutine Diff_MXM_2 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer           :: N, M, O_Of_D, ldy
  integer, OPTIONAL :: IEO

  REALTYPE, dimension(ldy,N) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: BlockSize, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_Offset, i_d_s

#if defined (PARALLEL_OPENMP)
  m_Block   = OPERATOR_NUMBER_THREADS(Aux_Diff)
  BlockSize = M/m_Block
#else
  BlockSize = OPERATOR_BLOCKSIZE(Aux_Diff)

  if (BlockSize == 0) BlockSize = M

  BlockSize = MIN(BlockSize, M)
  m_Block   = M/BlockSize
#endif

  m_Left    = M - m_Block*BlockSize

  i_Offset = (O_Of_D-1)*MXM_FULL_SIZE(N)
  i_d_s    = MXM_MATRIX_START(N) + i_Offset

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call MXM_Diff_2 (N, m_Vector, Aux_Diff(i_d_s), y(m_Start,1), dy(m_Start,1))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call MXM_Diff_2 (N, m_Vector, Aux_Diff(i_d_s), y(m_Start,1), dy(m_Start,1))
  endif

  CONTAINS

    Subroutine MXM_Diff_2 (N, M, D, y, dy)

    integer  :: N, M
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(N,N)   :: D
    REALTYPE, dimension(ldy,N) :: y, dy

    LD_A = SIZE(y , DIM=1) ; Row_A = M ; Col_A = N
    LD_B = SIZE(D , DIM=1)
    LD_C = SIZE(dy, DIM=1) ;             Col_C = N

    call GEMM ('n', 't', Row_A, Col_C, Col_A,  &
                         ONE,  y, LD_A, D, LD_B, ZERO, dy, LD_C)

    END Subroutine MXM_Diff_2 

  END Subroutine Diff_MXM_2 

END MODULE Differentiation_MXM_2
