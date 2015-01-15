#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_MXM_1
! DESCRIPTION  : Smooth the vectors in the Matrix y
!                using the plain Matrix-vector algorithm
!
! Note : Smoothing is done on the first  index.
!
!      N   = number of grid points
!      M   = number of vectors
!      ldy = leading dimension of y and of dy
! ----------------------------------------------------------------------

MODULE Smoothing_MXM_1
#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE

PRIVATE

PUBLIC :: Smooth_MXM_1

CONTAINS

  Subroutine Smooth_MXM_1 (N, M, y, ldy, Aux_Smth, IEO)

  integer           :: N, M, ldy
  integer, OPTIONAL :: IEO

  REALTYPE, dimension(ldy,M) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_d_s

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

  i_d_s = MXM_MATRIX_START(N)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call MXM_Smooth_1 (N, m_Vector, Aux_Smth(i_d_s), y(1,m_Start))
  enddo
!$OMP END PARALLEL DO
 
  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call MXM_Smooth_1 (N, m_Vector, Aux_Smth(i_d_s), y(1,m_Start))
  endif

  CONTAINS

    Subroutine MXM_Smooth_1 (N, M, S, y)

    integer  :: N, M
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(N,N)   :: S
    REALTYPE, dimension(ldy,M) :: y
    REALTYPE, dimension(N,M)   :: Sy

    LD_A = SIZE(S , DIM=1) ; Row_A = N ; Col_A = N
    LD_B = SIZE(Sy, DIM=1)
    LD_C = SIZE( y, DIM=1) ;             Col_C = M

    Sy = y(1:N,:)

    call GEMM ('n', 'n', Row_A, Col_C, Col_A,  &
                         ONE, S, LD_A, Sy, LD_B, ZERO, y, LD_C)

    END Subroutine MXM_Smooth_1 

  END Subroutine Smooth_MXM_1 

END MODULE Smoothing_MXM_1



! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_MXM_2
! DESCRIPTION  : Smooth the vectors in the Matrix y
!                using the plain Matrix-vector algorithm
!
! Note : Smoothing is done on the second index.
!
!      N   = number of grid points
!      M   = number of vectors
!      ldy = leading dimension of y and of dy
! ----------------------------------------------------------------------

MODULE Smoothing_MXM_2
#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE

PRIVATE

PUBLIC :: Smooth_MXM_2

CONTAINS

  Subroutine Smooth_MXM_2 (N, M, y, ldy, Aux_Smth, IEO)

  integer           :: N, M, ldy
  integer, OPTIONAL :: IEO

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_d_s

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

  i_d_s = MXM_MATRIX_START(N)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call MXM_Smooth_2 (N, m_Vector, Aux_Smth(i_d_s), y(m_Start,1))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call MXM_Smooth_2 (N, m_Vector, Aux_Smth(i_d_s), y(m_Start,1))
  endif

  CONTAINS

    Subroutine MXM_Smooth_2 (N, M, S, y)

    integer  :: N, M
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(N,N)   :: S
    REALTYPE, dimension(ldy,N) :: y
    REALTYPE, dimension(M,N)   :: Sy

    LD_A = SIZE(Sy, DIM=1) ; Row_A = M ; Col_A = N
    LD_B = SIZE(S , DIM=1)
    LD_C = SIZE( y, DIM=1) ;             Col_C = N

    Sy = y(1:M,:)

    call GEMM ('n', 't', Row_A, Col_C, Col_A,  &
                         ONE, Sy, LD_A, S, LD_B, ZERO, y, LD_C)

    END Subroutine MXM_Smooth_2 

  END Subroutine Smooth_MXM_2 

END MODULE Smoothing_MXM_2
