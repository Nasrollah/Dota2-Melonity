#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_EOD_1
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : Differentiates the vectors in the Matrix y
!                using the Even-Odd algorithm along the First Index
!               
! Note : 
!        In this version, y and dy are assumed to contain contiguous vectors.
!
!        In order to minimize temporary storage, vectors are processed
!          in groups of size BlockSize.
! 
!   N = Total number of grid points before Symmetry
!
!   M = number of vectors
!
! IEO = Even_Or_Odd Control
!       = -1 Odd     functions : call Diff_EOD_1_BOTH
!       =  0 General functions : call Diff_EOD_1_NONE
!       =  1 Even    functions : call Diff_EOD_1_BOTH
!
!  Routine: Diff_EOD_1_BOTH do differentiation on Even/ODD functions
!           Diff_EOD_1_NONE do differentiation on General  functions
!
! ----------------------------------------------------------------------

MODULE Differentiation_EOD_1
#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE

PRIVATE 

PUBLIC  :: Diff_EOD_1

CONTAINS

  Subroutine Diff_EOD_1      (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer  :: M, N, O_Of_D, ldy

  integer , OPTIONAL :: IEO
  logical            :: L_IEO

  REALTYPE, dimension(ldy,M)   :: y, dy
  REALTYPE, dimension(*)       :: Aux_Diff

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call   Diff_EOD_1_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)
    else
      call   Diff_EOD_1_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
    endif
  else
      call   Diff_EOD_1_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
  endif

  END Subroutine Diff_EOD_1      

!
! ----------------------------------------------------------------------
!
  Subroutine Diff_EOD_1_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)

  integer  :: N, M, O_Of_D, ldy
  integer  :: i

  REALTYPE, dimension(ldy,M)   :: y, dy
  REALTYPE, dimension(*)       :: Aux_Diff

  integer  :: Method, BlockSize, EOD_Code
  integer  :: N0, N1, N2
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_Offset, i_d_e, i_d_o, i_m_s
  REALTYPE :: Factor

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

     Method = OPERATOR_METHOD(Aux_Diff)

  SELECT CASE (Method)
    CASE (0)
      N0 = N   ; N1 = N+1 

    CASE (1,2)
      N0 = N-1 ; N1 = N   

  END SELECT

  i_Offset = (O_Of_D-1)*EOD_HALF_SIZE(N1)
  i_d_e    = EOD_MATRIX_START_EVEN(N1) + i_Offset
  i_d_o    = EOD_MATRIX_START_ODD(N1)  + i_Offset

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call EOD_Diff_1_NONE (N0, m_Vector, Aux_Diff(i_d_e), Aux_Diff(i_d_o),  &
                                        y(1,m_Start)   , dy(1,m_Start)  )
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call EOD_Diff_1_NONE (N0, m_Vector, Aux_Diff(i_d_e), Aux_Diff(i_d_o),  &
                                        y(1,m_Start)   , dy(1,m_Start)  )
  endif

   EOD_Code = OPERATOR_EOD_CODE(Aux_Diff)

  if (EOD_Code == 1) then
    i_m_s = EOD_METRIC_START(N1)

    call Apply_Mapping_1 (N-1, M, ldy, Aux_Diff(i_m_s), dy)
  endif

  CONTAINS

    Subroutine EOD_Diff_1_NONE (N, M, D_Even, D_Odd, y, dy)
                                    
    integer  :: N, M
    integer  :: i, j
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(0:ldy-1,M)   :: y, dy
    REALTYPE, dimension(0:N/2,0:N/2) :: D_Even, D_Odd
    REALTYPE, dimension(0:N/2,M)     :: y_Even, y_Odd, dy_Even, dy_Odd

    Factor = (-1)**(O_Of_D-1)

      N2 = N/2

    LD_A = SIZE( D_Even, DIM=1) ; Row_A = N/2+1 ; Col_A = N/2+1
    LD_B = SIZE( y_Even, DIM=1)
    LD_C = SIZE(dy_Even, DIM=1) ;                 Col_C = M

    SELECT CASE (Method)
      CASE (0)
        y_Even(0,:) = y(0,:) + y(0,:)
        y_Odd (0,:) = ZERO

      CASE (1,2)
        y_Even(0,:) = y(0,:) + y(N,:)
        y_Odd (0,:) = y(0,:) - y(N,:)

    END SELECT

    do j = 1,M
      do i = 1,N2
        y_Even(i,j) = y(i,j) + y(N-i,j)
        y_Odd (i,j) = y(i,j) - y(N-i,j)
      enddo
    enddo

    call GEMM ('n', 'n', Row_A, Col_C, Col_A, ONE, D_Even, LD_A, &
                              y_Even, LD_B, ZERO, dy_Even, LD_C)
    call GEMM ('n', 'n', Row_A, Col_C, Col_A, ONE, D_Odd , LD_A, &
                              y_Odd , LD_B, ZERO, dy_Odd , LD_C)

    SELECT CASE (Method)
      CASE (0)
        dy(0,:) =          dy_Even(0,:) + dy_Odd(0,:)

      CASE (1,2)
        dy(0,:) =          dy_Even(0,:) + dy_Odd(0,:)
        dy(N,:) = Factor*(-dy_Even(0,:) + dy_Odd(0,:))

    END SELECT 

    do j = 1,M
      do i = 1,N2
        dy(  i,j) =          dy_Even(i,j) + dy_Odd(i,j)
        dy(N-i,j) = Factor*(-dy_Even(i,j) + dy_Odd(i,j))
      enddo
    enddo

    END Subroutine EOD_Diff_1_NONE 
! ----------------------------------------------------------------------
    Subroutine Apply_Mapping_1   (N, M, ldy, Slope, dy)

    integer  :: N, M, ldy
    integer  :: j

    REALTYPE, dimension(0:ldy-1,M) :: dy
    REALTYPE, dimension(0:N)       :: Slope

!$OMP PARALLEL DO
    do j = 1,M
      dy(0:N,j) = dy(0:N,j)*Slope
    enddo
!$OMP END PARALLEL DO

    END Subroutine Apply_Mapping_1   

  END Subroutine Diff_EOD_1_NONE 

! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_EOD_1_BOTH
! ----------------------------------------------------------------------
  Subroutine Diff_EOD_1_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer  :: N, M, O_Of_D, ldy, IEO
  integer  :: i

  REALTYPE, dimension(ldy,M)   :: y, dy
  REALTYPE, dimension(*)       :: Aux_Diff

  integer  :: Method, BlockSize
  integer  :: N0, N1
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_Offset, i_d

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

     Method = OPERATOR_METHOD(Aux_Diff)

  SELECT CASE (Method)
    CASE (0)
      N0 = N   ; N1 = N+1 

    CASE (1,2)
      N0 = N-1 ; N1 = N   

  END SELECT

  i_Offset = (O_Of_D-1)*EOD_HALF_SIZE(N1)

  if (IEO == 1) then
    i_d = EOD_MATRIX_START_EVEN(N1) + i_Offset
  else
    i_d = EOD_MATRIX_START_ODD(N1)  + i_Offset
  endif

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call EOD_Diff_BOTH_1 (N0, m_Vector, Aux_Diff(i_d),                &
                                        y(1,m_Start) , dy(1,m_Start))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call EOD_Diff_BOTH_1 (N0, m_Vector, Aux_Diff(i_d),                &
                                        y(1,m_Start) , dy(1,m_Start))
  endif

  if (IEO == 1) then
    if (MOD(O_Of_D,2) /= 0) then
      if (Method    == 0) dy(1    ,:) = ZERO
      if (MOD(N0,2) == 0) dy(N/2+1,:) = ZERO
    endif
  else
    if (MOD(O_Of_D,2) == 0) then
      if (Method    == 0) dy(1    ,:) = ZERO
      if (MOD(N0,2) == 0) dy(N/2+1,:) = ZERO
    endif
  endif

  CONTAINS

    Subroutine EOD_Diff_BOTH_1 (N, M, D, y, dy)

    integer  :: N, M
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(ldy,M)       :: y, dy
    REALTYPE, dimension(0:N/2,0:N/2) :: D   

    LD_A = SIZE(D , DIM=1) ; Row_A = N/2+1 ; Col_A = N/2+1
    LD_B = SIZE(y , DIM=1)
    LD_C = SIZE(dy, DIM=1) ;                 Col_C = M

    call GEMM ('n', 'n', Row_A, Col_C, Col_A,  &
                         ONE , D, LD_A,  y, LD_B, ZERO, dy, LD_C)

    END Subroutine EOD_Diff_BOTH_1 

  END Subroutine Diff_EOD_1_BOTH 

END MODULE Differentiation_EOD_1

! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_EOD_2
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : Differentiates the vectors in the Matrix y
!                using the Even-Odd algorithm along the Second Index
!
! Note : 
!        In this version, y and dy are assumed to contain contiguous vectors.
!
!        In order to minimize temporary storage, vectors are processed
!          in groups of size BlockSize.
!
!   N = Total number of grid points before Symmetry
!
!   M = number of vectors
!
! IEO = Even_Or_Odd Control
!       = -1 Odd     functions : call Diff_EOD_2_BOTH
!       =  0 General functions : call Diff_EOD_2_NONE
!       =  1 Even    functions : call Diff_EOD_2_BOTH
!
!  Routine: Diff_EOD_2_BOTH do differentiation on Even/ODD functions
!           Diff_EOD_2_NONE do differentiation on General  functions
!
! ----------------------------------------------------------------------

MODULE Differentiation_EOD_2
#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE

PRIVATE 

PUBLIC  :: Diff_EOD_2

CONTAINS

  Subroutine Diff_EOD_2      (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer  :: M, N, O_Of_D, ldy

  integer , OPTIONAL :: IEO
  logical            :: L_IEO 

  REALTYPE, dimension(ldy,*)   :: y, dy
  REALTYPE, dimension(*)       :: Aux_Diff

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call   Diff_EOD_2_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)
    else
      call   Diff_EOD_2_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
    endif
  else
      call   Diff_EOD_2_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
  endif

  END Subroutine Diff_EOD_2      


! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_EOD_2_NONE
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : Differentiates the vectors in the Matrix f
!                using the Even-Odd algorithm
! ----------------------------------------------------------------------
  Subroutine Diff_EOD_2_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)

  integer  :: N, M, O_Of_D, ldy
  integer  :: i

  REALTYPE, dimension(ldy,*)   :: y, dy
  REALTYPE, dimension(*)       :: Aux_Diff

  integer  :: Method, BlockSize, EOD_Code
  integer  :: N0, N1
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_Offset, i_d_e, i_d_o, i_m_s
  REALTYPE :: Factor

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

     Method = OPERATOR_METHOD(Aux_Diff)

  SELECT CASE (Method)
    CASE (0)
      N0 = N   ; N1 = N+1 

    CASE (1,2)
      N0 = N-1 ; N1 = N   

  END SELECT

  i_Offset = (O_Of_D-1)*EOD_HALF_SIZE(N1)
  i_d_e    = EOD_MATRIX_START_EVEN(N1) + i_Offset
  i_d_o    = EOD_MATRIX_START_ODD(N1)  + i_Offset

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call EOD_Diff_2_NONE (N0, m_Vector, Aux_Diff(i_d_e), Aux_Diff(i_d_o),  &
                                        y(m_Start,1)   , dy(m_Start,1)  )
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call EOD_Diff_2_NONE (N0, m_Vector, Aux_Diff(i_d_e), Aux_Diff(i_d_o),  &
                                        y(m_Start,1)   , dy(m_Start,1)  )
  endif

   EOD_Code = OPERATOR_EOD_CODE(Aux_Diff)

  if (EOD_Code == 1) then
    i_m_s = EOD_METRIC_START(N1)

    call Apply_Mapping_2 (N-1, M, ldy, Aux_Diff(i_m_s), dy)
  endif

  CONTAINS

    Subroutine EOD_Diff_2_NONE (N, M, D_Even, D_Odd, y, dy)

    integer  :: N, M
    integer  :: i, j, N2
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(ldy,0:N)     :: y, dy
    REALTYPE, dimension(0:N/2,0:N/2) :: D_Even, D_Odd
    REALTYPE, dimension(0:N/2,M)     :: y_Even, y_Odd, dy_Even, dy_Odd

    Factor = (-1)**(O_Of_D-1)

      N2 = N/2

    LD_A = SIZE( D_Even, DIM=1) ; Row_A = N/2+1 ; Col_A = N/2+1
    LD_B = SIZE( y_Even, DIM=1)
    LD_C = SIZE(dy_Even, DIM=1) ;                 Col_C = M

    SELECT CASE (Method)
      CASE (0)
        y_Even(0,:) = y(1:M,0) + y(1:M,0)
        y_Odd (0,:) = ZERO

      CASE (1,2)
        y_Even(0,:) = y(1:M,0) + y(1:M,N)
        y_Odd (0,:) = y(1:M,0) - y(1:M,N)

    END SELECT

    do j = 1,M
      do i = 1,N2
        y_Even(i,j) = y(j,i) + y(j,N-i)
        y_Odd (i,j) = y(j,i) - y(j,N-i)
      enddo
    enddo

    call GEMM ('n', 'n', Row_A, Col_C, Col_A, ONE, D_Even, LD_A, &
                              y_Even, LD_B, ZERO, dy_Even, LD_C)
    call GEMM ('n', 'n', Row_A, Col_C, Col_A, ONE, D_Odd , LD_A, &
                              y_Odd , LD_B, ZERO, dy_Odd , LD_C)

    SELECT CASE (Method)
      CASE (0)
        dy(1:M,0) =          dy_Even(0,:) + dy_Odd(0,:)

      CASE (1,2)
        dy(1:M,0) =          dy_Even(0,:) + dy_Odd(0,:)
        dy(1:M,N) = Factor*(-dy_Even(0,:) + dy_Odd(0,:))

    END SELECT

    do j = 1,M
      do i = 1,N2
        dy(j,  i) =          dy_Even(i,j) + dy_Odd(i,j)
        dy(j,N-i) = Factor*(-dy_Even(i,j) + dy_Odd(i,j))
      enddo
    enddo

    END Subroutine EOD_Diff_2_NONE 
! ----------------------------------------------------------------------
    Subroutine Apply_Mapping_2   (N, M, ldy, Slope, dy)

    integer  :: N, M, ldy
    integer  :: j

    REALTYPE, dimension(ldy,0:N) :: dy
    REALTYPE, dimension(0:N)     :: Slope

!$OMP PARALLEL DO
    do j = 0,N
      dy(1:M,j) = dy(1:M,j)*Slope(j)
    enddo
!$OMP END PARALLEL DO

    END Subroutine Apply_Mapping_2   

  END Subroutine Diff_EOD_2_NONE 

! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_EOD_2_BOTH
! ----------------------------------------------------------------------
  Subroutine Diff_EOD_2_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer  :: N, M, O_Of_D, ldy, IEO
  integer  :: i

  REALTYPE, dimension(ldy,*)   :: y, dy
  REALTYPE, dimension(*)       :: Aux_Diff

  integer  :: Method, BlockSize
  integer  :: N0, N1
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_Offset, i_d

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

     Method = OPERATOR_METHOD(Aux_Diff)

  SELECT CASE (Method)
    CASE (0)
      N0 = N   ; N1 = N+1 

    CASE (1,2)
      N0 = N-1 ; N1 = N   

  END SELECT

  i_Offset = (O_Of_D-1)*EOD_HALF_SIZE(N1)

  if (IEO == 1) then
    i_d  = EOD_MATRIX_START_EVEN(N1) + i_Offset
  else
    i_d  = EOD_MATRIX_START_ODD(N1)  + i_Offset
  endif

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call EOD_Diff_BOTH_2 (N0, m_Vector, Aux_Diff(i_d),   &
                                        y(m_Start,1), dy(m_Start,1))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1

    call EOD_Diff_BOTH_2 (N0, m_Vector, Aux_Diff(i_d),   &
                                        y(m_Start,1), dy(m_Start,1))
  endif

  if (IEO == 1) then
    if (MOD(O_Of_D,2) /= 0) then
      if (Method    == 0) dy(1:M,1    ) = ZERO
      if (MOD(N0,2) == 0) dy(1:M,N/2+1) = ZERO
    endif
  else
    if (MOD(O_Of_D,2) == 0) then
      if (Method    == 0) dy(1:M,1    ) = ZERO
      if (MOD(N0,2) == 0) dy(1:M,N/2+1) = ZERO
    endif
  endif

  CONTAINS

    Subroutine EOD_Diff_BOTH_2 (N, M, D, y, dy)

    integer  :: N, M
    integer  :: Row_A, Col_A, Col_C, LD_A, LD_B, LD_C

    REALTYPE, dimension(ldy,0:N/2)   :: y, dy
    REALTYPE, dimension(0:N/2,0:N/2) :: D   

    LD_A = SIZE(y , DIM=1) ; Row_A = M ; Col_A = N/2+1
    LD_B = SIZE(D , DIM=1)
    LD_C = SIZE(dy, DIM=1) ;             Col_C = N/2+1

    call GEMM ('n', 't', Row_A, Col_C, Col_A,  &
                         ONE, y, LD_A, D, LD_B, ZERO, dy, LD_C)

    END Subroutine EOD_Diff_BOTH_2

  END Subroutine Diff_EOD_2_BOTH 

END MODULE Differentiation_EOD_2
