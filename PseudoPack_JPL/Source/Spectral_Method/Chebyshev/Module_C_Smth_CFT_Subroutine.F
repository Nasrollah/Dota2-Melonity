#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_CFT_1
! DESCRIPTION  : does the Chebyshev solution smoothing on m vectors
!                using the transform-recursion algorithm.
!                This is a "high-level" Subroutine that does
!                unpacking and blocking.
!
! ----------------------------------------------------------------------

MODULE Chebyshev_Smoothing_CFT_1
  USE PS_WorkSpace_Allocation, Work=>WorkSpace

#if defined (ESSL)
  USE ESSL_FFT
#endif
#if defined (VFFT)
  USE VFFT_FFT
#endif
#if defined (CRAY)
  USE CRAY_FFT
#endif

implicit NONE

PRIVATE

PUBLIC  :: C_Smooth_CFT_1

CONTAINS

  Subroutine C_Smooth_CFT_1      (N, M, y, ldy, Aux_Smth, IEO)

  integer            :: N, M, ldy

  integer , OPTIONAL :: IEO
  logical            :: L_IEO 

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call C_Smooth_CFT_1_BOTH (N, M, y, ldy, Aux_Smth, IEO)
    else
      call C_Smooth_CFT_1_NONE (N, M, y, ldy, Aux_Smth)
    endif
  else
      call C_Smooth_CFT_1_NONE (N, M, y, ldy, Aux_Smth)
  endif

  END Subroutine C_Smooth_CFT_1

! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_CFT_1_NONE
! ----------------------------------------------------------------------
  Subroutine C_Smooth_CFT_1_NONE (N, M, y, ldy, Aux_Smth)

  integer            :: N, M, ldy

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, N1, N_Pnts, i
  integer  :: m_Block, m_Left, m_OffSet, m_Vector, m_Start
  integer  :: i_w_s, i_s_s

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

  N1 = N ; N_Pnts = N

  i_w_s = CFT_WSAVE_START_FULL(N_Pnts)
  i_s_s = CFT_SIGMA_START_FULL(N_Pnts)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_OffSet = 0
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call CFT_Smooth_1_NONE (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_s_s), &
                                          y(1,m_Start))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined(ESSL)
    m_OffSet =   BlockSize-m_Left
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined(VFFT) || defined(CRAY)
           i =   m_Block
    m_OffSet = 0
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    call CFT_Smooth_1_NONE (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_s_s), &
                                          y(1,m_Start))

  endif

  CONTAINS

    Subroutine CFT_Smooth_1_NONE (N1, M, WSave, Sigma, y)
      
    integer  :: N1, M

    REALTYPE, dimension(*)         :: WSave
    REALTYPE, dimension(0:N1-1)    :: Sigma
    REALTYPE, dimension(0:ldy-1,M) :: y

    integer  :: N, i, j

#if defined(ESSL) 
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(M,0:N1-1)   :: w1
#endif

    N = N1-1 

#if defined(ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,1,ldy,  y,1,ldy, N2, M, Scale, WSave, wc_Size,  &
                                                    Work , wk_Size)
                  
!   if (Sigma(0) /= ONE) then
!     y(0,m_OffSet+1:M) = y(0,m_OffSet+1:M)*Sigma(0)
!   endif

    do j = m_OffSet+1,M
      y(1:N,j) = y(1:N,j)*Sigma(1:N)
    enddo

    call COSF (0,  y,1,ldy,  y,1,ldy, N2, M, Scale, WSave, wc_Size,  &
                                                    Work , wk_Size)
 
    call DeAllocate_WorkSpace (wk_Size)

#endif

#if defined(VFFT)

    ldw = SIZE(w1, DIM=1)

    do j = 1,M
      do i = 0,N
        w1(j,i) = y(i,j)
      enddo
    enddo

    call VCosT  (M, N+1, w1,  y, ldw, WSave)
 
!   if (Sigma(0) /= ONE) then
!     w1(:,0) = w1(:,0)*Sigma(0)
!   endif

    do i = 1,N
      w1(:,i) = w1(:,i)*Sigma(i)
    enddo

    call VCosT  (M, N+1, w1,  y, ldw, WSave)
 
    do i = 0,N
      do j = 1,M
        y(i,j) = w1(j,i)
      enddo
    enddo

#endif

#if defined(CRAY)

    call Cray_CFT_1 (N, M,  y,ldy,  y,ldy, WSave)

!   if (Sigma(0) /= ONE) then
!     y(0,:) = y(0,:)*Sigma(0)
!   endif

    do j = 1,M
      y(1:N,j) = y(1:N,j)*Sigma(1:N)
    enddo

    call Cray_CFT_1 (N, M,  y,ldy,  y,ldy, WSave)

#endif

    END Subroutine CFT_Smooth_1_NONE 

  END Subroutine C_Smooth_CFT_1_NONE 

! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_CFT_1_BOTH
! ----------------------------------------------------------------------
  Subroutine C_Smooth_CFT_1_BOTH (N, M, y, ldy, Aux_Smth, IEO) 

  integer            :: N, M, ldy, IEO

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, N1, N_Pnts, i
  integer  :: m_Block, m_Left, m_OffSet, m_Vector, m_Start
  integer  :: i_w_s, i_q_s, i_s_s

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

  N1 = (N+1)/2 ; N_Pnts = N

  i_w_s = CFT_WSAVE_START_HALF(N_Pnts)
  i_q_s = CFT_QWAVE_START_HALF(N_Pnts)
  i_s_s = CFT_SIGMA_START_HALF(N_Pnts)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_OffSet = 0
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    if (IEO == 1) then
      call CFT_Smooth_1_EVEN (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_q_s),&
                                            Aux_Smth(i_s_s), y(1,m_Start)  )
    else
      call CFT_Smooth_1_ODD  (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_q_s),&
                                            Aux_Smth(i_s_s), y(1,m_Start)  )
    endif

  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined(ESSL)
    m_OffSet =   BlockSize-m_Left
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined(VFFT) || defined(CRAY)
           i =   m_Block
    m_OffSet = 0
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    if (IEO == 1) then
      call CFT_Smooth_1_EVEN (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_q_s),&
                                            Aux_Smth(i_s_s), y(1,m_Start)  )
    else
      call CFT_Smooth_1_ODD  (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_q_s),&
                                            Aux_Smth(i_s_s), y(1,m_Start)  )
    endif

  endif

  CONTAINS

    Subroutine CFT_Smooth_1_EVEN (N1, M, WSave, QWave, Sigma, y) 

    integer  :: N1, M

    REALTYPE, dimension(*)          :: WSave, QWave
    REALTYPE, dimension(0:2*(N1-1)) :: Sigma
    REALTYPE, dimension(0:ldy-1,M)  :: y

    integer  :: N, i, j

#if defined(ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(M,0:N1-1)   :: w1
#endif

    N = N1-1 

#if defined(ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,1,ldy,  y,1,ldy, N2, M, Scale, WSave, wc_Size,  &
                                                    Work , wk_Size)
 
!   if (Sigma(0) /= ONE) then
!     y(0,m_OffSet+1:M) = y(0,m_OffSet+1:M)*Sigma(0)
!   endif

    do j = m_OffSet+1,M
      do i = 1,N
        y(i,j) = y(i,j)*Sigma(2*i)
      enddo
    enddo

    call COSF (0,  y,1,ldy,  y,1,ldy, N2, M, Scale, WSave, wc_Size,  &
                                                    Work , wk_Size)
 
    call DeAllocate_WorkSpace (wk_Size)

#endif

#if defined(VFFT)

    ldw = SIZE(w1, DIM=1)

    do j = 1,M
      do i = 0,N
        w1(j,i) = y(i,j)
      enddo
    enddo

    call VCosT  (M, N+1, w1,  y, ldw, WSave)
 
!   if (Sigma(0) /= ONE) then
!     w1(:,0) = w1(:,0)*Sigma(0)
!   endif

    do i = 1,N
      w1(:,i) = w1(:,i)*Sigma(2*i)
    enddo

    call VCosT  (M, N+1, w1,  y, ldw, WSave)
 
    do i = 0,N
      do j = 1,M
        y(i,j) = w1(j,i)
      enddo
    enddo

#endif

#if defined(CRAY)

    call Cray_CFT_1 (N, M,  y,ldy,  y,ldy, WSave)

!   if (Sigma(0) /= ONE) then
!     y(0,1:M) = y(0,1:M)*Sigma(0)
!   endif

    do j = 1,M
      do i = 1,N
        y(i,j) = y(i,j)*Sigma(2*i)
      enddo
    enddo

    call Cray_CFT_1 (N, M,  y,ldy,  y,ldy, WSave)

#endif

    END Subroutine CFT_Smooth_1_EVEN 
! ----------------------------------------------------------------------
    Subroutine CFT_Smooth_1_ODD  (N1, M, WSave, QWave, Sigma, y) 
      
    integer  :: N1, M

    REALTYPE, dimension(*)          :: WSave, QWave
    REALTYPE, dimension(0:2*(N1-1)) :: Sigma
    REALTYPE, dimension(0:ldy-1,M)  :: y

    integer  :: N, i, j

#if defined(ESSL)
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(M,0:N1-1)   :: w1
#endif

    N = N1-1 

#if defined(ESSL)

    call ESSL_QWAVE_FWRD_1 (1, N, M, y, ldy, y, ldy, QWave)

    do j = m_OffSet+1,M
      do i = 0,N-1
        y(i,j) = y(i,j)*Sigma(2*i+1)
      enddo
    enddo

    call ESSL_QWAVE_BWRD_1 (1, N, M, y, ldy, y, ldy, QWave)

#endif

#if defined(VFFT)

    ldw = SIZE(w1, DIM=1)

    do j = 1,M
      do i = 0,N
        w1(j,i) = y(i,j)
      enddo
    enddo

    call VCosQF (M, N  , w1, y, ldw, QWave)

    do i = 0,N-1
      w1(:,i) = w1(:,i)*Sigma(2*i+1)
    enddo

    call VCosQB (M, N  , w1, y, ldw, QWave)
 
    do i = 0,N
      do j = 1,M
        y(i,j) = w1(j,i)
      enddo
    enddo

#endif

#if defined(CRAY)

    call Cray_QWAVE_FWRD_1 (N, M,  y,ldy,  y,ldy, QWave, WSave)

    do j = 1,M
      do i = 0,N-1
        y(i,j) = y(i,j)*Sigma(2*i+1)
      enddo
    enddo

    call Cray_QWAVE_BWRD_1 (N, M,  y,ldy,  y,ldy, QWave, WSave)

#endif

    END Subroutine CFT_Smooth_1_ODD  

  END Subroutine C_Smooth_CFT_1_BOTH 

END MODULE Chebyshev_Smoothing_CFT_1


! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_CFT_2
! DESCRIPTION  : does the Chebyshev solution smoothing on m vectors
!                using the transform-recursion algorithm.
!                This is a "high-level" Subroutine that does
!                unpacking and blocking.
!
! ----------------------------------------------------------------------

MODULE Chebyshev_Smoothing_CFT_2
  USE PS_WorkSpace_Allocation, Work=>WorkSpace

#if defined (ESSL)
  USE ESSL_FFT
#endif
#if defined (VFFT)
  USE VFFT_FFT
#endif
#if defined (CRAY)
  USE CRAY_FFT
#endif

implicit NONE

PRIVATE

PUBLIC  :: C_Smooth_CFT_2

CONTAINS

  Subroutine C_Smooth_CFT_2      (N, M, y, ldy, Aux_Smth, IEO)

  integer            :: N, M, ldy

  integer , OPTIONAL :: IEO
  logical            :: L_IEO 

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call C_Smooth_CFT_2_BOTH (N, M, y, ldy, Aux_Smth, IEO)
    else
      call C_Smooth_CFT_2_NONE (N, M, y, ldy, Aux_Smth)
    endif
  else
      call C_Smooth_CFT_2_NONE (N, M, y, ldy, Aux_Smth)
  endif

  END Subroutine C_Smooth_CFT_2

! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_CFT_2_NONE
! ----------------------------------------------------------------------
  Subroutine C_Smooth_CFT_2_NONE (N, M, y, ldy, Aux_Smth)

  integer            :: N, M, ldy

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, N1, N_Pnts, i
  integer  :: m_Block, m_Left, m_OffSet, m_Vector, m_Start
  integer  :: i_w_s, i_s_s

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

  N1 = N ; N_Pnts = N

  i_w_s = CFT_WSAVE_START_FULL(N_Pnts)
  i_s_s = CFT_SIGMA_START_FULL(N_Pnts)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_OffSet = 0
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call CFT_Smooth_2_NONE (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_s_s), &
                                          y(m_Start,1))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined(ESSL)
    m_OffSet =   BlockSize-m_Left
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined(VFFT) || defined(CRAY)
           i =   m_Block
    m_OffSet = 0
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    call CFT_Smooth_2_NONE (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_s_s), &
                                          y(m_Start,1))

  endif

  CONTAINS

    Subroutine CFT_Smooth_2_NONE (N1, M, WSave, Sigma, y)

    integer  :: N1, M

    REALTYPE, dimension(*)          :: WSave
    REALTYPE, dimension(0:N1-1)     :: Sigma
    REALTYPE, dimension(ldy,0:N1-1) :: y

    integer  :: N, i

#if defined(ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    REALTYPE, dimension(ldy,0:N1-1) :: w1
#endif

    N = N1-1 

#if defined(ESSL)
 
    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 2)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,ldy,1,  y,ldy,1, N2, M, Scale, WSave, wc_Size,  &
                                                    Work , wk_Size)
 
!   if (Sigma(0) /= ONE) then
!     y(m_OffSet+1:M,0) = y(m_OffSet+1:M,0)*Sigma(0)
!   endif

    do i = 1,N
      y(m_OffSet+1:M,i) = y(m_OffSet+1:M,i)*Sigma(i)
    enddo

    call COSF (0,  y,ldy,1,  y,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                    Work , wk_Size)
 
    call DeAllocate_WorkSpace (wk_Size)

#endif

#if defined(VFFT)

    call VCosT  (M, N+1,  y, w1, ldy, WSave)
 
!   if (Sigma(0) /= ONE) then
!     y(1:M,0) = y(1:M,0)*Sigma(0)
!   endif

    do i = 1,N
      y(1:M,i) = y(1:M,i)*Sigma(i)
    enddo

    call VCosT  (M, N+1,  y, w1, ldy, WSave)
 
#endif

#if defined(CRAY)
 
    call Cray_CFT_2 (N, M,  y,ldy,  y,ldy, WSave)
 
!   if (Sigma(0) /= ONE) then
!     y(1:M,0) = y(1:M,0)*Sigma(0)
!   endif

    do i = 1,N
      y(1:M,i) = y(1:M,i)*Sigma(i)
    enddo

    call Cray_CFT_2 (N, M,  y,ldy,  y,ldy, WSave)

#endif
 
    END Subroutine CFT_Smooth_2_NONE 

  ENd Subroutine C_Smooth_CFT_2_NONE 

! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_CFT_2_BOTH
! ----------------------------------------------------------------------
  Subroutine C_Smooth_CFT_2_BOTH (N, M, y, ldy, Aux_Smth, IEO)

  integer            :: N, M, ldy, IEO

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, N1, N_Pnts, i
  integer  :: m_Block, m_Left, m_OffSet, m_Vector, m_Start
  integer  :: i_w_s, i_q_s, i_s_s

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

  N1 = (N+1)/2 ; N_Pnts = N

  i_w_s = CFT_WSAVE_START_HALF(N_Pnts)
  i_q_s = CFT_QWAVE_START_HALF(N_Pnts)
  i_s_s = CFT_SIGMA_START_HALF(N_Pnts)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_OffSet = 0
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    if (IEO == 1) then
      call CFT_Smooth_2_EVEN (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_q_s),&
                                            Aux_Smth(i_s_s), y(m_Start,1))
     else
      call CFT_Smooth_2_ODD  (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_q_s),&
                                            Aux_Smth(i_s_s), y(m_Start,1))
     endif

  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined(ESSL)
    m_OffSet =   BlockSize-m_Left
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined(VFFT) || defined(CRAY)
           i =   m_Block
    m_OffSet = 0
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    if (IEO == 1) then
      call CFT_Smooth_2_EVEN (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_q_s),&
                                            Aux_Smth(i_s_s), y(m_Start,1))
     else
      call CFT_Smooth_2_ODD  (N1, m_Vector, Aux_Smth(i_w_s), Aux_Smth(i_q_s),&
                                            Aux_Smth(i_s_s), y(m_Start,1))
     endif

  endif

  CONTAINS

    Subroutine CFT_Smooth_2_EVEN (N1, M, WSave, QWave, Sigma, y)

    integer  :: N1, M

    REALTYPE, dimension(*)          :: WSave, QWave
    REALTYPE, dimension(0:2*(N1-1)) :: Sigma
    REALTYPE, dimension(ldy,0:N1-1) :: y

    integer  :: N, i, j

#if defined(ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    REALTYPE, dimension(ldy,0:N1-1) :: w1
#endif

    N = N1-1 

#if defined(ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 2)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,ldy,1,  y,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                    Work , wk_Size)

!   if (Sigma(0) /= ONE) then
!     y(m_OffSet+1:M,0) = y(m_OffSet+1:M,0)*Sigma(0)
!   endif

    do i = 1,N
      y(m_OffSet+1:M,i) = y(m_OffSet+1:M,i)*Sigma(2*i)
    enddo

    call COSF (0,  y,ldy,1,  y,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                    Work , wk_Size)

    call DeAllocate_WorkSpace (wk_Size)

#endif

#if defined(VFFT)

    call VCosT  (M, N+1,  y, w1, ldy, WSave)

!   if (Sigma(0) /= ONE) then
!     y(1:M,0) = y(1:M,0)*Sigma(0)
!   endif

    do i = 1,N
      y(1:M,i) = y(1:M,i)*Sigma(2*i)
    enddo

    call VCosT  (M, N+1,  y, w1, ldy, WSave)

#endif

#if defined(CRAY)

    call Cray_CFT_2 (N, M,  y,ldy,  y,ldy, WSave)

!   if (Sigma(0) /= ONE) then
!     y(1:M,0) = y(1:M,0)*Sigma(0)
!   endif

    do i = 1,N
      y(1:M,i) = y(1:M,i)*Sigma(2*i)
    enddo

    call Cray_CFT_2 (N, M,  y,ldy,  y,ldy, WSave)

#endif

    END Subroutine CFT_Smooth_2_EVEN 
! ----------------------------------------------------------------------
    Subroutine CFT_Smooth_2_ODD  (N1, M, WSave, QWave, Sigma, y)

    integer  :: N1, M

    REALTYPE, dimension(*)          :: WSave, QWave
    REALTYPE, dimension(0:2*(N1-1)) :: Sigma
    REALTYPE, dimension(ldy,0:N1-1) :: y

    integer  :: N, i, j

#if defined(ESSL)
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    REALTYPE, dimension(ldy,0:N1-1) :: w1
#endif

    N = N1-1 

#if defined(ESSL)

    call ESSL_QWAVE_FWRD_2 (1, N, M, y, ldy, y, ldy, QWave)

    do i = 0,N-1
      y(m_OffSet+1:M,i) = y(m_OffSet+1:M,i)*Sigma(2*i+1)
    enddo

    call ESSL_QWAVE_BWRD_2 (1, N, M, y, ldy, y, ldy, QWave)

#endif

#if defined(VFFT)

    call VCosQF (M, N  , y, w1, ldy, QWave)

    do i = 0,N-1
      y(1:M,i) = y(1:M,i)*Sigma(2*i+1)
    enddo

    call VCosQB (M, N  , y, w1, ldy, QWave)

#endif

#if defined(CRAY)

    call Cray_QWAVE_FWRD_2 (N, M,  y,ldy,  y,ldy, QWave, WSave)

    do i = 0,N-1
      y(1:M,i) = y(1:M,i)*Sigma(2*i+1)
    enddo

    call Cray_QWAVE_BWRD_2 (N, M,  y,ldy,  y,ldy, QWave, WSave)

#endif

    END Subroutine CFT_Smooth_2_ODD  

  END Subroutine C_Smooth_CFT_2_BOTH 

END MODULE Chebyshev_Smoothing_CFT_2
