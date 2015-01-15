#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------
! FUNCTION NAME: F_Smooth_FFT_1
! DESCRIPTION  : does the Fourier solution smoothing on m vectors
!                using the transform-recursion algorithm.
!                This is a "high-level" Subroutine that does
!                unpacking and blocking.
!
! ----------------------------------------------------------------------

MODULE Fourier_Smoothing_FFT_1
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

PUBLIC  :: F_Smooth_FFT_1

CONTAINS

  Subroutine F_Smooth_FFT_1      (N, M, y, ldy, Aux_Smth, IEO)

  integer            :: N, M, ldy

  integer , OPTIONAL :: IEO
  logical            :: L_IEO

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call F_Smooth_FFT_1_BOTH (N, M, y, ldy, Aux_Smth, IEO)
    else
      call F_Smooth_FFT_1_NONE (N, M, y, ldy, Aux_Smth)
    endif
  else
      call F_Smooth_FFT_1_NONE (N, M, y, ldy, Aux_Smth)
  endif

  END Subroutine F_Smooth_FFT_1

! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_FFT_1_NONE
! ----------------------------------------------------------------------

  Subroutine F_Smooth_FFT_1_NONE (N, M, y, ldy, Aux_Smth)

  integer            :: N, M, ldy

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, N0, N1, i
  integer  :: m_Block, m_Left, m_OffSet, m_Vector, m_Start
  integer  :: i_fft, i_s_s

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

  N0 = N ; N1 = N

  i_fft = FFT_FWAVE_START_FULL(N1)
  i_s_s = FFT_SIGMA_START_FULL(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_OffSet = 0
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call FFT_Smooth_1_NONE (N0, m_Vector, Aux_Smth(i_fft), Aux_Smth(i_s_s), &
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

    call FFT_Smooth_1_NONE (N0, m_Vector, Aux_Smth(i_fft), Aux_Smth(i_s_s), &
                                          y(1,m_Start))
                                
  endif

  CONTAINS
 
    Subroutine FFT_Smooth_1_NONE (N0, M, WSave, Sigma, y)
      
    integer  :: N0, M

    REALTYPE, dimension(*)         :: WSave
    REALTYPE, dimension(0:N0/2)    :: Sigma
    REALTYPE, dimension(0:ldy-1,M) :: y

    integer  :: N, i, j

#if defined(ESSL) 
    integer  :: ISIGN
    integer  :: N2, N3, i_fft, i_bft
    integer  :: wf_Size, wb_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(0:N0+1,M)   :: w1
    REALTYPE, dimension(1)          :: w3
#endif

#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(M,0:N0-1)   :: w1
#endif

#if defined(CRAY)
    integer  :: ldw, ISIGN, Stride
    integer , dimension(19)         :: Factor

    REALTYPE, dimension(0:N0+2,M)   :: w1
    REALTYPE, dimension(0:N0-1,M,2) :: w2
#endif

    N  = N0

#if defined(ESSL)

    N2 = N ; Scale = SQRT(ONE/N) ; N3 = N/2+1

      i_fft = 1
      i_bft = i_fft + ESSL_FFT_Twiddle_Factor (N2)

    wf_Size = ESSL_FFT_Twiddle_Factor (N2)
    wb_Size = ESSL_FFT_Twiddle_Factor (N2)
    wk_Size = ESSL_FFT_Work           (N2)

    call   Allocate_WorkSpace (wk_Size)

    ISIGN =  1 

    call RCFT (0,  y, ldy, w1,  N3, N2, M, ISIGN, Scale,  &
                      WSave(i_fft), wf_Size, Work, wk_Size, w3, 0)

!   if (Sigma(0) /= ONE) then
!     w1(0,m_OffSet+1:M) = w1(0,m_OffSet+1:M)*Sigma(0)
!   endif

    do j = m_OffSet+1,M
      do i = 1,N/2
        w1(2*i  ,j) = w1(2*i  ,j)*Sigma(i)
        w1(2*i+1,j) = w1(2*i+1,j)*Sigma(i)
      enddo
    enddo

    ISIGN = -1

    call CRFT (0, w1,  N3,  y, ldy, N2, M, ISIGN, Scale,  &
                      WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

    call DeAllocate_WorkSpace (wk_Size)

#endif

#if defined(VFFT)

    ldw = SIZE(w1, DIM=1)

    do j = 1,M
      w1(j,0:N-1) = y(0:N-1,j)
    enddo

    call VRFFTF (M, N, w1,  y, ldw, WSave)
 
!   if (Sigma(0) /= ONE) then
!     w1(:,0) = w1(:,0)*Sigma(0)
!   endif

    do i = 1,N/2-1
      w1(:,2*i-1) = w1(:,2*i-1)*Sigma(i)
      w1(:,2*i  ) = w1(:,2*i  )*Sigma(i)
    enddo

    if (MOD(N,2) == 0) then
      w1(:,N-1) = w1(:,N-1)*Sigma(N/2)
    endif

    call VRFFTB (M, N, w1,  y, ldw, WSave)
 
    do j = 1,M
      y(0:N-1,j) = w1(j,0:N-1)
    enddo

#endif

#if defined(CRAY)

    ldw = SIZE(w1, DIM=1)

    w1(0:N-1,:) = y(0:N-1,:)

    Factor = INT(WSave(1:19))

    ISIGN  = -1 ; Stride = 1

    call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M, ISIGN)

!   if (Sigma(0) /= ONE) then
!     w1(0,:) = w1(0,:)*Sigma(0)
!   endif

    do j = 1,M
      do i = 1,N/2
        w1(2*i  ,j) = w1(2*i  ,j)*Sigma(i)
        w1(2*i+1,j) = w1(2*i+1,j)*Sigma(i)
      enddo
    enddo

    ISIGN  =  1

    call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M, ISIGN)

    y(0:N-1,:) = w1(0:N-1,:)

#endif

    END Subroutine FFT_Smooth_1_NONE 

  END Subroutine F_Smooth_FFT_1_NONE 

! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_FFT_1_BOTH
! ----------------------------------------------------------------------

  Subroutine F_Smooth_FFT_1_BOTH (N, M, y, ldy, Aux_Smth, IEO)

  integer            :: N, M, ldy, IEO

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, N0, N1, i
  integer  :: m_Block, m_Left, m_OffSet, m_Vector, m_Start
  integer  :: i_cos, i_sin, i_s_s

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

  N0 = N/2 ; N1 = 2*N0+1

  i_cos = FFT_CWAVE_START_HALF(N1)
  i_sin = FFT_SWAVE_START_HALF(N1)
  i_s_s = FFT_SIGMA_START_HALF(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_OffSet = 0
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    if (IEO == 1) then
      call FFT_Smooth_1_EVEN (N0, m_Vector, Aux_Smth(i_cos), Aux_Smth(i_s_s), &
                                            y(1,m_Start))
    else
      call FFT_Smooth_1_ODD  (N0, m_Vector, Aux_Smth(i_sin), Aux_Smth(i_s_s), &
                                            y(1,m_Start))
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
      call FFT_Smooth_1_EVEN (N0, m_Vector, Aux_Smth(i_cos), Aux_Smth(i_s_s), &
                                            y(1,m_Start))
    else
      call FFT_Smooth_1_ODD  (N0, m_Vector, Aux_Smth(i_sin), Aux_Smth(i_s_s), &
                                            y(1,m_Start))
    endif

  endif

  CONTAINS

    Subroutine FFT_Smooth_1_EVEN (N0, M, WSave, Sigma, y)
      
    integer  :: N0, M

    REALTYPE, dimension(*)         :: WSave
    REALTYPE, dimension(0:N0)      :: Sigma
    REALTYPE, dimension(0:ldy-1,M) :: y

    integer  :: N, i, j

#if defined(ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(M,0:N0)   :: w1
#endif

    N  = N0

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
      w1(j,0:N) = y(0:N,j)
    enddo

    call VCosT  (M, N+1, w1,  y, ldw, WSave)
 
!   if (Sigma(0) /= ONE) then
!     w1(:,0) = w1(:,0)*Sigma(0)
!   endif

    do i = 1,N
      w1(:,i) = w1(:,i)*Sigma(i)
    enddo

    call VCosT  (M, N+1, w1,  y, ldw, WSave)
 
    do j = 1,M
      y(0:N,j) = w1(j,0:N)
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

    END Subroutine FFT_Smooth_1_EVEN 
! ----------------------------------------------------------------------
    Subroutine FFT_Smooth_1_ODD  (N0, M, WSave, Sigma, y)
      
    integer  :: N0, M

    REALTYPE, dimension(*)         :: WSave
    REALTYPE, dimension(0:N0)      :: Sigma
    REALTYPE, dimension(0:ldy-1,M) :: y

    integer  :: N, i, j

#if defined(ESSL)
    integer  :: N2
    integer  :: ws_Size, wk_Size
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(M,0:N0)   :: w1
#endif

    N  = N0

#if defined(ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)
   
    ws_Size = ESSL_SFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call   Allocate_WorkSpace (wk_Size)

    call SINF (0,  y,1,ldy,  y,1,ldy, N2, M, Scale, WSave, ws_Size,  &
                                                    Work , wk_Size)

    do j = m_OffSet+1,M
      y(1:N-1,j) = y(1:N-1,j)*Sigma(1:N-1)
    enddo

    call SINF (0,  y,1,ldy,  y,1,ldy, N2, M, Scale, WSave, ws_Size,  &
                                                    Work , wk_Size)

    call DeAllocate_WorkSpace (wk_Size)

#endif

#if defined(VFFT)

    ldw = SIZE(w1, DIM=1)

    do j = 1,M
      w1(j,0:N-1) = y(0:N-1,j)
    enddo

    call VSinT  (M, N-1, w1(1,1),  y, ldw, WSave)

    do i = 1,N-1
      w1(:,i) = w1(:,i)*Sigma(i)
    enddo

    call VSinT  (M, N-1, w1(1,1),  y, ldw, WSave)
 
    do j = 1,M
      y(1:N-1,j) = w1(j,1:N-1)
    enddo
 
    y(0,:) = ZERO ; y(N,:) = ZERO

#endif

#if defined(CRAY)

    call Cray_SFT_1 (N, M,  y,ldy,  y,ldy, WSave)

    do j = 1,M
      y(1:N-1,j) = y(1:N-1,j)*Sigma(1:N-1)
    enddo

    call Cray_SFT_1 (N, M,  y,ldy,  y,ldy, WSave)

#endif

    END Subroutine FFT_Smooth_1_ODD  

  END Subroutine F_Smooth_FFT_1_BOTH 

END MODULE Fourier_Smoothing_FFT_1



! ----------------------------------------------------------------------
! FUNCTION NAME: F_Smooth_FFT_2
! DESCRIPTION  : does the Fourier solution smoothing on m vectors
!                using the transform-recursion algorithm.
!                This is a "high-level" Subroutine that does
!                unpacking and blocking.
!
! ----------------------------------------------------------------------

MODULE Fourier_Smoothing_FFT_2
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

PUBLIC  :: F_Smooth_FFT_2

CONTAINS

  Subroutine F_Smooth_FFT_2      (N, M, y, ldy, Aux_Smth, IEO)

  integer            :: N, M, ldy

  integer , OPTIONAL :: IEO
  logical            :: L_IEO 

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call F_Smooth_FFT_2_BOTH (N, M, y, ldy, Aux_Smth, IEO)
    else
      call F_Smooth_FFT_2_NONE (N, M, y, ldy, Aux_Smth)
    endif
  else
      call F_Smooth_FFT_2_NONE (N, M, y, ldy, Aux_Smth)
  endif

  END Subroutine F_Smooth_FFT_2

! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_FFT_2_NONE
! ----------------------------------------------------------------------

  Subroutine F_Smooth_FFT_2_NONE (N, M, y, ldy, Aux_Smth)

  integer            :: N, M, ldy

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, N0, N1, i
  integer  :: m_Block, m_Left, m_OffSet, m_Vector, m_Start
  integer  :: i_fft, i_s_s

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

  N0 = N ; N1 = N

  i_fft = FFT_FWAVE_START_FULL(N1)
  i_s_s = FFT_SIGMA_START_FULL(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_OffSet = 0
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call FFT_Smooth_2_NONE (N0, m_Vector, Aux_Smth(i_fft), Aux_Smth(i_s_s), &
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

    call FFT_Smooth_2_NONE (N0, m_Vector, Aux_Smth(i_fft), Aux_Smth(i_s_s), &
                                          y(m_Start,1))

  endif

  CONTAINS

    Subroutine FFT_Smooth_2_NONE (N0, M, WSave, Sigma, y)

    integer  :: N0, M

    REALTYPE, dimension(*)          :: WSave
    REALTYPE, dimension(0:N0/2)     :: Sigma
    REALTYPE, dimension(ldy,0:N0-1) :: y

    integer  :: N, i, j

#if defined(ESSL)
    integer  :: ldw, ISIGN
    integer  :: N2, N3, i_fft, i_bft
    integer  :: wf_Size, wb_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(0:N0+1,M)   :: w1
    REALTYPE, dimension(1)          :: w3
#endif

#if defined(VFFT)
    REALTYPE, dimension(ldy,0:N0-1) :: w1
#endif

#if defined(CRAY)
    integer  :: ldw, ISIGN, Stride
    integer , dimension(19)         :: Factor

    REALTYPE, dimension(0:N0+2,M)   :: w1
    REALTYPE, dimension(0:N0-1,M,2) :: w2
#endif

    N  = N0

#if defined(ESSL)

    ldw = SIZE(w1, DIM=1)

    do j = 1,M
      w1(0:N-1,j) = y(j,0:N-1)
    enddo

    N2 = N ; Scale = SQRT(ONE/N) ; N3 = N/2+1

      i_fft = 1
      i_bft = i_fft + ESSL_FFT_Twiddle_Factor (N2)

    wf_Size = ESSL_FFT_Twiddle_Factor (N2)
    wb_Size = ESSL_FFT_Twiddle_Factor (N2)
    wk_Size = ESSL_FFT_Work           (N2)

    call   Allocate_WorkSpace (wk_Size)

    ISIGN = 1 

    call RCFT (0, w1, ldw, w1,  N3, N2, M, ISIGN, Scale,  &
                      WSave(i_fft), wf_Size, Work, wk_Size, w3, 0)

!   if (Sigma(0) /= ONE) then
!     w1(0,m_OffSet+1:M) = w1(0,m_OffSet+1:M)*Sigma(0)
!   endif

    do j = m_OffSet+1,M
      do i = 1,N/2
        w1(2*i  ,j) = w1(2*i  ,j)*Sigma(i)
        w1(2*i+1,j) = w1(2*i+1,j)*Sigma(i)
      enddo
    enddo

    ISIGN = -1

    call CRFT (0, w1,  N3, w1, ldw, N2, M, ISIGN, Scale,  &
                      WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

    do j = 1,M
      y(j,0:N-1) = w1(0:N-1,j)
    enddo
          
    call DeAllocate_WorkSpace (wk_Size)

#endif

#if defined(VFFT)

    call VRFFTF (M, N,  y, w1, ldy, WSave)

!   if (Sigma(0) /= ONE) then
!     y(1:M,0) = y(1:M,0)*Sigma(0)
!   endif

    do i = 1,N/2-1
      y(1:M,2*i-1) = y(1:M,2*i-1)*Sigma(i)
      y(1:M,2*i  ) = y(1:M,2*i  )*Sigma(i)
    enddo

    if (MOD(N,2) == 0) then
      y(1:M,N-1) = y(1:M,N-1)*Sigma(N/2)
    endif

    call VRFFTB (M, N,  y, w1, ldy, WSave)
 
#endif

#if defined(CRAY)
 
    ldw = SIZE(w1, DIM=1)

    do j = 1,M
      w1(0:N-1,j) = y(j,0:N-1)
    enddo

    Factor = INT(WSave(1:19))

    ISIGN  = -1 ; Stride = 1

    call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M, ISIGN)

!   if (Sigma(0) /= ONE) then
!     w1(0,1:M) = w1(0,1:M)*Sigma(0)
!   endif

    do j = 1,M
      do i = 1,N/2
        w1(2*i  ,j) = w1(2*i  ,j)*Sigma(i)
        w1(2*i+1,j) = w1(2*i+1,j)*Sigma(i)
      enddo
    enddo

    ISIGN  =  1 

    call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  ISIGN)

    do j = 1,M
      y(j,0:N-1) = w1(0:N-1,j)
    enddo

#endif

    END Subroutine FFT_Smooth_2_NONE 

  END Subroutine F_Smooth_FFT_2_NONE

! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth_FFT_2_BOTH
! ----------------------------------------------------------------------

  Subroutine F_Smooth_FFT_2_BOTH (N, M, y, ldy, Aux_Smth, IEO) 

  integer            :: N, M, ldy, IEO

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: BlockSize, N0, N1, i
  integer  :: m_Block, m_Left, m_OffSet, m_Vector, m_Start
  integer  :: i_cos, i_sin, i_s_s

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

  N0 = N/2 ; N1 = 2*N0+1

  i_cos = FFT_CWAVE_START_HALF(N1)
  i_sin = FFT_SWAVE_START_HALF(N1)
  i_s_s = FFT_SIGMA_START_HALF(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_OffSet = 0
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    if (IEO == 1) then
      call FFT_Smooth_2_EVEN (N0, m_Vector, Aux_Smth(i_cos), Aux_Smth(i_s_s), &
                                            y(m_Start,1))
    else
      call FFT_Smooth_2_ODD  (N0, m_Vector, Aux_Smth(i_sin), Aux_Smth(i_s_s), &
                                            y(m_Start,1))
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
      call FFT_Smooth_2_EVEN (N0, m_Vector, Aux_Smth(i_cos), Aux_Smth(i_s_s), &
                                            y(m_Start,1))
    else
      call FFT_Smooth_2_ODD  (N0, m_Vector, Aux_Smth(i_sin), Aux_Smth(i_s_s), &
                                            y(m_Start,1))
    endif

  endif

  CONTAINS

    Subroutine FFT_Smooth_2_EVEN (N0, M, WSave, Sigma, y)

    integer  :: N0, M

    REALTYPE, dimension(*)          :: WSave
    REALTYPE, dimension(0:N0)       :: Sigma
    REALTYPE, dimension(ldy,0:N0)   :: y

    integer  :: N, i, j

#if defined(ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    REALTYPE, dimension(ldy,0:N0)   :: w1
#endif

    N  = N0

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

    call COSF (0,  y,ldy,1,  y,ldy,1, N2, M, Scale, WSave, wc_Size,  &
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

    END Subroutine FFT_Smooth_2_EVEN 
! ----------------------------------------------------------------------
    Subroutine FFT_Smooth_2_ODD  (N0, M, WSave, Sigma, y)

    integer  :: N0, M

    REALTYPE, dimension(*)          :: WSave
    REALTYPE, dimension(0:N0)       :: Sigma
    REALTYPE, dimension(ldy,0:N0)   :: y

    integer  :: N, i, j

#if defined(ESSL)
    integer  :: N2
    integer  :: ws_Size, wk_Size
    REALTYPE :: Scale
#endif

#if defined(VFFT)
    REALTYPE, dimension(ldy,0:N0)   :: w1
#endif

    N  = N0

#if defined(ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    ws_Size = ESSL_SFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 2)

    call   Allocate_WorkSpace (wk_Size)

    call SINF (0,  y,ldy,1,  y,ldy,1, N2, M, Scale, WSave, ws_Size,  &
                                                    Work , wk_Size)

    do i = 1,N-1
      y(m_OffSet+1:M,i) = y(m_OffSet+1:M,i)*Sigma(i)
    enddo

    call SINF (0,  y,ldy,1,  y,ldy,1, N2, M, Scale, WSave, ws_Size,  &
                                                    Work , wk_Size)

    call DeAllocate_WorkSpace (wk_Size)

#endif

#if defined(VFFT)

    call VSinT  (M, N-1,  y(1,1), w1, ldy, WSave)

    do i = 1,N-1
      y(1:M,i) = y(1:M,i)*Sigma(i)
    enddo

    call VSinT  (M, N-1,  y(1,1), w1, ldy, WSave)

    y(1:M,N) = ZERO

#endif

#if defined(CRAY)

    call Cray_SFT_2 (N, M,  y,ldy,  y,ldy, WSave)

    do i = 1,N-1
      y(1:M,i) = y(1:M,i)*Sigma(i)
    enddo

    call Cray_SFT_2 (N, M,  y,ldy,  y,ldy, WSave)

#endif

    END Subroutine FFT_Smooth_2_ODD

  END Subroutine F_Smooth_FFT_2_BOTH

END MODULE Fourier_Smoothing_FFT_2
