#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: F_Diff_FFT_1
! DESCRIPTION  : does the Fourier solution smoothing on m vectors
!                using the transform-recursion algorithm.
!                This is a "high-level" Subroutine that does
!                unpacking and blocking.
!
! ----------------------------------------------------------------------

MODULE Fourier_Differentiation_FFT_1
  USE PS_WorkSpace_Allocation, Work=>WorkSpace
  USE FFT_Auxillary_Subroutine

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

PUBLIC  :: F_Diff_FFT_1

CONTAINS

  Subroutine F_Diff_FFT_1      (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer            :: N, M, ldy, O_Of_D

  integer , OPTIONAL :: IEO
  logical            :: L_IEO 

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call F_Diff_FFT_1_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)
    else
      call F_Diff_FFT_1_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
    endif
  else
      call F_Diff_FFT_1_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
  endif

  END Subroutine F_Diff_FFT_1

! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_FFT_1_NONE
! ----------------------------------------------------------------------
  Subroutine F_Diff_FFT_1_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)

  integer            :: N, M, ldy, O_Of_D

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: Map, Smooth, BlockSize
  integer  :: N0, N1, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_fft, i_s_s, i_m_s

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

        Map = OPERATOR_MAP(Aux_Diff)
     Smooth = OPERATOR_SMOOTH_1(Aux_Diff)

  N0 = N ; N1 = N

  i_fft = FFT_FWAVE_START_FULL(N1)
  i_s_s = FFT_SIGMA_START_FULL(N1)
  i_m_s = FFT_METRIC_START_FULL(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call FFT_Diff_1_NONE (N0, m_Vector, Aux_Diff(i_fft), Aux_Diff(i_m_s),  &
                                        Aux_Diff(i_s_s),                   &
                                        y(1,m_Start)   , dy(1,m_Start))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined(ESSL)
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined(VFFT) || defined(CRAY)
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    call FFT_Diff_1_NONE (N0, m_Vector, Aux_Diff(i_fft), Aux_Diff(i_m_s),  &
                                        Aux_Diff(i_s_s),                   &
                                        y(1,m_Start)   , dy(1,m_Start))

  endif

  CONTAINS

    Subroutine FFT_Diff_1_NONE (N0, M, WSave, Metric, Sigma, y, dy)
      
    integer  :: N0, M

    REALTYPE, dimension(*)          :: WSave
    REALTYPE, dimension(0:N0-1,10)  :: Metric 
    REALTYPE, dimension(0:N0/2)     :: Sigma
    REALTYPE, dimension(0:ldy-1,M)  :: y, dy

    integer  :: N, i, j, k, Even_Or_Odd
    REALTYPE :: dy_Real, dy_Imag

    REALTYPE, dimension(0:N0/2,4)   :: b

#if defined(ESSL) 
    integer  :: N2, N3, ISIGN
    integer  :: i_fft, i_bft
    integer  :: wf_Size, wb_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(1:N0-1,M)   ::  a
    REALTYPE, dimension(0:N0-1,M)   :: dd
    REALTYPE, dimension(0:N0+1,M)   :: w1
    REALTYPE, dimension(1)          :: w3
#endif
#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(M,1:N0-1)   ::  a 
    REALTYPE, dimension(0:N0-1,M)   :: dd 
    REALTYPE, dimension(M,0:N0-1)   :: w1 
#endif
#if defined(CRAY)
    integer  :: ldw, Stride
    integer, dimension(19)          :: Factor

    REALTYPE, dimension(0:N0+1,M)   ::  a 
    REALTYPE, dimension(0:N0-1,M)   :: dd 
    REALTYPE, dimension(0:N0+2,M)   :: w1
    REALTYPE, dimension(0:N0-1,M,2) :: w2
#endif

#if defined(DEBUG)
    ASSERT(ldy>=N0,"F_Diff_1_None:ldy<N")
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N0 ; Even_Or_Odd =  0

    call FFT_D_Coefficients (Map, Smooth, Even_Or_Odd, O_Of_D, N/2, Sigma, b)
      
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

    ISIGN = -1

    w1(0,:) = ZERO ; w1(N,:) = ZERO

! -------------------------------------------------------------------
!    This case is for NO MAPPING
! -------------------------------------------------------------------
    if (Map == 0) then

      if (MOD(O_Of_D,2) == 1) then
        do j = 1,M
          do i = 1,N/2-1
                dy_Real = -b(i,O_Of_D)*w1(2*i+1,j)
                dy_Imag =  b(i,O_Of_D)*w1(2*i  ,j)

            w1(2*i  ,j) = dy_Real
            w1(2*i+1,j) = dy_Imag
          enddo
        enddo
      else
        do j = 1,M
          do i = 1,N/2-1
            w1(2*i  ,j) =  b(i,O_Of_D)*w1(2*i  ,j)
            w1(2*i+1,j) =  b(i,O_Of_D)*w1(2*i+1,j)
          enddo
        enddo
      endif

      call CRFT (0, w1,  N3, dy, ldy, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

    endif

! --------------------------------------------------------------------
!    Here starts the MAPPING CASE
! --------------------------------------------------------------------
    if (Map /= 0) then

      a(1:N-1,:) = w1(1:N-1,:)

! FIRST DERIVATIVE

      do j = 1,M
        do i = 1,N/2-1
          w1(2*i  ,j) = -b(i,1)*a(2*i+1,j)
          w1(2*i+1,j) =  b(i,1)*a(2*i  ,j)
        enddo
      enddo

      call CRFT (0, w1,  N3, dy, ldy, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N-1,j) =               dy(0:N-1,j)*Metric(0:N-1,k)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N-1,j) =               dy(0:N-1,j)*Metric(0:N-1,k)
        enddo
      endif

! SECOND DERIVATIVE

      do j = 1,M
        do i = 1,N/2-1
          w1(2*i  ,j) =  b(i,2)*a(2*i  ,j)
          w1(2*i+1,j) =  b(i,2)*a(2*i+1,j)
        enddo
      enddo

      call CRFT (0, w1,  N3, dy, ldy, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N-1,j) = dd(0:N-1,j) + dy(0:N-1,j)*Metric(0:N-1,k+1)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + dy(0:N-1,j)*Metric(0:N-1,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do j = 1,M
        do i = 1,N/2-1
          w1(2*i  ,j) = -b(i,3)*a(2*i+1,j)
          w1(2*i+1,j) =  b(i,3)*a(2*i  ,j)
        enddo
      enddo

      call CRFT (0, w1,  N3, dy, ldy, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N-1,j) = dd(0:N-1,j) + dy(0:N-1,j)*Metric(0:N-1,k+2)
        enddo
  
        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + dy(0:N-1,j)*Metric(0:N-1,k+2)
        enddo
      endif

! FOURTH DERIVATIVE

      do j = 1,M
        do i = 1,N/2-1
          w1(2*i  ,j) =  b(i,4)*a(2*i  ,j)
          w1(2*i+1,j) =  b(i,4)*a(2*i+1,j)
        enddo
      enddo

      call CRFT (0, w1,  N3, dy, ldy, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

        do j = 1,M
          dy(0:N-1,j) = dd(0:N-1,j) + dy(0:N-1,j)*Metric(0:N-1,k+3)
        enddo

    endif

        call DeAllocate_WorkSpace (wk_Size) ; RETURN

#endif

#if defined(VFFT)
  
    ldw = SIZE(w1, DIM=1)

    do i = 0,N-1
      w1(:,i) = y(i,:)
    enddo

    call VRFFTF (M, N, w1, dy, ldw, WSave)

! -------------------------------------------------------------------
!    This is the case for NO MAPPING
! -------------------------------------------------------------------
    if (Map == 0) then

      w1(:,0) = ZERO

      if (MOD(O_Of_D,2) == 1) then
        do i = 1,N/2-1
          do j = 1,M
                dy_Real = -b(i,O_Of_D)*w1(j,2*i  )
                dy_Imag =  b(i,O_Of_D)*w1(j,2*i-1)

            w1(j,2*i-1) = dy_Real
            w1(j,2*i  ) = dy_Imag
          enddo
        enddo
      else
        do i = 1,N/2-1
          w1(:,2*i-1) =  b(i,O_Of_D)*w1(:,2*i-1)
          w1(:,2*i  ) =  b(i,O_Of_D)*w1(:,2*i  )
        enddo
      endif

      if (Mod(N,2) == 0) w1(:,N-1) = ZERO

      call VRFFTB (M, N, w1, dy, ldw, WSave)

      do i = 0,N-1
        dy(i,:) = w1(:,i)
      enddo

    endif

! --------------------------------------------------------------------
!    Here starts the MAPPING CASE
! --------------------------------------------------------------------
    if (Map /= 0) then

      a(:,1:N-1) = w1(:,1:N-1)

! FIRST DERIVATIVE

      w1(:,0) = ZERO

      do i = 1,N/2-1
        w1(:,2*i-1) = -b(i,1)*a(:,2*i  )
        w1(:,2*i  ) =  b(i,1)*a(:,2*i-1)
      enddo

      if (Mod(N,2) == 0) w1(:,N-1) = ZERO
 
      call VRFFTB (M, N, w1, dy, ldw, WSave)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N-1,j) =               w1(j,0:N-1)*Metric(0:N-1,k) 
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N-1,j) =               w1(j,0:N-1)*Metric(0:N-1,k) 
        enddo
      endif

! SECOND DERIVATIVE

      w1(:,0) = ZERO

      do i = 1,N/2-1
        w1(:,2*i-1) =  b(i,2)*a(:,2*i-1)
        w1(:,2*i  ) =  b(i,2)*a(:,2*i  )
      enddo

      if (Mod(N,2) == 0) w1(:,N-1) = ZERO

      call VRFFTB (M, N, w1, dy, ldw, WSave)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N-1,j) = dd(0:N-1,j) + w1(j,0:N-1)*Metric(0:N-1,k+1)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + w1(j,0:N-1)*Metric(0:N-1,k+1) 
        enddo
      endif

! THIRD DERIVATIVE

      w1(:,0) = ZERO

      do i = 1,N/2-1
        w1(:,2*i-1) = -b(i,3)*a(:,2*i  )
        w1(:,2*i  ) =  b(i,3)*a(:,2*i-1)
      enddo

      if (Mod(N,2) == 0) w1(:,N-1) = ZERO

      call VRFFTB (M, N, w1, dy, ldw, WSave)

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N-1,j) = dd(0:N-1,j) + w1(j,0:N-1)*Metric(0:N-1,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + w1(j,0:N-1)*Metric(0:N-1,k+2) 
        enddo
      endif

! FOURTH DERIVATIVE

      w1(:,0) = ZERO

      do i = 1,N/2-1
        w1(:,2*i-1) =  b(i,4)*a(:,2*i-1)
        w1(:,2*i  ) =  b(i,4)*a(:,2*i  )
      enddo

      if (Mod(N,2) == 0) w1(:,N-1) = ZERO

      call VRFFTB (M, N, w1, dy, ldw, WSave)

        do j = 1,M
          dy(0:N-1,j) = dd(0:N-1,j) + w1(j,0:N-1)*Metric(0:N-1,k+3)
        enddo

    endif

#endif

#if defined(CRAY)

    ldw = SIZE(w1, DIM=1)

    w1(0:N-1,:) = y(0:N-1,:)

    Factor = INT(WSave(1:19))

    Stride = 1

    call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M, -1)

! -------------------------------------------------------------------
!    This case is for NO MAPPING
! -------------------------------------------------------------------
    if (Map == 0) then

      if (MOD(O_Of_D,2) == 1) then
        do j = 1,M
          do i = 1,N/2
                dy_Real = -b(i,O_Of_D)*w1(2*i+1,j)
                dy_Imag =  b(i,O_Of_D)*w1(2*i  ,j)

            w1(2*i  ,j) = dy_Real
            w1(2*i+1,j) = dy_Imag
          enddo
        enddo
      else
        do j = 1,M
          do i = 1,N/2
            w1(2*i  ,j) =  b(i,O_Of_D)*w1(2*i  ,j)
            w1(2*i+1,j) =  b(i,O_Of_D)*w1(2*i+1,j)
          enddo
        enddo
      endif

      w1(0:1,:) = ZERO ; w1(N+1,:) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

      do j = 1,M
        dy(0:N-1,j) =  w1(0:N-1,j)
      enddo

    endif

! --------------------------------------------------------------------
!    Here starts the MAPPING CASE
! --------------------------------------------------------------------
    if (Map /= 0) then

      a(:,:) = w1(0:N+1,:)

! FIRST DERIVATIVE

      do j = 1,M
        do i = 1,N/2
          w1(2*i  ,j) = -b(i,1)*a(2*i+1,j)
          w1(2*i+1,j) =  b(i,1)*a(2*i  ,j)
        enddo
      enddo

      w1(0:1,:) = ZERO ; w1(N+1,:) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N-1,j) =               w1(0:N-1,j)*Metric(0:N-1,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N-1,j) =               w1(0:N-1,j)*Metric(0:N-1,k)
        enddo
      endif

! SECOND DERIVATIVE

      do j = 1,M
        do i = 1,N/2
          w1(2*i  ,j) =  b(i,2)*a(2*i  ,j)
          w1(2*i+1,j) =  b(i,2)*a(2*i+1,j)
        enddo
      enddo

      w1(0:1,:) = ZERO ; w1(N+1,:) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N-1,j) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+1)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do j = 1,M
        do i = 1,N/2
          w1(2*i  ,j) = -b(i,3)*a(2*i+1,j)
          w1(2*i+1,j) =  b(i,3)*a(2*i  ,j)
        enddo
      enddo

      w1(0:1,:) = ZERO ; w1(N+1,:) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N-1,j) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+2)
        enddo
      endif

! FOURTH DERIVATIVE

      do j = 1,M
        do i = 1,N/2
          w1(2*i  ,j) =  b(i,4)*a(2*i  ,j)
          w1(2*i+1,j) =  b(i,4)*a(2*i+1,j)
        enddo
      enddo

      w1(0:1,:) = ZERO ; w1(N+1,:) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

        do j = 1,M
          dy(0:N-1,j) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+3)
        enddo

    endif
#endif

    END Subroutine FFT_Diff_1_NONE 

  END Subroutine F_Diff_FFT_1_NONE 

! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_FFT_1_BOTH
! ----------------------------------------------------------------------
  Subroutine F_Diff_FFT_1_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer            :: N, M, ldy, O_Of_D, IEO

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: Map, Smooth, BlockSize
  integer  :: N0, N1, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_cos, i_sin, i_s_s, i_m_s

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

        Map = OPERATOR_MAP(Aux_Diff)
     Smooth = OPERATOR_SMOOTH_1(Aux_Diff)

  N0 = N/2 ; N1 = 2*N0+1

  i_cos = FFT_CWAVE_START_HALF(N1)
  i_sin = FFT_SWAVE_START_HALF(N1)
  i_s_s = FFT_SIGMA_START_HALF(N1)
  i_m_s = FFT_METRIC_START_HALF(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    if (IEO == 1) then
      call FFT_Diff_1_EVEN (N0, m_Vector, Aux_Diff(i_cos), Aux_Diff(i_sin), &
                                          Aux_Diff(i_m_s), Aux_Diff(i_s_s), &
                                          y(1,m_Start)   , dy(1,m_Start))
    else
      call FFT_Diff_1_ODD  (N0, m_Vector, Aux_Diff(i_cos), Aux_Diff(i_sin), &
                                          Aux_Diff(i_m_s), Aux_Diff(i_s_s), &
                                          y(1,m_Start)   , dy(1,m_Start))
    endif
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined(ESSL)
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined(VFFT) || defined(CRAY)
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    if (IEO == 1) then
      call FFT_Diff_1_EVEN (N0, m_Vector, Aux_Diff(i_cos), Aux_Diff(i_sin), &
                                          Aux_Diff(i_m_s), Aux_Diff(i_s_s), &
                                          y(1,m_Start)   , dy(1,m_Start))
    else
      call FFT_Diff_1_ODD  (N0, m_Vector, Aux_Diff(i_cos), Aux_Diff(i_sin), &
                                          Aux_Diff(i_m_s), Aux_Diff(i_s_s), &
                                          y(1,m_Start)   , dy(1,m_Start))
    endif

  endif

  CONTAINS

    Subroutine FFT_Diff_1_EVEN (N0, M, CWave, SWave, Metric, Sigma, y, dy)

    integer  :: N0, M

    REALTYPE, dimension(*)          :: CWave, SWave
    REALTYPE, dimension(0:N0,10)    :: Metric
    REALTYPE, dimension(0:N0)       :: Sigma
    REALTYPE, dimension(0:ldy-1,M)  :: y, dy

    integer  :: N, i, j, k, Even_Or_Odd
    REALTYPE, dimension(0:N0,4)     :: b

#if defined(ESSL)
    integer  :: N2
    integer  :: wc_Size, ws_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(0:N0,M)     :: a, dd
#endif
#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(M,0:N0)     :: w1, a
    REALTYPE, dimension(0:N0,M)     :: dd
#endif
#if defined(CRAY)
    REALTYPE, dimension(0:N0,M)     :: a, dd
#endif

#if defined(DEBUG)
    ASSERT(ldy>=N0+1,"F_Diff_1_Even:ldy<N+1")
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N0 ; Even_Or_Odd =  1

    call FFT_D_Coefficients (Map, Smooth, Even_Or_Odd, O_Of_D, N  , Sigma, b)

#if defined(ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    ws_Size = ESSL_SFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,1,ldy, dy,1,ldy, N2, M, Scale, CWave, wc_Size,  &
                                                    Work , wk_Size)

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then
      do j = 1,M
        dy(0:N,j) = b(0:N,O_Of_D)*dy(0:N,j)
      enddo

      if (MOD(O_Of_D,2) == 0) then
        call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, CWave, wc_Size,  &
                                                        Work , wk_Size)
      else
        dy(N,:) = ZERO

        call SINF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, SWave, ws_Size,  &
                                                        Work , wk_Size)

        dy(0,:) = ZERO ; dy(N,:) = ZERO
      endif                                                         
    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then

      a = dy(0:N,1:M)

! First Derivative

      do j = 1,M
        dy(1:N-1,j) = b(1:N-1,1)*a(1:N-1,j)
      enddo

      dy(N,:) = ZERO

      call SINF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, SWave, ws_Size,  &
                                                      Work , wk_Size)

      dy(0,:) = ZERO ; dy(N,:) = ZERO

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             dy(0:N,j)*Metric(0:N,k)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N,j) =             dy(0:N,j)*Metric(0:N,k)
        enddo
      endif

! Second Derivative

      do j = 1,M
        dy(0:N,j) = b(0:N,2)*a(0:N,j)
      enddo

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, CWave, wc_Size,  &
                                                      Work , wk_Size)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+1)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do j = 1,M
        dy(0:N-1,j) = b(0:N-1,3)*a(0:N-1,j)
      enddo

      dy(N,:) = ZERO

      call SINF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, SWave, ws_Size,  &
                                                      Work , wk_Size)

      dy(0,:) = ZERO ; dy(N,:) = ZERO

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+2)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+2)
        enddo
      endif

! Fourth Derivative

      do j = 1,M
        dy(0:N,j) = b(0:N,4)*a(0:N,j)
      enddo

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, CWave, wc_Size,  &
                                                      Work , wk_Size)

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+3)
        enddo

    endif

        call DeAllocate_WorkSpace (wk_Size) ; RETURN

#endif

#if defined(VFFT)

    ldw = SIZE(w1, DIM=1)

    do i = 0,N
      w1(:,i) = y(i,:)
    enddo

    call VCosT  (M, N+1, w1     , dy, ldw, CWave)
 
!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then

      do i = 0,N
        w1(:,i) = b(i,O_Of_D)*w1(:,i)
      enddo

      if (MOD(O_Of_D,2) == 0) then
        call VCosT  (M, N+1, w1     , dy, ldw, CWave)

        do i = 0,N
          dy(i,:) = w1(:,i)
        enddo
      else
        call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

        dy(0,:) = ZERO ; dy(N,:) = ZERO

        do i = 1,N-1
          dy(i,:) = w1(:,i)
        enddo
      endif

    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then
      a = w1(1:M,0:N)

! First Derivative

      do i = 0,N
        w1(:,i) = b(i,1)*a(:,i)
      enddo

      call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

      w1(:,0) = ZERO ; w1(:,N) = ZERO

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             w1(j,0:N)*Metric(0:N,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) =             w1(j,0:N)*Metric(0:N,k)
        enddo
      endif

! Second Derivative

      do i = 0,N
        w1(:,i) = b(i,2)*a(:,i)
      enddo

      call VCosT  (M, N+1, w1     , dy, ldw, CWave)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+1)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do i = 0,N
        w1(:,i) = b(i,3)*a(:,i)
      enddo

      call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

      w1(:,0) = ZERO ; w1(:,N) = ZERO

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+2)
        enddo
      endif

! Fourth Derivative

      do i = 0,N
        w1(:,i) = b(i,4)*a(:,i)
      enddo

      call VCosT  (M, N+1, w1     , dy, ldw, CWave)

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+3)
        enddo

    endif

#endif

#if defined(CRAY)

    call Cray_CFT_1 (N, M,  y,ldy, dy,ldy, CWave)

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then

      do j = 1,M
        dy(0:N,j) = b(0:N,O_Of_D)*dy(0:N,j)
      enddo

      if (MOD(O_Of_D,2) == 0) then
        call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, CWave)
      else
        call Cray_SFT_1 (N, M, dy,ldy, dy,ldy, SWave)

        dy(N,:) = ZERO
      endif

    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then

      a = dy(0:N,1:M)

! First Derivative

      do j = 1,M
        dy(0:N-1,j) = b(0:N-1,1)*a(0:N-1,j)
      enddo

      call Cray_SFT_1 (N, M, dy,ldy, dy,ldy, SWave)

      dy(N,:) = ZERO

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             dy(0:N,j)*Metric(0:N,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) =             dy(0:N,j)*Metric(0:N,k)
        enddo
      endif

! Second Derivative

      do j = 1,M
        dy(0:N,j) = b(0:N,2)*a(0:N,j)
      enddo

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, CWave)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+1)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do j = 1,M
        dy(0:N-1,j) = b(0:N-1,3)*a(0:N-1,j)
      enddo

      call Cray_SFT_1 (N, M, dy,ldy, dy,ldy, SWave)

      dy(N,:) = ZERO

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+2)
        enddo
      endif

! Fourth Derivative

      do j = 1,M
        dy(0:N,j) = b(0:N,4)*a(0:N,j)
      enddo

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, CWave)

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+3)
        enddo

    endif

#endif

    END Subroutine FFT_Diff_1_EVEN 
! ----------------------------------------------------------------------
    Subroutine FFT_Diff_1_ODD  (N0, M, CWave, SWave, Metric, Sigma, y, dy)

    integer  :: N0, M

    REALTYPE, dimension(*)          :: CWave, SWave
    REALTYPE, dimension(0:N0,10)    :: Metric
    REALTYPE, dimension(0:N0)       :: Sigma
    REALTYPE, dimension(0:ldy-1,M)  :: y, dy

    integer  :: N, i, j, k, Even_Or_Odd
    REALTYPE, dimension(0:N0,4)     :: b

#if defined(ESSL)
    integer  :: N2
    integer  :: wc_Size, ws_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(0:N0,M)     :: a, dd
#endif
#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(M,0:N0)     :: w1, a
    REALTYPE, dimension(0:N0,M)     :: dd
#endif
#if defined(CRAY)
    REALTYPE, dimension(0:N0,M)     :: a, dd
#endif

#if defined(DEBUG)
    ASSERT(ldy>=N0+1,"F_Diff_1_Odd:ldy<N+1")
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N0 ; Even_Or_Odd = -1

    call FFT_D_Coefficients (Map, Smooth, Even_Or_Odd, O_Of_D, N  , Sigma, b)

#if defined(ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    ws_Size = ESSL_SFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call   Allocate_WorkSpace (wk_Size)

    call SINF (0,  y,1,ldy, dy,1,ldy, N2, M, Scale, SWave, ws_Size,  &
                                                    Work , wk_Size)

    dy(N,:) = ZERO

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then
      do j = 1,M
        dy(0:N,j) = b(0:N,O_Of_D)*dy(0:N,j)
      enddo

      if (MOD(O_Of_D,2) /= 0) then
        call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, CWave, wc_Size,  &
                                                        Work , wk_Size)
      else
        call SINF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, SWave, ws_Size,  &
                                                        Work , wk_Size)

        dy(0,:) = ZERO ; dy(N,:) = ZERO
      endif                                                         
    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then

      a = dy(0:N,1:M)

! First Derivative

      do j = 1,M
        dy(0:N,j) = b(0:N,1)*a(0:N,j)
      enddo

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, CWave, wc_Size,  &
                                                      Work , wk_Size)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             dy(0:N,j)*Metric(0:N,k)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N,j) =             dy(0:N,j)*Metric(0:N,k)
        enddo
      endif

! Second Derivative

      do j = 1,M
        dy(0:N,j) = b(0:N,2)*a(0:N,j)
      enddo

      call SINF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, SWave, ws_Size,  &
                                                      Work , wk_Size)

      dy(0,:) = ZERO ; dy(N,:) = ZERO

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+1)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do j = 1,M
        dy(0:N,j) = b(0:N,3)*a(0:N,j)
      enddo

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, CWave, wc_Size,  &
                                                      Work , wk_Size)

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+2)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+2)
        enddo
      endif

! Fourth Derivative

      do j = 1,M
        dy(0:N,j) = b(0:N,4)*a(0:N,j)
      enddo

      call SINF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, SWave, ws_Size,  &
                                                      Work , wk_Size)

      dy(0,:) = ZERO ; dy(N,:) = ZERO

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+3)
        enddo

    endif

        call DeAllocate_WorkSpace (wk_Size) ; RETURN

#endif

#if defined(VFFT)

    ldw = SIZE(w1, DIM=1)

    do i = 0,N
      w1(:,i) = y(i,:)
    enddo

    call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)
 
    w1(:,N) = ZERO

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then

      do i = 0,N
        w1(:,i) = b(i,O_Of_D)*w1(:,i)
      enddo

      if (MOD(O_Of_D,2) /= 0) then
        call VCosT  (M, N+1, w1     , dy, ldw, CWave)

        do i = 0,N
          dy(i,:) = w1(:,i)
        enddo
      else
        call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

        dy(0,:) = ZERO ; dy(N,:) = ZERO

        do i = 1,N-1
          dy(i,:) = w1(:,i)
        enddo
      endif

    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then
      a = w1(1:M,0:N)

! First Derivative

      do i = 0,N
        w1(:,i) = b(i,1)*a(:,i)
      enddo

      call VCosT  (M, N+1, w1     , dy, ldw, CWave)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             w1(j,0:N)*Metric(0:N,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) =             w1(j,0:N)*Metric(0:N,k)
        enddo
      endif

! Second Derivative

      do i = 0,N
        w1(:,i) = b(i,2)*a(:,i)
      enddo

      call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

      w1(:,0) = ZERO ; w1(:,N) = ZERO

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+1)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do i = 0,N
        w1(:,i) = b(i,3)*a(:,i)
      enddo

      call VCosT  (M, N+1, w1     , dy, ldw, CWave)

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+2)
        enddo
      endif

! Fourth Derivative

      do i = 0,N
        w1(:,i) = b(i,4)*a(:,i)
      enddo

      call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

      w1(:,0) = ZERO ; w1(:,N) = ZERO

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+3)
        enddo

    endif

#endif

#if defined(CRAY)

    call Cray_SFT_1 (N, M,  y,ldy, dy,ldy, SWave)

    dy(N,:) = ZERO

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then

      do j = 1,M
        dy(0:N,j) = b(0:N,O_Of_D)*dy(0:N,j)
      enddo

      if (MOD(O_Of_D,2) /= 0) then
        call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, CWave)
      else
        call Cray_SFT_1 (N, M, dy,ldy, dy,ldy, SWave)

        dy(N,:) = ZERO
      endif

    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then

      a = dy(0:N,1:M)

! First Derivative

      do j = 1,M
        dy(0:N,j) = b(0:N,1)*a(0:N,j)
      enddo

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, CWave)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             dy(0:N,j)*Metric(0:N,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) =             dy(0:N,j)*Metric(0:N,k)
        enddo
      endif

! Second Derivative

      do j = 1,M
        dy(0:N-1,j) = b(0:N-1,2)*a(0:N-1,j)
      enddo

      call Cray_SFT_1 (N, M, dy,ldy, dy,ldy, SWave)

      dy(N,:) = ZERO

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+1)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do j = 1,M
        dy(0:N,j) = b(0:N,3)*a(0:N,j)
      enddo

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, CWave)

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+2)
        enddo
      endif

! Fourth Derivative

      do j = 1,M
        dy(0:N-1,j) = b(0:N-1,4)*a(0:N-1,j)
      enddo

      call Cray_SFT_1 (N, M, dy,ldy, dy,ldy, SWave)

      dy(N,:) = ZERO

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+3)
        enddo

    endif

#endif

    END Subroutine FFT_Diff_1_ODD

  END Subroutine F_Diff_FFT_1_BOTH 

END MODULE Fourier_Differentiation_FFT_1



! ----------------------------------------------------------------------
! FUNCTION NAME: F_Diff_FFT_2
! DESCRIPTION  : does the Fourier solution smoothing on m vectors
!                using the transform-recursion algorithm.
!                This is a "high-level" Subroutine that does
!                unpacking and blocking.
!
! ----------------------------------------------------------------------

MODULE Fourier_Differentiation_FFT_2
  USE PS_WorkSpace_Allocation, Work=>WorkSpace
  USE FFT_Auxillary_Subroutine

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

PUBLIC  :: F_Diff_FFT_2

CONTAINS

  Subroutine F_Diff_FFT_2      (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer            :: N, M, ldy, O_Of_D

  integer , OPTIONAL :: IEO
  logical            :: L_IEO 

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /=  0) then
      call F_Diff_FFT_2_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)
    else
      call F_Diff_FFT_2_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
    endif
  else
      call F_Diff_FFT_2_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
  endif

  END Subroutine F_Diff_FFT_2

! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_FFT_2_NONE
! ----------------------------------------------------------------------
  Subroutine F_Diff_FFT_2_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)

  integer            :: N, M, ldy, O_Of_D

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: Map, Smooth, BlockSize
  integer  :: N0, N1, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_fft, i_s_s, i_m_s

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

        Map = OPERATOR_MAP(Aux_Diff)
     Smooth = OPERATOR_SMOOTH_1(Aux_Diff)

  N0 = N ; N1 = N

  i_fft = FFT_FWAVE_START_FULL(N1)
  i_s_s = FFT_SIGMA_START_FULL(N1)
  i_m_s = FFT_METRIC_START_FULL(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call FFT_Diff_2_NONE (N0, m_Vector, Aux_Diff(i_fft), Aux_Diff(i_m_s),  &
                                        Aux_Diff(i_s_s),                   &
                                        y(m_Start,1)   , dy(m_Start,1))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined(ESSL)
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined(VFFT) || defined(CRAY)
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    call FFT_Diff_2_NONE (N0, m_Vector, Aux_Diff(i_fft), Aux_Diff(i_m_s),  &
                                        Aux_Diff(i_s_s),                   &
                                        y(m_Start,1)   , dy(m_Start,1))

  endif

  CONTAINS

    Subroutine FFT_Diff_2_NONE (N0, M, WSave, Metric, Sigma, y, dy)

    integer  :: N0, M

    REALTYPE, dimension(*)          :: WSave
    REALTYPE, dimension(0:N0-1,10)  :: Metric
    REALTYPE, dimension(0:N0/2)     :: Sigma
    REALTYPE, dimension(ldy,0:N0-1) :: y, dy

    integer  :: N, i, j, k, Even_Or_Odd
    REALTYPE :: dy_Real, dy_Imag

    REALTYPE, dimension(0:N0/2,4)   :: b

#if defined(ESSL)
    integer  :: N2, N3, ISIGN, ldw
    integer  :: i_fft, i_bft
    integer  :: wf_Size, wb_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(1:N0-1,M)   ::  a
    REALTYPE, dimension(0:N0-1,M)   :: dd
    REALTYPE, dimension(0:N0+1,M)   :: w1
    REALTYPE, dimension(1)          :: w3
#endif
#if defined(VFFT)
    REALTYPE, dimension(M,1:N0-1)   ::  a
    REALTYPE, dimension(M,0:N0-1)   :: dd
    REALTYPE, dimension(ldy,0:N0-1) :: w1
#endif
#if defined(CRAY)
    integer  :: ldw, Stride
    integer, dimension(19)          :: Factor

    REALTYPE, dimension(0:N0+1,M)   ::  a
    REALTYPE, dimension(0:N0-1,M)   :: dd
    REALTYPE, dimension(0:N0+2,M)   :: w1
    REALTYPE, dimension(0:N0-1,M,2) :: w2
#endif

#if defined(DEBUG)
    ASSERT(ldy>=M,"F_Diff_2_None:ldy<M")
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N0 ; Even_Or_Odd =  0

    call FFT_D_Coefficients (Map, Smooth, Even_Or_Odd, O_Of_D, N/2, Sigma, b)

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

    ISIGN =  1 

    call RCFT (0, w1, ldw, w1,  N3, N2, M, ISIGN, Scale,  &
                  WSave(i_fft), wf_Size, Work, wk_Size, w3, 0)

    ISIGN = -1

! -------------------------------------------------------------------
!    This case is for NO MAPPING
! -------------------------------------------------------------------
    if (Map == 0) then

      w1(0,:) = ZERO ; w1(N,:) = ZERO

      if (MOD(O_Of_D,2) == 1) then
        do j = 1,M
          do i = 1,N/2-1
                dy_Real = -b(i,O_Of_D)*w1(2*i+1,j)
                dy_Imag =  b(i,O_Of_D)*w1(2*i  ,j)

            w1(2*i  ,j) = dy_Real
            w1(2*i+1,j) = dy_Imag
          enddo
        enddo
      else
        do j = 1,M
          do i = 1,N/2-1
            w1(2*i  ,j) =  b(i,O_Of_D)*w1(2*i  ,j)
            w1(2*i+1,j) =  b(i,O_Of_D)*w1(2*i+1,j)
          enddo
        enddo
      endif

      call CRFT (0, w1,  N3, w1, ldw, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

      do j = 1,M
        dy(j,0:N-1) = w1(0:N-1,j)
      enddo

    endif

! --------------------------------------------------------------------
!    Here starts the MAPPING CASE
! --------------------------------------------------------------------
    if (Map /= 0) then

      a = w1(1:N-1,:)

! FIRST DERIVATIVE

      w1(0,:) = ZERO ; w1(N,:) = ZERO

      do j = 1,M
        do i = 1,N/2-1
          w1(2*i  ,j) = -b(i,1)*a(2*i+1,j)
          w1(2*i+1,j) =  b(i,1)*a(2*i  ,j)
        enddo
      enddo

      call CRFT (0, w1,  N3, w1, ldw, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(j,0:N-1) =               w1(0:N-1,j)*Metric(0:N-1,k)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N-1,j) =               w1(0:N-1,j)*Metric(0:N-1,k)
        enddo
      endif

! SECOND DERIVATIVE

      w1(0,:) = ZERO ; w1(N,:) = ZERO

      do j = 1,M
        do i = 1,N/2-1
          w1(2*i  ,j) =  b(i,2)*a(2*i  ,j)
          w1(2*i+1,j) =  b(i,2)*a(2*i+1,j)
        enddo
      enddo

      call CRFT (0, w1,  N3, w1, ldw, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(j,0:N-1) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+1)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      w1(0,:) = ZERO ; w1(N,:) = ZERO

      do j = 1,M
        do i = 1,N/2-1
          w1(2*i  ,j) = -b(i,3)*a(2*i+1,j)
          w1(2*i+1,j) =  b(i,3)*a(2*i  ,j)
        enddo
      enddo

      call CRFT (0, w1,  N3, w1, ldw, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

      if (O_Of_D == 3) then
        do j = 1,M
          dy(j,0:N-1) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+2)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+2)
        enddo
      endif

! FOURTH DERIVATIVE

      w1(0,:) = ZERO ; w1(N,:) = ZERO

      do j = 1,M
        do i = 1,N/2-1
          w1(2*i  ,j) =  b(i,4)*a(2*i  ,j)
          w1(2*i+1,j) =  b(i,4)*a(2*i+1,j)
        enddo
      enddo

      call CRFT (0, w1,  N3, w1, ldw, N2, M, ISIGN, Scale,  &
                    WSave(i_bft), wb_Size, Work, wk_Size, w3, 0)

        do j = 1,M
          dy(j,0:N-1) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+3)
        enddo

    endif

        call DeAllocate_WorkSpace (wk_Size) ; RETURN

#endif

#if defined(VFFT)

    w1(1:M,:) = y(1:M,:)

    call VRFFTF (M, N, w1, dy, ldy, WSave)
 
! -------------------------------------------------------------------
!    This is the case for NO MAPPING
! -------------------------------------------------------------------
    if (Map == 0) then

      dy(1:M,0) = ZERO

      if (MOD(O_Of_D,2) == 1) then
        do i = 1,N/2-1
          dy(1:M,2*i-1) = -b(i,O_Of_D)*w1(1:M,2*i  )
          dy(1:M,2*i  ) =  b(i,O_Of_D)*w1(1:M,2*i-1)
        enddo
      else
        do i = 1,N/2-1
          dy(1:M,2*i-1) =  b(i,O_Of_D)*w1(1:M,2*i-1)
          dy(1:M,2*i  ) =  b(i,O_Of_D)*w1(1:M,2*i  )
        enddo
      endif

      if (Mod(N,2) == 0) dy(1:M,N-1) = ZERO

      call VRFFTB (M, N, dy, w1, ldy, WSave)

    endif

! --------------------------------------------------------------------
!    Here starts the MAPPING CASE
! --------------------------------------------------------------------
    if (Map /= 0) then

      a(:,1:N-1) = w1(1:M,1:N-1)

! FIRST DERIVATIVE

      dy(1:M,0) = ZERO

      do i = 1,N/2-1
        dy(1:M,2*i-1) = -b(i,1)*a(:,2*i  )
        dy(1:M,2*i  ) =  b(i,1)*a(:,2*i-1)
      enddo

      if (Mod(N,2) == 0) dy(1:M,N-1) = ZERO
 
      call VRFFTB (M, N, dy, w1, ldy, WSave)

      if (O_Of_D == 1) then
        do i = 0,N-1
          dy(1:M,i)   =               dy(1:M,i)  *Metric(i,k)
        enddo

        RETURN
      else
        do i = 0,N-1
          dd(1:M,i)   =               dy(1:M,i)  *Metric(i,k)
        enddo
      endif

! SECOND DERIVATIVE

      dy(1:M,0) = ZERO

      do i = 1,N/2-1
        dy(1:M,2*i-1) =  b(i,2)*a(:,2*i-1)
        dy(1:M,2*i  ) =  b(i,2)*a(:,2*i  )
      enddo

      if (Mod(N,2) == 0) dy(1:M,N-1) = ZERO

      call VRFFTB (M, N, dy, w1, ldy, WSave)

      if (O_Of_D == 2) then
        do i = 0,N-1
          dy(1:M,i)   = dd(1:M,i)   + dy(1:M,i)  *Metric(i,k+1)
        enddo

        RETURN
      else
        do i = 0,N-1
          dd(1:M,i)   = dd(1:M,i)   + dy(1:M,i)  *Metric(i,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      dy(1:M,0) = ZERO

      do i = 1,N/2-1
        dy(1:M,2*i-1) = -b(i,3)*a(:,2*i  )
        dy(1:M,2*i  ) =  b(i,3)*a(:,2*i-1)
      enddo

      if (Mod(N,2) == 0) dy(1:M,N-1) = ZERO

      call VRFFTB (M, N, dy, w1, ldy, WSave)

      if (O_Of_D == 3) then
        do i = 0,N-1
          dy(1:M,i)   = dd(1:M,i)   + dy(1:M,i)  *Metric(i,k+2)
        enddo

        RETURN
      else
        do i = 0,N-1
          dd(1:M,i)   = dd(1:M,i)   + dy(1:M,i)  *Metric(i,k+2)
        enddo
      endif

! FOURTH DERIVATIVE

      dy(1:M,0) = ZERO

      do i = 1,N/2-1
        dy(1:M,2*i-1) =  b(i,4)*a(:,2*i-1)
        dy(1:M,2*i  ) =  b(i,4)*a(:,2*i  )
      enddo

      if (Mod(N,2) == 0) dy(1:M,N-1) = ZERO

      call VRFFTB (M, N, dy, w1, ldy, WSave)

        do i = 0,N-1
          dy(1:M,i)   = dd(1:M,i)   + dy(1:M,i)  *Metric(i,k+3)
        enddo

    endif

#endif

#if defined(CRAY)
 
    do j = 1,M
      w1(0:N-1,j) = y(j,0:N-1)
    enddo

    ldw = SIZE(w1, DIM=1)

    Factor = INT(WSave(1:19))

    Stride = 1

    call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M, -1)

! -------------------------------------------------------------------
!    This case is for NO MAPPING
! -------------------------------------------------------------------
    if (Map == 0) then

      if (MOD(O_Of_D,2) == 1) then
        do j = 1,M
          do i = 1,N/2
                dy_Real = -b(i,O_Of_D)*w1(2*i+1,j)
                dy_Imag =  b(i,O_Of_D)*w1(2*i  ,j)

            w1(2*i  ,j) = dy_Real
            w1(2*i+1,j) = dy_Imag
          enddo
        enddo
      else
        do j = 1,M
          do i = 1,N/2
            w1(2*i  ,j) =  b(i,O_Of_D)*w1(2*i  ,j)
            w1(2*i+1,j) =  b(i,O_Of_D)*w1(2*i+1,j)
          enddo
        enddo
      endif

      w1(0:1,1:M) = ZERO ; w1(N+1,1:M) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

      do j = 1,M
        dy(j,0:N-1) = w1(0:N-1,j)
      enddo

    endif

! --------------------------------------------------------------------
!    Here starts the MAPPING CASE
! --------------------------------------------------------------------
    if (Map /= 0) then

      a = w1(0:N+1,:)

! FIRST DERIVATIVE

      do j = 1,M
        do i = 1,N/2
          w1(2*i  ,j) = -b(i,1)*a(2*i+1,j)
          w1(2*i+1,j) =  b(i,1)*a(2*i  ,j)
        enddo
      enddo

      w1(0:1,1:M) = ZERO ; w1(N+1,1:M) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(j,0:N-1) =               w1(0:N-1,j)*Metric(0:N-1,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N-1,j) =               w1(0:N-1,j)*Metric(0:N-1,k)
        enddo
      endif

! SECOND DERIVATIVE

      do j = 1,M
        do i = 1,N/2
          w1(2*i  ,j) =  b(i,2)*a(2*i  ,j)
          w1(2*i+1,j) =  b(i,2)*a(2*i+1,j)
        enddo
      enddo

      w1(0:1,1:M) = ZERO ; w1(N+1,1:M) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(j,0:N-1) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+1)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do j = 1,M
        do i = 1,N/2
          w1(2*i  ,j) = -b(i,3)*a(2*i+1,j)
          w1(2*i+1,j) =  b(i,3)*a(2*i  ,j)
        enddo
      enddo

      w1(0:1,1:M) = ZERO ; w1(N+1,1:M) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

      if (O_Of_D == 3) then
        do j = 1,M
          dy(j,0:N-1) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N-1,j) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+2)
        enddo
      endif

! FOURTH DERIVATIVE

      do j = 1,M
        do i = 1,N/2
          w1(2*i  ,j) =  b(i,4)*a(2*i  ,j)
          w1(2*i+1,j) =  b(i,4)*a(2*i+1,j)
        enddo
      enddo

      w1(0:1,1:M) = ZERO ; w1(N+1,1:M) = ZERO

      call rfftmlt (w1, w2, WSave(20), Factor, Stride, ldw, N, M,  1)

        do j = 1,M
          dy(j,0:N-1) = dd(0:N-1,j) + w1(0:N-1,j)*Metric(0:N-1,k+3)
        enddo

    endif

#endif

    END Subroutine FFT_Diff_2_NONE 

  END Subroutine F_Diff_FFT_2_NONE 

! ----------------------------------------------------------------------
! FUNCTION NAME: Diff_FFT_2_BOTH
! ----------------------------------------------------------------------
  Subroutine F_Diff_FFT_2_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer            :: N, M, ldy, O_Of_D, IEO

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: Map, Smooth, BlockSize
  integer  :: N0, N1, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_cos, i_sin, i_s_s, i_m_s

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

        Map = OPERATOR_MAP(Aux_Diff)
     Smooth = OPERATOR_SMOOTH_1(Aux_Diff)

  N0 = N/2 ; N1 = 2*N0+1

  i_cos = FFT_CWAVE_START_HALF(N1)
  i_sin = FFT_SWAVE_START_HALF(N1)
  i_s_s = FFT_SIGMA_START_HALF(N1)
  i_m_s = FFT_METRIC_START_HALF(N1)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    if (IEO == 1) then
      call FFT_Diff_2_EVEN (N0, m_Vector, Aux_Diff(i_cos), Aux_Diff(i_sin), &
                                          Aux_Diff(i_m_s), Aux_Diff(i_s_s), &
                                          y(m_Start,1)   , dy(m_Start,1))
    else
      call FFT_Diff_2_ODD  (N0, m_Vector, Aux_Diff(i_cos), Aux_Diff(i_sin), &
                                          Aux_Diff(i_m_s), Aux_Diff(i_s_s), &
                                          y(m_Start,1)   , dy(m_Start,1))
    endif
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined(ESSL)
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined(VFFT) || defined(CRAY)
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    if (IEO == 1) then
      call FFT_Diff_2_EVEN (N0, m_Vector, Aux_Diff(i_cos), Aux_Diff(i_sin), &
                                          Aux_Diff(i_m_s), Aux_Diff(i_s_s), &
                                          y(m_Start,1)   , dy(m_Start,1))
    else
      call FFT_Diff_2_ODD  (N0, m_Vector, Aux_Diff(i_cos), Aux_Diff(i_sin), &
                                          Aux_Diff(i_m_s), Aux_Diff(i_s_s), &
                                          y(m_Start,1)   , dy(m_Start,1))
    endif

  endif

  CONTAINS

    Subroutine FFT_Diff_2_EVEN (N0, M, CWave, SWave, Metric, Sigma, y, dy)

    integer  :: N0, M

    REALTYPE, dimension(*)          :: CWave, SWave
    REALTYPE, dimension(0:N0,10)    :: Metric
    REALTYPE, dimension(0:N0)       :: Sigma
    REALTYPE, dimension(ldy,0:N0)   :: y, dy

    integer  :: N, i, k, Even_Or_Odd
    REALTYPE, dimension(0:N0,4)     :: b

#if defined(ESSL)
    integer  :: N2
    integer  :: wc_Size, ws_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(M,0:N0)     :: a, dd
#endif
#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(ldy,0:N0)   :: w1
    REALTYPE, dimension(M,0:N0)     :: a, dd
#endif
#if defined(CRAY)
    REALTYPE, dimension(M,0:N0)     :: a, dd
#endif

#if defined(DEBUG)
    ASSERT(ldy>=M,"F_Diff_2_Even:ldy<M")
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N0 ; Even_Or_Odd =  1

    call FFT_D_Coefficients (Map, Smooth, Even_Or_Odd, O_Of_D, N  , Sigma, b)

#if defined(ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    ws_Size = ESSL_SFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,ldy,1, dy,ldy,1, N2, M, Scale, CWave, wc_Size,  &
                                                    Work , wk_Size)

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then
      do i = 0,N
        dy(1:M,i) = b(i,O_Of_D)*dy(1:M,i)
      enddo

      if (MOD(O_Of_D,2) == 0) then
        call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, CWave, wc_Size,  &
                                                        Work , wk_Size)
      else
        dy(1:M,N) = ZERO

        call SINF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, SWave, ws_Size,  &
                                                        Work , wk_Size)
        dy(1:M,0) = ZERO ; dy(1:M,N) = ZERO
      endif                                                         
    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then

      a = dy(1:M,0:N)

! First Derivative

      do i = 0,N-1
        dy(1:M,i) = b(i,1)*a(1:M,i)
      enddo

      dy(1:M,N) = ZERO

      call SINF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, SWave, ws_Size,  &
                                                      Work , wk_Size)

      dy(1:M,0) = ZERO ; dy(1:M,N) = ZERO

      if (O_Of_D == 1) then
        do i = 0,N
          dy(1:M,i) =             dy(1:M,i)*Metric(i,k)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do i = 0,N
          dd(1:M,i) =             dy(1:M,i)*Metric(i,k)
        enddo
      endif

! Second Derivative

      do i = 0,N
        dy(1:M,i) = b(i,2)*a(1:M,i)
      enddo

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, CWave, wc_Size,  &
                                                      Work , wk_Size)

      if (O_Of_D == 2) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+1)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do i = 0,N-1
        dy(1:M,i) = b(i,3)*a(1:M,i)
      enddo

      dy(1:M,N) = ZERO

      call SINF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, SWave, ws_Size,  &
                                                      Work , wk_Size)

      dy(1:M,0) = ZERO ; dy(1:M,N) = ZERO

      if (O_Of_D == 3) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+2)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+2)
        enddo
      endif

! Fourth Derivative

      do i = 0,N
        dy(1:M,i) = b(i,4)*a(1:M,i)
      enddo

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, CWave, wc_Size,  &
                                                      Work , wk_Size)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+3)
        enddo

    endif

        call DeAllocate_WorkSpace (wk_Size) ; RETURN

#endif

#if defined(VFFT)

    ldw = SIZE(w1, DIM=1)

    w1(1:M,:) = y(1:M,:)

    call VCosT  (M, N+1, w1     , dy, ldw, CWave)
 
!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then

      do i = 0,N
        w1(1:M,i) = b(i,O_Of_D)*w1(1:M,i)
      enddo

      if (MOD(O_Of_D,2) == 0) then
        call VCosT  (M, N+1, w1     , dy, ldw, CWave)

        dy(1:M,:) = w1(1:M,:)
      else
        call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

        dy(1:M,0) = ZERO ; dy(1:M,N) = ZERO

        dy(1:M,1:N-1) = w1(1:M,1:N-1)
      endif

    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then
      a = w1(1:M,0:N)

! First Derivative

      do i = 0,N
        w1(1:M,i) = b(i,1)*a(1:M,i)
      enddo

      call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

      w1(1:M,0) = ZERO ; w1(1:M,N) = ZERO

      if (O_Of_D == 1) then
        do i = 0,N
          dy(1:M,i) =             w1(1:M,i)*Metric(i,k)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) =             w1(1:M,i)*Metric(i,k)
        enddo
      endif

! Second Derivative

      do i = 0,N
        w1(1:M,i) = b(i,2)*a(1:M,i)
      enddo

      call VCosT  (M, N+1, w1     , dy, ldw, CWave)

      if (O_Of_D == 2) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+1)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do i = 0,N
        w1(1:M,i) = b(i,3)*a(1:M,i)
      enddo

      call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

      w1(1:M,0) = ZERO ; w1(1:M,N) = ZERO

      if (O_Of_D == 3) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+2)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+2)
        enddo
      endif

! Fourth Derivative

      do i = 0,N
        w1(1:M,i) = b(i,4)*a(1:M,i)
      enddo

      call VCosT  (M, N+1, w1     , dy, ldw, CWave)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+3)
        enddo

    endif

#endif

#if defined(CRAY)

    call Cray_CFT_2 (N, M,  y,ldy, dy,ldy, CWave)

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then

      do i = 0,N
        dy(1:M,i) = b(i,O_Of_D)*dy(1:M,i)
      enddo

      if (MOD(O_Of_D,2) == 0) then
        call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, CWave)
      else
        call Cray_SFT_2 (N, M, dy,ldy, dy,ldy, SWave)

        dy(1:M,N) = ZERO
      endif

    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then

      a = dy(1:M,0:N)

! First Derivative

      do i = 0,N-1
        dy(1:M,i) = b(i,1)*a(1:M,i)
      enddo

      call Cray_SFT_2 (N, M, dy,ldy, dy,ldy, SWave)

      dy(1:M,N) = ZERO

      if (O_Of_D == 1) then
        do i = 0,N
          dy(1:M,i) =             dy(1:M,i)*Metric(i,k)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) =             dy(1:M,i)*Metric(i,k)
        enddo
      endif

! Second Derivative

      do i = 0,N
        dy(1:M,i) = b(i,2)*a(1:M,i)
      enddo

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, CWave)

      if (O_Of_D == 2) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+1)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do i = 0,N-1
        dy(1:M,i) = b(i,3)*a(1:M,i)
      enddo

      call Cray_SFT_2 (N, M, dy,ldy, dy,ldy, SWave)

      dy(1:M,N) = ZERO

      if (O_Of_D == 3) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+2)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+2)
        enddo
      endif

! Fourth Derivative

      do i = 0,N
        dy(1:M,i) = b(i,4)*a(1:M,i)
      enddo

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, CWave)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+3)
        enddo

    endif

#endif

    END Subroutine FFT_Diff_2_EVEN
! ----------------------------------------------------------------------
    Subroutine FFT_Diff_2_ODD  (N0, M, CWave, SWave, Metric, Sigma, y, dy)

    integer  :: N0, M

    REALTYPE, dimension(*)          :: CWave, SWave
    REALTYPE, dimension(0:N0,10)    :: Metric
    REALTYPE, dimension(0:N0)       :: Sigma
    REALTYPE, dimension(ldy,0:N0)   :: y, dy

    integer  :: N, i, k, Even_Or_Odd
    REALTYPE, dimension(0:N0,4)     :: b

#if defined(ESSL)
    integer  :: N2
    integer  :: wc_Size, ws_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(M,0:N0)     :: a, dd
#endif
#if defined(VFFT)
    integer  :: ldw

    REALTYPE, dimension(ldy,0:N0)   :: w1
    REALTYPE, dimension(M,0:N0)     :: a, dd
#endif
#if defined(CRAY)
    REALTYPE, dimension(M,0:N0)     :: a, dd
#endif

#if defined(DEBUG)
    ASSERT(ldy>=M,"F_Diff_2_Odd:ldy<M")
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N0 ; Even_Or_Odd = -1

    call FFT_D_Coefficients (Map, Smooth, Even_Or_Odd, O_Of_D, N  , Sigma, b)

#if defined(ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    ws_Size = ESSL_SFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call   Allocate_WorkSpace (wk_Size)

    call SINF (0,  y,ldy,1, dy,ldy,1, N2, M, Scale, SWave, ws_Size,  &
                                                    Work , wk_Size)

    dy(1:M,N) = ZERO

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then
      do i = 0,N
        dy(1:M,i) = b(i,O_Of_D)*dy(1:M,i)
      enddo

      if (MOD(O_Of_D,2) /= 0) then
        call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, CWave, wc_Size,  &
                                                        Work , wk_Size)
      else
        call SINF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, SWave, ws_Size,  &
                                                        Work , wk_Size)
        dy(1:M,0) = ZERO ; dy(1:M,N) = ZERO
      endif                                                         
    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then

      a = dy(1:M,0:N)

! First Derivative

      do i = 0,N
        dy(1:M,i) = b(i,1)*a(1:M,i)
      enddo

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, CWave, wc_Size,  &
                                                      Work , wk_Size)

      if (O_Of_D == 1) then
        do i = 0,N
          dy(1:M,i) =             dy(1:M,i)*Metric(i,k)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do i = 0,N
          dd(1:M,i) =             dy(1:M,i)*Metric(i,k)
        enddo
      endif

! Second Derivative

      do i = 0,N
        dy(1:M,i) = b(i,2)*a(1:M,i)
      enddo

      call SINF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, SWave, ws_Size,  &
                                                      Work , wk_Size)

      dy(1:M,0) = ZERO ; dy(1:M,N) = ZERO

      if (O_Of_D == 2) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+1)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do i = 0,N
        dy(1:M,i) = b(i,3)*a(1:M,i)
      enddo

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, CWave, wc_Size,  &
                                                      Work , wk_Size)

      if (O_Of_D == 3) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+2)
        enddo

        call DeAllocate_WorkSpace (wk_Size) ; RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+2)
        enddo
      endif

! Fourth Derivative

      do i = 0,N
        dy(1:M,i) = b(i,4)*a(1:M,i)
      enddo

      call SINF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, SWave, ws_Size,  &
                                                      Work , wk_Size)

      dy(1:M,0) = ZERO ; dy(1:M,N) = ZERO

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+3)
        enddo

    endif

        call DeAllocate_WorkSpace (wk_Size) ; RETURN

#endif

#if defined(VFFT)

    ldw = SIZE(w1, DIM=1)

    w1(1:M,:) = y(1:M,:)

    call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)
 
    w1(1:M,N) = ZERO

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then

      do i = 0,N
        w1(1:M,i) = b(i,O_Of_D)*w1(1:M,i)
      enddo

      if (MOD(O_Of_D,2) /= 0) then
        call VCosT  (M, N+1, w1     , dy, ldw, CWave)

        dy(1:M,:) = w1(1:M,:)
      else
        call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

        dy(1:M,0) = ZERO ; dy(1:M,N) = ZERO

        dy(1:M,1:N-1) = w1(1:M,1:N-1)
      endif

    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then
      a = w1(1:M,0:N)

! First Derivative

      do i = 0,N
        w1(1:M,i) = b(i,1)*a(1:M,i)
      enddo

      call VCosT  (M, N+1, w1     , dy, ldw, CWave)

      if (O_Of_D == 1) then
        do i = 0,N
          dy(1:M,i) =             w1(1:M,i)*Metric(i,k)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) =             w1(1:M,i)*Metric(i,k)
        enddo
      endif

! Second Derivative

      do i = 0,N
        w1(1:M,i) = b(i,2)*a(1:M,i)
      enddo

      call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

      w1(1:M,0) = ZERO ; w1(1:M,N) = ZERO

      if (O_Of_D == 2) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+1)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do i = 0,N
        w1(1:M,i) = b(i,3)*a(1:M,i)
      enddo

      call VCosT  (M, N+1, w1     , dy, ldw, CWave)

      if (O_Of_D == 3) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+2)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+2)
        enddo
      endif

! Fourth Derivative

      do i = 0,N
        w1(1:M,i) = b(i,4)*a(1:M,i)
      enddo

      call VSinT  (M, N-1, w1(1,1), dy, ldw, SWave)

      w1(1:M,0) = ZERO ; w1(1:M,N) = ZERO

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+3)
        enddo

    endif

#endif

#if defined(CRAY)

    call Cray_SFT_2 (N, M,  y,ldy, dy,ldy, SWave)

    dy(1:M,N) = ZERO

!---------------------------
! No Mapping
!---------------------------
    if (Map == 0) then

      do i = 0,N
        dy(1:M,i) = b(i,O_Of_D)*dy(1:M,i)
      enddo

      if (MOD(O_Of_D,2) /= 0) then
        call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, CWave)
      else
        call Cray_SFT_2 (N, M, dy,ldy, dy,ldy, SWave)

        dy(1:M,N) = ZERO
      endif

    endif

!---------------------------
! With Mapping
!---------------------------
    if (Map /= 0) then

      a = dy(1:M,0:N)

! First Derivative

      do i = 0,N
        dy(1:M,i) = b(i,1)*a(1:M,i)
      enddo

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, CWave)

      if (O_Of_D == 1) then
        do i = 0,N
          dy(1:M,i) =             dy(1:M,i)*Metric(i,k)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) =             dy(1:M,i)*Metric(i,k)
        enddo
      endif

! Second Derivative

      do i = 0,N-1
        dy(1:M,i) = b(i,2)*a(1:M,i)
      enddo

      call Cray_SFT_2 (N, M, dy,ldy, dy,ldy, SWave)

      dy(1:M,N) = ZERO

      if (O_Of_D == 2) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+1)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+1)
        enddo
      endif

! THIRD DERIVATIVE

      do i = 0,N
        dy(1:M,i) = b(i,3)*a(1:M,i)
      enddo

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, CWave)

      if (O_Of_D == 3) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+2)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+2)
        enddo
      endif

! Fourth Derivative

      do i = 0,N-1
        dy(1:M,i) = b(i,4)*a(1:M,i)
      enddo

      call Cray_SFT_2 (N, M, dy,ldy, dy,ldy, SWave)

      dy(1:M,N) = ZERO

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+3)
        enddo

    endif

#endif

    END Subroutine FFT_Diff_2_ODD

  END Subroutine F_Diff_FFT_2_BOTH 

END MODULE Fourier_Differentiation_FFT_2
