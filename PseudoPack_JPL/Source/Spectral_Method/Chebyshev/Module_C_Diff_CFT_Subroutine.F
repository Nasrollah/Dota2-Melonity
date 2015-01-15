#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------
! FUNCTION NAME: C_Diff_CFT_1
! DESCRIPTION  : does the chebyshev derivative on m vectors
!                using the transform-recursion algorithm.
!                This is a "high-level" Subroutine that does
!                unpacking and blocking.
! ----------------------------------------------------------------------

MODULE Chebyshev_Differentiation_CFT_1
  USE PS_WorkSpace_Allocation, Work=>WorkSpace
  USE CFT_Auxillary_Subroutine

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

PUBLIC  :: C_Diff_CFT_1

CONTAINS

  Subroutine C_Diff_CFT_1      (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer            :: N, M, ldy, O_Of_D

  integer , OPTIONAL :: IEO
  logical            :: L_IEO 

  REALTYPE, dimension(ldy,M) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /= 0) then
      call C_Diff_CFT_1_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)
    else
      call C_Diff_CFT_1_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
    endif
  else
      call C_Diff_CFT_1_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
  endif

  END Subroutine C_Diff_CFT_1

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Diff_CFT_1_NONE
! ----------------------------------------------------------------------
  Subroutine C_Diff_CFT_1_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)

  integer            :: N, M, ldy, O_Of_D

  REALTYPE, dimension(ldy,M) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: Map, Smooth_1, Smooth_2, BlockSize
  integer  :: N1, N_Pnts, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_w_s, i_s_s, i_m_s

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
   Smooth_1 = OPERATOR_SMOOTH_1(Aux_Diff)
   Smooth_2 = OPERATOR_SMOOTH_2(Aux_Diff)

  N1 = N ; N_Pnts = N

  i_w_s = CFT_WSAVE_START_FULL(N_Pnts)
  i_s_s = CFT_SIGMA_START_FULL(N_Pnts)
  i_m_s = CFT_METRIC_START_FULL(N_Pnts)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call CFT_Diff_1_NONE (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                        Aux_Diff(i_s_s),                  &
                                        y(1,m_Start)   , dy(1,m_Start))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined (ESSL)
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined (VFFT) || defined (CRAY)
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    call CFT_Diff_1_NONE (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                        Aux_Diff(i_s_s),                  &
                                        y(1,m_Start)   , dy(1,m_Start))
  endif

  CONTAINS

    Subroutine CFT_Diff_1_NONE (N1, M, Metric, WSave, Sigma, y, dy)
      
    integer  :: N1, M

    REALTYPE, dimension(0:ldy-1,M)  :: y, dy
    REALTYPE, dimension(0:N1-1,10)  :: Metric
    REALTYPE, dimension(*)          :: WSave
    REALTYPE, dimension(0:N1-1)     :: Sigma

    integer  :: N, ldw, i, j, k

#if defined (ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(0:N1-1,M)   :: a, dd
#endif
#if defined (VFFT)
    REALTYPE, dimension(M,0:N1-1)   :: w1, w2, a
    REALTYPE, dimension(0:N1-1,M)   :: dd
#endif
#if defined (CRAY)
    REALTYPE, dimension(0:N1-1,M)   :: a, dd
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N1-1 

#if defined (ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size,  &
                                                    Work , wk_Size)

      if (Smooth_1 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        call CFT_1_Copy      (N, M,        dy, a )
        call CFT_1_Recursive (N, M,        a , dy)
      enddo
          
      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

                         call CFT_1_Copy      (N, M,        dy, a )
! FIRST derivative

                         call CFT_1_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 1) call CFT_1_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
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
     
! SECOND derivative

                         call CFT_1_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 2) call CFT_1_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
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

! THIRD derivative

                         call CFT_1_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 3) call CFT_1_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
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

! FOURTH derivative

                         call CFT_1_Recursive (N, M,        a , dy)
      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+3)
        enddo

    endif 

        call DeAllocate_WorkSpace (wk_Size) ; RETURN

#endif

#if defined (VFFT)

    ldw = SIZE(w1, DIM=1) 

    do j = 1,M
      w2(j,0:N) = y(0:N,j)
    enddo

    call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (Smooth_1 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        call CFT_2_Copy      (N, M,        w2, a )
        call CFT_2_Recursive (N, M,        a , w2)
      enddo

      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      do j = 1,M
        dy(0:N,j) = w2(j,0:N)
      enddo

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

                         call CFT_2_Copy      (N, M,        w2, a )
! FIRST derivative

                         call CFT_2_Recursive (N, M,        a , w2)
      if (O_Of_D   /= 1) call CFT_2_Copy      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             w2(j,0:N)*Metric(0:N,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) =             w2(j,0:N)*Metric(0:N,k)
        enddo
      endif

! SECOND Derivative

                         call CFT_2_Recursive (N, M,        a , w2)
      if (O_Of_D   /= 2) call CFT_2_Copy      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w2(j,0:N)*Metric(0:N,k+1)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w2(j,0:N)*Metric(0:N,k+1)
        enddo
      endif

! THIRD derivative

                         call CFT_2_Recursive (N, M,        a , w2)
      if (O_Of_D   /= 3) call CFT_2_Copy      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 3) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w2(j,0:N)*Metric(0:N,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w2(j,0:N)*Metric(0:N,k+2)
        enddo
      endif

! FOURTH Derivative

                         call CFT_2_Recursive (N, M,        a , w2)
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)
 
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w2(j,0:N)*Metric(0:N,k+3)
        enddo

    endif
#endif

#if defined (CRAY)

    call Cray_CFT_1 (N, M,  y,ldy, dy,ldy, WSave)

      if (Smooth_1 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        call CFT_1_Copy      (N, M,        dy, a )
        call CFT_1_Recursive (N, M,        a , dy)
      enddo

      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, WSave)

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

                         call CFT_1_Copy      (N, M,        dy, a )
! FIRST derivative

                         call CFT_1_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 1) call CFT_1_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, WSave)

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

! SECOND derivative

                         call CFT_1_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 2) call CFT_1_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, WSave)

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

! THIRD derivative

                         call CFT_1_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 3) call CFT_1_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, WSave)

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

! FOURTH derivative

                         call CFT_1_Recursive (N, M,        a , dy)
      if (Smooth_2 == 1) call CFT_1_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, WSave)

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+3)
        enddo

    endif
#endif

    END Subroutine CFT_Diff_1_NONE 

  END Subroutine C_Diff_CFT_1_NONE

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Diff_CFT_1_BOTH
! ----------------------------------------------------------------------
  Subroutine C_Diff_CFT_1_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer            :: N, M, ldy, O_Of_D, IEO

  REALTYPE, dimension(ldy,M) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: Map, Smooth_1, Smooth_2, BlockSize
  integer  :: N1, N_Pnts, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_w_s, i_q_s, i_s_s, i_m_s

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
   Smooth_1 = OPERATOR_SMOOTH_1(Aux_Diff)
   Smooth_2 = OPERATOR_SMOOTH_2(Aux_Diff)

  N1 = (N+1)/2 ; N_Pnts = N

  i_w_s = CFT_WSAVE_START_HALF(N_Pnts)
  i_q_s = CFT_QWAVE_START_HALF(N_Pnts)
  i_s_s = CFT_SIGMA_START_HALF(N_Pnts)
  i_m_s = CFT_METRIC_START_HALF(N_Pnts)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    if (IEO == 1) then
      call CFT_Diff_1_EVEN (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                          Aux_Diff(i_q_s), Aux_Diff(i_s_s), &
                                          y(1,m_Start)   , dy(1,m_Start))
    else
      call CFT_Diff_1_ODD  (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                          Aux_Diff(i_q_s), Aux_Diff(i_s_s), &
                                          y(1,m_Start)   , dy(1,m_Start))
    endif

  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined (ESSL)
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined (VFFT) || defined (CRAY)
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    if (IEO == 1) then
      call CFT_Diff_1_EVEN (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                          Aux_Diff(i_q_s), Aux_Diff(i_s_s), &
                                          y(1,m_Start)   , dy(1,m_Start))
    else
      call CFT_Diff_1_ODD  (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                          Aux_Diff(i_q_s), Aux_Diff(i_s_s), &
                                          y(1,m_Start)   , dy(1,m_Start))
    endif

  endif

  CONTAINS

    Subroutine CFT_Diff_1_EVEN (N1, M, Metric, WSave, QWave, Sigma, y, dy)
      
    integer  :: N1, M

    REALTYPE, dimension(0:N1-1,10)  :: Metric
    REALTYPE, dimension(*)          :: WSave, QWave
    REALTYPE, dimension(0:2*(N1-1)) :: Sigma
    REALTYPE, dimension(0:ldy-1,M)  :: y, dy

    integer  :: N, ldw, i, j, k

#if defined (ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(0:N1-1,M)   :: w1
    REALTYPE, dimension(0:N1-1,M)   :: a, dd
#endif
#if defined (VFFT)
    REALTYPE, dimension(M,0:N1-1)   :: w1, w2, a
    REALTYPE, dimension(0:N1-1,M)   :: dd
#endif
#if defined (CRAY)
    REALTYPE, dimension(0:N1-1,M)   :: w1
    REALTYPE, dimension(0:N1-1,M)   :: a, dd
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N1-1 

    ldw = SIZE(w1, DIM=1)

#if defined (ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size,  &
                                                    Work , wk_Size)
 
    call DeAllocate_WorkSpace (wk_Size)

        if (Smooth_1 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_1_Recursive_Even (N, M,        dy, w1)
        if (MOD(k,2) == 0) call CFT_1_Recursive_Odd  (N, M,        w1, dy)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

        call ESSL_QWAVE_BWRD_1 (1, N, M, w1, ldw, dy, ldy, QWave)
      else
        if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

        call   Allocate_WorkSpace (wk_Size)

        call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
                                                        Work , wk_Size)

        call DeAllocate_WorkSpace (wk_Size)
      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_1_Recursive_Even (N, M,        dy, w1)
      if (O_Of_D   /= 1) call CFT_1_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

      call ESSL_QWAVE_BWRD_1 (1, N, M, w1, ldw, dy, ldy, QWave)

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

! SECOND derivative

                         call CFT_1_Recursive_Odd  (N, M,        a , dy)
      if (O_Of_D   /= 2) call CFT_1_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

      call   Allocate_WorkSpace (wk_Size) 

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

      call DeAllocate_WorkSpace (wk_Size) 

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

! THIRD derivative

                         call CFT_1_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 3) call CFT_1_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

      call ESSL_QWAVE_BWRD_1 (1, N, M, w1, ldw, dy, ldy, QWave)

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

! FOURTH derivative

                         call CFT_1_Recursive_Odd  (N, M,        a , dy)
      if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

      call   Allocate_WorkSpace (wk_Size) 

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

      call DeAllocate_WorkSpace (wk_Size) 

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+3)
        enddo

    endif

#endif

#if defined (VFFT)

    do i = 0,N
      w1(:,i) = y(i,:)
    enddo

    call VCosT  (M, N+1, w1, w2, ldw, WSave)

        if (Smooth_1 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_2_Recursive_Even (N, M,        w1, w2)
        if (MOD(k,2) == 0) call CFT_2_Recursive_Odd  (N, M,        w2, w1)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w2)

        call VCosQB (M, N  , w2, w1, ldw, QWave)

        do j = 1,M
          dy(0:N-1,j) = w2(j,0:N-1)
        enddo

        dy(N,1:M) = ZERO

      else
        if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

        call VCosT  (M, N+1, w1, w2, ldw, WSave)

        do j = 1,M
          dy(0:N,j) = w1(j,0:N)
        enddo

      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_2_Recursive_Even (N, M,        w1, w2)
      if (O_Of_D   /= 1) call CFT_2_Copy_Odd       (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w2)

      call VCosQB (M, N  , w2, w1, ldw, QWave)

      w2(1:M,N) = ZERO

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             w2(j,0:N)*Metric(0:N,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) =             w2(j,0:N)*Metric(0:N,k)
        enddo
      endif

! SECOND derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , w1)
      if (O_Of_D   /= 2) call CFT_2_Copy_Even      (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

      call VCosT  (M, N+1, w1, w2, ldw, WSave)

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

! THIRD derivative

                         call CFT_2_Recursive_Even (N, M,        a , w2)
      if (O_Of_D   /= 3) call CFT_2_Copy_Odd       (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w2)

      call VCosQB (M, N  , w2, w1, ldw, QWave)

      w2(1:M,N) = ZERO

      if (O_Of_D == 3) then
       do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w2(j,0:N)*Metric(0:N,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w2(j,0:N)*Metric(0:N,k+2)
        enddo
      endif

! FOURTH derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , w1)
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

      call VCosT  (M, N+1, w1, w2, ldw, WSave)

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+3)
        enddo

    endif
#endif

#if defined (CRAY)

    call Cray_CFT_1 (N, M,  y,ldy, dy,ldy, WSave)

        if (Smooth_1 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_1_Recursive_Even (N, M,        dy, w1)
        if (MOD(k,2) == 0) call CFT_1_Recursive_Odd  (N, M,        w1, dy)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

        call Cray_QWAVE_BWRD_1 (N, M, w1,ldw, dy,ldy, QWave, WSave)
      else
        if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

        call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, WSave)
      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_1_Recursive_Even (N, M,        dy, w1)
      if (O_Of_D   /= 1) call CFT_1_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

      call Cray_QWAVE_BWRD_1 (N, M, w1,ldw, w1,ldw, QWave, WSave)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             w1(0:N,j)*Metric(0:N,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) =             w1(0:N,j)*Metric(0:N,k)
        enddo
      endif

! SECOND derivative

                         call CFT_1_Recursive_Odd  (N, M,        a , dy)
      if (O_Of_D   /= 2) call CFT_1_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, WSave)

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

! THIRD derivative

                         call CFT_1_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 3) call CFT_1_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

      call Cray_QWAVE_BWRD_1 (N, M, w1,ldw, w1,ldw, QWave, WSave)

      if (O_Of_D == 3) then
       do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(0:N,j)*Metric(0:N,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w1(0:N,j)*Metric(0:N,k+2)
        enddo
      endif

! FOURTH derivative

                         call CFT_1_Recursive_Odd  (N, M,        a , dy)
      if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, WSave)

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+3)
        enddo

    endif
#endif

    END Subroutine CFT_Diff_1_EVEN

! ----------------------------------------------------------------------
! FUNCTION NAME: CFT_Diff_1_ODD
! ----------------------------------------------------------------------
    Subroutine CFT_Diff_1_ODD  (N1, M, Metric, WSave, QWave, Sigma, y, dy)

    integer  :: N1, M

    REALTYPE, dimension(0:N1-1,10)  :: Metric
    REALTYPE, dimension(*)          :: WSave, QWave
    REALTYPE, dimension(0:2*(N1-1)) :: Sigma
    REALTYPE, dimension(0:ldy-1,M)  :: y, dy

    integer  :: N, ldw, i, j, k

#if defined (ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(0:N1-1,M)   :: w1
    REALTYPE, dimension(0:N1-1,M)   :: a, dd
#endif
#if defined (VFFT)
    REALTYPE, dimension(M,0:N1-1)   :: w1, w2, a
    REALTYPE, dimension(0:N1-1,M)   :: dd
#endif
#if defined (CRAY)
    REALTYPE, dimension(0:N1-1,M)   :: w1
    REALTYPE, dimension(0:N1-1,M)   :: a, dd
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N1-1 

    ldw = SIZE(w1, DIM=1)

#if defined (ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call ESSL_QWAVE_FWRD_1 (1, N, M,  y, ldy, w1, ldw, QWave)

        if (Smooth_1 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_1_Recursive_Odd  (N, M,        w1, dy)
        if (MOD(k,2) == 0) call CFT_1_Recursive_Even (N, M,        dy, w1)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

        call   Allocate_WorkSpace (wk_Size)

        call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
                                                        Work , wk_Size)

        call DeAllocate_WorkSpace (wk_Size)
      else
        if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

        call ESSL_QWAVE_BWRD_1 (1, N, M, w1, ldw, dy, ldy, QWave)
      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_1_Recursive_Odd  (N, M,        w1, dy)
      if (O_Of_D   /= 1) call CFT_1_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

      call   Allocate_WorkSpace (wk_Size)

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

      call DeAllocate_WorkSpace (wk_Size)

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

! SECOND derivative

                         call CFT_1_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 2) call CFT_1_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

      call ESSL_QWAVE_BWRD_1 (1, N, M, w1, ldw, dy, ldy, QWave)

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

! THIRD derivative

                         call CFT_1_Recursive_Odd  (N, M,        a , dy)
      if (O_Of_D   /= 3) call CFT_1_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

      call   Allocate_WorkSpace (wk_Size) 

      call COSF (0, dy,1,ldy, dy,1,ldy, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

      call DeAllocate_WorkSpace (wk_Size) 

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

! FOURTH derivative

                         call CFT_1_Recursive_Even (N, M,        a , w1)
      if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

      call ESSL_QWAVE_BWRD_1 (1, N, M, w1, ldw, dy, ldy, QWave)

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + dy(0:N,j)*Metric(0:N,k+3)
        enddo

    endif
#endif

#if defined (VFFT)

    do j = 1,M
      w1(j,0:N) = y(0:N,j)
    enddo

    call VCosQF (M, N  , w1, w2, ldw, QWave)

        if (Smooth_1 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_2_Recursive_Odd  (N, M,        w1, w2)
        if (MOD(k,2) == 0) call CFT_2_Recursive_Even (N, M,        w2, w1)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w2)

        call VCosT  (M, N+1, w2, w1, ldw, WSave)

        do j = 1,M
          dy(0:N,j) = w2(j,0:N)
        enddo

      else
        if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

        call VCosQB (M, N  , w1, w2, ldw, QWave)

        do j = 1,M
          dy(0:N-1,j) = w1(j,0:N)
        enddo

          dy(N,1:M) = ZERO
      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_2_Recursive_Odd  (N, M,        w1, w2)
      if (O_Of_D   /= 1) call CFT_2_Copy_Even      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 1) then
        do j = 1,M
          dy(0:N,j) =             w2(j,0:N)*Metric(0:N,k)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) =             w2(j,0:N)*Metric(0:N,k)
        enddo
      endif

! SECOND derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 2) call CFT_2_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call VCosQB (M, N  , w1, w2, ldw, QWave)

      w1(1:M,N) = ZERO

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

! THIRD derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , w2)
      if (O_Of_D   /= 3) call CFT_2_Copy_Even      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 3) then
       do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w2(j,0:N)*Metric(0:N,k+2)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w2(j,0:N)*Metric(0:N,k+2)
        enddo
      endif

! FOURTH derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call VCosQB (M, N  , w1, w2, ldw, QWave)

      w1(1:M,N) = ZERO

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(j,0:N)*Metric(0:N,k+3)
        enddo

    endif

#endif

#if defined (CRAY)

    call Cray_QWAVE_FWRD_1 (N, M,  y,ldy, w1,ldw, QWave, WSave)

        if (Smooth_1 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_1_Recursive_Odd  (N, M,        w1, dy)
        if (MOD(k,2) == 0) call CFT_1_Recursive_Even (N, M,        dy, w1)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

        call Cray_CFT_1 (N, M, dy,ldy, dy,ldy, WSave)
      else
        if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

        call Cray_QWAVE_BWRD_1 (N, M, w1,ldw, dy,ldy, QWave, WSave)
      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_1_Recursive_Odd  (N, M,        w1, dy)
      if (O_Of_D   /= 1) call CFT_1_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldw, WSave)

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

! SECOND derivative

                         call CFT_1_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 2) call CFT_1_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

      call Cray_QWAVE_BWRD_1 (N, M, w1,ldw, w1,ldw, QWave, WSave)

      if (O_Of_D == 2) then
        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(0:N,j)*Metric(0:N,k+1)
        enddo

        RETURN
      else
        do j = 1,M
          dd(0:N,j) = dd(0:N,j) + w1(0:N,j)*Metric(0:N,k+1)
        enddo
      endif

! THIRD derivative

                         call CFT_1_Recursive_Odd  (N, M,        a , dy)
      if (O_Of_D   /= 3) call CFT_1_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_1_Smooth_Even    (N, M, Sigma, dy)

      call Cray_CFT_1 (N, M, dy,ldy, dy,ldw, WSave)

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

! FOURTH derivative

                         call CFT_1_Recursive_Even (N, M,        a , w1)
      if (Smooth_2 == 1) call CFT_1_Smooth_Odd     (N, M, Sigma, w1)

      call Cray_QWAVE_BWRD_1 (N, M, w1,ldw, w1,ldw, QWave, WSave)

        do j = 1,M
          dy(0:N,j) = dd(0:N,j) + w1(0:N,j)*Metric(0:N,k+3)
        enddo

    endif
#endif

    END Subroutine CFT_Diff_1_ODD  

  END Subroutine C_Diff_CFT_1_BOTH

END MODULE Chebyshev_Differentiation_CFT_1


! ----------------------------------------------------------------------
! FUNCTION NAME: C_Diff_CFT_2
! DESCRIPTION  : does the chebyshev derivative on m vectors
!                using the transform-recursion algorithm.
!                This is a "high-level" Subroutine that does
!                unpacking and blocking.
! ----------------------------------------------------------------------

MODULE Chebyshev_Differentiation_CFT_2
  USE PS_WorkSpace_Allocation, Work=>WorkSpace
  USE CFT_Auxillary_Subroutine

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

PUBLIC  :: C_Diff_CFT_2

CONTAINS

  Subroutine C_Diff_CFT_2      (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer            :: N, M, ldy, O_Of_D

  integer , OPTIONAL :: IEO
  logical            :: L_IEO

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  L_IEO = PRESENT(IEO)

  if (L_IEO) then
    if (IEO /= 0) then
      call C_Diff_CFT_2_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)
    else
      call C_Diff_CFT_2_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
    endif
  else
      call C_Diff_CFT_2_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)
  endif

  END Subroutine C_Diff_CFT_2

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Diff_CFT_2_NONE
! ----------------------------------------------------------------------
  Subroutine C_Diff_CFT_2_NONE (N, M, O_Of_D, y, dy, ldy, Aux_Diff)

  integer            :: N, M, ldy, O_Of_D

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: Map, Smooth_1, Smooth_2, BlockSize
  integer  :: N1, N_Pnts, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_w_s, i_s_s, i_m_s

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
   Smooth_1 = OPERATOR_SMOOTH_1(Aux_Diff)
   Smooth_2 = OPERATOR_SMOOTH_2(Aux_Diff)

  N1 = N ; N_Pnts = N

  i_w_s = CFT_WSAVE_START_FULL(N_Pnts)
  i_s_s = CFT_SIGMA_START_FULL(N_Pnts)
  i_m_s = CFT_METRIC_START_FULL(N_Pnts)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    call CFT_Diff_2_NONE (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                        Aux_Diff(i_s_s),                  &
                                        y(m_Start,1)   , dy(m_Start,1))
  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined (ESSL)
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined (VFFT) || defined (CRAY)
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    call CFT_Diff_2_NONE (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                        Aux_Diff(i_s_s),                  &
                                        y(m_Start,1)   , dy(m_Start,1))
  endif

  CONTAINS

    Subroutine CFT_Diff_2_NONE (N1, M, Metric, WSave, Sigma, y, dy)

    integer  :: N1, M

    REALTYPE, dimension(0:N1-1,10)  :: Metric
    REALTYPE, dimension(*)          :: WSave
    REALTYPE, dimension(0:N1-1)     :: Sigma
    REALTYPE, dimension(ldy,0:N1-1) :: y, dy

    integer  :: N, ldw, i, k

#if defined (ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(M,0:N1-1)   :: a, dd
#endif
#if defined (VFFT)
    REALTYPE, dimension(M,0:N1-1)   :: w1, w2
    REALTYPE, dimension(M,0:N1-1)   :: a, dd
#endif
#if defined (CRAY)
    REALTYPE, dimension(M,0:N1-1)   :: a, dd
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N1-1

#if defined (ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 2)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size,  &
                                                    Work , wk_Size)
 
      if (Smooth_1 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        call CFT_2_Copy      (N, M,        dy, a )
        call CFT_2_Recursive (N, M,        a , dy)
      enddo

      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size,  &
                                                      Work , wk_Size)

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

                         call CFT_2_Copy      (N, M,        dy, a )
! FIRST derivative

                         call CFT_2_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 1) call CFT_2_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size,  &
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

! SECOND derivative

                         call CFT_2_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 2) call CFT_2_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size, &
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
     
! THIRD derivative

                         call CFT_2_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 3) call CFT_2_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size, &
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

! FOURTH derivative

                         call CFT_2_Recursive (N, M,        a , dy)
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+3)
        enddo

    endif

        call DeAllocate_WorkSpace (wk_Size) ; RETURN

#endif

#if defined (VFFT)

    ldw = SIZE(w1, DIM=1) 
      
    w2 = y(1:M,:)

    call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (Smooth_1 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        call CFT_2_Copy      (N, M,        w2, a )
        call CFT_2_Recursive (N, M,        a , w2)
      enddo

      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      dy(1:M,:) = w2

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

                         call CFT_2_Copy      (N, M,        w2, a )
! FIRST derivative

                         call CFT_2_Recursive (N, M,        a , w2)
      if (O_Of_D   /= 1) call CFT_2_Copy      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 1) then
        do i = 0,N
          dy(1:M,i) =             w2(1:M,i)*Metric(i,k)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) =             w2(1:M,i)*Metric(i,k)
        enddo
      endif

! SECOND Derivative

                         call CFT_2_Recursive (N, M,        a , w2)
      if (O_Of_D   /= 2) call CFT_2_Copy      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 2) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w2(1:M,i)*Metric(i,k+1)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + w2(1:M,i)*Metric(i,k+1)
        enddo
      endif

! THIRD derivative

                         call CFT_2_Recursive (N, M,        a , w2)
      if (O_Of_D   /= 3) call CFT_2_Copy      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 3) then
        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w2(1:M,i)*Metric(i,k+2)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + w2(1:M,i)*Metric(i,k+2)
        enddo
      endif

! FOURTH Derivative

                         call CFT_2_Recursive (N, M,        a , w2)
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w2(1:M,i)*Metric(i,k+3)
        enddo

    endif

#endif

#if defined (CRAY)
 
    call Cray_CFT_2 (N, M,  y,ldy, dy,ldy, WSave)

      if (Smooth_1 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        call CFT_2_Copy      (N, M,        dy, a )
        call CFT_2_Recursive (N, M,        a , dy)
      enddo

      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, WSave)

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

                         call CFT_2_Copy      (N, M,        dy, a )
! FIRST derivative

                         call CFT_2_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 1) call CFT_2_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, WSave)

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

! SECOND derivative

                         call CFT_2_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 2) call CFT_2_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, WSave)

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
     
! THIRD derivative

                         call CFT_2_Recursive (N, M,        a , dy)
      if (O_Of_D   /= 3) call CFT_2_Copy      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, WSave)

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

! FOURTH derivative

                         call CFT_2_Recursive (N, M,        a , dy)
      if (Smooth_2 == 1) call CFT_2_Smooth    (N, M, Sigma, dy)

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, WSave)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+3)
        enddo

    endif

#endif
 
    END Subroutine CFT_Diff_2_NONE 

  END Subroutine C_Diff_CFT_2_NONE

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Diff_CFT_2_BOTH
! ----------------------------------------------------------------------
  Subroutine C_Diff_CFT_2_BOTH (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

  integer            :: N, M, ldy, O_Of_D, IEO

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: Map, Smooth_1, Smooth_2, BlockSize
  integer  :: N1, N_Pnts, i
  integer  :: m_Block, m_Left, m_Vector, m_Start
  integer  :: i_w_s, i_q_s, i_s_s, i_m_s

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
   Smooth_1 = OPERATOR_SMOOTH_1(Aux_Diff)
   Smooth_2 = OPERATOR_SMOOTH_2(Aux_Diff)

  N1 = (N+1)/2 ; N_Pnts = N

  i_w_s = CFT_WSAVE_START_HALF(N_Pnts)
  i_q_s = CFT_QWAVE_START_HALF(N_Pnts)
  i_s_s = CFT_SIGMA_START_HALF(N_Pnts)
  i_m_s = CFT_METRIC_START_HALF(N_Pnts)

!$OMP PARALLEL DO
  do i = 0, m_Block-1
    m_Vector =   BlockSize
    m_Start  = i*BlockSize+1

    if (IEO == 1) then
      call CFT_Diff_2_EVEN (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                          Aux_Diff(i_q_s), Aux_Diff(i_s_s), &
                                          y(m_Start,1)   , dy(m_Start,1))
    else
      call CFT_Diff_2_ODD  (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                          Aux_Diff(i_q_s), Aux_Diff(i_s_s), &
                                          y(m_Start,1)   , dy(m_Start,1))
    endif

  enddo
!$OMP END PARALLEL DO

  if (m_Left >  0) then

#if defined (ESSL)
    m_Vector =   BlockSize
    m_Start  = M-BlockSize+1
#endif

#if defined (VFFT) || defined (CRAY)
           i =   m_Block
    m_Vector =   m_Left
    m_Start  = i*BlockSize+1
#endif

    if (IEO == 1) then
      call CFT_Diff_2_EVEN (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                          Aux_Diff(i_q_s), Aux_Diff(i_s_s), &
                                          y(m_Start,1)   , dy(m_Start,1))
    else
      call CFT_Diff_2_ODD  (N1, m_Vector, Aux_Diff(i_m_s), Aux_Diff(i_w_s), &
                                          Aux_Diff(i_q_s), Aux_Diff(i_s_s), &
                                          y(m_Start,1)   , dy(m_Start,1))
    endif

  endif

  CONTAINS

! ----------------------------------------------------------------------
! FUNCTION NAME: CFT_Diff_2_EVEN
! ----------------------------------------------------------------------
    Subroutine CFT_Diff_2_EVEN (N1, M, Metric, WSave, QWave, Sigma, y, dy)

    integer  :: N1, M

    REALTYPE, dimension(0:N1-1,10)  :: Metric
    REALTYPE, dimension(*)          :: WSave, QWave
    REALTYPE, dimension(0:2*(N1-1)) :: Sigma
    REALTYPE, dimension(ldy,0:N1-1) :: y, dy

    integer  :: N, ldw, i, j, k

#if defined (ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(M,0:N1-1)   :: w1
    REALTYPE, dimension(M,0:N1-1)   :: a, dd
#endif
#if defined (VFFT)
    REALTYPE, dimension(M,0:N1-1)   :: w1, w2
    REALTYPE, dimension(M,0:N1-1)   :: a, dd
#endif
#if defined (CRAY)
    REALTYPE, dimension(M,0:N1-1)   :: w1
    REALTYPE, dimension(M,0:N1-1)   :: a, dd
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N1-1 

    ldw = SIZE(w1, DIM=1)

#if defined (ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 2)

    call   Allocate_WorkSpace (wk_Size)

    call COSF (0,  y,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size,  &
                                                    Work , wk_Size)

    call DeAllocate_WorkSpace (wk_Size) 

        if (Smooth_1 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_2_Recursive_Even (N, M,        dy, w1)
        if (MOD(k,2) == 0) call CFT_2_Recursive_Odd  (N, M,        w1, dy)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

        call ESSL_QWAVE_BWRD_2 (1, N, M, w1, ldw, dy, ldy, QWave)
      else
        if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

        call   Allocate_WorkSpace (wk_Size)

        call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                        Work , wk_Size)

        call DeAllocate_WorkSpace (wk_Size) 
      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_2_Recursive_Even (N, M,        dy, w1)
      if (O_Of_D   /= 1) call CFT_2_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call ESSL_QWAVE_BWRD_2 (1, N, M, w1, ldw, dy, ldy, QWave)

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

! SECOND derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , dy)
      if (O_Of_D   /= 2) call CFT_2_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

      call   Allocate_WorkSpace (wk_Size)

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

      call DeAllocate_WorkSpace (wk_Size) 

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

! THIRD derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 3) call CFT_2_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call ESSL_QWAVE_BWRD_2 (1, N, M, w1, ldw, dy, ldy, QWave)

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

! FOURTH derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , dy)
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

      call   Allocate_WorkSpace (wk_Size) 

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

      call DeAllocate_WorkSpace (wk_Size) 

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+3)
        enddo

    endif

#endif

#if defined (VFFT)

    w1 = y(1:M,:)

    call VCosT  (M, N+1, w1, w2, ldw, WSave)

        if (Smooth_1 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_2_Recursive_Even (N, M,        w1, w2)
        if (MOD(k,2) == 0) call CFT_2_Recursive_Odd  (N, M,        w2, w1)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w2)

        call VCosQB (M, N  , w2, w1, ldw, QWave)

        do i = 0,N-1
          dy(1:M,i) = w2(1:M,i)
        enddo

        dy(1:M,N) = ZERO

      else
        if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

        call VCosT  (M, N+1, w1, w2, ldw, WSave)

        do i = 0,N
          dy(1:M,i) = w1(1:M,i)
        enddo

      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_2_Recursive_Even (N, M,        w1, w2)
      if (O_Of_D   /= 1) call CFT_2_Copy_Odd       (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w2)

      call VCosQB (M, N  , w2, w1, ldw, QWave)

      w2(1:M,N) = ZERO

      if (O_Of_D == 1) then
        do i = 0,N
          dy(1:M,i) =             w2(1:M,i)*Metric(i,k)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) =             w2(1:M,i)*Metric(i,k)
        enddo
      endif

! SECOND derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , w1)
      if (O_Of_D   /= 2) call CFT_2_Copy_Even      (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

      call VCosT  (M, N+1, w1, w2, ldw, WSave)

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

! THIRD derivative

                         call CFT_2_Recursive_Even (N, M,        a , w2)
      if (O_Of_D   /= 3) call CFT_2_Copy_Odd       (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w2)

      call VCosQB (M, N  , w2, w1, ldw, QWave)

      w2(1:M,N) = ZERO

      if (O_Of_D == 3) then
       do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w2(1:M,i)*Metric(i,k+2)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + w2(1:M,i)*Metric(i,k+2)
        enddo
      endif

! FOURTH derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , w1)
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

      call VCosT  (M, N+1, w1, w2, ldw, WSave)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+3)
        enddo

    endif
#endif

#if defined (CRAY)

    call Cray_CFT_2 (N, M,  y,ldy, dy,ldy, WSave)

        if (Smooth_1 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_2_Recursive_Even (N, M,        dy, w1)
        if (MOD(k,2) == 0) call CFT_2_Recursive_Odd  (N, M,        w1, dy)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

        call Cray_QWAVE_BWRD_2 (N, M, w1,ldw, dy,ldy, QWave, WSave)
      else
        if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

        call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, WSave)
      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_2_Recursive_Even (N, M,        dy, w1)
      if (O_Of_D   /= 1) call CFT_2_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call Cray_QWAVE_BWRD_2 (N, M, w1,ldw, w1,ldw, QWave, WSave)

      dy(1:M,N) = ZERO

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

! SECOND derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , w1)
      if (O_Of_D   /= 2) call CFT_2_Copy_Even      (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

      call Cray_CFT_2 (N, M, w1,ldw, w1,ldw, WSave)

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

! THIRD derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 3) call CFT_2_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call Cray_QWAVE_BWRD_2 (N, M, w1,ldw, w1,ldw, QWave, WSave)

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

! FOURTH derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , w1)
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w1)

      call Cray_CFT_2 (N, M, w1,ldw, w1,ldw, WSave)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+3)
        enddo

    endif
#endif
 
    END Subroutine CFT_Diff_2_EVEN 

! ----------------------------------------------------------------------
! FUNCTION NAME: CFT_Diff_2_ODD 
! ----------------------------------------------------------------------
    Subroutine CFT_Diff_2_ODD  (N1, M, Metric, WSave, QWave, Sigma, y, dy)

    integer  :: N1, M

    REALTYPE, dimension(0:N1-1,10)  :: Metric
    REALTYPE, dimension(*)          :: WSave, QWave
    REALTYPE, dimension(0:2*(N1-1)) :: Sigma
    REALTYPE, dimension(ldy,0:N1-1) :: y, dy

    integer  :: N, ldw, i, j, k

#if defined (ESSL)
    integer  :: N2
    integer  :: wc_Size, wk_Size
    REALTYPE :: Scale

    REALTYPE, dimension(ldy,0:N1-1) :: w1
    REALTYPE, dimension(M,0:N1-1)   :: a, dd
#endif
#if defined (VFFT)
    REALTYPE, dimension(M,0:N1-1)   :: w1, w2
    REALTYPE, dimension(M,0:N1-1)   :: a, dd
#endif
#if defined (CRAY)
    REALTYPE, dimension(M,0:N1-1)   :: w1
    REALTYPE, dimension(M,0:N1-1)   :: a, dd
#endif

    if (Map /= 0) then
      if (O_Of_D == 1) k = 1
      if (O_Of_D == 2) k = 2
      if (O_Of_D == 3) k = 4
      if (O_Of_D == 4) k = 7
    endif

    N = N1-1 

    ldw = SIZE(w1, DIM=1)

#if defined (ESSL)

    N2 = 2*N ; Scale = SQRT(TWO/N)

    wc_Size = ESSL_CFT_Twiddle_Factor (N2)
    wk_Size = ESSL_CFT_Work           (N2, M, 1)

    call ESSL_QWAVE_FWRD_2 (1, N, M,  y, ldy, w1, ldw, QWave)

        if (Smooth_1 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_2_Recursive_Odd  (N, M,        w1, dy)
        if (MOD(k,2) == 0) call CFT_2_Recursive_Even (N, M,        dy, w1)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

        call   Allocate_WorkSpace (wk_Size)

        call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                        Work , wk_Size)

        call DeAllocate_WorkSpace (wk_Size) 
      else
        if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

        call ESSL_QWAVE_BWRD_2 (1, N, M, w1, ldw, dy, ldy, QWave)
      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_2_Recursive_Odd  (N, M,        w1, dy)
      if (O_Of_D   /= 1) call CFT_2_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

      call   Allocate_WorkSpace (wk_Size)

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

      call DeAllocate_WorkSpace (wk_Size) 

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

! SECOND derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 2) call CFT_2_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call ESSL_QWAVE_BWRD_2 (1, N, M, w1, ldw, dy, ldy, QWave)

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

! THIRD derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , dy)
      if (O_Of_D   /= 3) call CFT_2_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

      call   Allocate_WorkSpace (wk_Size)

      call COSF (0, dy,ldy,1, dy,ldy,1, N2, M, Scale, WSave, wc_Size, &
                                                      Work , wk_Size)

      call DeAllocate_WorkSpace (wk_Size) 

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

! FOURTH derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call ESSL_QWAVE_BWRD_2 (1, N, M, w1, ldw, dy, ldy, QWave)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + dy(1:M,i)*Metric(i,k+3)
        enddo

    endif
#endif

#if defined (VFFT)

    w1 = y(1:M,:)

    call VCosQF (M, N  , w1, w2, ldw, QWave)

        if (Smooth_1 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_2_Recursive_Odd  (N, M,        w1, w2)
        if (MOD(k,2) == 0) call CFT_2_Recursive_Even (N, M,        w2, w1)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w2)

        call VCosT  (M, N+1, w2, w1, ldw, WSave)

        do i = 0,N
          dy(1:M,i) = w2(1:M,i)
        enddo

      else
        if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

        call VCosQB (M, N  , w1, w2, ldw, QWave)

        do i = 0,N
          dy(1:M,i) = w1(1:M,i)
        enddo

      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_2_Recursive_Odd  (N, M,        w1, w2)
      if (O_Of_D   /= 1) call CFT_2_Copy_Even      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 1) then
        do i = 0,N
          dy(1:M,i) =             w2(1:M,i)*Metric(i,k)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) =             w2(1:M,i)*Metric(i,k)
        enddo
      endif

! SECOND derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 2) call CFT_2_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call VCosQB (M, N  , w1, w2, ldw, QWave)

      w1(1:M,N) = ZERO

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

! THIRD derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , w2)
      if (O_Of_D   /= 3) call CFT_2_Copy_Even      (N, M,        w2, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, w2)

      call VCosT  (M, N+1, w2, w1, ldw, WSave)

      if (O_Of_D == 3) then
       do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w2(1:M,i)*Metric(i,k+2)
        enddo

        RETURN
      else
        do i = 0,N
          dd(1:M,i) = dd(1:M,i) + w2(1:M,i)*Metric(i,k+2)
        enddo
      endif

! FOURTH derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call VCosQB (M, N  , w1, w2, ldw, QWave)

      w1(1:M,N) = ZERO

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+3)
        enddo

    endif

#endif

#if defined (CRAY)

    call Cray_QWAVE_FWRD_2 (N, M,  y,ldy, w1,ldw, QWave, WSave)

        if (Smooth_1 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

! -----------------------------------------------------------
! NO MAPPING CASE
! -----------------------------------------------------------
    if (Map == 0) then

      do k = 1,O_Of_D
        if (MOD(k,2) == 1) call CFT_2_Recursive_Odd  (N, M,        w1, dy)
        if (MOD(k,2) == 0) call CFT_2_Recursive_Even (N, M,        dy, w1)
      enddo

      if (MOD(O_Of_D,2) == 1) then
        if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

        call Cray_CFT_2 (N, M, dy,ldy, dy,ldy, WSave)
      else
        if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

        call Cray_QWAVE_BWRD_2 (N, M, w1,ldw, dy,ldy, QWave, WSave)
      endif

    endif

! -----------------------------------------------------------
! Here starts the case WITH MAPPING
! -----------------------------------------------------------
    if (Map /= 0) then

! FIRST derivative

                         call CFT_2_Recursive_Odd  (N, M,        w1, dy)
      if (O_Of_D   /= 1) call CFT_2_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldw, WSave)

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

! SECOND derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (O_Of_D   /= 2) call CFT_2_Copy_Odd       (N, M,        w1, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call Cray_QWAVE_BWRD_2 (N, M, w1,ldw, w1,ldw, QWave, WSave)

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

! THIRD derivative

                         call CFT_2_Recursive_Odd  (N, M,        a , dy)
      if (O_Of_D   /= 3) call CFT_2_Copy_Even      (N, M,        dy, a )
      if (Smooth_2 == 1) call CFT_2_Smooth_Even    (N, M, Sigma, dy)

      call Cray_CFT_2 (N, M, dy,ldy, dy,ldw, WSave)

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

! FOURTH derivative

                         call CFT_2_Recursive_Even (N, M,        a , w1)
      if (Smooth_2 == 1) call CFT_2_Smooth_Odd     (N, M, Sigma, w1)

      call Cray_QWAVE_BWRD_2 (N, M, w1,ldw, w1,ldw, QWave, WSave)

        do i = 0,N
          dy(1:M,i) = dd(1:M,i) + w1(1:M,i)*Metric(i,k+3)
        enddo

    endif
#endif
 
    END Subroutine CFT_Diff_2_ODD  

  END Subroutine C_Diff_CFT_2_BOTH

END MODULE Chebyshev_Differentiation_CFT_2
