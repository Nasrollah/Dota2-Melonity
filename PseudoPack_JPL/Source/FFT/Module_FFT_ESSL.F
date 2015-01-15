#if defined (ESSL)
  
#include "pseudopack.h"
  
! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! ----------------------------------------------------------------------
!
!  N+1 = number of gridpoints in a vector
!  M   = number of vectors to be transformed
!  A   = array of vectors to be transformed
!  B   = destination array
!
!  Note :
!   N         must be an EVEN number.
!   A, B      can be the same array.
!   SinCos    is filled out by ESSL_QVE_Init with vector lengths N.
!   C         work array of size (N+2)*M
!
! Date    : 8-10-99
! Author  : Wai Sun Don
! Comment : Fortran 90 ized.
!
! ----------------------------------------------------------------------
  
  
MODULE ESSL_FFT
  USE Memory_Allocation_Status

implicit NONE
  
PRIVATE

PUBLIC  :: ESSL_QWAVE_FWRD_1 
PUBLIC  :: ESSL_QWAVE_BWRD_1 
PUBLIC  :: ESSL_QWAVE_FWRD_2 
PUBLIC  :: ESSL_QWAVE_BWRD_2 

PUBLIC  :: ESSL_FFT_Init 
PUBLIC  :: ESSL_COS_Init 
PUBLIC  :: ESSL_SIN_Init 
PUBLIC  :: ESSL_CFT_Init 
PUBLIC  :: ESSL_QWV_Init 

PUBLIC  :: ESSL_FFT_Twiddle_Factor
PUBLIC  :: ESSL_CFT_Twiddle_Factor
PUBLIC  :: ESSL_SFT_Twiddle_Factor

PUBLIC  :: ESSL_FFT_Work
PUBLIC  :: ESSL_CFT_Work

PUBLIC  :: ESSL_Adjust_N 

CONTAINS
  
! ----------------------------------------------------------------------
!                  <--  Along the First Index  -->
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_QWAVE_FWRD_1
! DESCRIPTION  : Forward Quarter-Wave Cosine transform of an Odd sequenc
! ----------------------------------------------------------------------

  Subroutine ESSL_QWAVE_FWRD_1 (Init, N, M, A, LD_A, B, LD_B, SinCos)
  
  integer  :: Init, N, M, LD_A, LD_B
  
  REALTYPE, dimension(0:N)      :: SinCos
  REALTYPE, dimension(1)        :: w3
  
  REALTYPE, dimension(0:LD_A-1,M) :: A
  REALTYPE, dimension(0:LD_B-1,M) :: B
  REALTYPE, dimension(0:N+1,M)    :: C
  
  integer  :: wf_Size, wk_Size
  integer  :: N2, ldB, ldC, ISIGN, j, k
  REALTYPE :: Scale, A_Even, A_Odd
  
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  REALTYPE, dimension(:), ALLOCATABLE :: WSave
  
#ifdef DEBUG
  ASSERT(MOD(N,2)==0,"E_QWAVE_F_1:N2")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2  = N/2  
  ldB = N+2  
  ldC = N2+1  
  
  wf_Size = ESSL_FFT_Twiddle_Factor (N)
  wk_Size = ESSL_FFT_Work           (N)
  
#if 0
  call   Allocate_WSave (WSave, wf_Size)
  call   Allocate_Work  (Work , wk_Size)
#else
! Allocate (WSave(wf_Size), Work(wk_Size))

  if (.NOT. ALLOCATED(WSave)) then
    ALLOCATE (WSave(wf_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wf_Size
    Subroutine_Name  = '  Allocate_WSave'

    call Memory_Error_Check (0)
#endif
  endif

  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(wk_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

#endif
  
  ISIGN =  1 ; Scale = SQRT(HALF/N)  
  
  if (Init == 1) then  
    call RCFT (1, C,ldB, C,ldC, N, M, ISIGN, Scale, WSave, wf_Size, &
                  Work , wk_Size, w3, 0)
  endif  
  
  do j = 1,M  
    do k = 1,N2-1  
        A_Even = A(k,j)+A(N-k,j)  
        A_Odd  = A(k,j)-A(N-k,j)  
  
      C(k  ,j) = SinCos(N-k)*A_Odd +SinCos(k)*A_Even  
      C(N-k,j) = SinCos(N-k)*A_Even-SinCos(k)*A_Odd  
    enddo  
  enddo  
  
  k = N2  
  do j = 1,M  
    C(k,j) = TWO*SinCos(k)*A(k,j)  
    C(0,j) = A(0,j)  
    C(N,j) = A(0,j)  
  enddo  
  
    call RCFT (0, C,ldB, C,ldC, N, M, ISIGN, Scale, WSave, wf_Size,  &
                  Work , wk_Size, w3, 0)
  
#if 0
  call DeAllocate_WSave (WSave, wf_Size)
  call DeAllocate_Work  (Work , wk_Size)
#else
! DeAllocate (Work, WSave)

  if (ALLOCATED(WSave)) then
    DEALLOCATE (WSave, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wf_Size
    Subroutine_Name  = 'DeAllocate_WSave'

    call Memory_Error_Check (1)
#endif
  endif

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = 'DeAllocate_Work'

    call Memory_Error_Check (1)
#endif
  endif

#endif
  
    B(0,:) = C(0,:)  
  
  do j = 1,M  
    do k = 2,N-1,2  
      B(k-1,j) = C(k,j)-C(k+1,j)  
      B(k  ,j) = C(k,j)+C(k+1,j)  
    enddo  
  enddo  
  
    B(N-1,:) = C(N,:)  
  
  END Subroutine ESSL_QWAVE_FWRD_1
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_QWAVE_BWRD_1
! DESCRIPTION  : Backward Quarter-Wave Cosine transform of an Odd sequen
! ----------------------------------------------------------------------
  
  Subroutine ESSL_QWAVE_BWRD_1 (Init, N, M, A, LD_A, B, LD_B, SinCos)
  
  integer  :: Init, N, M, LD_A, LD_B
  
  REALTYPE, dimension(0:N)      :: SinCos
  REALTYPE, dimension(1)        :: w3
  
  REALTYPE, dimension(0:LD_A-1,M) :: A
  REALTYPE, dimension(0:LD_B-1,M) :: B
  REALTYPE, dimension(0:N+1,M)    :: C
  
  integer  :: wf_Size, wk_Size
  integer  :: N2, ldB, ldC, ISIGN, j, k
  REALTYPE :: Scale, C_Even, C_Odd
  
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  REALTYPE, dimension(:), ALLOCATABLE :: WSave
  
#ifdef DEBUG
  ASSERT(MOD(N,2)==0,"E_QWAVE_B_1:N2")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2  = N/2  
  ldB = N2+1  
  ldC = N+2  
  
  wf_Size = ESSL_FFT_Twiddle_Factor (N)
  wk_Size = ESSL_FFT_Work           (N)
  
#if 0
  call   Allocate_WSave (WSave, wf_Size)
  call   Allocate_Work  (Work , wk_Size)
#else
! Allocate (WSave(wf_Size), Work(wk_Size))

  if (.NOT. ALLOCATED(WSave)) then
    ALLOCATE (WSave(wf_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wf_Size
    Subroutine_Name  = '  Allocate_WSave'

    call Memory_Error_Check (0)
#endif
  endif

  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(wk_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

#endif
  
  ISIGN = -1 ; Scale = HALF*SQRT(HALF/N)  
  
  if (Init == 1) then  
    call CRFT (1, C,ldB, C,ldC, N, M, ISIGN, Scale, WSave, wf_Size, &
                  Work , wk_Size, w3, 0)
  endif  
  
  do j = 1,M  
    C(0,j) = TWO*A(0,j)  
    C(1,j) = ZERO  
    do k = 2,N-1,2  
      C(k  ,j) = A(k-1,j)+A(k,j)  
      C(k+1,j) = A(k-1,j)-A(k,j)  
    enddo  
    C(N  ,j) = TWO*A(N-1,j)  
    C(N+1,j) = ZERO  
  enddo  
  
    call CRFT (0, C,ldB, C,ldC, N, M, ISIGN, Scale, WSave, wf_Size, &
                  Work , wk_Size, w3, 0)
  
#if 0
  call DeAllocate_Work  (Work , wk_Size)
  call DeAllocate_WSave (WSave, wf_Size)
#else
! DeAllocate (Work, WSave)

  if (ALLOCATED(WSave)) then
    DEALLOCATE (WSave, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wf_Size
    Subroutine_Name  = 'DeAllocate_WSave'

    call Memory_Error_Check (1)
#endif
  endif

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = 'DeAllocate_Work'

    call Memory_Error_Check (1)
#endif
  endif

#endif
  
    B(0,:) = TWO*C(0,:)  
  
  do j = 1,M  
    do k = 1,N2-1  
        C_Even = C(k,j)+C(N-k,j)  
        C_Odd  = C(k,j)-C(N-k,j)  
  
      B(k  ,j) = SinCos(N-k)*C_Even-SinCos(k)*C_Odd  
      B(N-k,j) = SinCos(N-k)*C_Odd +SinCos(k)*C_Even  
    enddo  
  enddo  
  
    B(N2,:) = TWO*SinCos(N2)*C(N2,:)  
    B(N, :) = ZERO  
  
  END Subroutine ESSL_QWAVE_BWRD_1
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_QWAVE_FWRD_2
! DESCRIPTION  : Forward Quarter-Wave Cosine transform of an Odd sequenc
! ----------------------------------------------------------------------
  
  Subroutine ESSL_QWAVE_FWRD_2 (Init, N, M, A, LD_A, B, LD_B, SinCos)
  
  integer  :: Init, N, M, LD_A, LD_B
  
  REALTYPE, dimension(0:N)      :: SinCos
  REALTYPE, dimension(1)        :: w3
  
  REALTYPE, dimension(LD_A,0:N) :: A
  REALTYPE, dimension(LD_B,0:N) :: B
  REALTYPE, dimension(0:N+1,M)  :: C
  
  integer  :: wf_Size, wk_Size
  integer  :: N2, ldB, ldC, ISIGN, j, k
  REALTYPE :: Scale, A_Even, A_Odd
  
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  REALTYPE, dimension(:), ALLOCATABLE :: WSave
  
#ifdef DEBUG
  ASSERT(MOD(N,2)==0,"E_QWAVE_F_2:N2")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2  = N/2  
  ldB = N+2  
  ldC = N2+1  
  
  wf_Size = ESSL_FFT_Twiddle_Factor (N)
  wk_Size = ESSL_FFT_Work           (N)
  
#if 0
  call   Allocate_WSave (WSave, wf_Size)
  call   Allocate_Work  (Work , wk_Size)
#else
! Allocate (WSave(wf_Size), Work(wk_Size))

  if (.NOT. ALLOCATED(WSave)) then
    ALLOCATE (WSave(wf_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wf_Size
    Subroutine_Name  = '  Allocate_WSave'

    call Memory_Error_Check (0)
#endif
  endif

  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(wk_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

#endif
  
  ISIGN =  1 ; Scale = SQRT(HALF/N)  
  
  if (Init == 1) then  
    call RCFT (1, C,ldB, C,ldC, N, M, ISIGN, Scale, WSave, wf_Size, &
                  Work , wk_Size, w3, 0)
  endif  
  
  do j = 1,M  
    do k = 1,N2-1  
        A_Even = A(j,k)+A(j,N-k)  
        A_Odd  = A(j,k)-A(j,N-k)  
  
      C(k  ,j) = SinCos(N-k)*A_Odd +SinCos(k)*A_Even  
      C(N-k,j) = SinCos(N-k)*A_Even-SinCos(k)*A_Odd  
    enddo  
  enddo  
  
  k = N2  
  do j = 1,M  
    C(k,j) = TWO*SinCos(k)*A(j,k)  
    C(0,j) = A(j,0)  
    C(N,j) = A(j,0)  
  enddo  
  
    call RCFT (0, C,ldB, C,ldC, N, M, ISIGN, Scale, WSave, wf_Size, &
                  Work , wk_Size, w3, 0)
  
#if 0
  call DeAllocate_Work  (Work , wk_Size)
  call DeAllocate_WSave (WSave, wf_Size)
#else
! DeAllocate (Work, WSave)

  if (ALLOCATED(WSave)) then
    DEALLOCATE (WSave, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wf_Size
    Subroutine_Name  = 'DeAllocate_WSave'

    call Memory_Error_Check (1)
#endif
  endif

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = 'DeAllocate_Work'

    call Memory_Error_Check (1)
#endif
  endif

#endif
  
    B(1:M,0) = C(0,:)  
  
  do j = 1,M  
    do k = 2,N-1,2  
      B(j,k-1) = C(k,j)-C(k+1,j)  
      B(j,k  ) = C(k,j)+C(k+1,j)  
    enddo  
  enddo  
  
    B(1:M,N-1) = C(N,:)  
  
  END Subroutine ESSL_QWAVE_FWRD_2
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_QWAVE_BWRD_2
! DESCRIPTION  : Backward Quarter-Wave Cosine transform of an Odd sequen
! ----------------------------------------------------------------------
  
  Subroutine ESSL_QWAVE_BWRD_2 (Init, N, M, A, LD_A, B, LD_B, SinCos)
  
  integer  :: Init, N, M, LD_A, LD_B
  
  REALTYPE, dimension(0:N)      :: SinCos
  REALTYPE, dimension(1)        :: w3
  
  REALTYPE, dimension(LD_A,0:N) :: A
  REALTYPE, dimension(LD_B,0:N) :: B
  REALTYPE, dimension(0:N+1,M)  :: C
  
  integer  :: wf_Size, wk_Size
  integer  :: N2, ldB, ldC, ISIGN, j, k
  REALTYPE :: Scale, C_Even, C_Odd
  
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  REALTYPE, dimension(:), ALLOCATABLE :: WSave
  
#ifdef DEBUG
  ASSERT(MOD(N,2)==0,"E_QWAVE_B_2:N2")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2  = N/2  
  ldB = N2+1  
  ldC = N+2  
  
  wf_Size = ESSL_FFT_Twiddle_Factor (N)
  wk_Size = ESSL_FFT_Work           (N)
  
#if 0
  call   Allocate_WSave (WSave, wf_Size)
  call   Allocate_Work  (Work , wk_Size)
#else
! Allocate (WSave(wf_Size), Work(wk_Size))

  if (.NOT. ALLOCATED(WSave)) then
    ALLOCATE (WSave(wf_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wf_Size
    Subroutine_Name  = '  Allocate_WSave'

    call Memory_Error_Check (0)
#endif
  endif

  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(wk_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

#endif
  
  ISIGN = -1 ; Scale = HALF*SQRT(HALF/N)  
  
  if (Init == 1) then  
    call CRFT (1, C,ldB, C,ldC, N, M, ISIGN, Scale, WSave, wf_Size, &
                  Work , wk_Size, w3, 0)
  endif  
  
  do j = 1,M  
    C(0,j) = TWO*A(j,0)  
    C(1,j) = ZERO  
    do k = 2,N-1,2  
      C(k  ,j) = A(j,k-1)+A(j,k)  
      C(k+1,j) = A(j,k-1)-A(j,k)  
    enddo  
    C(N,j) = TWO*A(j,N-1)  
  enddo  
  
    call CRFT (0, C,ldB, C,ldC, N, M, ISIGN, Scale, WSave, wf_Size, &
                  Work , wk_Size, w3, 0)
  
#if 0
  call DeAllocate_Work  (Work , wk_Size)
  call DeAllocate_WSave (WSave, wf_Size)
#else
! DeAllocate (Work, WSave)

  if (ALLOCATED(WSave)) then
    DEALLOCATE (WSave, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wf_Size
    Subroutine_Name  = 'DeAllocate_WSave'

    call Memory_Error_Check (1)
#endif
  endif

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = 'DeAllocate_Work'

    call Memory_Error_Check (1)
#endif
  endif

#endif
  
    B(1:M,0) = TWO*C(0,:)  
  
  do j = 1,M  
    do k = 1,N2-1  
        C_Even = C(k,j)+C(N-k,j)  
        C_Odd  = C(k,j)-C(N-k,j)  
  
      B(j,k  ) = SinCos(N-k)*C_Even-SinCos(k)*C_Odd  
      B(j,N-k) = SinCos(N-k)*C_Odd +SinCos(k)*C_Even  
    enddo  
  enddo  
  
    B(1:M,N2) = TWO*SinCos(N2)*C(N2,:)  
    B(1:M,N ) = ZERO  
  
  END Subroutine ESSL_QWAVE_BWRD_2
  
! ----------------------------------------------------------------------
!                  <--  Initiailization Routines  -->
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_FFT_Init
! DESCRIPTION  : fills out the WSave structure for the FFT Transform
! ----------------------------------------------------------------------
  
  Subroutine ESSL_FFT_Init (N0, ldy, BlockSize, Index, WSave, ISIGN)
  
  integer  :: N0, ldy, BlockSize, Index, ISIGN  
  
  REALTYPE, dimension(1) :: y, w3
  REALTYPE, dimension(*) :: WSave
  
  integer  :: wf_Size, wk_Size
  integer  :: N2, M, N3, ldw
  REALTYPE :: Scale
  
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  
  N2 = N0  ; M = BlockSize ; Scale = SQRT(ONE/N0)
  N3 = N0/2+1  
  
  wf_Size = ESSL_FFT_Twiddle_Factor (N2)
  wk_Size = ESSL_FFT_Work           (N2)
  
#if 0
  call   Allocate_Work  (Work , wk_Size)
#else
! Allocate (Work(wk_Size))

  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(wk_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

#endif
  
  if (Index == 1) then  
    if (ISIGN ==  1) call RCFT (1, y, ldy, y, N3, N2, M, ISIGN, Scale, &
                                   WSave, wf_Size, Work , wk_Size, w3, 0)
    if (ISIGN == -1) call CRFT (1, y, N3, y, ldy, N2, M, ISIGN, Scale, &
                                   WSave, wf_Size, Work , wk_Size, w3, 0)
  endif  
  
  if (Index == 2) then  
    ldw = N0+2  
  
    if (ISIGN ==  1) call RCFT (1, y, ldw, y, N3, N2, M, ISIGN, Scale, &
                                   WSave, wf_Size, Work , wk_Size, w3, 0)
    if (ISIGN == -1) call CRFT (1, y, N3, y, ldw, N2, M, ISIGN, Scale, &
                                   WSave, wf_Size, Work , wk_Size, w3, 0)
  endif  
  
#if 0
  call DeAllocate_Work  (Work , wk_Size)
#else
! DeAllocate (Work)

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = 'DeAllocate_Work'

    call Memory_Error_Check (1)
#endif
  endif

#endif
  
  END Subroutine ESSL_FFT_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_COS_Init
! DESCRIPTION  : fills out the WSave structure for the Cosine transform
! ----------------------------------------------------------------------
  
  Subroutine ESSL_COS_Init (N0, ldy, BlockSize, Index, CWave)
  
  integer  :: N0, ldy, BlockSize, Index  
  
  REALTYPE, dimension(1) :: y
  REALTYPE, dimension(*) :: CWave
  
  integer  :: wc_Size, wk_Size
  integer  :: N2, M
  REALTYPE :: Scale
  
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  
  N2 = 2*N0 ; M = BlockSize ; Scale = SQRT(TWO/N0)  
  
  wc_Size = ESSL_CFT_Twiddle_Factor (N2)
  wk_Size = ESSL_CFT_Work           (N2, M, Index)
  
#if 0
  call   Allocate_Work  (Work , wk_Size)
#else
! Allocate (Work(wk_Size))

  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(wk_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

#endif
  
  if (Index == 1) then  
    call COSF (1, y,1,ldy, y,1,ldy, N2, M, Scale, CWave, wc_Size, Work , wk_Size)
  else  
    call COSF (1, y,ldy,1, y,ldy,1, N2, M, Scale, CWave, wc_Size, Work , wk_Size)
  endif  
  
#if 0
  call DeAllocate_Work  (Work , wk_Size)
#else
! DeAllocate (Work)

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = 'DeAllocate_Work'

    call Memory_Error_Check (1)
#endif
  endif

#endif
  
  END Subroutine ESSL_COS_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_SIN_Init
! DESCRIPTION  : fills out the WSave structure for the Sine Transform
! ----------------------------------------------------------------------
  
  Subroutine ESSL_SIN_Init (N0, ldy, BlockSize, Index, SWave)
  
  integer  :: N0, ldy, BlockSize, Index  
  
  REALTYPE, dimension(1) :: y
  REALTYPE, dimension(*) :: SWave
  
  integer  :: ws_Size, wk_Size
  integer  :: N2, M
  REALTYPE :: Scale
  
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  
  N2 = 2*N0 ; M = BlockSize ; Scale = SQRT(TWO/N0)  
  
  ws_Size = ESSL_SFT_Twiddle_Factor (N2)
  wk_Size = ESSL_CFT_Work           (N2, M, Index)
  
#if 0
  call   Allocate_Work  (Work , wk_Size)
#else
! Allocate (Work(wk_Size))

  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(wk_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

#endif
  
  if (Index == 1) then  
    call SINF (1, y,1,ldy, y,1,ldy, N2, M, Scale, SWave, ws_Size, Work , wk_Size)
  else  
    call SINF (1, y,ldy,1, y,ldy,1, N2, M, Scale, SWave, ws_Size, Work , wk_Size)
  endif  
  
#if 0
  call DeAllocate_Work  (Work , wk_Size)
#else
! DeAllocate (Work)

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = 'DeAllocate_Work'

    call Memory_Error_Check (1)
#endif
  endif

#endif
  
  END Subroutine ESSL_SIN_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_CFT_Init
! DESCRIPTION  : fills out the WSave structure for the CFT
! ----------------------------------------------------------------------
  
  Subroutine ESSL_CFT_Init (N0, ldy, BlockSize, Index, WSave)
  
  integer  :: N0, ldy, BlockSize, Index
  
  REALTYPE, dimension(1) :: y
  REALTYPE, dimension(*) :: WSave
  
  integer  :: wc_Size, wk_Size
  integer  :: N2, M
  REALTYPE :: Scale
  
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  
  N2 = 2*N0 ; M = BlockSize ; Scale = SQRT(TWO/N0)  
  
  wc_Size = ESSL_CFT_Twiddle_Factor (N2)
  wk_Size = ESSL_CFT_Work           (N2, M, Index)
  
#if 0
  call   Allocate_Work  (Work , wk_Size)
#else
! Allocate (Work(wk_Size))

  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(wk_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

#endif
  
  if (Index == 1) then  
    call COSF (1, y,1,ldy, y,1,ldy, N2, M, Scale, WSave, wc_Size, Work , wk_Size)
  else  
    call COSF (1, y,ldy,1, y,ldy,1, N2, M, Scale, WSave, wc_Size, Work , wk_Size)
  endif  
  
#if 0
  call DeAllocate_Work  (Work , wk_Size)
#else
! DeAllocate (Work)

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = 'DeAllocate_Work'

    call Memory_Error_Check (1)
#endif
  endif

#endif
  
  END Subroutine ESSL_CFT_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_QWV_Init
! DESCRIPTION  : fills out the WSave structure for the Quarter-Wave tran
! ----------------------------------------------------------------------
  
  Subroutine ESSL_QWV_Init (N_Pnts, QWave)  
  
  integer                  :: N_Pnts
  REALTYPE                 :: PiN2
  REALTYPE, dimension(0:*) :: QWave
  
  integer  :: i, N, N2
  
  if (N_Pnts <= 2) RETURN  
  
    N  = N_Pnts-1  
   N2  = N/2  
  PiN2 = HALF*PI/N  
  
  do i = 0,N2  
    QWave(  i) = SIN(i*PiN2)  
    QWave(N-i) = COS(i*PiN2)  
  enddo  
  
  if (MOD(N,2) == 0) QWave(N2) = HALF*SQRT(TWO)  
  
  END Subroutine ESSL_QWV_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_FFT_Twiddle_Factor
! DESCRIPTION  : Twiddle_Factor Size for RCFT, CRFT transforms
! ----------------------------------------------------------------------
  
  Integer FUNCTION ESSL_FFT_Twiddle_Factor (N)
  
  integer           :: N, Twiddle_Factor
  integer           :: ESSL_FFT, ESSL_FFT_16384, ESSL_FFT_4096
  
#if defined (SINGLE_PRECISION)
  ESSL_FFT       (n) = 25000
  ESSL_FFT_16384 (n) = 20000 + 0.82*n

  if (N <= 16384) Twiddle_Factor = ESSL_FFT       (N)
  if (N >  16384) Twiddle_Factor = ESSL_FFT_16384 (N)
#else
  ESSL_FFT       (n) = 22000
  ESSL_FFT_4096  (n) = 20000 + 0.82*(2*n)

  if (N <= 4096 ) Twiddle_Factor = ESSL_FFT       (N)
  if (N >  4096 ) Twiddle_Factor = ESSL_FFT_4096  (N)
#endif

  ESSL_FFT_Twiddle_Factor = Twiddle_Factor

  END FUNCTION ESSL_FFT_Twiddle_Factor
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_CFT_Twiddle_Factor
! DESCRIPTION  : Twiddle_Factor Size for COSF transforms
! ----------------------------------------------------------------------
  
  Integer FUNCTION ESSL_CFT_Twiddle_Factor (N)

  integer           :: N, Twiddle_Factor
  integer           :: ESSL_CFT, ESSL_CFT_16384

#if defined (SINGLE_PRECISION)
  ESSL_CFT       (n) = 40000
  ESSL_CFT_16384 (n) = 20000 + 0.3*n
#else
  ESSL_CFT       (n) = 35000
  ESSL_CFT_16384 (n) = 20000 + 0.3*(2*n)
#endif
  
  if (N <= 16384) Twiddle_Factor = ESSL_CFT       (N)
  if (N >  16384) Twiddle_Factor = ESSL_CFT_16384 (N)

  ESSL_CFT_Twiddle_Factor = Twiddle_Factor
  
  END FUNCTION ESSL_CFT_Twiddle_Factor
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_SFT_Twiddle_Factor
! DESCRIPTION  : Twiddle_Factor Size for SINF transforms
! ----------------------------------------------------------------------
  
  Integer FUNCTION ESSL_SFT_Twiddle_Factor (N)
  
  integer           :: N, Twiddle_Factor
  integer           :: ESSL_SFT, ESSL_SFT_16384
  
#if defined (SINGLE_PRECISION)
  ESSL_SFT       (n) = 60000
  ESSL_SFT_16384 (n) = 20000 + 0.3*n
#else
  ESSL_SFT       (n) = 50000
  ESSL_SFT_16384 (n) = 20000 + 0.3*(2*n)
#endif

  if (N <= 16384) Twiddle_Factor = ESSL_SFT       (N)
  if (N >  16384) Twiddle_Factor = ESSL_SFT_16384 (N)

  ESSL_SFT_Twiddle_Factor = Twiddle_Factor

  END FUNCTION ESSL_SFT_Twiddle_Factor
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_FFT_Work
! DESCRIPTION  : WorkSpace Size for RCFT, CRFT transforms
! ----------------------------------------------------------------------
  
  Integer FUNCTION ESSL_FFT_Work (N)

  integer           :: N, Work
  integer           :: ESSL_FFT, ESSL_FFT_4096, ESSL_FFT_13684

#if defined (SINGLE_PRECISION)
  ESSL_FFT       (n)   = 20000
  ESSL_FFT_16384 (n)   = 20000 + 0.57*n

  if (N <= 16384) Work = ESSL_FFT       (N)
  if (N >  16384) Work = ESSL_FFT_16384 (N)
#else
  ESSL_FFT       (n)   = 20000
  ESSL_FFT_4096  (n)   = 20000 + 0.57*(2*n)
  
  if (N <= 4096 ) Work = ESSL_FFT       (N)
  if (N >  4096 ) Work = ESSL_FFT_4096  (N)
#endif

  ESSL_FFT_Work = Work

  END FUNCTION ESSL_FFT_Work 
  
! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_CFT_Work
! DESCRIPTION  : WorkSpace Size for COSF, SINF transforms
! ----------------------------------------------------------------------
  
  Integer FUNCTION ESSL_CFT_Work (N, M, Index)

  integer           :: N, Work
  integer           :: M, Index
  integer           :: ESSL_CFT, ESSL_CFT_16384, ESSL_CFT_ADD

#if defined (SINGLE_PRECISION)
  ESSL_CFT       (n,m) = 25000
  ESSL_CFT_16384 (n,m) = 20000 + 0.32*n
  ESSL_CFT_ADD   (n,m) = MIN(128,m)*(257+(2*(n-1))/4)
#else
  ESSL_CFT       (n,m) = 20000
  ESSL_CFT_16384 (n,m) = 20000 + 0.32*(2*n)
  ESSL_CFT_ADD   (n,m) = MIN(128,m)*(257+(2*(n-1))/2)
#endif
  
  if (N <= 16384) Work = ESSL_CFT       (N,M)
  if (N >  16384) Work = ESSL_CFT_16384 (N,M)

  if (Index == 2) then  ! Transpose Algorithm and N >= 252
    if (N >= 252) Work = Work + ESSL_CFT_ADD (N,M)
  endif

  ESSL_CFT_Work = Work

  END FUNCTION ESSL_CFT_Work 

! ----------------------------------------------------------------------
! FUNCTION NAME: ESSL_Adjust_N
! DESCRIPTION  : Adjust vector lengths to conform to ESSL FFT requiremen
! ----------------------------------------------------------------------
  
  Subroutine ESSL_Adjust_N (Method, NN, M, Index)
  
  integer                   :: Method, NN, M, Index, N
  REALTYPE, dimension(NN,M) :: f
  Character(LEN=8)          :: S2030  
  
  integer  :: wc_Size, wk_Size
  integer  :: N2
  REALTYPE :: Scale
  
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  REALTYPE, dimension(:), ALLOCATABLE :: WSave
  
  External ENoTrm  
  
  call EInfo  (0)  
  call ErrSav (2030, S2030)  
  call ErrSet (2030, 0, -1, 1, ENoTrm, 0)  
  
  N = NN-1; N2 = 2*N; Scale = SQRT(TWO/N)  
  
  wc_Size = ESSL_CFT_Twiddle_Factor (N2)
  wk_Size = ESSL_CFT_Work           (N2, M, Index)
  
#if 0
  call   Allocate_WSave (WSave, wc_Size)
  call   Allocate_Work  (Work , wk_Size)
#else
! Allocate (WSave(wc_Size), Work(wk_Size))

  if (.NOT. ALLOCATED(WSave)) then
    ALLOCATE (WSave(wc_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wc_Size
    Subroutine_Name  = '  Allocate_WSave'

    call Memory_Error_Check (0)
#endif
  endif

  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(wk_Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

#endif
  
  call COSF (1, f,1,NN, f,1,NN, N2,1, Scale, WSave, wc_Size, Work , wk_Size, *1)  
  
  goto 10
  
   1 CONTINUE  
  
  if (Method == 0) NN = N2  
  if (Method == 1) NN = N2/2+1  
  
  10 CONTINUE
  
  call ErrStr (2030, S2030)  
  
#if 0
  call DeAllocate_Work  (Work , wk_Size)
  call DeAllocate_WSave (WSave, wc_Size)
#else
! DeAllocate (Work, WSave)

  if (ALLOCATED(WSave)) then
    DEALLOCATE (WSave, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wc_Size
    Subroutine_Name  = 'DeAllocate_WSave'

    call Memory_Error_Check (1)
#endif
  endif

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = wk_Size
    Subroutine_Name  = 'DeAllocate_Work'

    call Memory_Error_Check (1)
#endif
  endif

#endif
  
  END Subroutine ESSL_Adjust_N
#if 0
!
!=======================================================================
!
  Subroutine Allocate_WSave (WSave, Size)
  
  integer  :: Size
  REALTYPE, dimension(:), ALLOCATABLE :: WSave
  
  if (.NOT. ALLOCATED(WSave)) then
    ALLOCATE (WSave(Size), STAT=M_Error)
  
#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = Size
    Subroutine_Name  = '  Allocate_WSave'
 
    call Memory_Error_Check (0)
#endif
  endif
  
  END Subroutine Allocate_WSave
!
!=======================================================================
!
  Subroutine DeAllocate_WSave (WSave, Size)
 
  integer  :: Size
  REALTYPE, dimension(:), ALLOCATABLE :: WSave
 
  if (ALLOCATED(WSave)) then
    DEALLOCATE (WSave, STAT=M_Error)
 
#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = Size
    Subroutine_Name  = 'DeAllocate_WSave'

    call Memory_Error_Check (1)
#endif
  endif
  
  END Subroutine DeAllocate_WSave
!
!=======================================================================
!
  Subroutine Allocate_Work (Work, Size)
  
  integer  :: Size
  REALTYPE, dimension(:), ALLOCATABLE :: Work
  
  if (.NOT. ALLOCATED(Work)) then
    ALLOCATE (Work(Size), STAT=M_Error)
  
#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = Size
    Subroutine_Name  = '  Allocate_Work'

    call Memory_Error_Check (0)
#endif
  endif

  END Subroutine Allocate_Work
!
!=======================================================================
!
  Subroutine DeAllocate_Work (Work, Size)

  integer  :: Size
  REALTYPE, dimension(:), ALLOCATABLE :: Work

  if (ALLOCATED(Work)) then
    DEALLOCATE (Work, STAT=M_Error)
  
#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = Size
    Subroutine_Name  = 'DeAllocate_Work'
  
    call Memory_Error_Check (1)
#endif
  endif
  
  END Subroutine DeAllocate_Work
#endif
  
END MODULE ESSL_FFT

#else

MODULE ESSL_FFT

implicit NONE

END MODULE ESSL_FFT

#endif
