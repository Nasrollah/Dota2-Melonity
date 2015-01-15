#if defined (CRAY)
  
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
!   SinCos    is filled out by Cray_QWAVE_Init with vector lengths N.
!   WSave     is filled out by Cray_CFT_Init   with vector lengths N.
!   WorkSpace is of size (N+3)*M*3.
!               Padded WorkSpace so that its have odd vector lengths to
!               memory contention for vectors loads and stores.
!
!  For Sine   Transform, only the first N points are used k=0,..,N-1
!  For Cosine, Quarter wave Transform, all N+1 points are used k=0,..,N
!
! Date    : 8-10-99
! Author  : Wai Sun Don
! Comment : Fortran 90 ized.
!
! ----------------------------------------------------------------------
  
  
MODULE CRAY_FFT
  
#if defined (LOCAL_CRAY_FFT)
  USE CRAY_FFT_JM
#endif

implicit NONE

PRIVATE
  
PUBLIC  :: Cray_SFT_1 
PUBLIC  :: Cray_CFT_1 
PUBLIC  :: Cray_QWAVE_FWRD_1 
PUBLIC  :: Cray_QWAVE_BWRD_1 
  
PUBLIC  :: Cray_SFT_2 
PUBLIC  :: Cray_CFT_2 
PUBLIC  :: Cray_QWAVE_FWRD_2 
PUBLIC  :: Cray_QWAVE_BWRD_2 
  
PUBLIC  :: Cray_FFT_Init 
PUBLIC  :: Cray_COS_Init 
PUBLIC  :: Cray_SIN_Init 
PUBLIC  :: Cray_CFT_Init 
PUBLIC  :: Cray_QWV_Init 
  
PUBLIC  :: FFTFAX_wrap 
  
CONTAINS
  
! ----------------------------------------------------------------------
!                  <--  Along the First Index  -->
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_SFT_1
! DESCRIPTION  : The Cooley-Lewis-Welch algorithm for the Sine transform
! ----------------------------------------------------------------------
  
  Subroutine Cray_SFT_1 (N, M, A, LD_A, B, LD_B, WSave)
  
  integer  :: N, M, LD_A, LD_B

  REALTYPE, dimension(0:LD_A-1,M) :: A
  REALTYPE, dimension(0:LD_B-1,M) :: B
  REALTYPE, dimension(*)          :: WSave
  REALTYPE, dimension(0:N+2,M,3)  :: WorkSpace

  integer  :: N2, N_Skip, M_Skip, i, j
  REALTYPE :: Scale

  integer , dimension(19) :: Factor

#ifdef DEBUG
  ASSERT(Mod(N,2)==0,"C_SFT_1:N not Even")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2   = N/2  
  do j = 1,M  
    WorkSpace(0,j,1) = TWO*A(1,j)  
    WorkSpace(1,j,1) = ZERO  
    do i = 2,N-2,2  
      WorkSpace(i  ,j,1) =  A(i+1,j)-A(i-1,j)  
      WorkSpace(i+1,j,1) = -A(i  ,j)  
    enddo  
    WorkSpace(N  ,j,1) = -TWO*A(N-1,j)  
    WorkSpace(N+1,j,1) = ZERO  
  enddo  
  
  N_Skip = 1  
  M_Skip = N + 3  
  
  Factor = INT(WSave(N2+1:N2+19))

  call rfftmlt (WorkSpace(0,1,1), WorkSpace(0,1,2), &
                WSave(N2+20), Factor, N_Skip, M_Skip, N, M, 1)
  
  Scale = ONE/SQRT(TWO*N)  
  
  do j = 1,M  
    do i = 1,N2  
      WorkSpace(i,j,2) =     HALF*(WorkSpace(i,j,1) - WorkSpace(N-i,j,1))
      WorkSpace(i,j,3) = WSave(i)*(WorkSpace(i,j,1) + WorkSpace(N-i,j,1))
    enddo  
  enddo  
  
  do j = 1,M  
    do i = 1,N2-1  
      B(  i,j) = Scale*(WorkSpace(i,j,3) + WorkSpace(i,j,2))  
      B(N-i,j) = Scale*(WorkSpace(i,j,3) - WorkSpace(i,j,2))  
    enddo  
  enddo  
  
  i = N2  
  do j = 1,M  
    B(i,j) = Scale*(WorkSpace(i,j,3))  
    B(0,j) = ZERO  
  enddo  
  
  END Subroutine Cray_SFT_1
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_CFT_1
! DESCRIPTION  : The Cooley-Lewis-Welch algorithm for the cosine transform
! ----------------------------------------------------------------------
  
  Subroutine Cray_CFT_1 (N, M, A, LD_A, B, LD_B, WSave)
  
  integer  :: N, M, LD_A, LD_B

  REALTYPE, dimension(0:LD_A-1,M) :: A
  REALTYPE, dimension(0:LD_B-1,M) :: B
  REALTYPE, dimension(*)          :: WSave
  REALTYPE, dimension(0:N+2,M,3)  :: WorkSpace

  integer  :: N2, N_Skip, M_Skip, i, j
  REALTYPE :: Scale, A_Temp

  integer , dimension(19) :: Factor

#ifdef DEBUG
  ASSERT(Mod(N,2)==0,"C_CFT_1:N not Even")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2   = N/2  
  do j = 1,M  
    WorkSpace(0,j,1) = A(0,j)  
    WorkSpace(1,j,1) = ZERO  
    do i = 2,N-2,2  
      WorkSpace(i  ,j,1) = A(i  ,j)  
      WorkSpace(i+1,j,1) = A(i+1,j)-A(i-1,j)  
    enddo  
    WorkSpace(N  ,j,1) = A(N,j)  
    WorkSpace(N+1,j,1) = ZERO  
  enddo  
  
  N_Skip = 1  
  M_Skip = N + 3  
  
  Factor = INT(WSave(N2+1:N2+19))

  call rfftmlt (WorkSpace(0,1,1), WorkSpace(0,1,2), &
                WSave(N2+20), Factor, N_Skip, M_Skip, N, M, 1)
  
  do j = 1,M  
    A_Temp = ZERO  
    do i = 1,N,2  
      A_Temp = A_Temp + A(i,j)  
    enddo  
    WorkSpace(N+2,j,1) = TWO*A_Temp  
  enddo  
  
  Scale = ONE/SQRT(TWO*N)  
  
  do j = 1,M  
    do i = 1,N2  
      WorkSpace(i,j,2) =     HALF*(WorkSpace(i,j,1) + WorkSpace(N-i,j,1))
      WorkSpace(i,j,3) = WSave(i)*(WorkSpace(i,j,1) - WorkSpace(N-i,j,1))
    enddo  
  enddo  
  
  do j = 1,M  
    do i = 1,N2-1  
      B(  i,j) = Scale*(WorkSpace(i,j,2) + WorkSpace(i,j,3))  
      B(N-i,j) = Scale*(WorkSpace(i,j,2) - WorkSpace(i,j,3))  
    enddo  
  enddo  
  
  i = N2  
  do j = 1,M  
    B(i,j) = Scale*(WorkSpace(i,j,2) + WorkSpace(i  ,j,3))  
    B(0,j) = Scale*(WorkSpace(0,j,1) + WorkSpace(N+2,j,1))  
    B(N,j) = Scale*(WorkSpace(0,j,1) - WorkSpace(N+2,j,1))  
  enddo  
  
  END Subroutine Cray_CFT_1
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_QWAVE_FWRD_1
! DESCRIPTION  : Forward Quarter-Wave Cosine transform of an Odd sequence 
! ----------------------------------------------------------------------
  
  Subroutine Cray_QWAVE_FWRD_1 (N, M, A, LD_A, B, LD_B, SinCos, WSave)
  
  integer  :: N, M, LD_A, LD_B

  REALTYPE, dimension(0:LD_A-1,M) :: A
  REALTYPE, dimension(0:LD_B-1,M) :: B
  REALTYPE, dimension(0:N)        :: SinCos
  REALTYPE, dimension(*)          :: WSave
  REALTYPE, dimension(0:N+2,M,3)  :: WorkSpace

  integer  :: N2, N_Skip, M_Skip, j, k
  REALTYPE :: Scale, A_Even, A_Odd

  integer , dimension(19) :: Factor

#ifdef DEBUG
  ASSERT(MOD(N,2)==0,"C_QWAVE_F_1:N not Even")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2 = N/2  
  do j = 1,M  
    do k = 1,N2-1  
      A_Even = A(k,j)+A(N-k,j)  
      A_Odd  = A(k,j)-A(N-k,j)  
  
      WorkSpace(k  ,j,1) = SinCos(N-k)*A_Odd +SinCos(k)*A_Even  
      WorkSpace(N-k,j,1) = SinCos(N-k)*A_Even-SinCos(k)*A_Odd  
    enddo  
  enddo  
  
  k = N2  
  do j = 1,M  
    WorkSpace(k,j,1) = TWO*SinCos(k)*A(k,j)  
    WorkSpace(0,j,1) = A(0,j)  
    WorkSpace(N,j,1) = A(0,j)  
  enddo  
  
  N_Skip = 1  
  M_Skip = N + 3  
  
  Factor = INT(WSave(N2+1:N2+19))

  call rfftmlt (WorkSpace(0,1,1), WorkSpace(0,1,2), &
                WSave(N2+20), Factor, N_Skip, M_Skip, N, M, -1)
  
  Scale = HALF*SQRT(TWO*N)  
  
  do j = 1,M  
    B(0,j) = Scale*WorkSpace(0,j,1)  
    do k = 2,N-1,2  
      B(k-1,j) = Scale*(WorkSpace(k,j,1)-WorkSpace(k+1,j,1))  
      B(k  ,j) = Scale*(WorkSpace(k,j,1)+WorkSpace(k+1,j,1))  
    enddo  
    B(N-1,j) = Scale*WorkSpace(N,j,1)  
  enddo  
  
  END Subroutine Cray_QWAVE_FWRD_1
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_QWAVE_BWRD_1
! DESCRIPTION  : Backward Quarter-Wave Cosine transform of an Odd sequence
! ----------------------------------------------------------------------
  
  Subroutine Cray_QWAVE_BWRD_1 (N, M, A, LD_A, B, LD_B, SinCos, WSave)
  
  integer  :: N, M, LD_A, LD_B

  REALTYPE, dimension(0:LD_A-1,M) :: A
  REALTYPE, dimension(0:LD_B-1,M) :: B
  REALTYPE, dimension(0:N)        :: SinCos
  REALTYPE, dimension(*)          :: WSave
  REALTYPE, dimension(0:N+2,M,3)  :: WorkSpace

  integer  :: N2, N_Skip, M_Skip, j, k
  REALTYPE :: Scale, B_Even, B_Odd

  integer , dimension(19) :: Factor

#ifdef DEBUG
  ASSERT(MOD(N,2)==0,"C_QWAVE_B_1:N not Even")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2 = N/2  
  do j = 1,M  
    WorkSpace(0,j,1) = TWO*A(0,j)  
    WorkSpace(1,j,1) = ZERO  
    do k = 2,N-1,2  
      WorkSpace(k  ,j,1) = A(k-1,j)+A(k,j)  
      WorkSpace(k+1,j,1) = A(k-1,j)-A(k,j)  
    enddo  
    WorkSpace(N,j,1) = TWO*A(N-1,j)  
  enddo  
  
  N_Skip = 1  
  M_Skip = N + 3  
  
  Factor = INT(WSave(N2+1:N2+19))

  call rfftmlt (WorkSpace(0,1,1), WorkSpace(0,1,2), &
                WSave(N2+20), Factor, N_Skip, M_Skip, N, M,  1)
  
  Scale = HALF/SQRT(TWO*N)  
  
  do j = 1,M  
    B(0 ,j) = Scale*TWO*WorkSpace(0,j,1)  
  enddo  
  
  do j = 1,M  
    do k = 1,N2-1  
        B_Even = Scale*(SinCos(k)+SinCos(N-k))  
        B_Odd  = Scale*(SinCos(k)-SinCos(N-k))  
  
      B(k  ,j) = WorkSpace(N-k,j,1)*B_Even-WorkSpace(k,j,1)*B_Odd  
      B(N-k,j) = WorkSpace(N-k,j,1)*B_Odd +WorkSpace(k,j,1)*B_Even  
    enddo  
  enddo  
  
  k = N2  
  do j = 1,M  
    B(k,j) = Scale*TWO*SinCos(k)*WorkSpace(k,j,1)  
    B(N,j) = ZERO  
  enddo  
  
  END Subroutine Cray_QWAVE_BWRD_1
  
! ----------------------------------------------------------------------
!                  <--  Along the Second Index  -->
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_SFT_2
! DESCRIPTION  : The Cooley-Lewis-Welch algorithm for the Sine   transform
! ----------------------------------------------------------------------
  
  Subroutine Cray_SFT_2 (N, M, A, LD_A, B, LD_B, WSave)
  
  integer  :: N, M, LD_A, LD_B

  REALTYPE, dimension(LD_A,0:N)   :: A
  REALTYPE, dimension(LD_B,0:N)   :: B
  REALTYPE, dimension(*)          :: WSave
  REALTYPE, dimension(0:N+2,M,3)  :: WorkSpace

  integer  :: N2, N_Skip, M_Skip, i, j
  REALTYPE :: Scale

  integer , dimension(19) :: Factor

#ifdef DEBUG
  ASSERT(Mod(N,2)==0,"C_SFT_2:N not Even")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2   = N/2  
  do j = 1,M  
    WorkSpace(0,j,1) = TWO*A(j,1)  
    WorkSpace(1,j,1) = ZERO  
    do i = 2,N-2,2  
      WorkSpace(i  ,j,1) =  A(j,i+1)-A(j,i-1)  
      WorkSpace(i+1,j,1) = -A(j,i  )  
    enddo  
    WorkSpace(N  ,j,1) = -TWO*A(j,N-1)  
    WorkSpace(N+1,j,1) = ZERO  
  enddo  
  
  N_Skip = 1  
  M_Skip = N + 3  
  
  Factor = INT(WSave(N2+1:N2+19))

  call rfftmlt (WorkSpace(0,1,1), WorkSpace(0,1,2), &
                WSave(N2+20), Factor, N_Skip, M_Skip, N, M, 1)
  
  Scale = ONE/SQRT(TWO*N)  
  
  do j = 1,M  
    do i = 1,N2  
      WorkSpace(i,j,2) =     HALF*(WorkSpace(i,j,1) - WorkSpace(N-i,j,1))
      WorkSpace(i,j,3) = WSave(i)*(WorkSpace(i,j,1) + WorkSpace(N-i,j,1))
    enddo  
  enddo  
  
  do i = 1,N2-1  
    do j = 1,M  
      B(j,  i) = Scale*(WorkSpace(i,j,3) + WorkSpace(i,j,2))  
      B(j,N-i) = Scale*(WorkSpace(i,j,3) - WorkSpace(i,j,2))  
    enddo  
  enddo  
  
  i = N2  
  do j = 1,M  
    B(j,i) = Scale*(WorkSpace(i,j,3))  
    B(j,0) = ZERO  
  enddo  
  
  END Subroutine Cray_SFT_2
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_CFT_2
! DESCRIPTION  : The Cooley-Lewis-Welch algorithm for the cosine transform
! ----------------------------------------------------------------------
  
  Subroutine Cray_CFT_2 (N, M, A, LD_A, B, LD_B, WSave)
  
  integer  :: N, M, LD_A, LD_B

  REALTYPE, dimension(LD_A,0:N)   :: A
  REALTYPE, dimension(LD_B,0:N)   :: B
  REALTYPE, dimension(*)          :: WSave
  REALTYPE, dimension(0:N+2,M,3)  :: WorkSpace

  integer  :: N2, N_Skip, M_Skip, i, j
  REALTYPE :: Scale, A_Temp

  integer , dimension(19) :: Factor

#ifdef DEBUG
  ASSERT(Mod(N,2)==0,"C_CFT_2:N not Even")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2   = N/2  
  do j = 1,M  
    WorkSpace(0,j,1) = A(j,0)  
    WorkSpace(1,j,1) = ZERO  
    do i = 2,N-2,2  
      WorkSpace(i  ,j,1) = A(j,i)  
      WorkSpace(i+1,j,1) = A(j,i+1)-A(j,i-1)  
    enddo  
    WorkSpace(N  ,j,1) = A(j,N)  
    WorkSpace(N+1,j,1) = ZERO  
  enddo  
  
  N_Skip = 1  
  M_Skip = N + 3  
  
  Factor = INT(WSave(N2+1:N2+19))

  call rfftmlt (WorkSpace(0,1,1), WorkSpace(0,1,2), &
                WSave(N2+20), Factor, N_Skip, M_Skip, N, M, 1)
  
  do j = 1,M  
    A_Temp = ZERO  
    do i = 1,N,2  
      A_Temp = A_Temp + A(j,i)  
    enddo  
    WorkSpace(N+2,j,1) = TWO*A_Temp  
  enddo  
  
  Scale = ONE/SQRT(TWO*N)  
  
  do j = 1,M  
    do i = 1,N2  
      WorkSpace(i,j,2) =     HALF*(WorkSpace(i,j,1) + WorkSpace(N-i,j,1))
      WorkSpace(i,j,3) = WSave(i)*(WorkSpace(i,j,1) - WorkSpace(N-i,j,1))
    enddo  
  enddo  
  
  do j = 1,M  
    do i = 1,N2-1  
      B(j,  i) = Scale*(WorkSpace(i,j,2) + WorkSpace(i,j,3))  
      B(j,N-i) = Scale*(WorkSpace(i,j,2) - WorkSpace(i,j,3))  
    enddo  
  enddo  
  
  i = N2  
  do j = 1,M  
    B(j,i) = Scale*(WorkSpace(i,j,2) + WorkSpace(i  ,j,3))  
    B(j,0) = Scale*(WorkSpace(0,j,1) + WorkSpace(N+2,j,1))  
    B(j,N) = Scale*(WorkSpace(0,j,1) - WorkSpace(N+2,j,1))  
  enddo  
  
  END Subroutine Cray_CFT_2
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_QWAVE_FWRD_2
! DESCRIPTION  : Forward Quarter-Wave Cosine transform of an Odd sequence
! ----------------------------------------------------------------------
  
  Subroutine Cray_QWAVE_FWRD_2 (N, M, A, LD_A, B, LD_B, SinCos, WSave)
  
  integer  :: N, M, LD_A, LD_B

  REALTYPE, dimension(LD_A,0:N)   :: A
  REALTYPE, dimension(LD_B,0:N)   :: B
  REALTYPE, dimension(0:N)        :: SinCos
  REALTYPE, dimension(*)          :: WSave
  REALTYPE, dimension(0:N+2,M,3)  :: WorkSpace

  integer  :: N2, N_Skip, M_Skip, j, k
  REALTYPE :: Scale, A_Even, A_Odd

  integer , dimension(19) :: Factor

#ifdef DEBUG
  ASSERT(MOD(N,2)==0,"C_QWAVE_F_2:N not Even")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2 = N/2  
  do j = 1,M  
    do k = 1,N2-1  
      A_Even = A(j,k)+A(j,N-k)  
      A_Odd  = A(j,k)-A(j,N-k)  
  
      WorkSpace(k  ,j,1) = SinCos(N-k)*A_Odd +SinCos(k)*A_Even  
      WorkSpace(N-k,j,1) = SinCos(N-k)*A_Even-SinCos(k)*A_Odd  
    enddo  
  enddo  
  
  k = N2  
  do j = 1,M  
    WorkSpace(k,j,1) = TWO*SinCos(k)*A(j,k)  
    WorkSpace(0,j,1) = A(j,0)  
    WorkSpace(N,j,1) = A(j,0)  
  enddo  
  
  N_Skip = 1  
  M_Skip = N + 3  
  
  Factor = INT(WSave(N2+1:N2+19))

  call rfftmlt (WorkSpace(0,1,1), WorkSpace(0,1,2), &
                WSave(N2+20), Factor, N_Skip, M_Skip, N, M, -1)
  
  Scale = HALF*SQRT(TWO*N)  
  
  do j = 1,M  
    B(j,0) = Scale*WorkSpace(0,j,1)  
    do k = 2,N-1,2  
      B(j,k-1) = Scale*(WorkSpace(k,j,1)-WorkSpace(k+1,j,1))  
      B(j,k  ) = Scale*(WorkSpace(k,j,1)+WorkSpace(k+1,j,1))  
    enddo  
    B(j,N-1) = Scale*WorkSpace(N,j,1)  
  enddo  
  
  END Subroutine Cray_QWAVE_FWRD_2
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_QWAVE_BWRD_2
! DESCRIPTION  : Backward Quarter-Wave Cosine transform of an Odd sequence
! ----------------------------------------------------------------------
  
  Subroutine Cray_QWAVE_BWRD_2 (N, M, A, LD_A, B, LD_B, SinCos, WSave)
  
  integer  :: N, M, LD_A, LD_B

  REALTYPE, dimension(LD_A,0:N)   :: A
  REALTYPE, dimension(LD_B,0:N)   :: B
  REALTYPE, dimension(0:N)        :: SinCos
  REALTYPE, dimension(*)          :: WSave
  REALTYPE, dimension(0:N+2,M,3)  :: WorkSpace

  integer  :: N2, N_Skip, M_Skip, j, k
  REALTYPE :: Scale, B_Even, B_Odd

  integer , dimension(19) :: Factor

#ifdef DEBUG
  ASSERT(MOD(N,2)==0,"C_QWAVE_B_1:N not Even")  
#endif
  
  if ((M <= 0) .OR. (N <= 1)) RETURN  
  
  N2 = N/2  
  do j = 1,M  
    WorkSpace(0,j,1) = TWO*A(j,0)  
    WorkSpace(1,j,1) = ZERO  
    do k = 2,N-1,2  
      WorkSpace(k  ,j,1) = A(j,k-1)+A(j,k)  
      WorkSpace(k+1,j,1) = A(j,k-1)-A(j,k)  
    enddo  
    WorkSpace(N,j,1) = TWO*A(j,N-1)  
  enddo  
  
  N_Skip = 1  
  M_Skip = N + 3  
  
  Factor = INT(WSave(N2+1:N2+19))

  call rfftmlt (WorkSpace(0,1,1), WorkSpace(0,1,2), &
                WSave(N2+20), Factor, N_Skip, M_Skip, N, M,  1)
  
  Scale = HALF/SQRT(TWO*N)  
  
  do j = 1,M  
    B(j,0) = Scale*TWO*WorkSpace(0,j,1)  
  enddo  
  
  do j = 1,M  
    do k = 1,N2-1  
        B_Even = Scale*(SinCos(k)+SinCos(N-k))  
        B_Odd  = Scale*(SinCos(k)-SinCos(N-k))  
  
      B(j,k  ) = WorkSpace(N-k,j,1)*B_Even-WorkSpace(k,j,1)*B_Odd  
      B(j,N-k) = WorkSpace(N-k,j,1)*B_Odd +WorkSpace(k,j,1)*B_Even  
    enddo  
  enddo  
  
  k = N2  
  do j = 1,M  
    B(j,k) = Scale*TWO*SinCos(k)*WorkSpace(k,j,1)  
    B(j,N) = ZERO  
  enddo  
  
  END Subroutine Cray_QWAVE_BWRD_2
  
! ----------------------------------------------------------------------
!                  <--  Initiailization Routines  -->
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_FFT_Init
! DESCRIPTION  : fills out the WSave structure for the cosine transform
! ----------------------------------------------------------------------
  
  Subroutine Cray_FFT_Init (N, WSave)  
  
  integer                  :: N
  REALTYPE, dimension(*)   :: WSave
  
  integer , dimension(19) :: Factor

  if (N <= 1) RETURN  
  
  call FFTFAX_wrap (N, Factor, WSave(20))  
  
  WSave(1:19) = Factor
  
  END Subroutine Cray_FFT_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_COS_Init
! DESCRIPTION  : fills out the WSave structure for the cosine transform
! ----------------------------------------------------------------------
  
  Subroutine Cray_COS_Init (N, WSave)  
  
  integer                  :: N
  REALTYPE, dimension(*)   :: WSave

  integer  :: N2, i
  REALTYPE :: PiN

  integer , dimension(19) :: Factor

  if (N <= 1) RETURN  
  
   N2 = N/2  
  PiN = PI/N  
  
  do i = 1,N2  
    WSave(i) = QUARTER/SIN(i*PiN)  
  enddo  
  
  if (MOD(N,2) == 0) WSave(N2) = QUARTER  
  
  call FFTFAX_wrap (N, Factor, WSave(N2+20))  
  
  WSave(N2+1:N2+19) = Factor
  
  END Subroutine Cray_COS_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_SIN_Init
! DESCRIPTION  : fills out the WSave structure for the cosine transform
! ----------------------------------------------------------------------
  
  Subroutine Cray_SIN_Init (N, WSave)  
  
  integer                  :: N
  REALTYPE, dimension(*)   :: WSave
  
  integer  :: N2, i
  REALTYPE :: PiN

  integer , dimension(19) :: Factor

  if (N <= 1) RETURN  
  
   N2 = N/2  
  PiN = PI/N  
  
  do i = 1,N2  
    WSave(i) = QUARTER/SIN(i*PiN)  
  enddo  
  
  if (MOD(N,2) == 0) WSave(N2) = QUARTER  
  
  call FFTFAX_wrap (N, Factor, WSave(N2+20))  
  
  WSave(N2+1:N2+19) = Factor
  
  END Subroutine Cray_SIN_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_CFT_Init
! DESCRIPTION  : fills out the WSave structure for the cosine transform
! ----------------------------------------------------------------------
  
  Subroutine Cray_CFT_Init (N_Pnts, WSave)  
  
  integer                  :: N_Pnts, N
  REALTYPE, dimension(*)   :: WSave
  
  integer  :: N2, i
  REALTYPE :: PiN

  integer , dimension(19) :: Factor

  if (N_Pnts <= 2) RETURN  
  
    N = N_Pnts-1  
   N2 = N/2  
  PiN = PI/N  
  
  do i = 1,N2  
    WSave(i) = QUARTER/SIN(i*PiN)  
  enddo  
  
  if (MOD(N,2) == 0) WSave(N2) = QUARTER  
  
  call FFTFAX_wrap (N, Factor, WSave(N2+20))  
  
  WSave(N2+1:N2+19) = Factor
  
  END Subroutine Cray_CFT_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: Cray_QWV_Init
! DESCRIPTION  : fills out the WSave structure for the Quarter-Wave tranform
! ----------------------------------------------------------------------
  
  Subroutine Cray_QWV_Init (N_Pnts, QWave)  
  
  integer                  :: N_Pnts, N
  REALTYPE, dimension(0:*) :: QWave
  
  integer  :: N2, i
  REALTYPE :: PiN2

  if (N_Pnts <= 2) RETURN  
  
    N  = N_Pnts-1  
   N2  = N/2  
  PiN2 = HALF*PI/N  
  
  do i = 0,N2  
    QWave(  i) = SIN(i*PiN2)  
    QWave(N-i) = COS(i*PiN2)  
  enddo  
  
  if (MOD(N,2) == 0) QWave(N2) = HALF*SQRT(TWO)  
  
  END Subroutine Cray_QWV_Init
  
! ----------------------------------------------------------------------
! FUNCTION NAME: FFTFAX_wrap
! DESCRIPTION  : wrapper for FFTFAX that checks for failure
! ----------------------------------------------------------------------
  
  Subroutine FFTFAX_wrap (N, Factors, Trigs)  
  
  integer                  :: N
  integer , dimension(19)  :: Factors
  REALTYPE, dimension(2*N) :: Trigs
  
  call FFTFAX (N, Factors, Trigs)  
  
  if (Factors(1) == -99) then  
    write (6 ,*) 'FFTFAX : Invalid Transform length N= ',N  
    write (99,*) 'FFTFAX : Invalid Transform length N= ',N  
    STOP  
  endif  
  
  END Subroutine FFTFAX_wrap
  
END MODULE CRAY_FFT
  
#else
  
MODULE CRAY_FFT
  
implicit NONE

END MODULE CRAY_FFT
  
#endif
