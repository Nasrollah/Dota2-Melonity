#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Diff_EOD
! DESCRIPTION  : sets up the aux Diff structure for the (possibly Mapped)
!                Legendre Derivative using the plain EOD algorithm.
!                Also computes the x that functions should be
!                sampled at.
!
!  Point_Type = 1 Lobatto Point
!             = 2 Radau   Point
!             = 3 Gauss   Point
!
!  Index      = 1 if Differentiating along the first index,
!               2 if Differentiating along the second index.
!
!  x_Left     = Left  Physical Domain Point
!  x_Right    = Right Physical Domain Point
!
!  N_Pnts     = length of the vectors
!  M          = total number of vectors.
!  ldy        = the leading dimension of the array to be processed. (NOT USED)
!  BlockSize  = number of vectors to be processed at one time.
!                 if BlockSize=0, BlockSize=M  (NOT USED HERE)
!  Symmetry   = 0 if no symmetry assumed
!               1 if data set is symmetric/anti-symmetric about mid-points.
!
! For Smoothing the Derivative...
!
!     call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)
!
!  Filter_Choice = Choice of Filtering function
!                    1. Exponential Filter (Mode_CutOff, Omega, Order)
!                    2. Lanczos Filter (2nd order)
!                    3. Raised Cosine Filter (2nd order)
!                    4. 8th order Sharpened Rasied Cosine Filter
!  Mode_CutOff   = Cut-Off Frequency for the Exponential Filter
!                    if equal N_Pnts-1, Filter_Choice=0
!  Omega         = constant in which the value of the Nth mode = Exp(-Omega)
!  Order         = Order of exponential filter = Exp(-a(k/N)**Order)
!  Sigma         = Filter function
!  Smooth_1      = Smooth the Flux BEFORE Recursion
!                    0 = NO
!                    1 = Yes
!  Smooth_2      = Smooth the Flux AFTER  Recursion
!                    0 = NO
!                    1 = Yes
!
! ----------------------------------------------------------------------


MODULE Chebyshev_Differentiation_EOD
  USE Chebyshev_Grid_Point

  USE EOD_Auxillary_Subroutine

  USE Chebyshev_GL_EOD
  USE Chebyshev_Ga_EOD

  USE Grid_Generation, ONLY: Grid_And_Mapping

  USE Differentiation_EOD_1, C_Diff_EOD_1=>Diff_EOD_1
  USE Differentiation_EOD_2, C_Diff_EOD_2=>Diff_EOD_2

implicit NONE

PRIVATE

PUBLIC  :: C_Setup_Aux_Diff_EOD
PUBLIC  :: C_Deriv_Matrix_EOD
PUBLIC  :: C_Diff_EOD_1
PUBLIC  :: C_Diff_EOD_2

CONTAINS

  Subroutine C_Setup_Aux_Diff_EOD (Index, Method, Point_Type,          &
                                   Max_Order_Of_D,                     &
                                   Algorithm, Symmetry,                &
                                   N_Pnts, M, ldy, BlockSize,          &
                                   x, x_Left,x_Right,                  &
                                   Map, Manual, alpha, beta,           &
                                   Filter_Choice, Mode_CutOff,         &
                                   Omega, Order, Smooth_1, Smooth_2,   &
                                   Aux_Diff, EOD_Code)

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry
  integer  :: N_Pnts, M, ldy, BlockSize, Map, Manual
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  integer  :: EOD_Code
  integer  :: O_Of_D
  REALTYPE :: x_Left,x_Right
  REALTYPE :: alpha, beta, Omega, Order

  REALTYPE, dimension(:) :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: N, N0, N1, N2
  integer  :: i_m_s, i_e_s, i_o_s, i_Offset

  N  = N_Pnts-1
  N0 = N       
  N1 = N+1
  N2 = N/2

  if (Symmetry == 1) N0 = N2

  i_m_s = EOD_METRIC_START(N1)

  call Grid_And_Mapping (Method, Point_Type,  &
                         N, N0, Map, Manual, alpha, beta,  &
                         x_Left, x_Right, x, Aux_Diff(i_m_s))

  i_e_s = EOD_MATRIX_START_EVEN(N1)
  i_o_s = EOD_MATRIX_START_ODD(N1)

  call C_Deriv_Matrix_EOD (Max_Order_Of_D, Point_Type,  &
                           Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2,  & 
                           Aux_Diff(i_e_s), Aux_Diff(i_o_s))

  Aux_Diff(_I_MATRIX_START)     = i_e_s
  Aux_Diff(_I_ODD_MATRIX_START) = i_o_s
  Aux_Diff(_I_SIGMA_START)      = 0
  Aux_Diff(_I_METRIC_START)     = i_m_s

  if (Symmetry == 0) then
    do O_Of_D = 1,Max_Order_Of_D

      i_Offset = (O_Of_D-1)*EOD_HALF_SIZE(N1)
      i_e_s    = EOD_MATRIX_START_EVEN(N1) + i_Offset
      i_o_s    = EOD_MATRIX_START_ODD(N1)  + i_Offset

      call Adjust_EOD_Matrix (N2, Aux_Diff(i_e_s), Aux_Diff(i_o_s))
    enddo
  endif

  EOD_Code = 0

  if (Map /= 0) then
    if (Symmetry == 0) call Slope_Check (N, Aux_Diff(i_m_s), EOD_Code)

    if (EOD_Code == 0) then
      i_e_s = EOD_MATRIX_START_EVEN(N1)
      i_o_s = EOD_MATRIX_START_ODD(N1)

      call Apply_Mapping_EOD (Map, N0, N2, Max_Order_Of_D,  &
                              Aux_Diff(i_m_s),              &
                              Aux_Diff(i_e_s), Aux_Diff(i_o_s)) 
    else
      if (Symmetry == 1) then
        write ( 6,100) ; write (99,100) ; STOP
      endif
      if (Max_Order_Of_D >  1) then
        write ( 6,110) ; write (99,110) ; STOP
      endif
    endif
  endif

 100 format (1x,'FATAL ERROR : EOD Algorithm '/  &
             1x,'  Non-Anti-Symmetric mapping metric (EOD_Code=1) '/ &
             1x,'  with (Symmetry=1) for Even-Odd Algorithm (Algorithm=1)')
 110 format (1x,'FATAL ERROR : EOD Algorithm '/  &
             1x,'  Non-Anti-Symmetric mapping metric (EOD_Code=1) '/ &
             1x,'  with Max_Order_Of_Derivative > 1! '// &
             1x,'  This is inefficient for Higher Derivatives.'/ &
             1x,'  Use MXM Algorithm instead') 

  END Subroutine C_Setup_Aux_Diff_EOD 

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Deriv_Matrix_EOD
! AUTHOR   : Wai-Sun Don and Alex Solomonoff
! AUTHOR   : Wai-Sun Don and Bruno Costa  (Rio)
! DESCRIPTION  : computes the elements of the two matrices used in the 
!        Even-Odd algorithm for the Chebyshev Smoothing. Does some
!        crufty Stuff to avoid using O(N^2) temporary storage.
!
!  N+1 = number of grid points
! ----------------------------------------------------------------------

  Subroutine C_Deriv_Matrix_EOD (Max_Order_Of_D, Point_Type,  &
                                 Filter_Choice, N, Mode_CutOff, &
                                 Omega, Order, Smooth_1, Smooth_2,  &
                                 D_Even, D_Odd) 

  integer  :: Point_Type
  integer  :: Max_Order_Of_D, N, O_Of_D
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N/2,0:N/2,Max_Order_Of_D) :: D_Even, D_Odd
  REALTYPE, dimension(0:N/2,0:N/2)                :: S_Even, S_Odd

#if defined (DEBUG)
  ASSERT(N>0, "cheb_der_m_eo:N<=0")
#endif

  SELECT CASE (Point_Type)
    CASE (1)
      call C_GL_EOD_D (N, D_Even, D_Odd)

#ifdef DEBUG
    CASE (2)
      write(6 ,*) 'FATAL ERROR! Chebyshev Radau pt is not allowed for EOD'
      write(99,*) 'FATAL ERROR! Chebyshev Radau pt is not allowed for EOD'
      STOP
#endif

    CASE (3)
      call C_Ga_EOD_D (N, D_Even, D_Odd)

  END SELECT

  if (Max_Order_Of_D >  1)  &
    call C_EOD_High_Order (Point_Type, N, Max_Order_Of_D, D_Even, D_Odd)

  if (Filter_Choice == 0) then
    Smooth_1 = 0 ; Smooth_2 = 0
  else
    SELECT CASE (Point_Type)
      CASE (1)
        call C_GL_EOD_S_Even (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                              S_Even)

        call C_GL_EOD_S_Odd  (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                              S_Odd )

      CASE (3)
        call C_Ga_EOD_S_Even (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                              S_Even)

        call C_Ga_EOD_S_Odd  (Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                              S_Odd )

    END SELECT

    do O_Of_D = 1,Max_Order_Of_D
      call C_EOD_SDS (O_Of_D, N, Smooth_1, Smooth_2,  &
                      S_Even, S_Odd, D_Even(0,0,O_Of_D), D_Odd(0,0,O_Of_D))
    enddo
  endif

  END Subroutine C_Deriv_Matrix_EOD 
!
! ----------------------------------------------------------------------
!
  Subroutine C_EOD_SDS (O_Of_D, N, Smooth_1, Smooth_2,   &
                        S_Even, S_Odd, D_Even, D_Odd)

  integer  :: O_Of_D, N, Smooth_1, Smooth_2, N2

  REALTYPE, dimension(0:N/2,0:N/2) :: D_Even, D_Odd
  REALTYPE, dimension(0:N/2,0:N/2) :: S_Even, S_Odd

  N2 = N/2

  if (Smooth_1 == 1) then
      call C_EOD_DS (D_Even, S_Even)
      call C_EOD_DS (D_Odd , S_Odd )
  endif

  if (Smooth_2 == 1) then
    if (MOD(O_Of_D,2) == 0) then
      call C_EOD_SD (S_Even, D_Even)
      call C_EOD_SD (S_Odd , D_Odd )
    else
      call C_EOD_SD (S_Odd , D_Even)
      call C_EOD_SD (S_Even, D_Odd )
    endif
  endif

  CONTAINS

    Subroutine C_EOD_DS (D, S)

    integer                        :: i, j, k
    REALTYPE, dimension(0:N2,0:N2) :: D, S
    REALTYPE, dimension(0:N2)      :: Tmp

    do k = 0,N2

      do j = 0,N2
        Tmp(j) = ZERO
        do i = 0,N2
          Tmp(j) = Tmp(j) + D(k,i)*S(i,j)
        enddo
      enddo

      D(k,:) = Tmp
    enddo

    END Subroutine C_EOD_DS 
!
! ----------------------------------------------------------------------
!
    Subroutine C_EOD_SD (S, D)

    integer                        :: i, j, k
    REALTYPE, dimension(0:N2,0:N2) :: S, D
    REALTYPE, dimension(0:N2)      :: Tmp

    do j = 0,N2

      do k = 0,N2
        Tmp(k) = ZERO
        do i = 0,N2
          Tmp(k) = Tmp(k) + S(k,i)*D(i,j)
        enddo
      enddo

      D(:,j) = Tmp
    enddo

    END Subroutine C_EOD_SD 

  END Subroutine C_EOD_SDS
!
! -------------------------------------------------------------------------
!
  Subroutine C_EOD_High_Order (Point_Type, N, O_Of_D, D_Even, D_Odd)

  integer  :: Point_Type, N, M, O_Of_D, IGauss
  integer  :: j, k, l, N0, N2

  REALTYPE :: PiN, D_kk, D_kj, D_kl, x_kj, x_kl

  REALTYPE, dimension(0:N/2,0:N/2,O_Of_D) :: D_Even, D_Odd
  REALTYPE, dimension(-N:2*N)             :: Sin2
  REALTYPE, dimension(0:N/2)              :: x
    
#if defined (DEBUG)
  if (O_Of_D == 1) RETURN
#endif

  SELECT CASE (Point_Type)
    CASE (1)
      call C_GL_Point_Half (N, x)

      IGauss = 0

         PiN = PI/(2*N)

      Sin2(0) = ZERO
      do k = 1,N-1
        Sin2(k) = ONE/SIN(PiN*k)
      enddo
      Sin2(N) = ONE

      do k = 1,N
        Sin2(N+k) =  Sin2(N-k)
        Sin2( -k) = -Sin2(  k)
    enddo

    CASE (3)
      call C_Ga_Point_Half (N, x)

      IGauss = 1

         PiN = PI/(2*(N+1))

      do k = -N,-1
        Sin2(k) = ONE/SIN(PiN*k)
      enddo
      Sin2(0) = ZERO
      do k = 1,2*N + 1
        Sin2(k) = ONE/SIN(PiN*k)
      enddo

    END SELECT

                     N2 = N/2
                     N0 = N2
  if (MOD(N,2) == 0) N0 = N2-1

  do M = 2,O_Of_D

    do j = 0,N2
      do k = 0,N0
 
        if (k /= j) then
          l    = N-j

          x_kj = -HALF*Sin2(k+j+IGauss)*Sin2(k-j)
          x_kl = -HALF*Sin2(k+l+IGauss)*Sin2(k-l)

          D_kk =  HALF*(D_Even(k,k,M-1)+D_Odd(k,k,M-1))
          D_kj = x_kj*HALF*(D_Even(k,j,M-1)+D_Odd(k,j,M-1))
          D_kl = x_kl*HALF*(D_Even(k,j,M-1)-D_Odd(k,j,M-1))

          D_Even(k,j,M) = M*( D_kk*D_Even(k,j,1) - (D_kj + D_kl) )
          D_Odd (k,j,M) = M*( D_kk*D_Odd (k,j,1) - (D_kj - D_kl) )

        endif

      enddo
    enddo

    if (MOD(N,2) == 0) then
      k = N2
      do j = 0,N0
 
        l    = N-j

        x_kj = -HALF*Sin2(k+j+IGauss)*Sin2(k-j)
        x_kl = -HALF*Sin2(k+l+IGauss)*Sin2(k-l)

        D_kk =            D_Even(k,k,M-1)
        D_kj = x_kj*HALF*(D_Even(k,j,M-1)+D_Odd(k,j,M-1))
        D_kl = x_kl*HALF*(D_Even(k,j,M-1)-D_Odd(k,j,M-1))

        D_Even(k,j,M) = M*(D_kk*D_Even(k,j,1) - (D_kj + D_kl))
        D_Odd (k,j,M) = M*(D_kk*D_Odd (k,j,1) - (D_kj - D_kl))

      enddo
    endif

#if 1
    call Diagonal_Element_By_Sum_EOD (N, M, x, D_Even, D_Odd)
#else
    call Diagonal_Element_By_MXM_EOD (N, M,    D_Even, D_Odd)
#endif

    if (MOD(N,2) == 0) then
        D_Odd (:,N2,M) = ZERO

      if (MOD(M,2) == 0) then
        D_Odd (N2,:,M) = ZERO
      else
        D_Even(N2,:,M) = ZERO
      endif
    endif

  enddo

  END Subroutine C_EOD_High_Order 

END MODULE Chebyshev_Differentiation_EOD
