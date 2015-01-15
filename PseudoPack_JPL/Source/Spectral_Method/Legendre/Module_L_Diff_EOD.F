#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Diff_EOD
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! DESCRIPTION  : sets up the aux Diff structure for the (possibly Mapped)
!                Legendre Derivative using the Even-Odd algorithm.
!                Also computes the x that functions should be
!                sampled at.
!
!  Index      = 1 if Differentiating along the first index, 
!               2 if Differentiating along the second index. 
!
!  x_Left     = Left  Physical Domain Point
!  x_Right    = Right Physical Domain Point
!
!  N_Pnts     = length of the vectors
!  M          = total number of vectors.
!  ldy        = the leading dimension of the array to be processed. 
!  BlockSize  = number of vectors to be processed at one time.
!               if BlockSize=0, BlockSize=M  
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
!  if Map = -1, the mapping function is not an odd function.
!
! ----------------------------------------------------------------------

MODULE Legendre_Differentiation_EOD
  USE Legendre_Grid_Point
  USE Filtering

  USE EOD_Auxillary_Subroutine

  USE Grid_Generation, ONLY: Grid_And_Mapping

  USE Differentiation_EOD_1, L_Diff_EOD_1=>Diff_EOD_1
  USE Differentiation_EOD_2, L_Diff_EOD_2=>Diff_EOD_2

implicit NONE

PRIVATE

PUBLIC  :: L_Setup_Aux_Diff_EOD
PUBLIC  :: L_Deriv_Matrix_EOD
PUBLIC  :: L_Diff_EOD_1
PUBLIC  :: L_Diff_EOD_2

CONTAINS

  Subroutine L_Setup_Aux_Diff_EOD (Index, Method, Point_Type,          &
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

  N = N_Pnts-1 ; N0 = N ; N1 = N+1 ; N2 = N/2

  if (Symmetry == 1) N0 = N2

  i_m_s = EOD_METRIC_START(N1)

  call Grid_And_Mapping (Method, Point_Type,  &
                         N, N0, Map, Manual, alpha, beta,  &
                         x_Left, x_Right, x, Aux_Diff(i_m_s))

  i_e_s = EOD_MATRIX_START_EVEN(N1)
  i_o_s = EOD_MATRIX_START_ODD(N1)

  call L_Deriv_Matrix_EOD (Max_Order_Of_D, Point_Type,  &
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
        write ( 6,1000) ; write (99,1000) ; STOP
      endif
      if (Max_Order_Of_D >  1) then
        write ( 6,1100) ; write (99,1100) ; STOP
      endif
    endif
  endif

 1000 format (1x,'FATAL ERROR : '/  &
              1x,'  Non-Anti-Symmetric mapping metric (EOD_Code=1) '/  &
              1x,'  with (Symmetry=1) for Even-Odd Algorithm (Algorithm=1)')
 1100 format (1x,'FATAL ERROR : '/  &
              1x,'  Non-Anti-Symmetric mapping metric (EOD_Code=1) '/  &
              1x,'  with Max_Order_Of_Derivative > 1! Change to MXM')

  END Subroutine L_Setup_Aux_Diff_EOD 

! ----------------------------------------------------------------------
! FUNCTION NAME: L_Deriv_Matrix_EOD
! DESCRIPTION  : computes the elements of the two matrices used in the
!                Even-Odd algorithm for the Legendre  Smoothing. Does some
!                crufty Stuff to avoid using O(N^2) temporary storage.
!
!  N+1 = number of grid points
! ----------------------------------------------------------------------
  Subroutine L_Deriv_Matrix_EOD (Max_Order_Of_D, Point_Type,  &
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
      ASSERT(N>0, "ledg_der_m_eo:N<=0")
#endif

  call L_EOD_D (Point_Type, N, D_Even, D_Odd)

  if (Max_Order_Of_D >  1)  &
    call L_EOD_High_Order (Point_Type, N, Max_Order_Of_D, D_Even, D_Odd)

  if (Filter_Choice == 0) then
    Smooth_1 = 0 ; Smooth_2 = 0
  else
    call L_EOD_S_BOTH (Point_Type,   &
                       Filter_Choice, N, Mode_CutOff, Omega, Order,  &
                       S_Even, S_Odd)

    do O_Of_D = 1,Max_Order_Of_D
      call L_EOD_SDS (O_Of_D, N, Smooth_1, Smooth_2,  &
                      S_Even, S_Odd, D_Even(0,0,O_Of_D), D_Odd(0,0,O_Of_D))
    enddo
  endif

  END Subroutine L_Deriv_Matrix_EOD 
!
! ----------------------------------------------------------------------
!
  Subroutine L_EOD_D (Point_Type, N, D_Even, D_Odd)

  integer  :: Point_Type, N
  REALTYPE :: Product, D_kj, D_kl

  REALTYPE, dimension(0:N)         :: x, c
  REALTYPE, dimension(0:N/2,0:N/2) :: D_Even, D_Odd

  integer  :: j, k, l, N0, N2

  SELECT CASE (Point_Type)
    CASE (1)
      call L_GL_Point_Full (N, x)

#ifdef DEBUG
    CASE (2)
      write(6 ,*) 'FATAL ERROR! Legendre  Radau pt is not allowed for EOD'
      write(99,*) 'FATAL ERROR! Legendre  Radau pt is not allowed for EOD'
      STOP
#endif

    CASE (3)
      call L_Ga_Point_Full (N, x)

  END SELECT

                     N2 = N/2

                     N0 = N2
  if (MOD(N,2) == 0) N0 = N2-1

  do k = 0,N

    Product = ZERO
    do l = 0,N
      if (l /= k) Product = Product + LOG(ABS(x(k)-x(l)))
    enddo
    c(k) = Product

  enddo

!     Calculate off-diagonal elements.
  do j = 0,N2
    do k = 0,N2
      if (k /= j) then
        l = N-j

        D_kj = (-1)**(k+j)*EXP(c(k)-c(j))/(x(k)-x(j))
        D_kl = (-1)**(k+l)*EXP(c(k)-c(l))/(x(k)-x(l))

        D_Even(k,j) = D_kj + D_kl
        D_Odd (k,j) = D_kj - D_kl
      endif
    enddo
  enddo

  if (MOD(N,2) == 0) then
    D_Even(N2,N2) = ZERO

    D_Even(:,N2) = HALF*D_Even(:,N2)
    D_Odd (:,N2) = ZERO
  endif

! Diagonal:
  call Diagonal_Element_By_Sum_EOD (N, 1, x, D_Even, D_Odd)

  if (MOD(N,2) == 0) then
    D_Even(N2,:) = ZERO
  endif

  END Subroutine L_EOD_D
!
! ----------------------------------------------------------------------
!
  Subroutine L_EOD_S_BOTH (Point_Type,  &
                           Filter_choice, N, Mode_CutOff, Omega, Order, &
                           S_Even, S_Odd)

  integer  :: Point_Type, N, Filter_choice, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N)         :: Sigma
  REALTYPE, dimension(0:N/2,0:N/2) :: S_Even, S_Odd

  integer  :: N2

  N2 = N/2

  call Filter (Filter_choice, N, Mode_CutOff, Omega, Order, Sigma)

  call Smoothing_Matrix (2, Point_Type, 0, N, N2, N2, Sigma, S_Even)
  call Smoothing_Matrix (2, Point_Type, 1, N, N2, N2, Sigma, S_Odd )

  END Subroutine L_EOD_S_BOTH

!
! ----------------------------------------------------------------------
!
  Subroutine L_EOD_SDS (O_Of_D, N, Smooth_1, Smooth_2,   &
                        S_Even, S_Odd, D_Even, D_Odd)

  integer  :: O_Of_D, N, Smooth_1, Smooth_2

  REALTYPE, dimension(0:N/2,0:N/2) :: D_Even, D_Odd
  REALTYPE, dimension(0:N/2,0:N/2) :: S_Even, S_Odd

  integer  :: N2

  N2 = N/2

  if (Smooth_1 == 1) then
      call L_EOD_DS (D_Even, S_Even)
      call L_EOD_DS (D_Odd , S_Odd )
  endif

  if (Smooth_2 == 1) then
    if (MOD(O_Of_D,2) == 0) then
      call L_EOD_SD (S_Even, D_Even)
      call L_EOD_SD (S_Odd , D_Odd )
    else
      call L_EOD_SD (S_Odd , D_Even)
      call L_EOD_SD (S_Even, D_Odd )
    endif
  endif

  CONTAINS

    Subroutine L_EOD_DS (D, S)

    REALTYPE, dimension(0:N2,0:N2) :: D, S
    REALTYPE, dimension(0:N2)      :: Tmp

    integer  :: i, j, k

    do k = 0,N2

      do j = 0,N2
        Tmp(j) = ZERO
        do i = 0,N2
          Tmp(j) = Tmp(j) + D(k,i)*S(i,j)
        enddo
      enddo

      D(k,:) = Tmp
    enddo

    END Subroutine L_EOD_DS
!
! ----------------------------------------------------------------------
!
    Subroutine L_EOD_SD (S, D)

    REALTYPE, dimension(0:N2,0:N2) :: S, D
    REALTYPE, dimension(0:N2)      :: Tmp

    integer  :: i, j, k

    do j = 0,N2

      do k = 0,N2
        Tmp(k) = ZERO
        do i = 0,N2
          Tmp(k) = Tmp(k) + S(k,i)*D(i,j)
        enddo
      enddo

      D(:,j) = Tmp
    enddo

    END Subroutine L_EOD_SD

  END Subroutine L_EOD_SDS

#if 0
! ----------------------------------------------------------------------
! FUNCTION NAME: L_Smooth_Matrix
! DESCRIPTION  : Construct the Legendre Filter Matrix corresponding
!                to the array of Filter coefficients Sigma for all points type.
!
!                Point_Type = 1, Lobatto
!                Point_Type = 2, Radua
!                Point_Type = 3, Gauss
!
!                if S_Type = 0, make S_Even
!                   S_Type = 1, make S_Odd
!                   S_Type = 2, make regular S.
!
!                S is of size (0:N0,0:N1)
!                N+1 is total number of grid points.
!                Sigma is the filter of size (0:N)
!
! ----------------------------------------------------------------------

  Subroutine L_Smooth_Matrix (Point_Type, S_Type, N, N0, N1, Sigma, S)

  integer  :: Point_Type, N, N0, N1, S_Type
  integer  :: Kind, Nfix, IERR
  logical  :: Symmetry, Even, Odd
  REALTYPE :: gamma, Ratio_1, Ratio_2

  REALTYPE, dimension(0:N0,0:N1) :: S
  REALTYPE, dimension(0:N)       :: Sigma

  REALTYPE, dimension(2)   :: qxfix
  REALTYPE, dimension(0:N) :: x, w, Work, Scale

  REALTYPE, dimension(0:N) :: L_Prev, L_Now, L_Next

  integer  :: i, j, k, N2

  Kind = 1 ! Legendre, w(x) = 1, (-1,1)
! Kind = 2 ! Chebyshev of the first kind, w(x) = 1/sqrt(1 - x*x), (-1,1)
! Kind = 3 ! Chebyshev of the second Kind, w(x) = sqrt(1 - x*x),  (-1,1)
! Kind = 4 ! Hermite, w(x) = exp(-x*x), (-inf,inf)
! Kind = 5 ! Jacobi, w(x) = (1-x)**alpha*(1+x)**beta, (-1,1), alpha, beta > -1
! Kind = 6 ! generalized Laguerre, w(x) = exp(-x)*x**alpha, (0, inf), alpha>-1

  Even     =       (S_Type == 0)
  Odd      =       (S_Type == 1)
  Symmetry = .NOT. (S_Type == 2)

  SELECT CASE (Point_Type)
    CASE (1)
      Nfix =  2 ; qxfix(1) = -ONE ; qxfix(2) =  ONE

    CASE (2)
      Nfix =  1 ; qxfix(1) = -ONE

      Symmetry = .FALSE.

    CASE (3)
      Nfix =  0

  END SELECT

  call GaussQ (Kind, N+1, ZERO, ZERO, Nfix, qxfix, work, x, w, IERR)

  if (IERR /= 0) then
    write (6 ,*) "L_Smooth_Matrix: Gauss Quadrature failed. Quitting."
    write (99,*) "L_Smooth_Matrix: Gauss Quadrature failed. Quitting."
    STOP
  endif

  if (Point_Type /= 2) then
    N2 = N/2

    if (MOD(N,2) == 0) x(N2) = ZERO

    do k = 0,N2
      x(N-k) = -x(k)
      w(N-k) =  w(k)
    enddo
  endif

  x = -x

  do k = 0,N
    gamma = ONE/(k+HALF) ; Scale(k) = Sigma(k)/gamma
  enddo

  if (Point_Type == 1) then
    gamma = TWO/N        ; Scale(N) = Sigma(N)/gamma
  endif

  L_Prev = ONE ; L_Now = x ; S = ZERO

  do i = 0,N
    if ((.NOT. Symmetry) .OR. (MOD(i,2) == S_Type)) then
      do j = 0,N1
        S(:,j) = S(:,j) + Scale(i)*L_Prev(0:N0)*L_Prev(j)
      enddo
    endif

    Ratio_1 = (2*i+3)/(i+TWO)
    Ratio_2 = (  i+1)/(i+TWO)

    L_Next = Ratio_1*x*L_Now - Ratio_2*L_Prev
    L_Prev = L_Now
    L_Now  = L_Next
  enddo

  do j = 0,N1
    S(:,j) = w(j)*S(:,j)
  enddo

  if (Symmetry) then
    S(:,0:N0) = TWO*S(:,0:N0)

    if (MOD(N,2) == 0) then
      if (Even) then
        S(:,N0) = HALF*S(:,N0)
      else
        S(:,N0) = ZERO
     endif
   endif

  endif

  END Subroutine L_Smooth_Matrix
#endif

!
!--------------------------------------------------------------------------
!
  Subroutine L_EOD_High_Order (Point_Type, N, Order, D_Even, D_Odd)

  integer  :: Point_Type, N, M, Order
  REALTYPE :: D_kk, D_kj, D_kl, x_kj, x_kl

  REALTYPE, dimension(0:N/2,0:N/2,Order) :: D_Even, D_Odd
  REALTYPE, dimension(0:N/2)             :: x

  integer  :: j, k, N0, N2

#if defined (DEBUG)
  if (Order == 1) RETURN
#endif

  SELECT CASE (Point_Type)
    CASE (1)
      call L_GL_Point_Half (N, x)

    CASE (3)
      call L_Ga_Point_Half (N, x)

  END SELECT

                     N2 = N/2

                     N0 = N2
  if (MOD(N,2) == 0) N0 = N2-1

  do M = 2,Order

    do j = 0,N2
      do k = 0,N0

        if (k /= j) then

          x_kj = ONE/(x(k)-x(j))
          x_kl = ONE/(x(k)+x(j))

          D_kk =      HALF*(D_Even(k,k,M-1)+D_Odd(k,k,M-1))
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

          x_kj = ONE/(x(k)-x(j))
          x_kl = ONE/(x(k)+x(j))

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

  END Subroutine L_EOD_High_Order 

END MODULE Legendre_Differentiation_EOD
