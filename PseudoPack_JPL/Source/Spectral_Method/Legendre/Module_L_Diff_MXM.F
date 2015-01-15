#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Diff_MXM
! DESCRIPTION  : sets up the aux Diff structure for the (possibly Mapped)
!                Legendre Derivative using the plain MXM algorithm.
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

MODULE Legendre_Differentiation_MXM
  USE Legendre_Grid_Point

  USE Legendre_GL_MXM
  USE Legendre_GR_MXM
  USE Legendre_Ga_MXM
  USE MXM_Auxillary_Subroutine

  USE Grid_Generation, ONLY: Grid_And_Mapping

  USE Differentiation_MXM_1, L_Diff_MXM_1=>Diff_MXM_1
  USE Differentiation_MXM_2, L_Diff_MXM_2=>Diff_MXM_2

implicit NONE

PRIVATE

PUBLIC  :: L_Setup_Aux_Diff_MXM
PUBLIC  :: L_Deriv_Matrix_MXM
PUBLIC  :: L_Diff_MXM_1
PUBLIC  :: L_Diff_MXM_2

CONTAINS

  Subroutine L_Setup_Aux_Diff_MXM (Index, Method, Point_Type,          &
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
  REALTYPE :: x_Left, x_Right
  REALTYPE :: alpha, beta, Omega, Order

  REALTYPE, dimension(:) :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: N, N0, N1
  integer  :: i_m_s, i_d_s

  N = N_Pnts-1 ; N0 = N ; N1 = N+1

  i_m_s = MXM_METRIC_START(N1)

  call Grid_And_Mapping (Method, Point_Type,  &
                         N, N0, Map, Manual, alpha, beta,  &
                         x_Left, x_Right, x, Aux_Diff(i_m_s))

  i_d_s = MXM_MATRIX_START(N1)

  call L_Deriv_Matrix_MXM (Max_Order_Of_D, Point_Type,   &
                           Filter_Choice, N, Mode_CutOff,   &
                           Omega, Order, Smooth_1, Smooth_2,  &
                           Aux_Diff(i_d_s))

  if (Map /= 0) call Apply_Mapping (Map, N, Max_Order_Of_D,  &
                                    Aux_Diff(i_m_s), Aux_Diff(i_d_s))

  Aux_Diff(_I_MATRIX_START)     = i_d_s
  Aux_Diff(_I_ODD_MATRIX_START) = i_d_s
  Aux_Diff(_I_SIGMA_START)      = 0
  Aux_Diff(_I_METRIC_START)     = i_m_s

  EOD_Code = 0

  END Subroutine L_Setup_Aux_Diff_MXM

! ----------------------------------------------------------------------
! FUNCTION NAME: L_Deriv_Matrix
! DESCRIPTION  : computes the elements of the (unMapped) Legendre  First
!                Derivative Matrix using trig identities and flipping.
!                with option for smoothing.
!
!  N+1 = number of grid points
! ----------------------------------------------------------------------
  Subroutine L_Deriv_Matrix_MXM (Max_Order_Of_D, Point_Type,    &
                                 Filter_Choice, N, Mode_CutOff, &
                                 Omega, Order, Smooth_1, Smooth_2, D)

  integer  :: Max_Order_Of_D, O_Of_D, N
  integer  :: Point_Type
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N,0:N,Max_Order_Of_D) :: D

#if defined (DEBUG)
  ASSERT(N>0,"L_Derv_MXM:N<=0")
  ASSERT(Point_Type<=3,"L_Derv_MXM:Point_Type>3")
  ASSERT(Point_Type>=1,"L_Derv_MXM:Point_Type<1")
#endif

  SELECT CASE (Point_Type)
    CASE (1)
      call L_GL_MXM_D (N, D)

    CASE (2)
      call L_GR_MXM_D (N, D)

    CASE (3)
      call L_Ga_MXM_D (N, D)

  END SELECT

  if (Max_Order_Of_D >  1)  &
    call L_MXM_High_Order (Point_Type, N, Max_Order_Of_D, D)

  if (Filter_Choice == 0) then
    Smooth_1 = 0 ; Smooth_2 = 0
  else
    do O_Of_D = 1,Max_Order_Of_D
      call L_MXM_SDS (D(0,0,O_Of_D))
    enddo
  endif

  CONTAINS

    Subroutine L_MXM_SDS (D)

    REALTYPE, dimension(0:N,0:N) :: D

    SELECT CASE (Point_Type)
      CASE (1)
        call L_GL_MXM_SDS (O_Of_D, Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2, D)

      CASE (2)
        call L_GR_MXM_SDS (O_Of_D, Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2, D)

      CASE (3)
        call L_Ga_MXM_SDS (O_Of_D, Filter_Choice, N, Mode_CutOff,  &
                           Omega, Order, Smooth_1, Smooth_2, D)

    END SELECT

    END Subroutine L_MXM_SDS

  END Subroutine L_Deriv_Matrix_MXM
!
! ----------------------------------------------------------------------
!
  Subroutine L_MXM_High_Order (Point_Type, N, Order, D)

  integer  :: Point_Type, N, Order, M, N2
  logical  :: Symmetry

  REALTYPE, dimension(0:N)           :: x
  REALTYPE, dimension(0:N,0:N,Order) :: D

  integer  :: j, k

#if defined (DEBUG)
  if (Order == 1) RETURN
#endif

  SELECT CASE (Point_Type)
    CASE (1)
      call L_GL_Point_Full (N, x)

      Symmetry = .TRUE.

    CASE (2)
      call L_GR_Point_Full (N, x)

      Symmetry = .FALSE.

    CASE (3)
      call L_Ga_Point_Full (N, x)

      Symmetry = .TRUE.
  END SELECT

  if (Symmetry) then ; N2 = N/2 ; else ; N2 = N ; endif

  do M = 2, Order

    do j = 0,N
      do k = 0,N2
        if (k == j) CYCLE

        D(k,j,M) = M*(D(k,k,M-1)*D(k,j,1)-D(k,j,M-1)/(x(k)-x(j)))
      enddo
    enddo

#if 1
    call Diagonal_Element_By_Sum_MXM (Symmetry, N, M, D)
#else
    call Diagonal_Element_By_MXM_MXM (Symmetry, N, M, D)
#endif
  enddo

    call Centro_Symmetry_Property    (Symmetry, -1, N, Order, D)

  END Subroutine L_MXM_High_Order  

END MODULE Legendre_Differentiation_MXM
