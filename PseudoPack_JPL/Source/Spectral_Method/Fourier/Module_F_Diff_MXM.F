#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! FUNCTION NAME: F_Setup_Aux_Diff_MXM
! DESCRIPTION  : sets up the aux Diff structure for the (possibly Mapped)
!                Fourier Differentiation using the plain MXM algorithm.
!                Also computes the x that functions should be
!                sampled at.
! 
!  Index      = 1 if Differentiating along the first index, 
!               2 if Differentiating along the second index.
!
!  N_Pnts     = length of the vectors
!  M          = total number of vectors.
!  ldy        = the leading dimension of the array to be processed. (NOT USED)
!  BlockSize  = number of vectors to be processed at one time.
!               if BlockSize=0, BlockSize=M  (NOT USED HERE)
!  Symmetry   = 0 if no symmetry assumed 
!               1 if data set is symmetric/anti-symmetric about mid-points.
!
!  Point_Type = 1 Lobatto Point
!             = 2 Radau   Point
!             = 3 Gauss   Point
!
! For Smoothing the Differentiation...
!
!     call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)
!
!  Filter_Choice = Choice of Filtering function
!                    1. Exponential Filter (Mode_CutOff, Omega, Order)
!                    2. Lanczos Filter (2nd order)
!                    3. Raised Cosine Filter (2nd order)
!                    4. 8th order Sharpened Rasied Cosine Filter
!  Mode_CutOff   = Cut-Off Frequency for the Exponential Filter
!                  if equal N_Pnts-1, Filter_Choice=0
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


MODULE Fourier_Differentiation_MXM
  USE MXM_Auxillary_Subroutine

  USE Fourier_GL_MXM
  USE Fourier_GR_MXM

  USE Grid_Generation, ONLY: Grid_And_Mapping

  USE Differentiation_MXM_1, F_Diff_MXM_1=>Diff_MXM_1
  USE Differentiation_MXM_2, F_Diff_MXM_2=>Diff_MXM_2

implicit NONE

PRIVATE

PUBLIC  :: F_Setup_Aux_Diff_MXM
PUBLIC  :: F_Deriv_Matrix_MXM
PUBLIC  :: F_Diff_MXM_1
PUBLIC  :: F_Diff_MXM_2

CONTAINS

  Subroutine F_Setup_Aux_Diff_MXM (Index, Method, Point_Type,          &
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

#if defined (DEBUG)
  ASSERT(N_Pnts> 0,"F_Der_M:N<=0")
  ASSERT((Point_Type<=3),"F_Set_a_d_m:Point_Type>3")
  ASSERT((Point_Type>=1),"F_Set_a_d_m:Point_Type<1")

  if ((Point_Type == 1) .AND. (MOD(N_Pnts,2) /= 0)) then
    write (6 ,*) "F_Set_a_d_m:Point_Type=1,N_Pnts=Odd"
    write (99,*) "F_Set_a_d_m:Point_Type=1,N_Pnts=Odd"
    STOP
  endif
      
  if ((Point_Type == 2) .AND. (MOD(N_Pnts,2) /= 1)) then
    write (6 ,*) "F_Set_a_d_m:Point_Type=2,N_Pnts=Even"
    write (99,*) "F_Set_a_d_m:Point_Type=2,N_Pnts=Even"
    STOP
  endif

  if ((Point_Type == 3) .AND. (MOD(N_Pnts,2) /= 0)) then
    write (6 ,*) "F_Set_a_d_m:Point_Type=1,N_Pnts=Odd"
    write (99,*) "F_Set_a_d_m:Point_Type=1,N_Pnts=Odd"
    STOP
  endif
      
#endif

  N = N_Pnts   ; N0 = N-1 ; N1 = N

  i_m_s = MXM_METRIC_START(N1)

  call Grid_And_Mapping (Method, Point_Type,                  &
                         N, N0, Map, Manual, alpha, beta,     &
                         x_Left, x_Right, x, Aux_Diff(i_m_s))

  i_d_s = MXM_MATRIX_START(N1)

  call F_Deriv_Matrix_MXM (Max_Order_of_D, Point_Type,        &
                           Filter_Choice, N, Mode_CutOff,     &
                           Omega, Order, Smooth_1, Smooth_2,  &
                           Aux_Diff(i_d_s)) 

  if (Map /= 0) call Apply_Mapping (Map, N0, Max_Order_of_D,  &
                                    Aux_Diff(i_m_s), Aux_Diff(i_d_s))

  Aux_Diff(_I_MATRIX_START)     = i_d_s
  Aux_Diff(_I_ODD_MATRIX_START) = i_d_s
  Aux_Diff(_I_SIGMA_START)      = 0
  Aux_Diff(_I_METRIC_START)     = i_m_s

  EOD_Code = 0

  END Subroutine F_Setup_Aux_Diff_MXM 

! ----------------------------------------------------------------------
! FUNCTION NAME: F_Deriv_Matrix
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa  (Rio)
! DESCRIPTION  : computes the elements of the (unMapped) Fourier First
!                Differentiation Matrix using trig identities and flipping.
!                with option for smoothing.
!
!  N = number of grid points
! ----------------------------------------------------------------------
  Subroutine F_Deriv_Matrix_MXM (Max_Order_Of_D, Point_Type,    &
                                 Filter_Choice, N, Mode_CutOff, &
                                 Omega, Order, Smooth_1, Smooth_2, D)

  integer  :: Max_Order_Of_D, N
  integer  :: Point_Type
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(0:N-1,0:N-1,Max_Order_Of_D) :: D

#if defined(DEBUG)
  ASSERT(N>0,"F_Derv_MXM:N<=0")
  ASSERT((Point_Type>=1),"F_Derv_MXM:Point_Type<1")
  ASSERT((Point_Type<=3),"F_Derv_MXM:Point_Type>3")
#endif

  if (Filter_Choice == 0) then
    SELECT CASE (Point_TYPE)
      CASE (1,3)
        call F_GL_MXM_D (N, Max_Order_of_D, D)

      CASE (2)
        call F_GR_MXM_D (N, Max_Order_of_D, D)

    END SELECT

    Smooth_1 = 0 ; Smooth_2 = 0
  else
    SELECT CASE (Point_TYPE)
      CASE (1,3)
        call F_GL_MXM_DS (Max_Order_Of_D, Filter_Choice, N,  Mode_CutOff,  &
                                          Omega, Order, D)

      CASE (2)
        call F_GR_MXM_DS (Max_Order_Of_D, Filter_Choice, N,  Mode_CutOff,  &
                                          Omega, Order, D)

    END SELECT
  endif

  END Subroutine F_Deriv_Matrix_MXM

END MODULE Fourier_Differentiation_MXM
