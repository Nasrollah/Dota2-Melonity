#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Diff_EOD
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa  (Rio)
! DESCRIPTION  : Sets up the Aux_Diff structure for the (possibly Mapped)
!                Fourier Differentiation using the Even-Odd algorithm.
!                Also computes the x that functions should be sampled at.
! ----------------------------------------------------------------------

MODULE Fourier_Differentiation_EOD
  USE Fourier_GL_EOD
  USE EOD_Auxillary_Subroutine

  USE Grid_Generation, ONLY: Grid_And_Mapping

  USE Differentiation_EOD_1, F_Diff_EOD_1=>Diff_EOD_1
  USE Differentiation_EOD_2, F_Diff_EOD_2=>Diff_EOD_2

implicit NONE

PRIVATE

PUBLIC  :: F_Setup_Aux_Diff_EOD
PUBLIC  :: F_Deriv_Matrix_EOD
PUBLIC  :: F_Diff_EOD_1
PUBLIC  :: F_Diff_EOD_2

CONTAINS

  Subroutine F_Setup_Aux_Diff_EOD (Index, Method, Point_Type,          &
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
  integer  :: O_Of_D
  integer  :: EOD_Code
  REALTYPE :: x_Left, x_Right
  REALTYPE :: alpha, beta, Omega, Order

  REALTYPE, dimension(:) :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: N, N0, N1, N2
  integer  :: i_m_s, i_e_s, i_o_s, i_Offset

  N  = N_Pnts

  if (Symmetry == 0) then
    N0 = N-1 ; N1 = N+1 ; N2 = N/2
  else
    N0 = N/2 ; N1 = N+1 ; N2 = N/2
  endif

  i_m_s = EOD_METRIC_START(N1)

  call Grid_And_Mapping (Method, Point_Type,  &
                         N, N0, Map, Manual, alpha, beta,  &
                         x_Left, x_Right, x, Aux_Diff(i_m_s))

  i_e_s = EOD_MATRIX_START_EVEN(N1)
  i_o_s = EOD_MATRIX_START_ODD(N1)

  Aux_Diff(_I_MATRIX_START)     = i_e_s
  Aux_Diff(_I_ODD_MATRIX_START) = i_o_s
  Aux_Diff(_I_SIGMA_START)      = 0
  Aux_Diff(_I_METRIC_START)     = i_m_s

  SELECT CASE (Point_Type)
    CASE (1,3)
      call F_Deriv_Matrix_EOD (Max_Order_Of_D, Filter_Choice, N, Mode_CutOff, &
                               Omega, Order, Smooth_1, Smooth_2, &
                               Aux_Diff(i_e_s), Aux_Diff(i_o_s))

   CASE (2)
     write(6 ,*) 'No EOD for the Fourier Gauss Radau points!'
     write(99,*) 'No EOD for the Fourier Gauss Radau points!'
     STOP

  END SELECT

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

      call Apply_Mapping_EOD (Map, N0, N2, Max_Order_Of_D,   &
                              Aux_Diff(i_m_s),               &
                              Aux_Diff(i_e_s), Aux_Diff(i_o_s))
    else
      if (Symmetry == 1) then
        write (6 ,1000) ; write (99,1000) ; STOP
      endif
    endif
  endif

 1000 format (1x,'FATAL ERROR : '/  &
              1x,'  Non-Anti-Symmetry mapping metric (EOD_Code=1) '/  &
              1x,'  with (Symmetry=1) for Even-Odd Algorithm (Algorithm=1)')

  END Subroutine F_Setup_Aux_Diff_EOD 

! ----------------------------------------------------------------------
! FUNCTION NAME: F_Deriv_Matrix_EOD
! DESCRIPTION  : computes the elements of the two matrices used in the 
!                Even-Odd algorithm for the Fourier Smoothing. 
! ----------------------------------------------------------------------
  Subroutine F_Deriv_Matrix_EOD (Max_Order_Of_D, Filter_Choice,     &
                                 N, Mode_CutOff,                    &
                                 Omega, Order, Smooth_1, Smooth_2,  &
                                 D_Even, D_Odd)

  integer  :: Max_Order_Of_D
  integer  :: Filter_Choice, N, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  REALTYPE, dimension(*) :: D_Even, D_Odd

#if defined(DEBUG)
  ASSERT(N>0,"F_DER_M_EOD:N<=0")
  ASSERT(Max_Order_Of_D>0,"F_DER_M_EOD:Max_Order_Of_D<=0")
  ASSERT(Filter_Choice>=0,"F_DER_M_EOD:Filter_Choice<0")
#endif

  if (Filter_Choice == 0) then 
    call F_GL_EOD_D  (N, Max_Order_Of_D, D_Even, D_Odd)

    Smooth_1 = 0 ; Smooth_2 = 0
  else
    call F_GL_EOD_DS (Max_Order_Of_D, Filter_Choice, N, Mode_CutOff,  &
                                      Omega, Order, D_Even, D_Odd)
  endif

  END Subroutine F_Deriv_Matrix_EOD 

END MODULE Fourier_Differentiation_EOD
