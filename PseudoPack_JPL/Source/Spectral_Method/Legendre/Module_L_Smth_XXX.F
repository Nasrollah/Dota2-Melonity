#include "pseudopack.h"


MODULE Legendre_Smoothing
  USE Filtering
  USE EOD_Auxillary_Subroutine

  USE Smoothing_MXM_1, L_Smooth_MXM_1=>Smooth_MXM_1
  USE Smoothing_MXM_2, L_Smooth_MXM_2=>Smooth_MXM_2

  USE Smoothing_EOD_1, L_Smooth_EOD_1=>Smooth_EOD_1
  USE Smoothing_EOD_2, L_Smooth_EOD_2=>Smooth_EOD_2

implicit NONE

PRIVATE

PUBLIC  :: L_Setup_Aux_Smth_MXM
PUBLIC  :: L_Smooth_MXM_1
PUBLIC  :: L_Smooth_MXM_2

PUBLIC  :: L_Setup_Aux_Smth_EOD
PUBLIC  :: L_Smooth_EOD_1
PUBLIC  :: L_Smooth_EOD_2

CONTAINS

! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------

  Subroutine L_Setup_Aux_Smth_MXM (Index, Method, Point_Type,          &
                                   Algorithm, Symmetry,                &
                                   N_Pnts, M, ldy, BlockSize,          &
                                   Filter_Choice, Mode_CutOff,         &
                                   Omega, Order, Aux_Smth)

  integer  :: Index, Method, Point_Type, Algorithm, Symmetry
  integer  :: N_Pnts, M, ldy, BlockSize
  integer  :: Filter_Choice, Mode_CutOff
  REALTYPE :: Omega, Order
  integer  :: N, N1
  integer  :: i_d_s

  REALTYPE, dimension(*) :: Aux_Smth

  N = N_Pnts-1 ; N1 = N+1 

  i_d_s = MXM_MATRIX_START(N1)

  call L_Smooth_Matrix_MXM (Aux_Smth(i_d_s))
  
  CONTAINS

    Subroutine L_Smooth_Matrix_MXM (S)

    REALTYPE, dimension(0:N)     :: Sigma
    REALTYPE, dimension(0:N,0:N) :: S

    call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)

    call Smoothing_Matrix (2, Point_Type, 2, N, N, N, Sigma, S)

    END Subroutine L_Smooth_Matrix_MXM 

  END Subroutine L_Setup_Aux_Smth_MXM 
!
! ----------------------------------------------------------------------
!
  Subroutine L_Setup_Aux_Smth_EOD (Index, Method, Point_Type,          &
                                   Algorithm, Symmetry,                &
                                   N_Pnts, M, ldy, BlockSize,          &
                                   Filter_Choice, Mode_CutOff,         &
                                   Omega, Order, Aux_Smth)

  integer  :: Index, Method, Point_Type, Algorithm, Symmetry
  integer  :: N_Pnts, M, ldy, BlockSize
  integer  :: Filter_Choice, Mode_CutOff
  REALTYPE :: Omega, Order

  REALTYPE, dimension(*) :: Aux_Smth

  integer  :: N, N1, N2
  integer  :: i_e_s, i_o_s

  N = N_Pnts-1 ; N1 = N+1 ; N2 = N/2

  i_e_s = EOD_MATRIX_START_EVEN(N1)
  i_o_s = EOD_MATRIX_START_EVEN(N1) + EOD_HALF_SIZE(N1)

  SELECT CASE (Point_Type)
    CASE (1,3)
      call L_EOD_S_BOTH (Aux_Smth(i_e_s), Aux_Smth(i_o_s))

#ifdef DEBUG
    CASE (2)
      write(6 ,*) 'FATAL ERROR! Legendre  Radau pt is not allowed for EOD'
      write(99,*) 'FATAL ERROR! Legendre  Radau pt is not allowed for EOD'
      STOP
#endif

  END SELECT

  if (Symmetry == 0)  &
    call Adjust_EOD_Matrix (N2, Aux_Smth(i_e_s), Aux_Smth(i_o_s))

  CONTAINS

    Subroutine L_EOD_S_BOTH (S_Even, S_Odd)

    REALTYPE, dimension(0:N)       :: Sigma
    REALTYPE, dimension(0:N2,0:N2) :: S_Even, S_Odd

    call Filter (Filter_choice, N, Mode_CutOff, Omega, Order, Sigma)

    call Smoothing_Matrix (2, Point_Type, 0, N, N2, N2, Sigma, S_Even)
    call Smoothing_Matrix (2, Point_Type, 1, N, N2, N2, Sigma, S_Odd )

    END Subroutine L_EOD_S_BOTH

  END Subroutine L_Setup_Aux_Smth_EOD 

END MODULE Legendre_Smoothing
