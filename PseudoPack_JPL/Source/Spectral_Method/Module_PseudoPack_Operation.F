#include "pseudopack.h"


MODULE PseudoPack_Operation
  USE PseudoPack_Differentiation, ONLY: Diff
  USE PseudoPack_Smoothing      , ONLY: Smooth

implicit NONE

INTERFACE PS_Diff
  MODULE PROCEDURE PS_Diff_1D_G
  MODULE PROCEDURE PS_Diff_2D_G
  MODULE PROCEDURE PS_Diff_1D
  MODULE PROCEDURE PS_Diff_2D
END INTERFACE 

INTERFACE PS_Smooth
  MODULE PROCEDURE PS_Smooth_1D_G
  MODULE PROCEDURE PS_Smooth_2D_G
  MODULE PROCEDURE PS_Smooth_1D
  MODULE PROCEDURE PS_Smooth_2D
END INTERFACE 

PRIVATE

PUBLIC  :: PS_Diff
PUBLIC  :: PS_Smooth

CONTAINS

  Subroutine PS_Diff_1D_G (LDY, D, y, dy, N_Vector, Order, Even_Or_Odd)
  
  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Order
  integer , OPTIONAL :: Even_Or_Odd
  integer            :: O_Of_D, Vector

  REALTYPE, dimension(*)     :: D
  REALTYPE, dimension(*)     :: y, dy

  Vector = 1 
  O_Of_D = 1 ; if (PRESENT(Order)   ) O_Of_D = Order

  call Diff (Vector, O_Of_D, y, dy, LDY, D, Even_Or_Odd)

  END Subroutine PS_Diff_1D_G
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Diff_2D_G (LDY, D, y, dy, N_Vector, Order, Even_Or_Odd)
  
  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Order
  integer , OPTIONAL :: Even_Or_Odd
  integer            :: O_Of_D, Vector

  REALTYPE, dimension(*)     :: D
  REALTYPE, dimension(LDY,*) :: y, dy

  Vector = 1 ; if (PRESENT(N_Vector)) Vector = N_Vector
  O_Of_D = 1 ; if (PRESENT(Order)   ) O_Of_D = Order

  call Diff (Vector, O_Of_D, y, dy, LDY, D, Even_Or_Odd)

  END Subroutine PS_Diff_2D_G
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Diff_1D   (D, y, dy, N_Vector, Order, Even_Or_Odd)
  
  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Order
  integer , OPTIONAL :: Even_Or_Odd
  integer            :: O_Of_D, Vector

  REALTYPE, dimension(*)     :: D
  REALTYPE, dimension(:)     :: y, dy

  LDY    = SIZE(y, DIM=1) 
  Vector = 1

  O_Of_D = 1 ; if (PRESENT(Order)   ) O_Of_D = Order

  call Diff (Vector, O_Of_D, y, dy, LDY, D, Even_Or_Odd)

  END Subroutine PS_Diff_1D
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Diff_2D    (D, y, dy, N_Vector, Order, Even_Or_Odd)
  
  integer            :: LDY, Index
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Order
  integer , OPTIONAL :: Even_Or_Odd
  integer            :: O_Of_D, Vector

  REALTYPE, dimension(*)     :: D
  REALTYPE, dimension(:,:)   :: y, dy

  Index  = D(_I_INDEX)

                  LDY    = SIZE(y, DIM=1) 
  if (Index == 1) Vector = SIZE(y, DIM=2)
  if (Index == 2) Vector = SIZE(y, DIM=1)

  if (PRESENT(N_Vector)) Vector = N_Vector

  O_Of_D = 1 ; if (PRESENT(Order)   ) O_Of_D = Order

  call Diff (Vector, O_Of_D, y, dy, LDY, D, Even_Or_Odd)

  END Subroutine PS_Diff_2D
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Smooth_1D_G (LDY, S, y, N_Vector, Even_Or_Odd)
  
  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Even_Or_Odd
  integer            :: Vector

  REALTYPE, dimension(*)     :: S
  REALTYPE, dimension(*)     :: y

  Vector = 1 

  call Smooth (Vector, y, LDY, S, Even_Or_Odd)

  END Subroutine PS_Smooth_1D_G
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Smooth_2D_G (LDY, S, y, N_Vector, Even_Or_Odd)
  
  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Even_Or_Odd
  integer            :: Vector

  REALTYPE, dimension(*)     :: S
  REALTYPE, dimension(LDY,*) :: y

  Vector = 1 ; if (PRESENT(N_Vector)) Vector = N_Vector

  call Smooth (Vector, y, LDY, S, Even_Or_Odd)

  END Subroutine PS_Smooth_2D_G
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Smooth_1D   (S, y, N_Vector, Even_Or_Odd)
  
  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Even_Or_Odd
  integer            :: Vector

  REALTYPE, dimension(*)     :: S
  REALTYPE, dimension(:)     :: y

  LDY    = SIZE(y, DIM=1) 
  Vector = 1

  call Smooth (Vector, y, LDY, S, Even_Or_Odd)

  END Subroutine PS_Smooth_1D
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Smooth_2D    (S, y, N_Vector, Even_Or_Odd)
  
  integer            :: LDY, Index
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Even_Or_Odd
  integer            :: Vector

  REALTYPE, dimension(*)     :: S
  REALTYPE, dimension(:,:)   :: y

  Index  = S(_I_INDEX)

                  LDY    = SIZE(y, DIM=1) 
  if (Index == 1) Vector = SIZE(y, DIM=2)
  if (Index == 2) Vector = SIZE(y, DIM=1)

  if (PRESENT(N_Vector)) Vector = N_Vector

  call Smooth (Vector, y, LDY, S, Even_Or_Odd)

  END Subroutine PS_Smooth_2D

END MODULE PseudoPack_Operation
