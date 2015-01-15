#include "pseudopack.h"


MODULE MPI_Differentiation

  USE Processor
  USE PseudoPack_Utility
  USE MPI_Data_Distribution
#if defined (PARALLEL_MPI)
# if defined (PESSL) || defined (SCALAPACK)
  USE MPI_Transpose_PESSL, PS_P_Transpose => P_Transpose
# endif
#endif
  USE PseudoPack_Operation

implicit NONE

INTERFACE PS_MPI_Diff
  MODULE PROCEDURE PS_MPI_Diff_1D_G
  MODULE PROCEDURE PS_MPI_Diff_2D_G
  MODULE PROCEDURE PS_MPI_Diff_1D
  MODULE PROCEDURE PS_MPI_Diff_2D
END INTERFACE

PRIVATE

PUBLIC  :: PS_MPI_Adjust_Grid_Index
PUBLIC  :: PS_MPI_Diff

CONTAINS
!
!=======================================================================
!
  Subroutine PS_MPI_Adjust_Grid_Index (Property, Grid_Index)

  TYPE (PS_Property)   :: Property
  TYPE (PS_Grid_Index) :: Grid_Index

#if defined (PARALLEL_MPI)
  integer  :: N, M

  if ((N_Processor == 1) .OR. (Property%Index /= 2)) RETURN

  N = Grid_Index%N
  M = Grid_Index%M

  Grid_Index%N   = M
  Grid_Index%M   = N
  Grid_Index%LDY = M
#endif

  END Subroutine PS_MPI_Adjust_Grid_Index
!
!=======================================================================
!
  Subroutine PS_MPI_Diff_1D_G (LDY, D, f, D_f, N_Vector, Order, Even_Or_Odd)

  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Order
  integer , OPTIONAL :: Even_Or_Odd

  REALTYPE, dimension(*)                :: D
  REALTYPE, dimension(*)                :: f, D_f

  call PS_Diff (LDY, D, f, D_f, N_Vector, Order, Even_Or_Odd)

  END Subroutine PS_MPI_Diff_1D_G 
!
!=======================================================================
!
  Subroutine PS_MPI_Diff_2D_G (LDY, D, f, D_f, N_Vector, Order, Even_Or_Odd)

  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Order
  integer , OPTIONAL :: Even_Or_Odd

  REALTYPE, dimension(*)                :: D
  REALTYPE, dimension(LDY,*)            :: f, D_f

#if defined (PARALLEL_MPI)
  REALTYPE, dimension(:,:), ALLOCATABLE :: g, u, D_u

  integer  :: N, n_Remain, n_0, n_Star
  integer  :: M, m_Remain, m_0, m_Star
  integer  :: NV
#endif

  SELECT CASE (OPERATOR_RANK(D))
    CASE DEFAULT
      call PS_Diff (LDY, D, f, D_f, N_Vector, Order, Even_Or_Odd)

    CASE (2)
#if defined (PARALLEL_MPI)
      if (N_Processor == 1) then
        call PS_Diff (LDY, D, f, D_f, N_Vector, Order, Even_Or_Odd)

        RETURN
      endif

      N = OPERATOR_N_VECTORS(D)

      call PS_Get_Dimension (D, M)

      call PS_MPI_Data_Distribution (M, m_Remain, m_0, m_Star)
      call PS_MPI_Data_Distribution (N, n_Remain, n_0, n_Star)

      NV = n_Star ; if (PRESENT(N_Vector)) NV = MIN(N_Vector, n_Star)

      ALLOCATE (g(1:N,1:m_Star), u(1:M,1:n_Star), D_u(1:M,1:n_Star))

      g = f(1:N,1:m_Star)

      call PS_P_Transpose (N, M, g, u)

      call PS_Diff (D, u, D_u, NV, Order, Even_Or_Odd)

      call PS_P_Transpose (M, N, D_u, g)

      NV = N ; if (PRESENT(N_Vector)) NV = N_Vector

      D_f(1:NV,1:m_Star) = g(1:NV,:)

      DEALLOCATE (D_u, u, g)
#else
      call PS_Diff (LDY, D, f, D_f, N_Vector, Order, Even_Or_Odd)
#endif

  END SELECT

  END Subroutine PS_MPI_Diff_2D_G
!
!=======================================================================
!
  Subroutine PS_MPI_Diff_1D   (D, f, D_f, N_Vector, Order, Even_Or_Odd)

  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Order
  integer , OPTIONAL :: Even_Or_Odd

  REALTYPE, dimension(*)                :: D
  REALTYPE, dimension(:)                :: f, D_f

  call PS_Diff (D, f, D_f, N_Vector, Order, Even_Or_Odd)

  END Subroutine PS_MPI_Diff_1D 
!
!=======================================================================
!
  Subroutine PS_MPI_Diff_2D   (D, f, D_f, N_Vector, Order, Even_Or_Odd)

  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Order
  integer , OPTIONAL :: Even_Or_Odd

  REALTYPE, dimension(*)                :: D
  REALTYPE, dimension(:,:)              :: f, D_f

#if defined (PARALLEL_MPI)
  REALTYPE, dimension(:,:), ALLOCATABLE :: g, u, D_u

  integer  :: N, n_Remain, n_0, n_Star
  integer  :: M, m_Remain, m_0, m_Star
  integer  :: NV
#endif

  SELECT CASE (OPERATOR_RANK(D))
    CASE DEFAULT
      call PS_Diff (D, f, D_f, N_Vector, Order, Even_Or_Odd)

    CASE (2)
#if defined (PARALLEL_MPI)
      if (N_Processor == 1) then
        call PS_Diff (D, f, D_f, N_Vector, Order, Even_Or_Odd)

        RETURN
      endif

      N = OPERATOR_N_VECTORS(D)

      call PS_Get_Dimension (D, M)

      call PS_MPI_Data_Distribution (M, m_Remain, m_0, m_Star)
      call PS_MPI_Data_Distribution (N, n_Remain, n_0, n_Star)

      NV = n_Star ; if (PRESENT(N_Vector)) NV = MIN(N_Vector, n_Star)

      ALLOCATE (g(1:N,1:m_Star), u(1:M,1:n_Star), D_u(1:M,1:n_Star))

      g = f(1:N,1:m_Star)

      call PS_P_Transpose (N, M, g, u)

      call PS_Diff        (D, u, D_u, NV, Order, Even_Or_Odd)

      call PS_P_Transpose (M, N, D_u, g)

      NV = N ; if (PRESENT(N_Vector)) NV = N_Vector

      D_f(1:NV,1:m_Star) = g(1:NV,:)

      DEALLOCATE (D_u, u, g)
#else
      call PS_Diff (D, f, D_f, N_Vector, Order, Even_Or_Odd)
#endif

  END SELECT

  END Subroutine PS_MPI_Diff_2D

END MODULE MPI_Differentiation
