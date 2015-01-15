#include "pseudopack.h"


MODULE MPI_Smoothing

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

INTERFACE PS_MPI_Smooth
  MODULE PROCEDURE PS_MPI_Smooth_1D_G
  MODULE PROCEDURE PS_MPI_Smooth_2D_G
  MODULE PROCEDURE PS_MPI_Smooth_1D
  MODULE PROCEDURE PS_MPI_Smooth_2D
END INTERFACE

PRIVATE

PUBLIC  :: PS_MPI_Smooth

CONTAINS
!
!=======================================================================
!
  Subroutine PS_MPI_Smooth_1D_G (LDY, S, f, N_Vector, Even_Or_Odd)

  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Even_Or_Odd

  REALTYPE, dimension(*)                :: S
  REALTYPE, dimension(*)                :: f

  call PS_Smooth (LDY, S, f, N_Vector, Even_Or_Odd)

  END Subroutine PS_MPI_Smooth_1D_G 
!
!=======================================================================
!
  Subroutine PS_MPI_Smooth_2D_G (LDY, S, f, N_Vector, Even_Or_Odd)

  integer            :: LDY
  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Even_Or_Odd

  REALTYPE, dimension(*)                :: S
  REALTYPE, dimension(LDY,*)            :: f

#if defined (PARALLEL_MPI)
  REALTYPE, dimension(:,:), ALLOCATABLE :: g, u

  integer  :: N, n_Remain, n_0, n_Star
  integer  :: M, m_Remain, m_0, m_Star
  integer  :: NV
#endif

  SELECT CASE (OPERATOR_INDEX(S))
    CASE DEFAULT
      call PS_Smooth (LDY, S, f, N_Vector, Even_Or_Odd)

    CASE (2)
#if defined (PARALLEL_MPI)
      if (N_Processor == 1) then
        call PS_Smooth (LDY, S, f, N_Vector, Even_Or_Odd)

        RETURN
      endif

      N = OPERATOR_N_VECTORS(S)

      call PS_Get_Dimension (S, M)

      call PS_MPI_Data_Distribution (M, m_Remain, m_0, m_Star)
      call PS_MPI_Data_Distribution (N, n_Remain, n_0, n_Star)

      NV = n_Star ; if (PRESENT(N_Vector)) NV = MIN(N_Vector, n_Star)

      ALLOCATE (g(1:N,1:m_Star), u(1:M,1:n_Star))

      g = f(1:N,1:m_Star)

      call PS_P_Transpose (N, M, g, u)

      call PS_Smooth (S, u, NV, Even_Or_Odd)

      call PS_P_Transpose (M, N, u, g)

      NV = N ; if (PRESENT(N_Vector)) NV = N_Vector

      f(1:NV,1:m_Star) = g(1:NV,:)

      DEALLOCATE (u, g)
#else
      call PS_Smooth (LDY, S, f, N_Vector, Even_Or_Odd)
#endif

  END SELECT

  END Subroutine PS_MPI_Smooth_2D_G
!
!=======================================================================
!
  Subroutine PS_MPI_Smooth_1D   (S, f, N_Vector, Even_Or_Odd)

  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Even_Or_Odd

  REALTYPE, dimension(*)                :: S
  REALTYPE, dimension(:)                :: f

  call PS_Smooth (S, f, N_Vector, Even_Or_Odd)

  END Subroutine PS_MPI_Smooth_1D 
!
!=======================================================================
!
  Subroutine PS_MPI_Smooth_2D   (S, f, N_Vector, Even_Or_Odd)

  integer , OPTIONAL :: N_Vector
  integer , OPTIONAL :: Even_Or_Odd

  REALTYPE, dimension(*)                :: S
  REALTYPE, dimension(:,:)              :: f

#if defined (PARALLEL_MPI)
  REALTYPE, dimension(:,:), ALLOCATABLE :: g, u

  integer  :: N, n_Remain, n_0, n_Star
  integer  :: M, m_Remain, m_0, m_Star
  integer  :: NV
#endif

  SELECT CASE (OPERATOR_INDEX(S))
    CASE DEFAULT
      call PS_Smooth (S, f, N_Vector, Even_Or_Odd)

    CASE (2)
#if defined (PARALLEL_MPI)
      if (N_Processor == 1) then
        call PS_Smooth (S, f, N_Vector, Even_Or_Odd)

        RETURN
      endif

      N = OPERATOR_N_VECTORS(S)

      call PS_Get_Dimension (S, M)

      call PS_MPI_Data_Distribution (M, m_Remain, m_0, m_Star)
      call PS_MPI_Data_Distribution (N, n_Remain, n_0, n_Star)

      NV = n_Star ; if (PRESENT(N_Vector)) NV = MIN(N_Vector, n_Star)

      ALLOCATE (g(1:N,1:m_Star), u(1:M,1:n_Star))

      g = f(1:N,1:m_Star)

      call PS_P_Transpose (N, M, g, u)

      call PS_Smooth (S, u, NV, Even_Or_Odd)

      call PS_P_Transpose (M, N, u, g)

      NV = N ; if (PRESENT(N_Vector)) NV = N_Vector

      f(1:NV,1:m_Star) = g(1:NV,:)

      DEALLOCATE (u, g)
#else
      call PS_Smooth (S, f, N_Vector, Even_Or_Odd)
#endif

  END SELECT

  END Subroutine PS_MPI_Smooth_2D

END MODULE MPI_Smoothing
