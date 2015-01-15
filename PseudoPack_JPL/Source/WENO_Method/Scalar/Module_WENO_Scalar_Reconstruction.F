#include "pseudopack.h"


Module WENO_Scalar_Reconstruction

  USE WENO_Option

  USE WENO_Scalar_Compact
  USE WENO_Scalar_Short
  USE WENO_Scalar_Long

implicit NONE

PRIVATE
PUBLIC  :: PS_WENO_Scalar_Reconstruction

CONTAINS

  Subroutine PS_WENO_Scalar_Reconstruction (Order, N0, N5, N2, N3, F, F_Half, &
                                                   Unroll_Loop, Epsilon, Power)

  integer                    :: Order, N0, N5, N2, N3

  REALTYPE, dimension(N0:N5) :: F
  REALTYPE, dimension(N0:N5) :: F_Half

  REALTYPE, OPTIONAL         :: Epsilon
  REALTYPE, OPTIONAL         :: Power

  integer , OPTIONAL         :: Unroll_Loop
  integer                    :: Unroll_Level

                            Unroll_Level = PS_WENO_Option%Unroll_Loop
  if (PRESENT(Unroll_Loop)) Unroll_Level = Unroll_Loop

  if ((Unroll_Level <  0) .OR. (Unroll_Level >  2)) Unroll_Level = 0
  
  SELECT CASE (Unroll_Level)
    CASE DEFAULT
      call PS_WENO_Scalar_Compact (Order, N0, N5, N2, N3, F, F_Half, Epsilon, Power)
  
    CASE (1)
      call PS_WENO_Scalar_Short   (Order, N0, N5, N2, N3, F, F_Half, Epsilon, Power)

    CASE (2)
      call PS_WENO_Scalar_Long    (Order, N0, N5, N2, N3, F, F_Half, Epsilon, Power)

  END SELECT

  END Subroutine PS_WENO_Scalar_Reconstruction 

END Module WENO_Scalar_Reconstruction
