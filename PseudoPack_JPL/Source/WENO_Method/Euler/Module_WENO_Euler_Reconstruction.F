#include "pseudopack.h"


Module WENO_Euler_Reconstruction

  USE WENO_Euler_Compact
  USE WENO_Euler_Short
  USE WENO_Euler_Long
  USE WENO_Option

implicit NONE

PRIVATE
PUBLIC  :: PS_WENO_Euler_Reconstruction

CONTAINS

  Subroutine PS_WENO_Euler_Reconstruction (Order, N0, N5, N2, N3, F, F_Half, &
                                           Unroll_Loop, Epsilon, Power, Weight_Type)

  integer                    :: Order, N0, N5, N2, N3

  REALTYPE, dimension(N0:N5,-(Order-1)/2:(Order-1)/2,2) :: F
  REALTYPE, dimension(N0:N5                           ) :: F_Half

  REALTYPE, OPTIONAL         :: Epsilon
  REALTYPE, OPTIONAL         :: Power

  integer , OPTIONAL         :: Unroll_Loop
  integer                    :: Unroll_Level

  integer , OPTIONAL         :: Weight_Type

                            Unroll_Level = PS_WENO_Option%Unroll_Loop
  if (PRESENT(Unroll_Loop)) Unroll_Level = Unroll_Loop

  if ((Unroll_Level <  0) .OR. (Unroll_Level >  2)) Unroll_Level = 0
  
  SELECT CASE (Unroll_Level)
    CASE DEFAULT
      call PS_WENO_Euler_Compact (Order, N0, N5, N2, N3, F, F_Half, & 
                                         Epsilon, Power, Weight_Type)
  
    CASE (1)
      call PS_WENO_Euler_Short   (Order, N0, N5, N2, N3, F, F_Half, &
                                         Epsilon, Power, Weight_Type)

    CASE (2)
      call PS_WENO_Euler_Long    (Order, N0, N5, N2, N3, F, F_Half, &
                                         Epsilon, Power, Weight_Type)

  END SELECT

  END Subroutine PS_WENO_Euler_Reconstruction 

END Module WENO_Euler_Reconstruction
