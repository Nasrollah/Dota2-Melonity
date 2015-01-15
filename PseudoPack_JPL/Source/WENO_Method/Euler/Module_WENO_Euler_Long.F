#include "pseudopack.h"


Module WENO_Euler_Long

  USE WENO_Option

  USE WENO_Euler_Long_JS
  USE WENO_Euler_Long_HAP
  USE WENO_Euler_Long_Zico

implicit NONE

PRIVATE
PUBLIC  :: PS_WENO_Euler_Long

CONTAINS

  Subroutine PS_WENO_Euler_Long (Order, N0, N5, i0, i1, F, F_Half,  &
                                           WENO_Epsilon, WENO_Power, WENO_Type)

  integer                    :: Order, N0, N5, i0, i1

  REALTYPE, dimension(N0:N5,-(Order-1)/2:(Order-1)/2,2) :: F
  REALTYPE, dimension(N0:N5                           ) :: F_Half

  REALTYPE, OPTIONAL         :: WENO_Epsilon
  REALTYPE, OPTIONAL         :: WENO_Power
  integer , OPTIONAL         :: WENO_Type
  integer                    :: WENO_Method

  REALTYPE                   :: Power
  logical                    :: Integer_Power

  REALTYPE                   :: p_r
  integer                    :: p_i

                          WENO_Method = PS_WENO_Option%Weight_Type
  IF (PRESENT(WENO_Type)) WENO_Method = WENO_Type

                           Power      = PS_WENO_Option%Power
  IF (PRESENT(WENO_Power)) Power      = WENO_Power
  
  Integer_Power = ABS(Power - INT(Power)) < SQRT(ABS(EPSILON(ONE))) 

  p_i = INT(Power)
  p_r =     Power 

  IF (Integer_Power) then
    SELECT CASE (WENO_Method)
      CASE DEFAULT
        call PS_WENO_Euler_Long_JS   (Order, N0, N5, i0, i1, F, F_Half,  &
                                             p_i, WENO_Epsilon)

      CASE (1)
        call PS_WENO_Euler_Long_HAP  (Order, N0, N5, i0, i1, F, F_Half,  &
                                             p_i, WENO_Epsilon)

      CASE (2)
        call PS_WENO_Euler_Long_Zico (Order, N0, N5, i0, i1, F, F_Half,  &
                                             p_i, WENO_Epsilon)

    END SELECT
  ELSE
    SELECT CASE (WENO_Method)
      CASE DEFAULT
        call PS_WENO_Euler_Long_JS   (Order, N0, N5, i0, i1, F, F_Half,  &
                                             p_r, WENO_Epsilon)

      CASE (1)
        call PS_WENO_Euler_Long_HAP  (Order, N0, N5, i0, i1, F, F_Half,  &
                                             p_r, WENO_Epsilon)

      CASE (2)
        call PS_WENO_Euler_Long_Zico (Order, N0, N5, i0, i1, F, F_Half,  &
                                             p_r, WENO_Epsilon)

    END SELECT
  ENDIF
 
  END Subroutine PS_WENO_Euler_Long 

END Module WENO_Euler_Long
