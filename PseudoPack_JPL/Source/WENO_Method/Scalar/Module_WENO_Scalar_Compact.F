#include "pseudopack.h"


Module WENO_Scalar_Compact

  USE WENO_Option

  USE WENO_Scalar_Compact_JS
  USE WENO_Scalar_Compact_HAP
  USE WENO_Scalar_Compact_Zico

implicit NONE

PRIVATE
PUBLIC  :: PS_WENO_Scalar_Compact

CONTAINS

  Subroutine PS_WENO_Scalar_Compact (Order, N0, N5, i0, i1, F, F_Half, &
                                     WENO_Epsilon, WENO_Power, WENO_Type)

  integer                    :: Order, N0, N5, i0, i1

  REALTYPE, dimension(N0:N5) :: F, F_Half

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
        call PS_WENO_Scalar_Compact_JS   (Order, N0, N5, i0, i1, F, F_Half,  &
                                                 p_i, WENO_Epsilon)

      CASE (1)
        call PS_WENO_Scalar_Compact_HAP  (Order, N0, N5, i0, i1, F, F_Half,  &
                                                 p_i, WENO_Epsilon)

      CASE (2)
        call PS_WENO_Scalar_Compact_Zico (Order, N0, N5, i0, i1, F, F_Half,  &
                                                 p_i, WENO_Epsilon)

    END SELECT
  ELSE
    SELECT CASE (WENO_Method)
      CASE DEFAULT
        call PS_WENO_Scalar_Compact_JS   (Order, N0, N5, i0, i1, F, F_Half,  &
                                                 p_r, WENO_Epsilon)

      CASE (1)
        call PS_WENO_Scalar_Compact_HAP  (Order, N0, N5, i0, i1, F, F_Half,  &
                                                 p_r, WENO_Epsilon)

      CASE (2)
        call PS_WENO_Scalar_Compact_Zico (Order, N0, N5, i0, i1, F, F_Half,  &
                                                 p_r, WENO_Epsilon)

    END SELECT
  ENDIF

  END SUBROUTINE PS_WENO_Scalar_Compact

END Module WENO_Scalar_Compact
