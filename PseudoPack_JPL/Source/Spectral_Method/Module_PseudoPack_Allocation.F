#include "pseudopack.h"


MODULE PseudoPack_Allocation
  USE PseudoPack_Type
  USE PseudoPack_Storage_Check

implicit NONE

INTERFACE PS_Get_Operator_Size
  MODULE PROCEDURE PS_Get_Operator_Size_0
  MODULE PROCEDURE PS_Get_Operator_Size_1
  MODULE PROCEDURE PS_Get_Operator_Size_2
END INTERFACE

PRIVATE

PUBLIC  :: PS_Property
PUBLIC  :: PS_Filtering

PUBLIC  :: PS_Get_Operator_Size
PUBLIC  :: PS_Get_D_Operator_Size
PUBLIC  :: PS_Get_S_Operator_Size

CONTAINS

  Subroutine PS_Get_D_Operator_Size (N, M_Operator, Property, Error)

  TYPE (PS_Property), OPTIONAL :: Property
  logical           , OPTIONAL :: Error

  integer  :: N, M_Operator
  integer  :: Method, Algorithm, Symmetry, Smooth
  logical  :: L_Property

  L_Property = PRESENT(Property)

  if (L_Property) then
    Method    = Property%Method
    Algorithm = Property%Algorithm
    Symmetry  = Property%Symmetry
  else
    Method = 0 ; Algorithm = 0 ; Symmetry = 0
  endif

  Smooth = 0

                     call Get_Aux_Diff_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  END Subroutine PS_Get_D_Operator_Size
!
!------------------------------------------------------------------------
!
  Subroutine PS_Get_S_Operator_Size (N, M_Operator, S_Smooth, Property, Error)

  integer           , OPTIONAL :: S_Smooth
  TYPE (PS_Property), OPTIONAL :: Property
  logical           , OPTIONAL :: Error

  integer  :: N, M_Operator
  integer  :: Method, Algorithm, Symmetry, Smooth
  logical  :: L_Property

  L_Property = PRESENT(Property)

  if (L_Property) then
    Method    = Property%Method
    Algorithm = Property%Algorithm
    Symmetry  = Property%Symmetry
  else
    Method = 0 ; Algorithm = 0 ; Symmetry = 0
  endif

  Smooth = 1 ; if (PRESENT(S_Smooth)) Smooth = S_Smooth

                     call Get_Aux_Smth_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  END Subroutine PS_Get_S_Operator_Size
!
!------------------------------------------------------------------------
!
  Subroutine PS_Get_Operator_Size_0 (D_Or_S, N, M_Operator, S_Smooth, &
                                          Property, Error)

  integer                      :: S_Smooth
  TYPE (PS_Property), OPTIONAL :: Property
  logical           , OPTIONAL :: Error

  character :: D_Or_S

  integer  :: N, M_Operator
  integer  :: Method, Algorithm, Symmetry, Smooth
  logical  :: L_Property

  L_Property = PRESENT(Property)

  if (L_Property) then
    Method    = Property%Method
    Algorithm = Property%Algorithm
    Symmetry  = Property%Symmetry
  else
    Method = 0 ; Algorithm = 0 ; Symmetry = 0
  endif

  Smooth = S_Smooth

  if (D_Or_S == 'd') D_Or_S = 'D'
  if (D_Or_S == 's') D_Or_S = 'S'

  if (D_Or_S == 'D') call Get_Aux_Diff_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  if (D_Or_S == 'S') call Get_Aux_Smth_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  END Subroutine PS_Get_Operator_Size_0 
!
!------------------------------------------------------------------------
!
  Subroutine PS_Get_Operator_Size_1 (D_Or_S, N, M_Operator, Property, Error)

  TYPE (PS_Property), OPTIONAL :: Property
  logical           , OPTIONAL :: Error

  character :: D_Or_S

  integer  :: N, M_Operator
  integer  :: Method, Algorithm, Symmetry, Smooth
  logical  :: L_Property

  L_Property = PRESENT(Property)

  if (L_Property) then
    Method    = Property%Method
    Algorithm = Property%Algorithm
    Symmetry  = Property%Symmetry
  else
    Method = 0 ; Algorithm = 0 ; Symmetry = 0
  endif

  Smooth = 1

  if (D_Or_S == 'd') D_Or_S = 'D'
  if (D_Or_S == 's') D_Or_S = 'S'

  if (D_Or_S == 'D') call Get_Aux_Diff_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  if (D_Or_S == 'S') call Get_Aux_Smth_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  END Subroutine PS_Get_Operator_Size_1 
!
!------------------------------------------------------------------------
!
  Subroutine PS_Get_Operator_Size_2 (D_Or_S, N, M_Operator, Filtering, & 
                                          Property, Error)

  TYPE (PS_Filtering)           :: Filtering
  TYPE (PS_Property) , OPTIONAL :: Property
  logical            , OPTIONAL :: Error

  character :: D_Or_S

  integer  :: N, M_Operator
  integer  :: Method, Algorithm, Symmetry, Smooth
  logical  :: L_Property

  L_Property = PRESENT(Property)

  if (L_Property) then
    Method    = Property%Method
    Algorithm = Property%Algorithm
    Symmetry  = Property%Symmetry
  else
    Method = 0 ; Algorithm = 0 ; Symmetry = 0
  endif

  Smooth = Filtering%Filter_Choice

  if (D_Or_S == 'd') D_Or_S = 'D'
  if (D_Or_S == 's') D_Or_S = 'S'

  if (D_Or_S == 'D') call Get_Aux_Diff_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  if (D_Or_S == 'S') call Get_Aux_Smth_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  END Subroutine PS_Get_Operator_Size_2 
!
!------------------------------------------------------------------------
!
  Subroutine PS_Get_Operator_Size_3 (D_Or_S, Method, Algorithm, Symmetry, &
                                          N, M_Operator, S_Smooth, Error)

  integer            , OPTIONAL :: S_Smooth
  logical            , OPTIONAL :: Error

  character :: D_Or_S

  integer  :: N, M_Operator
  integer  :: Method, Algorithm, Symmetry, Smooth

  if (D_Or_S == 'd') D_Or_S = 'D'
  if (D_Or_S == 's') D_Or_S = 'S'

  Smooth = 1 ; if (PRESENT(S_Smooth)) Smooth = S_Smooth

  if (D_Or_S == 'D') call Get_Aux_Diff_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  if (D_Or_S == 'S') call Get_Aux_Smth_Size (Method, Algorithm, Symmetry, &
                                             Smooth, N, M_Operator, Error)

  END Subroutine PS_Get_Operator_Size_3
!
!------------------------------------------------------------------------
!
  Subroutine Get_Aux_Diff_Size (Method, Algorithm, Symmetry,  &
                                S_Smooth, N, M_Aux_Diff, Error)

  logical            , OPTIONAL :: Error

  integer  :: Method, Algorithm, Symmetry, S_Smooth, N, M_Aux_Diff
  integer  :: M_Required

  SELECT CASE (Method)
    CASE (0)                    ! Fourier
      SELECT CASE (Algorithm)
        CASE (0)                !   MXM
            M_Aux_Diff =   MXM_AUX_FULL_SIZE(N)

        CASE (1)                !   EOD
          if (Symmetry == 0) then
            M_Aux_Diff =   EOD_AUX_FULL_SIZE(N+1)
          else
            M_Aux_Diff =   EOD_AUX_HALF_SIZE(N+1)
          endif

        CASE (2)                !   FFT
          if (Symmetry == 0) then
            M_Aux_Diff = F_FFT_AUX_DIFF_FULL(N)
          else
            M_Aux_Diff = F_FFT_AUX_DIFF_HALF(N)
          endif

      END SELECT

    CASE (1)                    ! Chebyshev
      SELECT CASE (Algorithm)
        CASE (0)                !   MXM
            M_Aux_Diff =   MXM_AUX_FULL_SIZE(N)

        CASE (1)                !   EOD
          if (Symmetry == 0) then
            M_Aux_Diff =   EOD_AUX_FULL_SIZE(N)
          else
            M_Aux_Diff =   EOD_AUX_HALF_SIZE(N)
          endif

        CASE (2)                !   CFT
          if (Symmetry == 0) then
            M_Aux_Diff = C_CFT_AUX_DIFF_FULL(N)
          else
            M_Aux_Diff = C_CFT_AUX_DIFF_HALF(N)
          endif

      END SELECT

    CASE (2)                    ! Legendre
      SELECT CASE (Algorithm)
        CASE (0)                !   MXM
            M_Aux_Diff =   MXM_AUX_FULL_SIZE(N)

        CASE (1)                !   EOD
          if (Symmetry == 0) then
            M_Aux_Diff =   EOD_AUX_FULL_SIZE(N)
          else
            M_Aux_Diff =   EOD_AUX_HALF_SIZE(N)
          endif

        CASE (2)                
            M_Aux_Diff = 0

      END SELECT

  END SELECT

  call Compute_Aux_Storage ('D', Method, Algorithm, Symmetry, S_Smooth, N, &
                                 M_Required, M_Aux_Diff, Error)

  END Subroutine Get_Aux_Diff_Size 
!
!------------------------------------------------------------------------
!
  Subroutine Get_Aux_Smth_Size (Method, Algorithm, Symmetry,  &
                                S_Smooth, N, M_Aux_Smth, Error)

  logical            , OPTIONAL :: Error

  integer  :: Method, Algorithm, Symmetry, S_Smooth, N, M_Aux_Smth
  integer  :: M_Required

  SELECT CASE (Method)
    CASE (0)                    ! Fourier
      SELECT CASE (Algorithm)
        CASE (0)                !   MXM
            M_Aux_Smth =   MXM_AUX_FULL_SIZE(N)

        CASE (1)                !   EOD
          if (Symmetry == 0) then
            M_Aux_Smth =   EOD_AUX_FULL_SIZE(N+1)
          else
            M_Aux_Smth =   EOD_AUX_HALF_SIZE(N+1)
          endif

        CASE (2)                !   FFT
          if (Symmetry == 0) then
            M_Aux_Smth = F_FFT_AUX_SMTH_FULL(N)
          else
            M_Aux_Smth = F_FFT_AUX_SMTH_HALF(N)
          endif

      END SELECT

    CASE (1)                    ! Chebyshev
      SELECT CASE (Algorithm)
        CASE (0)                !   MXM
            M_Aux_Smth =   MXM_AUX_FULL_SIZE(N)

        CASE (1)                !   EOD
          if (Symmetry == 0) then
            M_Aux_Smth =   EOD_AUX_FULL_SIZE(N)
          else
            M_Aux_Smth =   EOD_AUX_HALF_SIZE(N)
          endif

        CASE (2)                !   CFT
          if (Symmetry == 0) then
            M_Aux_Smth = C_CFT_AUX_SMTH_FULL(N)
          else
            M_Aux_Smth = C_CFT_AUX_SMTH_HALF(N)
          endif

      END SELECT

    CASE (2)                    ! Legendre
      SELECT CASE (Algorithm)
        CASE (0)                !   MXM
            M_Aux_Smth =   MXM_AUX_FULL_SIZE(N)

        CASE (1)                !   EOD
          if (Symmetry == 0) then
            M_Aux_Smth =   EOD_AUX_FULL_SIZE(N)
          else
            M_Aux_Smth =   EOD_AUX_HALF_SIZE(N)
          endif

        CASE (2)                !  UNSUPPORTED
            M_Aux_Smth = 0

      END SELECT

  END SELECT

  call Compute_Aux_Storage ('S', Method, Algorithm, Symmetry, S_Smooth, N, &
                                 M_Required, M_Aux_Smth, Error)

  END Subroutine Get_Aux_Smth_Size 

END MODULE PseudoPack_Allocation
