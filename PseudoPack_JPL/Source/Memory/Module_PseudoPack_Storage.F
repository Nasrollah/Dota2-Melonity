#include "pseudopack.h"


MODULE PseudoPack_Storage_Check
  USE PseudoPack_Aux_Storage

Implicit NONE

PRIVATE

PUBLIC  :: Aux_Storage_Check
PUBLIC  :: Compute_Aux_Storage

CONTAINS

! ----------------------------------------------------------------------
! FUNCTION NAME: Aux_Storage_Check
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! DESCRIPTION  : Check for Storage Allocation for Aux_Diff and Aux_Smooth 
!
!  Index       = 0 (Ignored)
!                1 if check is call for algorithm along the FIRST  index
!                2 if check is call for algorithm along the SECOND index
!  D_Or_S      = 'D'-Differentiation, 'S'-Smoothing
!  Method      = 0-Fourier, 1-Chebyshev, 2-Legendre
!  Algorithm   = 0-MXM, 1-EOD, 2-Transform  
!  Symmetry    = 0-No, 1-Yes
!  S_Smooth    = 0-No, 1-Yes  for the D_Or_S='S'
!  N           = Total number of grid points to be differentiated/Smoothed.
!  M_Specified = Amount of Storage Specified  
!  M_Required  = Amount of Storage Required   (Output)
!  M_Error     = Error Return Code (Output)
!                  0 - No Error
!                  1 - Error in storage allocation
!  lid_Error   = Ouput file unit for Insufficient storage allocation and 
!                  Input data Error
!  lid_Info    = Ouput file unit for Minimum storage alloaction 
!                  lid_Info = 0 would disable the output.
!------------------------------------------------------------------------

  Subroutine Aux_Storage_Check (Index, D_Or_S, Method, Algorithm,  &
                                Symmetry, S_Smooth, N,             &
                                M_Specified, M_Required, M_Error,  &
                                lid_Error, lid_Info)

  character         :: D_Or_S     
  integer           :: Index, Method, Algorithm,          &
                       Symmetry, S_Smooth, N,             &
                       M_Specified, M_Required, M_Error,  &
                       lid_Error, lid_Info

  character(LEN=8)  :: Direction
  character(LEN=15) :: Operation

  if (D_Or_S == 'd') D_Or_S = 'D'
  if (D_Or_S == 's') D_Or_S = 'S'

  if (D_Or_S == 'D') Operation = 'Differentiation'
  if (D_Or_S == 'S') Operation = 'Smoothing      '

                  Direction = 'N/A     '
  if (Index == 1) Direction = '(FIRST) '
  if (Index == 2) Direction = '(SECOND)'

#if defined (DEBUG)
  call Storage_Input_Check (D_Or_S, Method, Algorithm, Symmetry, &
                                    S_Smooth, N, 6, 99, Index)
#endif

  call Get_Aux_XXXX_Size (D_Or_S, Method, Algorithm, Symmetry, &
                                  S_Smooth, N, M_Required)
 
                                 M_Error = 0
  if (M_Specified <  M_Required) M_Error = 1

  if ((lid_Error /= 0) .AND. (M_Error == 1)) then
    if (M_Specified /= 0) write (lid_Error,100) Operation, Direction
    if (M_Specified == 0) write (lid_Error,200) Operation, Direction

    write (lid_Error,120) M_Specified, Algorithm, N, M_Required
  endif

  if (lid_Info == 0) RETURN

    write (lid_Info ,200) Operation, Direction
    write (lid_Info ,120) M_Specified, Algorithm, N, M_Required

 100 format (1x/1x,'-------'/  &
                1x,'Error : Insufficient Storage Allocation for ',  &
               a15,' in ',a8,' Index'/1x,'-------')
 200 format (1x/1x,'-------'/  &
                1x,'Info  : Minimum      Storage Allocation for ',  &
               a15,' in ',a8,' Index'/1x,'-------')
 120 format (9x,'              Specified = ',1x,1x,5x,i10/  &
             9x,'Algorithm, N, Required  = ',i1,1x,i5,i10)

  END Subroutine Aux_Storage_Check 

! ----------------------------------------------------------------------
! FUNCTION NAME: Compute_Aux_Storage
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! DESCRIPTION  : Check for Storage Allocation for Aux_Diff and Aux_Smooth 
!
!  D_Or_S      = 'D'-Differentiation, 'S'-Smoothing
!  Method      = 0-Fourier, 1-Chebyshev, 2-Legendre
!  Algorithm   = 0-MXM, 1-EOD, 2-Transform  
!  Symmetry    = 0-No, 1-Yes
!  S_Smooth    = 0-No, 1-Yes  for the D_Or_S='S' only. Ignored otherwise
!  N           = Total number of grid points to be differentiated/Smoothed.
!  M_Required  = Amount of Storage Required   (Output)
!  M_Specified = Amount of Storage Specified  ( Input) (OPTIONAL)
!  M_Error     = Error Flag (Logical)         (Output) (OPTIONAL)
!------------------------------------------------------------------------
 
  Subroutine Compute_Aux_Storage (D_Or_S, Method, Algorithm, Symmetry, &
                                          S_Smooth, N,                 &
                                          M_Required, M_Specified, M_Error)

  character          :: D_Or_S     
  integer            :: Method, Algorithm, Symmetry, S_Smooth, N, M_Required 

  logical , OPTIONAL :: M_Error
  integer , OPTIONAL :: M_Specified
  logical            :: L_Error, L_Specified

  if (D_Or_S == 'd') D_Or_S = 'D'
  if (D_Or_S == 's') D_Or_S = 'S'

  L_Error     = PRESENT(M_Error    )
  L_Specified = PRESENT(M_Specified)
  
#if defined (DEBUG)
  call Storage_Input_Check (D_Or_S, Method, Algorithm, Symmetry, &
                                    S_Smooth, N, 6, 99)
#endif

  call Get_Aux_XXXX_Size (D_Or_S, Method, Algorithm, Symmetry, &
                                  S_Smooth, N, M_Required)

  if (L_Specified) then
    if (M_Required >  M_Specified) then
      write (6 ,100) D_Or_S, Method, Algorithm, Symmetry, S_Smooth, N,  &
                     M_Specified, M_Required
      write (99,100) D_Or_S, Method, Algorithm, Symmetry, S_Smooth, N,  &
                     M_Specified, M_Required

      if (L_Error) then ; M_Error = .TRUE. ; else ; STOP ; endif
    endif
  endif

 100 format (1x/1x,'FATAL Error : Inconsistent Memory Specification'/ &
                1x,'              D_Or_S           : ',a1 / &
                1x,'              Method           : ',i5 / &
                1x,'              Algorithm        : ',i5 / &
                1x,'              Symmetry         : ',i5 / &
                1x,'              S_Smooth         : ',i5 / &
                1x,'              N                : ',i5 / &
                1x,'              Memory Specified : ',i15/ &
                1x,'              Memory Required  : ',i15/)

  END Subroutine Compute_Aux_Storage 

#if defined (DEBUG)
!
!------------------------------------------------------------------------
!
  Subroutine Storage_Input_Check (D_Or_S, Method, Algorithm, Symmetry, &
                                  S_Smooth, N, lid_Info, lid_Error, Index)

  character         :: D_Or_S
  integer           :: Method, Algorithm, Symmetry, S_Smooth, N
  integer           :: lid_Info, lid_Error
  integer, OPTIONAL :: Index
  logical           :: IO_Error

  if (D_Or_S == 'd') D_Or_S = 'D'
  if (D_Or_S == 's') D_Or_S = 'S'

  IO_Error = .FALSE.

  if (PRESENT(Index)) then
    if ((Index <  0) .OR. (Index >  2)) then
      write (lid_Info ,100) Index   ; write (lid_Error,100) Index 
      IO_Error = .TRUE.
    endif
  endif

  if ((D_Or_S /= 'D') .AND. (D_Or_S /= 'S')) then
    write (lid_Info ,110) D_Or_S    ; write (lid_Error,110) D_Or_S 
    IO_Error = .TRUE.
  endif

  if ((Method <  0) .OR. (Method >  2)) then
    write (lid_Info ,120) Method    ; write (lid_Error,120) Method 
    IO_Error = .TRUE.
  endif
       
  if ((Algorithm <  0) .OR. (Algorithm >  2)) then
    write (lid_Info ,130) Algorithm ; write (lid_Error,130) Algorithm 
    IO_Error = .TRUE.
  endif

  if ((Method == 2) .AND. (Algorithm == 2)) then
    write (lid_Info ,140)           ; write (lid_Error,140) 
    IO_Error = .TRUE.
  endif

  if (N <= 0) then
    write (lid_Info ,150) N         ; write (lid_Error,150) N 
    IO_Error = .TRUE.
  endif

  if ((Symmetry /= 0) .AND. (Algorithm == 0)) then
    write (lid_Info ,160)           ; write (lid_Error,160) 
    IO_Error = .TRUE.
  endif

  if ((D_Or_S == 'S') .AND. (S_Smooth <  0)) then
    write (lid_Info ,170)           ; write (lid_Error,170)
    IO_Error = .TRUE.
  endif

  if (IO_Error) then
    write (lid_Info ,180)           ; write (lid_Error,180)
    STOP
  endif

 100 format (1x,'STOP! Invalid choice of Index     (0, 1,2) = ',i2)
 110 format (1x,'STOP! Invalid choice of D_Or_S    (D,S)    = ',a1)
 120 format (1x,'STOP! Invalid choice of Method    (0,1,2)  = ',i2)
 130 format (1x,'STOP! Invalid choice of Algorithm (0,1,2)  = ',i2)
 140 format (1x,'STOP! Legendre Method (2) have no Transform Algorithm (2)')
 150 format (1x,'STOP! N<=0')
 160 format (1x,'STOP! Requested MXM Algorithm with Symmetry.'/  &
             1x,'      Switch to EOD Algorithm with Symmetry.')
 170 format (1x,'STOP! Invalid choice of D_Or_S = S and S_Smooth < 0')
 180 format (1x,'STOP! Error in Aux_Storage Input data')
     
  END Subroutine Storage_Input_Check 
#endif

END MODULE PseudoPack_Storage_Check
