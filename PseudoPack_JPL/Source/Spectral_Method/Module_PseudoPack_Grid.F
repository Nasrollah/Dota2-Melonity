#include "pseudopack.h"


MODULE PseudoPack_Grid
  USE PseudoPack_Type

  USE Grid_Generation, ONLY: Grid_And_Mapping

implicit NONE

INTERFACE PS_Collocation_Grid
  MODULE PROCEDURE PS_Collocation_Grid_1
  MODULE PROCEDURE PS_Collocation_Grid_2
  MODULE PROCEDURE PS_Collocation_Grid_3
END INTERFACE

PRIVATE

PUBLIC  :: PS_Property  
PUBLIC  :: PS_Grid_Index
PUBLIC  :: PS_Domain    
PUBLIC  :: PS_Mapping   

PUBLIC  :: PS_Collocation_Grid

CONTAINS

  Subroutine PS_Collocation_Grid_1 (x, N_Pnts, x0, x1, Map_F,        &
                                       Method, Point_Type, Symmetry, &
                                       Map, Manual, alpha, beta,     &
                                       Full_x, Error)

  REALTYPE, dimension(:) :: x

  integer ,   INTENT(IN) :: N_Pnts

  integer , OPTIONAL :: Method, Point_Type, Symmetry
  integer , OPTIONAL :: Map, Manual, Map_F
  REALTYPE, OPTIONAL :: alpha, beta
  REALTYPE           :: x0, x1

  logical , OPTIONAL :: Full_x
  logical , OPTIONAL :: Error

  logical  :: L_Full_x

  integer  :: I_Method, I_Point_Type, I_Symmetry
  integer  :: I_Map, I_Manual, I_Map_F
  REALTYPE :: R_alpha, R_beta, R_x0, R_x1
  integer  :: N, N0

  L_Full_x     = PRESENT(Full_x  )

                           I_Method     = Default_Property%Method     
  if (PRESENT(Method    )) I_Method     = Method
                           I_Point_Type = Default_Property%Point_Type 
  if (PRESENT(Point_Type)) I_Point_Type = Point_Type
                           I_Symmetry   = Default_Property%Symmetry   
  if (PRESENT(Symmetry  )) I_Symmetry   = Symmetry
    
  SELECT CASE (I_Method)
    CASE (0)
!     R_x0    = ZERO  ; if (PRESENT(x0   )) R_x0    = x0
!     R_x1    = TWOPI ; if (PRESENT(x1   )) R_x1    = x1
      I_Map_F = 0     ; if (PRESENT(Map_F)) I_Map_F = Map_F
 
    CASE (1,2)
!     R_x0    = -ONE  ; if (PRESENT(x0   )) R_x0    = x0
!     R_x1    =  ONE  ; if (PRESENT(x1   )) R_x1    = x1
      I_Map_F = 0     ; if (PRESENT(Map_F)) I_Map_F = Map_F

  END SELECT
    
  R_x0 = x0 ; R_x1 = x1

  I_Map    = Default_Mapping%Map    ; if (PRESENT(Map   )) I_Map    = Map
  I_Manual = Default_Mapping%Manual ; if (PRESENT(Manual)) I_Manual = Manual
  R_alpha  = Default_Mapping%alpha  ; if (PRESENT(alpha )) R_alpha  = alpha
  R_beta   = Default_Mapping%beta   ; if (PRESENT(beta  )) R_beta   = beta

  if (I_Method == 0) I_Manual = I_Map_F
  
  SELECT CASE (I_Method)
    CASE (0)
      N = N_Pnts   

      if ((I_Point_Type /= 2) .AND. (MOD(N,2) /= 0)) N = N-1
      if ((I_Point_Type == 2) .AND. (MOD(N,2) /= 1)) N = N-1

                     N0 = N-1 ; if (I_Symmetry == 1) N0 = N/2
      
    CASE (1,2)
      N = N_Pnts-1 ; N0 = N   ; if (I_Symmetry == 1) N0 = N/2

  END SELECT

  if (L_Full_x) then
    if (Full_x .AND. (I_Symmetry == 1)) x(1) = -99
  endif

#if defined (DEBUG)
  call Size_x_Check (x, N_Pnts, N0, I_Symmetry, L_Full_x, Full_x, Error)
#endif

  call Grid_And_Mapping (I_Method, I_Point_Type, N, N0,  &
                         I_Map, I_Manual, R_alpha, R_beta, &
                         R_x0, R_x1, x)

  END Subroutine PS_Collocation_Grid_1
!
! ----------------------------------------------------------------------
!
  Subroutine PS_Collocation_Grid_2 (x, Grid_Index, Property, Domain, &
                                       Mapping, Full_x, Error)

  REALTYPE, dimension(:)  :: x

  TYPE (PS_Grid_Index)           :: Grid_Index
  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping

  logical             , OPTIONAL :: Full_x
  logical             , OPTIONAL :: Error

  Logical  :: L_Property, L_Domain, L_Mapping, L_Full_x

  integer  :: Method, Point_Type, Symmetry
  integer  :: N, Map, Manual, Map_F
  REALTYPE :: alpha, beta, x0, x1
  integer  :: N_Pnts, N0

  L_Property   = PRESENT(Property)
  L_Domain     = PRESENT(Domain  )
  L_Mapping    = PRESENT(Mapping )
  L_Full_x     = PRESENT(Full_x  )

  if (L_Property) then
    Method         =         Property%Method
    Point_Type     =         Property%Point_Type
    Symmetry       =         Property%Symmetry
  else
    Method         = Default_Property%Method     
    Point_Type     = Default_Property%Point_Type 
    Symmetry       = Default_Property%Symmetry   
  endif

    N_Pnts         =         Grid_Index%N
    
  if (L_Domain) then
    x0             =         Domain%x0
    x1             =         Domain%x1

    Map_F          =         Domain%Map_F
  else
    SELECT CASE (Method)
      CASE (0)
        x0 = ZERO ; x1 = TWO*PI ; Map_F = 0
 
      CASE (1,2)
        x0 = -ONE ; x1 = ONE    ; Map_F = 0
    END SELECT
  endif
    
  if (L_Mapping) then
    Map            =         Mapping%Map
    Manual         =         Mapping%Manual
    alpha          =         Mapping%alpha
    beta           =         Mapping%beta
  else
    Map            = Default_Mapping%Map    
    Manual         = Default_Mapping%Manual 
    alpha          = Default_Mapping%alpha  
    beta           = Default_Mapping%beta   
  endif
    
  if (Method == 0) Manual = Map_F

  SELECT CASE (Method)
    CASE (0)
      N = N_Pnts   

      if ((Point_Type /= 2) .AND. (MOD(N,2) /= 0)) N = N-1
      if ((Point_Type == 2) .AND. (MOD(N,2) /= 1)) N = N-1

                     N0 = N-1 ; if (Symmetry == 1) N0 = N/2
      
    CASE (1,2)
      N = N_Pnts-1 ; N0 = N   ; if (Symmetry == 1) N0 = N/2

  END SELECT

  if (L_Full_x) then
    if (Full_x .AND. (Symmetry == 1)) x(1) = -99
  endif

#if defined (DEBUG)
  call Size_x_Check (x, N_Pnts, N0, Symmetry, L_Full_x, Full_x, Error)
#endif

  call Grid_And_Mapping (Method, Point_Type, N, N0, Map, Manual, alpha, beta, &
                         x0, x1, x)

  END Subroutine PS_Collocation_Grid_2
!
! ----------------------------------------------------------------------
!
  Subroutine PS_Collocation_Grid_3 (x, N_Pnts, Property, Domain, &
                                       Mapping, Full_x, Error)

  REALTYPE, dimension(:)  :: x

  integer                 :: N_Pnts

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping

  logical             , OPTIONAL :: Full_x
  logical             , OPTIONAL :: Error

  Logical  :: L_Property, L_Domain, L_Mapping, L_Full_x

  integer  :: Method, Point_Type, Symmetry
  integer  :: N, Map, Manual, Map_F
  REALTYPE :: alpha, beta, x0, x1
  integer  :: N0

  L_Property   = PRESENT(Property  )
  L_Domain     = PRESENT(Domain    )
  L_Mapping    = PRESENT(Mapping   )
  L_Full_x     = PRESENT(Full_x    )

  if (L_Property) then
    Method         =         Property%Method
    Point_Type     =         Property%Point_Type
    Symmetry       =         Property%Symmetry
  else
    Method         = Default_Property%Method     
    Point_Type     = Default_Property%Point_Type 
    Symmetry       = Default_Property%Symmetry   
  endif

  if (L_Domain) then
    x0             =         Domain%x0
    x1             =         Domain%x1

    Map_F          =         Domain%Map_F
  else
    SELECT CASE (Method)
      CASE (0)
        x0 = ZERO ; x1 = TWO*PI ; Map_F = 0
 
      CASE (1,2)
        x0 = -ONE ; x1 = ONE    ; Map_F = 0
    END SELECT
  endif
    
  if (L_Mapping) then
    Map            =         Mapping%Map
    Manual         =         Mapping%Manual
    alpha          =         Mapping%alpha
    beta           =         Mapping%beta
  else
    Map            = Default_Mapping%Map    
    Manual         = Default_Mapping%Manual 
    alpha          = Default_Mapping%alpha  
    beta           = Default_Mapping%beta   
  endif
    
  if (Method == 0) Manual = Map_F

  SELECT CASE (Method)
    CASE (0)
      N = N_Pnts   

      if ((Point_Type /= 2) .AND. (MOD(N,2) /= 0)) N = N-1
      if ((Point_Type == 2) .AND. (MOD(N,2) /= 1)) N = N-1

                     N0 = N-1 ; if (Symmetry == 1) N0 = N/2
      
    CASE (1,2)
      N = N_Pnts-1 ; N0 = N   ; if (Symmetry == 1) N0 = N/2

  END SELECT

  if (L_Full_x) then
    if (Full_x .AND. (Symmetry == 1)) x(1) = -99
  endif

#if defined (DEBUG)
  call Size_x_Check (x, N_Pnts, N0, Symmetry, L_Full_x, Full_x, Error)
#endif

  call Grid_And_Mapping (Method, Point_Type, N, N0, Map, Manual, alpha, beta, &
                         x0, x1, x)

  END Subroutine PS_Collocation_Grid_3
!
! ----------------------------------------------------------------------
!
#if defined (DEBUG)
  Subroutine Size_x_Check (x, N_Pnts, N0, Symmetry, L_Full_x, Full_x, Error)

  REALTYPE, dimension(:) :: x
  integer                :: N_Pnts, N0
  integer                :: Symmetry
  logical                :: L_Full_x
  logical , OPTIONAL     :: Full_x, Error

  if (PRESENT(Error)) Error = .FALSE.

  if (Symmetry /= 1) then
    if (SIZE(x) <  N_Pnts) then
      write (6,100) N_Pnts, SIZE(x) ; write (99,100) N_Pnts, SIZE(x)
 
      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  else
    if (L_Full_x) then
      if ((Size(x) <  N_Pnts) .AND. Full_x) then
        write (6,101) N_Pnts, SIZE(x) ; write (99,101) N_Pnts, SIZE(x)
 
        if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
      endif
    else
      if (Size(x) <  N0+1) then
        write (6,102) N0+1, SIZE(x) ; write (99,102) N0+1, SIZE(x)
 
        if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
      endif
    endif
  endif

 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Collocation_Point '/  &
             15x,'Insufficient storage was allocated to array x'/ &
             15x,'  for Non-Symmetric case and Full array'/ &
             15x,'Required  Size for x = ',i4/ &
             15x,'Actual    Size for x = ',i4/ &
             1x,60('*'))
 101 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Collocation_Point '/  &
             15x,'Insufficient storage was allocated to array x'/ &
             15x,'  for Symmetric case and Full array'/ &
             15x,'Required  Size for x = ',i4/ &
             15x,'Actual    Size for x = ',i4/ &
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Collocation_Point '/  &
             15x,'Insufficient storage was allocated to array x'/ &
             15x,'  for Symmetric case and Half array'/ &
             15x,'Required  Size for x = ',i4/ &
             15x,'Actual    Size for x = ',i4/ &
             1x,60('*'))

  END Subroutine Size_x_Check 
#endif

END MODULE PseudoPack_Grid
