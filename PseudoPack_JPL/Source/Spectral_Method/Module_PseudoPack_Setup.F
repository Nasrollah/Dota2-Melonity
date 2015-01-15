#include "pseudopack.h"


MODULE PseudoPack_Setup
  USE PseudoPack_Type

  USE PseudoPack_Differentiation, ONLY: Setup_Aux_Diff
  USE PseudoPack_Smoothing      , ONLY: Setup_Aux_Smth

  USE PS_IO_Unit

implicit NONE

INTERFACE PS_Setup_Operator
  MODULE PROCEDURE PS_Setup_Operator_X
  MODULE PROCEDURE PS_Setup_Operator_D
  MODULE PROCEDURE PS_Setup_Operator_S
END INTERFACE 

PRIVATE

PUBLIC  :: PS_Property  
PUBLIC  :: PS_Grid_Index
PUBLIC  :: PS_Domain    
PUBLIC  :: PS_Mapping   
PUBLIC  :: PS_Filtering 

PUBLIC  :: PS_Setup
PUBLIC  :: PS_Setup_Property
PUBLIC  :: PS_Setup_Grid_Index
PUBLIC  :: PS_Setup_Domain
PUBLIC  :: PS_Setup_Mapping
PUBLIC  :: PS_Setup_Filtering

PUBLIC  :: PS_Setup_Operator

CONTAINS

  Subroutine PS_Setup (Property  , Index, Method, Point_Type,                 &
                                     Max_Order_Of_D, Algorithm, Symmetry,     &
                       Grid_Index, N, M, LDY, BlockSize,                      &
                       Domain    , x0, x1, Map_F,                             &
                       Mapping   , Map, Manual, alpha, beta,                  &
                       Filtering , Filter_Choice, Mode_CutOff, Omega,         &
                                     Order, Smooth_1, Smooth_2                )
                       
  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Grid_Index), OPTIONAL :: Grid_Index
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping
  TYPE (PS_Filtering) , OPTIONAL :: Filtering

  integer , OPTIONAL :: Index, Method, Point_Type, Max_Order_Of_D,  &
                          Algorithm, Symmetry

  integer , OPTIONAL :: N, M, LDY, BlockSize

  REALTYPE, OPTIONAL :: x0, x1
  integer , OPTIONAL :: Map_F
  
  integer , OPTIONAL :: Map, Manual
  REALTYPE, OPTIONAL :: alpha, beta

  integer , OPTIONAL :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE, OPTIONAL :: Omega, Order

  if (PRESENT(Property  ))  &
    call PS_Setup_Property   (Property  , Index, Method, Point_Type,        &
                                            Max_Order_Of_D, Algorithm, Symmetry)

  if (PRESENT(Grid_Index))  &
    call PS_Setup_Grid_Index (Grid_Index, N, M, LDY, BlockSize)

  if (PRESENT(Domain    ))  &
    call PS_Setup_Domain     (Domain    , x0, x1, Map_F)

  if (PRESENT(Mapping   ))  &
    call PS_Setup_Mapping    (Mapping   , Map, Manual, alpha, beta)

  if (PRESENT(Filtering ))  &
    call PS_Setup_Filtering  (Filtering , Filter_Choice, Mode_CutOff, Omega, &
                                            Order, Smooth_1, Smooth_2)

#if defined (ESSL)
  if (Property%Algorithm /= 2)       RETURN
  if (PRESENT(M) .AND. PRESENT(LDY)) RETURN

  if (.NOT. PRESENT(M)) then
    write (lid6,100) ; write (lid99,100) 
  endif

  if (.NOT. PRESENT(LDY)) then
    write (lid6,110) ; write (lid99,110) 
  endif
 
 100 format (1x,60('*')/ &
             1x,'FATAL ERROR : PS_Setup ' / &
             3x,'Transform Algorithm in ESSL Library required that'/&
             3x,'the number of vectors M be provided in PS_Setup_Grid_Index'/&
             1x,60('*')/)
 110 format (1x,60('*')/ &
             1x,'FATAL ERROR : PS_Setup ' / &
             3x,'Transform Algorithm in ESSL Library required that'/&
             3x,'the leading dimension LDY be provided in PS_Setup_Grid_Index'/&
             1x,60('*')/)

  STOP
#endif

  END Subroutine PS_Setup
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Setup_Property (Property, Index, Method, Point_Type,  &
                                          Max_Order_Of_D, Algorithm, Symmetry)

  TYPE (PS_Property) ::   Property
  logical            :: L_Property

  integer , OPTIONAL :: Index, Method, Point_Type, Max_Order_Of_D,  &
                          Algorithm, Symmetry
  
    Property = Default_Property

  L_Property = PRESENT(Index     ) .OR. PRESENT(Method        ) .OR. &
               PRESENT(Point_Type) .OR. PRESENT(Max_Order_Of_D) .OR. &
               PRESENT(Algorithm ) .OR. PRESENT(Symmetry      )

  if (.NOT. L_Property) RETURN

  if (PRESENT(Index)         ) Property%Index          = Index
  if (PRESENT(Method)        ) Property%Method         = Method
  if (PRESENT(Point_Type)    ) Property%Point_Type     = Point_Type
  if (PRESENT(Max_Order_Of_D)) Property%Max_Order_Of_D = Max_Order_Of_D
  if (PRESENT(Algorithm)     ) Property%Algorithm      = Algorithm
  if (PRESENT(Symmetry)      ) Property%Symmetry       = Symmetry

  END Subroutine PS_Setup_Property 
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Setup_Grid_Index (Grid_Index, N, M, LDY, BlockSize)

  TYPE (PS_Grid_Index) ::   Grid_Index
  logical              :: L_Grid_Index

  integer , OPTIONAL :: N, M, LDY, BlockSize     

    Grid_Index = Default_Grid_Index

  L_Grid_Index = PRESENT(N) .OR. PRESENT(M) .OR. PRESENT(LDY)       &
                                            .OR. PRESENT(BlockSize)

  if (.NOT. L_Grid_Index) RETURN

  if (PRESENT(N        )) Grid_Index%N              = N     
  if (PRESENT(M        )) Grid_Index%M              = M     
  if (PRESENT(LDY)      ) Grid_Index%LDY            = LDY
  if (PRESENT(BlockSize)) Grid_Index%BlockSize      = BlockSize

  if ((.NOT. PRESENT(LDY)) .AND. PRESENT(N)) Grid_Index%LDY = N

#if defined (ESSL)
  if (PRESENT(M) .AND. PRESENT(LDY)) RETURN

  if (.NOT. PRESENT(M)) then
    write (lid99,*) 'CAUTION : Transform Algorithm in ESSL Library required that'
    write (lid99,*) '    M must be provided in PS_Setup_Grid_Index'
  endif

  if (.NOT. PRESENT(LDY)) then
    write (lid99,*) 'CAUTION : Transform Algorithm in ESSL Library required that'
    write (lid99,*) '  LDY must be provided in PS_Setup_Grid_Index'
  endif
#endif

  END Subroutine PS_Setup_Grid_Index 
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Setup_Domain (Domain, x0, x1, Map_F)

  TYPE (PS_Domain) ::   Domain
  logical          :: L_Domain

  REALTYPE, OPTIONAL :: x0, x1
  integer , OPTIONAL :: Map_F

    Domain = Default_Domain

  L_Domain = PRESENT(x0) .OR. PRESENT(x1) .OR. PRESENT(Map_F)

  if (.NOT. L_Domain) RETURN

  if (PRESENT(x0)   ) Domain%x0    = x0
  if (PRESENT(x1)   ) Domain%x1    = x1
  if (PRESENT(Map_F)) Domain%Map_F = Map_F

  END Subroutine PS_Setup_Domain 
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Setup_Mapping (Mapping, Map, Manual, alpha, beta)

  TYPE (PS_Mapping) ::   Mapping
  logical           :: L_Mapping

  integer , OPTIONAL :: Map, Manual
  REALTYPE, OPTIONAL :: alpha, beta

    Mapping = Default_Mapping

  L_Mapping = PRESENT(Map  ) .OR. PRESENT(Manual) .OR. &
              PRESENT(alpha) .OR. PRESENT(beta)

  if (.NOT. L_Mapping) RETURN

  if (PRESENT(Map)   ) Mapping%Map    = Map
  if (PRESENT(Manual)) Mapping%Manual = Manual
  if (PRESENT(alpha) ) Mapping%alpha  = alpha
  if (PRESENT(beta)  ) Mapping%beta   = beta

  END Subroutine PS_Setup_Mapping 
!
!-----------------------------------------------------------------------------
!
  Subroutine PS_Setup_Filtering (Filtering, Filter_Choice, Mode_CutOff, Omega, &
                                            Order, Smooth_1, Smooth_2)

  TYPE (PS_Filtering) ::   Filtering
  logical             :: L_Filtering

  integer , OPTIONAL :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE, OPTIONAL :: Omega, Order
  
    Filtering = Default_Filtering 

  L_Filtering = PRESENT(Filter_Choice) .OR. PRESENT(Mode_CutOff) .OR. &
                PRESENT(Omega        ) .OR. PRESENT(Order      ) .OR. &
                PRESENT(Smooth_1     ) .OR. PRESENT(Smooth_2   )

  if (.NOT. L_Filtering) RETURN

  if (PRESENT(Filter_Choice)) Filtering%Filter_Choice = Filter_Choice
  if (PRESENT(Mode_CutOff)  ) Filtering%Mode_CutOff   = Mode_CutOff
  if (PRESENT(Omega)        ) Filtering%Omega         = Omega
  if (PRESENT(Order)        ) Filtering%Order         = Order
  if (PRESENT(Smooth_1)     ) Filtering%Smooth_1      = Smooth_1
  if (PRESENT(Smooth_2)     ) Filtering%Smooth_2      = Smooth_2

  END Subroutine PS_Setup_Filtering 

!-----------------------------------------------------------------------------
!  PS_Operator Character(1)      Input
!    'D' or 'd'  for Differentiation Operator
!    'S' or 's'  for Smoothing       Operator
!
!  N_Pnts      integer           Input
!    Total number of grid points used
!  
!  Operator    Array(*)          Output
!  D           Array(*)          Output
!  S           Array(*)          Output
!
!  x           Array(*) Optional Output
! 
!  Property    Derived Data Type of PS_Property  , Optional, Input
!  Grid_Index  Derived Data Type of PS_Grid_Index, Optional, Input
!  Domain      Derived Data Type of PS_Domain    , Optional, Input
!  Mapping     Derived Data Type of PS_Mapping   , Optional, Input
!  Filtering   Derived Data Type of PS_Filtering , Optional, Input
!-----------------------------------------------------------------------------

  Subroutine PS_Setup_Operator_X (PS_Operator, N_Pnts, Operator, x,  &
                                  Property, Grid_Index, Domain,      &
                                  Mapping, Filtering)

  character(LEN=1) :: PS_Operator
  integer          :: N_Pnts

  REALTYPE, dimension(*)           :: Operator
  REALTYPE, dimension(:), OPTIONAL :: x

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Grid_Index), OPTIONAL :: Grid_Index
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping
  TYPE (PS_Filtering) , OPTIONAL :: Filtering

  logical  :: L_x

  SELECT CASE (PS_Operator)
    CASE ('d','D')

      L_x = PRESENT(x)

      if (.NOT. L_x) then
        write (lid6 ,100) ; write (lid99,100) ; STOP
      endif

      call PS_Setup_Operator_D (N_Pnts, Operator, x,           &
                                Property, Grid_Index, Domain,  &
                                Mapping, Filtering)

    CASE ('s','S')
      call PS_Setup_Operator_S (N_Pnts, Operator,              &
                                Property, Grid_Index, Filtering)

  END SELECT

 100 format (1x,60('*')/ &
             1x,'FATAL ERROR : PS_Setup_Operator ' / &
             3x,'The 1D array X must also be provided in order to setup'/&
             3x,'the differentiation Operator!'/&
             1x,60('*')/)

  END Subroutine PS_Setup_Operator_X
!
! ---------------------------------------------------------------------------
!
  Subroutine PS_Setup_Operator_D (N_Pnts, D, x,                    &
                                  Property, Grid_Index, Domain,    &
                                  Mapping, Filtering)
                               
  integer          :: N_Pnts

  REALTYPE, dimension(*) :: D
  REALTYPE, dimension(:) :: x

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Grid_Index), OPTIONAL :: Grid_Index
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping
  TYPE (PS_Filtering) , OPTIONAL :: Filtering

  logical  :: L_Property, L_Grid_Index, L_Domain, L_Mapping, L_Filtering

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry
  integer  :: N, M, ldy, BlockSize, Map, Manual, Map_F
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: x0, x1, alpha, beta, Omega, Order

  call Set_Default_Data ('D',                                                 &
                         N_Pnts    , L_Property, L_Grid_Index,                &
                                     L_Domain, L_Mapping, L_Filtering,        &
                         Property  , Index, Method, Point_Type,               &
                                       Max_Order_Of_D, Algorithm, Symmetry,   &
                         Grid_Index, N, M, LDY, BlockSize,                    &
                         Domain    , x0, x1, Map_F,                           &
                         Mapping   , Map, Manual, alpha, beta,                &
                         Filtering , Filter_Choice, Mode_CutOff, Omega, Order,&
                                       Smooth_1, Smooth_2                     )

  call Setup_Aux_Diff (Index, Method, Point_Type,            &
                       Max_Order_Of_D, Algorithm, Symmetry,  &
                       N, M, ldy, BlockSize,                 &
                       x, x0, x1, Map_F,                     &
                       Map, Manual, alpha, beta,             &
                       Filter_Choice, Mode_CutOff,           &
                       Omega, Order, Smooth_1, Smooth_2,     &
                       D)

  if (L_Property) then
    Property   = PS_Property   (Index, Method, Point_Type, Max_Order_Of_D, &
                                       Algorithm, Symmetry)
  endif

  if (L_Grid_Index) then
    Grid_Index = PS_Grid_Index (N, M, LDY, BlockSize)
  endif

  if (L_Domain) then
    Domain     = PS_Domain     (x0, x1, Map_F)
  endif

  if (L_Mapping) then
    Mapping    = PS_Mapping    (Map, Manual, alpha, beta)
  endif

  if (L_Filtering) then
    Filtering  = PS_Filtering  (Filter_Choice, Mode_CutOff, Omega, Order,  &
                                               Smooth_1, Smooth_2)
  endif

  END Subroutine PS_Setup_Operator_D
!
! ---------------------------------------------------------------------------
!
  Subroutine PS_Setup_Operator_S (N_Pnts, S,                        &
                                  Property, Grid_Index, Filtering)

  integer          :: N_Pnts

  REALTYPE, dimension(*) :: S

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Grid_Index), OPTIONAL :: Grid_Index
  TYPE (PS_Filtering) , OPTIONAL :: Filtering

  logical  :: L_Property, L_Grid_Index, L_Domain, L_Mapping, L_Filtering

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry
  integer  :: N, M, ldy, BlockSize
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: Omega, Order

  call Set_Default_Data ('S',                                               &
                         N_Pnts    , L_Property, L_Grid_Index,              &
                                     L_Domain, L_Mapping, L_Filtering,      &
                         Property  , Index, Method, Point_Type,             &
                                       Max_Order_Of_D, Algorithm, Symmetry, &
                         Grid_Index, N, M, LDY, BlockSize,                  &
                         Filtering=Filtering, Filter_Choice=Filter_Choice,  &
                                              Mode_CutOff  =Mode_CutOff  ,  &
                                              Omega        =Omega        ,  &
                                              Order        =Order        ,  &
                                              Smooth_1     =Smooth_1     ,  &
                                              Smooth_2     =Smooth_2        )

  Max_Order_Of_D = 0 ; Smooth_1 = 1 ; Smooth_2 = 1

  call Setup_Aux_Smth (Index, Method, Point_Type, Algorithm, Symmetry, &
                       N, M, ldy, BlockSize,                           &
                       Filter_Choice, Mode_CutOff, Omega, Order,       &
                       S)

  if (L_Property) then
    Property   = PS_Property   (Index, Method, Point_Type, Max_Order_Of_D, &
                                       Algorithm, Symmetry)
  endif

  if (L_Grid_Index) then
    Grid_Index = PS_Grid_Index (N, M, LDY, BlockSize)
  endif

  if (L_Filtering) then
    Filtering  = PS_Filtering  (Filter_Choice, Mode_CutOff, Omega, Order,  &
                                               Smooth_1, Smooth_2)
  endif

  END Subroutine PS_Setup_Operator_S
!
! --------------------------------------------------------------------------
!
  Subroutine Set_Default_Data (PS_Operator,                                  &
                               N_Pnts    , L_Property, L_Grid_Index,         &
                                           L_Domain, L_Mapping, L_Filtering, &
                               Property  , Index, Method, Point_Type,        &
                                           Max_Order_Of_D, Algorithm,        &
                                           Symmetry,                         &
                               Grid_Index, N, M, LDY, BlockSize,             & 
                               Domain    , x0, x1, Map_F,                    &
                               Mapping   , Map, Manual, alpha, beta,         &
                               Filtering , Filter_Choice, Mode_CutOff,       &
                                           Omega, Order, Smooth_1, Smooth_2)

  character(LEN=1) :: PS_Operator
  integer          :: N_Pnts

  logical  :: L_Property, L_Grid_Index, L_Domain, L_Mapping, L_Filtering

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Grid_Index), OPTIONAL :: Grid_Index
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping
  TYPE (PS_Filtering) , OPTIONAL :: Filtering

  integer , OPTIONAL :: Index, Method, Point_Type, Max_Order_Of_D,  &
                          Algorithm, Symmetry

  integer , OPTIONAL :: N, M, LDY, BlockSize

  REALTYPE, OPTIONAL :: x0, x1
  integer , OPTIONAL :: Map_F
  
  integer , OPTIONAL :: Map, Manual
  REALTYPE, OPTIONAL :: alpha, beta

  integer , OPTIONAL :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE, OPTIONAL :: Omega, Order

  L_Property   = PRESENT(Property  )
 
  L_Grid_Index = PRESENT(Grid_Index)
  L_Domain     = PRESENT(Domain    )
  L_Mapping    = PRESENT(Mapping   )
  L_Filtering  = PRESENT(Filtering )

  if (L_Property) then
    Index          =         Property%Index
    Method         =         Property%Method
    Point_Type     =         Property%Point_Type
    Max_Order_Of_D =         Property%Max_Order_Of_D
    Algorithm      =         Property%Algorithm
    Symmetry       =         Property%Symmetry
  else
    Index          = Default_Property%Index
    Method         = Default_Property%Method
    Point_Type     = Default_Property%Point_Type
    Max_Order_Of_D = Default_Property%Max_Order_Of_D
    Algorithm      = Default_Property%Algorithm
    Symmetry       = Default_Property%Symmetry
  endif

  if (L_Grid_Index) then
    N              =         Grid_Index%N
    M              =         Grid_Index%M
    LDY            =         Grid_Index%LDY
    BlockSize      =         Grid_Index%BlockSize
  else
    N              =                    N_Pnts
    M              = Default_Grid_Index%M
    LDY            =                    N_Pnts
    BlockSize      = Default_Grid_Index%BlockSize

    if (Index == 2) LDY = M
  endif
    
  if ((Method == 0) .AND. (Point_Type == 2) .AND. (MOD(N,2) /= 1)) N = N-1

  if ((PS_Operator == 'D') .OR. (PS_Operator == 'd')) then
    if (L_Domain) then
      x0           =         Domain%x0
      x1           =         Domain%x1

      Map_F        =         Domain%Map_F
    else
      SELECT CASE (Method)
        CASE (0)
          x0 = ZERO ; x1 = TWO*PI ; Map_F = 0

        CASE (1,2)
          x0 = -ONE ; x1 = ONE    ; Map_F = 0

      END SELECT
    endif
    
    if (L_Mapping) then
      Map          =         Mapping%Map
      Manual       =         Mapping%Manual
      alpha        =         Mapping%alpha
      beta         =         Mapping%beta
    else
      Map          = Default_Mapping%Map
      Manual       = Default_Mapping%Manual
      alpha        = Default_Mapping%alpha
      beta         = Default_Mapping%beta
    endif
  endif
    
  if (L_Filtering) then
    Filter_Choice  =         Filtering%Filter_Choice
    Mode_CutOff    =         Filtering%Mode_CutOff
    Omega          =         Filtering%Omega
    Order          =         Filtering%Order
    Smooth_1       =         Filtering%Smooth_1
    Smooth_2       =         Filtering%Smooth_2
  else
    Filter_Choice  = Default_Filtering%Filter_Choice
    Mode_CutOff    = Default_Filtering%Mode_CutOff
    Omega          = Default_Filtering%Omega
    Order          = Default_Filtering%Order
    Smooth_1       = Default_Filtering%Smooth_1
    Smooth_2       = Default_Filtering%Smooth_2
  endif

  END Subroutine Set_Default_Data 

END MODULE PseudoPack_Setup
