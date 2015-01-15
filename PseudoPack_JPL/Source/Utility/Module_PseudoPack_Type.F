#include "pseudopack.h"


MODULE PseudoPack_Type

implicit NONE

  TYPE PS_Property
    integer           :: Index, Method, Point_Type, Max_Order_Of_D 
    integer           :: Algorithm, Symmetry

  END TYPE PS_Property

  TYPE PS_Grid_Index
    integer           :: N, M, LDY, BlockSize

  END TYPE PS_Grid_Index

  TYPE PS_Domain
    REALTYPE          :: x0, x1
    integer           :: Map_F

  END TYPE PS_Domain

  TYPE PS_Mapping
    integer           :: Map, Manual
    REALTYPE          :: alpha, beta  

  END TYPE PS_Mapping

  TYPE PS_Filtering
    integer           :: Filter_Choice, Mode_CutOff
    REALTYPE          :: Omega, Order
    integer           :: Smooth_1, Smooth_2

  END TYPE PS_Filtering
    
  TYPE PS_Parallel
    integer           :: MPI_Distributed
  END TYPE PS_Parallel

  TYPE PS_Scheme
    TYPE (PS_Property)   :: Property
    TYPE (PS_Grid_Index) :: Grid_Index
    TYPE (PS_Domain)     :: Domain
    TYPE (PS_Mapping)    :: Mapping
    TYPE (PS_Filtering)  :: Filtering_D
    TYPE (PS_Filtering)  :: Filtering_S
    TYPE (PS_Filtering)  :: Filtering_R
    TYPE (PS_Parallel)   :: Parallel

  END TYPE PS_Scheme

#if defined (DEFAULT_FOURIER)
TYPE (PS_Property)   :: Default_Property   = PS_Property   (1, 0, 1, 1, 0, 0)
TYPE (PS_Grid_Index) :: Default_Grid_Index = PS_Grid_Index (32, 1, 32, 0)
TYPE (PS_Domain)     :: Default_Domain     = PS_Domain     (ZERO, TWO*PI, 0)
#else
TYPE (PS_Property)   :: Default_Property   = PS_Property   (1, 1, 1, 1, 0, 0)
TYPE (PS_Grid_Index) :: Default_Grid_Index = PS_Grid_Index (33, 1, 33, 0)
TYPE (PS_Domain)     :: Default_Domain     = PS_Domain     (-ONE, ONE,    0)
#endif
TYPE (PS_Mapping)    :: Default_Mapping    = PS_Mapping    (0, 0, ZERO, ZERO)
TYPE (PS_Filtering)  :: Default_Filtering  = PS_Filtering  (0, 0, ZERO, ZERO, &
                                                               1, 1)
PUBLIC 

PUBLIC  :: PS_Property
PUBLIC  :: PS_Grid_Index
PUBLIC  :: PS_Domain
PUBLIC  :: PS_Mapping
PUBLIC  :: PS_Filtering

PUBLIC  :: Default_Property
PUBLIC  :: Default_Grid_Index
PUBLIC  :: Default_Domain
PUBLIC  :: Default_Mapping
PUBLIC  :: Default_Filtering

END MODULE PseudoPack_Type
