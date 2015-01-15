#include "pseudopack.h"


MODULE PseudoPack_WENO_Type

implicit NONE

  TYPE PS_WENO_Property
    integer  :: Index, Distributed, Spatial_Order, Method, Symmetry
    integer  :: Cell_Center, Direction
    REALTYPE :: ISign

  END TYPE PS_WENO_Property

  TYPE PS_WENO_Grid_Index
    integer  :: NN_g, N_g, N0_g, N1_g, N2_g, N3_g, N4_g, N5_g
    integer  :: NN  , N  , N0  , N1  , N2  , N3  , N4  , N5  
    integer  :: n_Remain, n_0, n_Star
    integer  :: M_s

  END TYPE PS_WENO_Grid_Index

  TYPE PS_WENO_Domain
    REALTYPE          :: x0, x1
    integer           :: Map_F

  END TYPE PS_WENO_Domain

  TYPE PS_WENO_Mapping
    integer           :: Map, Manual
    REALTYPE          :: alpha, beta  

  END TYPE PS_WENO_Mapping

  TYPE PS_WENO_Scheme
    TYPE (PS_WENO_Property)   :: Property
    TYPE (PS_WENO_Grid_Index) :: Grid_Index
    TYPE (PS_WENO_Domain)     :: Domain
    TYPE (PS_WENO_Mapping)    :: Mapping

  END TYPE PS_WENO_Scheme

TYPE (PS_WENO_Property)   :: &
  Default_WENO_Property   = PS_WENO_Property   (1, 0, 5, 1, 0, 1, 1, ONE)
TYPE (PS_WENO_Grid_Index) :: &
  Default_WENO_Grid_Index = PS_WENO_Grid_Index (64, 64, -3, -1, 0, 64, 65, 67, &
                                                64, 64, -3, -1, 0, 64, 65, 67, &
                                                0, 0, 65,                      &
                                                3                              )
TYPE (PS_WENO_Domain)     :: &
  Default_WENO_Domain     = PS_WENO_Domain     (ZERO, ONE, 0)
TYPE (PS_WENO_Mapping)    :: &
  Default_WENO_Mapping    = PS_WENO_Mapping    (0, 0, ZERO, ZERO)

PUBLIC 

PUBLIC  :: PS_WENO_Property
PUBLIC  :: PS_WENO_Grid_Index
PUBLIC  :: PS_WENO_Domain
PUBLIC  :: PS_WENO_Mapping

PUBLIC  :: Default_WENO_Property
PUBLIC  :: Default_WENO_Grid_Index
PUBLIC  :: Default_WENO_Domain
PUBLIC  :: Default_WENO_Mapping

END MODULE PseudoPack_WENO_Type
