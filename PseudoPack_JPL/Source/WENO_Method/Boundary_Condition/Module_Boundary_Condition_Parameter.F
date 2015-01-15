#include "pseudopack.h"


MODULE Boundary_Condition_Parameter

IMPLICIT NONE

integer :: BC_On_Off_x, BC_Symmetry_x,  BC_Type_Left  , BC_Type_Right
integer :: BC_On_Off_y, BC_Symmetry_y,  BC_Type_Bottom, BC_Type_Top
integer :: BC_On_Off_z, BC_Symmetry_z,  BC_Type_Down  , BC_Type_Up

integer , dimension(2)   :: BC_Type_x
integer , dimension(2)   :: BC_Type_y
integer , dimension(2)   :: BC_Type_z

integer , dimension(3)   :: BC_On_Off, BC_Method, BC_Order, BC_Symmetry
logical , dimension(3)   :: BC_Distributed
integer , dimension(3,2) :: BC_Type

PUBLIC 

END MODULE Boundary_Condition_Parameter
