  call PS_WENO_BC_Setup   (BC_On_Off_x, Distributed_x, WENO_Order_x,         &
                                      Method_x, Symmetry_x, BC_Symmetry_x, &
                                      BC_Type_x,                           &
                           BC_On_Off_y, Distributed_y, WENO_Order_y,         &
                                      Method_y, Symmetry_y, BC_Symmetry_y, &
                                      BC_Type_y,                           &
                           BC_On_Off_z, Distributed_z, WENO_Order_z,         &
                                      Method_z, Symmetry_z, BC_Symmetry_z, &
                                      BC_Type_z,                           &
                           BC_On_Off, BC_Distributed, BC_Order, BC_Method,   &
                           BC_Symmetry, BC_Type)
