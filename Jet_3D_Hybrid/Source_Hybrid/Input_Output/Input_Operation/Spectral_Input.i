  call PS_Input (Index_x, Method_x, Point_Type_x, Max_Order_Of_D_x,  &
                   Algorithm_x, Symmetry_x,                          &
                 N_x, BlockSize_x, Left_x, Right_x, Map_F_x,         &
                   Map_x, Manual_x, alpha_x, Angle_x,                &
                 Filter_Choice_D_x, Mode_CutOff_D_x,                 &
                   Omega_D_x, Order_D_x, Smooth_1_x, Smooth_2_x,     &
                 Filter_Choice_S_x, Mode_CutOff_S_x,                 &
                   Omega_S_x, Order_S_x,                             &
                 lid)

  call PS_Input (Index_y, Method_y, Point_Type_y, Max_Order_Of_D_y,  &
                   Algorithm_y, Symmetry_y,                          &
                 N_y, BlockSize_y, Left_y, Right_y, Map_F_y,         &
                   Map_y, Manual_y, alpha_y, Angle_y,                &
                 Filter_Choice_D_y, Mode_CutOff_D_y,                 &
                   Omega_D_y, Order_D_y, Smooth_1_y, Smooth_2_y,     &
                 Filter_Choice_S_y, Mode_CutOff_S_y,                 &
                   Omega_S_y, Order_S_y,                             &
                 lid)

  call PS_Input (Index_z, Method_z, Point_Type_z, Max_Order_Of_D_z,  &
                   Algorithm_z, Symmetry_z,                          &
                 N_z, BlockSize_z, Left_z, Right_z, Map_F_z,         &
                   Map_z, Manual_z, alpha_z, Angle_z,                &
                 Filter_Choice_D_z, Mode_CutOff_D_z,                 &
                   Omega_D_z, Order_D_z, Smooth_1_z, Smooth_2_z,     &
                 Filter_Choice_S_z, Mode_CutOff_S_z,                 &
                   Omega_S_z, Order_S_z,                             &
                 lid)

  read (lid,'(//)')
  read (lid,*) RK_Smooth(1:5)
  read (lid,*) Local_Limit
  read (lid,*) Filtering_1
  read (lid,*) Order_1
