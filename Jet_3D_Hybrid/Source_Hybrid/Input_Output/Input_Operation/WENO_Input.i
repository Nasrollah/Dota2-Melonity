  read (lid,'(//)')
  read (lid,*) Index_x
  read (lid,*) Method_x
  read (lid,*) Symmetry_x
  read (lid,'(//)')
  read (lid,*) N_x 
  read (lid,'(//)')
  read (lid,*) Left_x
  read (lid,*) Right_x
#if 0
  read (lid,*)
  read (lid,*) WENO_Order_x        
  read (lid,*) WENO_Weight_Type_x
  read (lid,*) WENO_Epsilon_x
  read (lid,*) WENO_Power_x
  read (lid,*) WENO_Lambda_Factor_x
  read (lid,*) WENO_Sound_Speed_Min_x
  read (lid,*) WENO_Unroll_Loop_x
#endif

  read (lid,'(//)')
  read (lid,*) Index_y
  read (lid,*) Method_y
  read (lid,*) Symmetry_y
  read (lid,'(//)')
  read (lid,*) N_y 
  read (lid,'(//)')
  read (lid,*) Left_y
  read (lid,*) Right_y
#if 0
  read (lid,*)
  read (lid,*) WENO_Order_y        
  read (lid,*) WENO_Weight_Type_y
  read (lid,*) WENO_Epsilon_y
  read (lid,*) WENO_Power_y
  read (lid,*) WENO_Lambda_Factor_y
  read (lid,*) WENO_Sound_Speed_Min_y
  read (lid,*) WENO_Unroll_Loop_y
#endif

  read (lid,'(//)')
  read (lid,*) Index_z
  read (lid,*) Method_z
  read (lid,*) Symmetry_z
  read (lid,'(//)')
  read (lid,*) N_z 
  read (lid,'(//)')
  read (lid,*) Left_z
  read (lid,*) Right_z
#if 0
  read (lid,*)
  read (lid,*) WENO_Order_z        
  read (lid,*) WENO_Weight_Type_z
  read (lid,*) WENO_Epsilon_z
  read (lid,*) WENO_Power_z
  read (lid,*) WENO_Lambda_Factor_z
  read (lid,*) WENO_Sound_Speed_Min_z
  read (lid,*) WENO_Unroll_Loop_z
#endif

  read (lid,'(//)')
  read (lid,*) WENO_Order
  read (lid,*) WENO_Weight_Type
  read (lid,*) WENO_Power
  read (lid,*) WENO_Epsilon
  read (lid,*) WENO_Lambda_Factor
  read (lid,*) WENO_Sound_Speed_Min
  read (lid,*) WENO_Unroll_Loop

     WENO_Order_x = WENO_Order
     WENO_Order_y = WENO_Order
     WENO_Order_z = WENO_Order

#include "Adjust_WENO_Option.i"

  read (lid,'(//)')
  read (lid,*) CD_Order_x
  read (lid,*)
  read (lid,*) Filter_Choice_x
  read (lid,*) Filter_Order_x
  read (lid,*) Filter_Viscosity_x

  read (lid,'(//)')
  read (lid,*) CD_Order_y
  read (lid,*)
  read (lid,*) Filter_Choice_y
  read (lid,*) Filter_Order_y
  read (lid,*) Filter_Viscosity_y

  read (lid,'(//)')
  read (lid,*) CD_Order_z
  read (lid,*)
  read (lid,*) Filter_Choice_z
  read (lid,*) Filter_Order_z
  read (lid,*) Filter_Viscosity_z

  read (lid,'(//)')
  read (lid,*)   Viscous_Dissipation,    mu_0
  read (lid,*) NS_VD_Order_x, NS_VD_Order_y, NS_VD_Order_z
  read (lid,*)
  read (lid,*)   Thermal_Conduction , kappa_0
  read (lid,*) NS_TC_Order_x, NS_TC_Order_y, NS_TC_Order_z
  read (lid,*)
  read (lid,*) Molecular_Diffusion  ,  zeta_0 
  read (lid,*) NS_MD_Order_x, NS_MD_Order_y, NS_MD_Order_z

#include "Adjust_Spatial_Order.i"
