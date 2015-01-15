#if defined (NAVIER_STOKES)

if (Navier_Stokes) then
  if (Thermal_Conduction) then
    call PS_CD_Operator_Setup (NS_TC_Order_x, CD_TC_x, dx, Order_Of_D=2)
    call PS_CD_Operator_Setup (NS_TC_Order_y, CD_TC_y, dy, Order_Of_D=2)
    call PS_CD_Operator_Setup (NS_TC_Order_z, CD_TC_z, dz, Order_Of_D=2)
  endif

  if (Viscous_Dissipation) then
    call PS_CD_Operator_Setup (NS_VD_Order_x, CD_VD_x, dx, Order_Of_D=2)
    call PS_CD_Operator_Setup (NS_VD_Order_y, CD_VD_y, dy, Order_Of_D=2)
    call PS_CD_Operator_Setup (NS_VD_Order_z, CD_VD_z, dz, Order_Of_D=2)
  endif

  if (Molecular_Diffusion) then
    call PS_CD_Operator_Setup (NS_MD_Order_x, CD_MD_x, dx, Order_Of_D=2)
    call PS_CD_Operator_Setup (NS_MD_Order_y, CD_MD_y, dy, Order_Of_D=2)
    call PS_CD_Operator_Setup (NS_MD_Order_z, CD_MD_z, dz, Order_Of_D=2)
  endif
endif

#endif
