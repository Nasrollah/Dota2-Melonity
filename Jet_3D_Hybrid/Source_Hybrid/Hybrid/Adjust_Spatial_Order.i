#if defined (HYBRID)
  Central_Order_x = CD_Order_x
  Central_Order_y = CD_Order_y
  Central_Order_z = CD_Order_z

  Spatial_Order_x = MAX(WENO_Order_x, CD_Order_x)
  Spatial_Order_y = MAX(WENO_Order_y, CD_Order_y)
  Spatial_Order_z = MAX(WENO_Order_z, CD_Order_z)

  if (Filter_Choice_x /= 0) Spatial_Order_x = MAX(Spatial_Order_x, Filter_Order_x)
  if (Filter_Choice_y /= 0) Spatial_Order_y = MAX(Spatial_Order_y, Filter_Order_y)
  if (Filter_Choice_z /= 0) Spatial_Order_z = MAX(Spatial_Order_z, Filter_Order_z)
#else
  Spatial_Order_x =      WENO_Order_x
  Spatial_Order_y =      WENO_Order_y
  Spatial_Order_z =      WENO_Order_z
#endif

#if defined (NAVIER_STOKES)
  if (Navier_Stokes) then
    NS_Spatial_Order_x = 2*MAX(NS_TC_Order_x, NS_VD_Order_x, NS_MD_Order_x)
    NS_Spatial_Order_y = 2*MAX(NS_TC_Order_y, NS_VD_Order_y, NS_MD_Order_y)
    NS_Spatial_Order_z = 2*MAX(NS_TC_Order_z, NS_VD_Order_z, NS_MD_Order_z)

       Spatial_Order_x = MAX(Spatial_Order_x, NS_Spatial_Order_x)
       Spatial_Order_y = MAX(Spatial_Order_y, NS_Spatial_Order_y)
       Spatial_Order_z = MAX(Spatial_Order_z, NS_Spatial_Order_z)
  endif
#endif

