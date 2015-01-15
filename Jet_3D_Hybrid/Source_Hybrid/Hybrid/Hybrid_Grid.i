  call PS_WENO_Option_Setup (Order           = WENO_Order,             &
                             Epsilon         = WENO_Epsilon,           &
                             Power           = WENO_Power,             &
                             Lambda_Factor   = WENO_Lambda_Factor,     &
                             Sound_Speed_Min = WENO_Sound_Speed_Min,   &
                             Unroll_Loop     = WENO_Unroll_Loop,       &
                             Weight_Type     = WENO_Weight_Type        )

  M_s_x = (Spatial_Order_x+1)/2
  M_s_y = (Spatial_Order_y+1)/2
  M_s_z = (Spatial_Order_z+1)/2

  call PS_MPI_WENO_Grid_Index   (N_x, Index_x, Distributed_x,                  &
                                 Method_x, Symmetry_x, Spatial_Order_x,        &
                                 NN_g, N_g, N0_g, N1_g, N2_g, N3_g, N4_g, N5_g,&
                                 NN  , N  , N0  , N1  , N2  , N3  , N4  , N5  ,&
                                       n_Remain, n_0, n_Star,                  &
 
                                 N_y, Index_y, Distributed_y,                  &
                                 Method_y, Symmetry_y, Spatial_Order_y,        &
                                 MM_g, M_g, M0_g, M1_g, M2_g, M3_g, M4_g, M5_g,&
                                 MM  , M  , M0  , M1  , M2  , M3  , M4  , M5  ,&
                                       m_Remain, m_0, m_Star,                  &
 
                                 N_z, Index_z, Distributed_z,                  &
                                 Method_z, Symmetry_z, Spatial_Order_z,        &
                                 KK_g, K_g, K0_g, K1_g, K2_g, K3_g, K4_g, K5_g,&
                                 KK  , K  , K0  , K1  , K2  , K3  , K4  , K5  ,&
                                       k_Remain, k_0, k_Star,                  &

                                 Topology%MPI_Communicator)

#if defined (NAVIER_STOKES)
  if (Navier_Stokes) then
    call PS_Grid_Index_Update_NS  (N0, N2, N3, N5,  &
                                   M0, M2, M3, M5,  &
                                   K0, K2, K3, K5  )
  endif
#endif
