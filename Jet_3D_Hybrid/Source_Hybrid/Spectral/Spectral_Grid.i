  M_s_x = 0
  M_s_y = 0
  M_s_z = 0

  call PS_MPI_Spectral_Grid_Index (N_x, Index_x, Distributed_x,                  &
                                   Method_x, Symmetry_x, Point_Type_x,           &
                                   NN_g, N_g, N0_g, N1_g, N2_g, N3_g, N4_g, N5_g,&
                                   NN  , N  , N0  , N1  , N2  , N3  , N4  , N5  ,&
                                         n_Remain, n_0, n_Star,                  &
 
                                   N_y, Index_y, Distributed_y,                  &
                                   Method_y, Symmetry_y, Point_Type_y,           &
                                   MM_g, M_g, M0_g, M1_g, M2_g, M3_g, M4_g, M5_g,&
                                   MM  , M  , M0  , M1  , M2  , M3  , M4  , M5  ,&
                                         m_Remain, m_0, m_Star,                  &
   
                                   N_z, Index_z, Distributed_z,                  &
                                   Method_z, Symmetry_z, Point_Type_z,           &
                                   KK_g, K_g, K0_g, K1_g, K2_g, K3_g, K4_g, K5_g,&
                                   KK  , K  , K0  , K1  , K2  , K3  , K4  , K5  ,&
                                         k_Remain, k_0, k_Star,                  &

                                 Topology%MPI_Communicator)

