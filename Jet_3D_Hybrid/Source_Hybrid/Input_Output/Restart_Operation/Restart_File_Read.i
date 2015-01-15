if (Restart_Distributed_Read == 0) then
  call PS_Read_Restart_File (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                             M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                             K0_g, K5_g, K0, K5, k_0, M_s_z, NV,    &
                             Restart, Step, Time, dt, Q, Restart_Data_File)
else
  call PS_Read_Restart_File (Restart, Step, Time, dt, Q, Restart_Data_File)
endif

