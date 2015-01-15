  Subroutine MR_Flag_Plot (Filename, MR_Flag, Flag_WENO)

  Character(LEN=*)                                             :: Filename
  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,1:3) :: MR_Flag

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)           :: Flag
  integer , dimension(N0:N5,M0:M5,K0:K5,1:3), OPTIONAL :: Flag_WENO

  Character(LEN=16), dimension(3) :: Variable_Name
  Character(LEN=256)              :: FilenamA

  if (PRESENT(Flag_WENO)) then
    TecPlot_Header = .TRUE. ; FilenamA = TRIM(Filename)//'Flag_WENO'

    Variable_Name(1) = 'MR_Flag_x'
    Variable_Name(2) = 'MR_Flag_y'
    Variable_Name(3) = 'MR_Flag_z'

    Flag = Flag_WENO

    call PS_MPI_Ghost_Cell_Exchange (Index_x, Distributed_x,                 &
                                       N0, N2, N3, N5,                       &
                                     Index_y, Distributed_y,                 &
                                       M0, M2, M3, M5,                       &
                                     Index_z, Distributed_z,                 &
                                       K0, K2, K3, K5,                       &
                                     3, Flag, Topology%MPI_Communicator)

    call PS_Graphics (Graphic_Format, Data_Format,                           &
                      N0_g, N5_g, N0, N5, n_0, M_s_x,                        &
                      M0_g, M5_g, M0, M5, m_0, M_s_y,                        &
                      K0_g, K5_g, K0, K5, k_0, M_s_z,                        &
                      3 , x, y, z, Flag, Step, Time, dt,                     &
                      FilenamA, Variable_Name, TecPlot_Header, TecPlot_Title)
  endif

    TecPlot_Header = .TRUE. ; FilenamA = TRIM(Filename)//'Flag'

    Variable_Name(1) = 'MR_Flag_x'
    Variable_Name(2) = 'MR_Flag_y'
    Variable_Name(3) = 'MR_Flag_z'

    Flag = MR_Flag(N0:N5,M0:M5,K0:K5,:)

    call PS_MPI_Ghost_Cell_Exchange (Index_x, Distributed_x,                 &
                                       N0, N2, N3, N5,                       &
                                     Index_y, Distributed_y,                 &
                                       M0, M2, M3, M5,                       &
                                     Index_z, Distributed_z,                 &
                                       K0, K2, K3, K5,                       &
                                     3, Flag, Topology%MPI_Communicator)

    call PS_Graphics (Graphic_Format, Data_Format,                           &
                      N0_g, N5_g, N0, N5, n_0, M_s_x,                        &
                      M0_g, M5_g, M0, M5, m_0, M_s_y,                        &
                      K0_g, K5_g, K0, K5, k_0, M_s_z,                        &
                      3 , x, y, z, Flag, Step, Time, dt,                     &
                      FilenamA, Variable_Name, TecPlot_Header, TecPlot_Title)

    TecPlot_Header = .TRUE. ; FilenamA = TRIM(Filename)//'Scale'
  
    Variable_Name(1) = 'MR_scale_x'
    Variable_Name(2) = 'MR_scale_y'
    Variable_Name(3) = 'MR_scale_z'
  
    Flag = MR_Scale(N0:N5,M0:M5,K0:K5,:)
  
    call PS_MPI_Ghost_Cell_Exchange (Index_x, Distributed_x,                 &
                                       N0, N2, N3, N5,                       &
                                     Index_y, Distributed_y,                 &
                                       M0, M2, M3, M5,                       &
                                     Index_z, Distributed_z,                 &
                                       K0, K2, K3, K5,                       &
                                     3, Flag, Topology%MPI_Communicator)

    call PS_Graphics (Graphic_Format, Data_Format,                           &
                      N0_g, N5_g, N0, N5, n_0, M_s_x,                        &
                      M0_g, M5_g, M0, M5, m_0, M_s_y,                        &
                      K0_g, K5_g, K0, K5, k_0, M_s_z,                        &
                      3 , x, y, z, Flag, Step, Time, dt,                     &
                      FilenamA, Variable_Name, TecPlot_Header, TecPlot_Title)

    TecPlot_Header = .TRUE. ; FilenamA = TRIM(Filename)//'Q'

    Variable_Name(1) = 'MR_Q'

    Flag(:,:,:,1) = MR_Q(N0:N5,M0:M5,K0:K5,1)

    call PS_MPI_Ghost_Cell_Exchange (Index_x, Distributed_x,                 &
                                       N0, N2, N3, N5,                       &
                                     Index_y, Distributed_y,                 &
                                       M0, M2, M3, M5,                       &
                                     Index_z, Distributed_z,                 &
                                       K0, K2, K3, K5,                       &
                                     3, Flag, Topology%MPI_Communicator)

    call PS_Graphics (Graphic_Format, Data_Format,                           &
                      N0_g, N5_g, N0, N5, n_0, M_s_x,                        &
                      M0_g, M5_g, M0, M5, m_0, M_s_y,                        &
                      K0_g, K5_g, K0, K5, k_0, M_s_z,                        &
                      1 , x, y, z, Flag, Step, Time, dt,                     &
                      FilenamA, Variable_Name, TecPlot_Header, TecPlot_Title)

  End Subroutine MR_Flag_Plot
