  if (MR_Activation == 1) then

! DON : Comment out if use the input
    MR_Index = 1

    call PS_MR_Variable (Index_x, Distributed_x,                     &
                            N0,    N2,    N3,    N5,                 &
                         Index_y, Distributed_y,                     &
                            M0,    M2,    M3,    M5,                 &
                         Index_z, Distributed_z,                     &
                            K0,    K2,    K3,    K5,       NV,    Q, &
                         MR_N0, MR_N2, MR_N3, MR_N5,                 &
                         MR_M0, MR_M2, MR_M3, MR_M5,                 &
                         MR_K0, MR_K2, MR_K3, MR_K5, MR_Index, MR_Q, &
                         Topology%MPI_Communicator)

    call PS_MR_Analysis (MR_Order, MR_N0, MR_N5, MR_N2, MR_N3,             &
                         MR_Order, MR_M0, MR_M5, MR_M2, MR_M3,             &
                         MR_Order, MR_K0, MR_K5, MR_K2, MR_K3,             &
                          1, 1, 1, MR_Q(:,:,:,:), MR_Flag(:,:,:,1),        &
                         MR_Tolerance, Scale=MR_Scale(:,:,:,1) )

    call PS_MR_Analysis (MR_Order, MR_N0, MR_N5, MR_N2, MR_N3,             &
                         MR_Order, MR_M0, MR_M5, MR_M2, MR_M3,             &
                         MR_Order, MR_K0, MR_K5, MR_K2, MR_K3,             &
                          1, 2, 1, MR_Q(:,:,:,:), MR_Flag(:,:,:,2),        &
                         MR_Tolerance, Scale=MR_Scale(:,:,:,2) )

    call PS_MR_Analysis (MR_Order, MR_N0, MR_N5, MR_N2, MR_N3,             &
                         MR_Order, MR_M0, MR_M5, MR_M2, MR_M3,             &
                         MR_Order, MR_K0, MR_K5, MR_K2, MR_K3,             &
                          1, 3, 1, MR_Q(:,:,:,:), MR_Flag(:,:,:,3),        &
                         MR_Tolerance, Scale=MR_Scale(:,:,:,3) )

  endif
