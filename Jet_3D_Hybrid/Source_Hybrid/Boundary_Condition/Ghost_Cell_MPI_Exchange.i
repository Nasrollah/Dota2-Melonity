#if defined (PARALLEL_MPI)

  M_s = MIN(N3-N2+1,M_s_x)

  call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                           MPI_Comm_World, MPI_Error_Status  )

  do n = M_s_0, M_s_x
    call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q,  1, &
                                        MPI_Communicator)
    call PS_MPI_Exchange_Boundary_Data (Index_x, Distributed_x, M_s_x, Q, -1, &
                                        MPI_Communicator)
  enddo

  M_s = MIN(M3-M2+1,M_s_y)

  call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                           MPI_Comm_World, MPI_Error_Status  )

  do n = M_s_0, M_s_y
    call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q,  1, &
                                        MPI_Communicator)
    call PS_MPI_Exchange_Boundary_Data (Index_y, Distributed_y, M_s_y, Q, -1, &
                                        MPI_Communicator)
  enddo

  M_s = MIN(K3-K2+1,M_s_z)

  call MPI_ALLREDUCE (M_s, M_s_0, 1, MPI_INTEGER, MPI_MIN,  &
                           MPI_Comm_World, MPI_Error_Status  )

  do n = M_s_0, M_s_z
    call PS_MPI_Exchange_Boundary_Data (Index_z, Distributed_z, M_s_z, Q,  1, &
                                        MPI_Communicator)
    call PS_MPI_Exchange_Boundary_Data (Index_z, Distributed_z, M_s_z, Q, -1, &
                                        MPI_Communicator)
  enddo

#endif
