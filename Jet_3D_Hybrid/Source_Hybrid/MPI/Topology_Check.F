  Subroutine Topology_Check_1D (N2, N3, M_s_x)

  integer  :: N2, N3, M_s_x
 
  if (N3-N2+1 < M_s_x) then
    write (6,*) 'I_Am, N2, N3, M_s_x :', I_Am, N2, N3, M_s_x

    call PS_STOP ('Topology', ' ','1st Dim : N3-N2+1 < M_s_x')
  endif

  END Subroutine Topology_Check_1D 
! 
! ================================================================
!
  Subroutine Topology_Check_2D (N2, N3, M_s_x,  &
                                M2, M3, M_s_y   )

  integer  :: N2, N3, M_s_x
  integer  :: M2, M3, M_s_y
 
  if (N3-N2+1 < M_s_x) then
    write (6,*) 'I_Am, N2, N3, M_s_x :', I_Am, N2, N3, M_s_x

    call PS_STOP ('Topology', ' ','1st Dim : N3-N2+1 < M_s_x')
  endif

  if (M3-M2+1 < M_s_y) then
    write (6,*) 'I_Am, M2, M3, M_s_y :', I_Am, M2, M3, M_s_y

    call PS_STOP ('Topology', ' ','2st Dim : M3-M2+1 < M_s_y')
  endif

  END Subroutine Topology_Check_2D 
! 
! ================================================================
!
  Subroutine Topology_Check_3D (N2, N3, M_s_x,  &
                                M2, M3, M_s_y,  &
                                K2, K3, M_s_z   )

  integer  :: N2, N3, M_s_x
  integer  :: M2, M3, M_s_y
  integer  :: K2, K3, M_s_z
 
  if (N3-N2+1 <  M_s_x) then
    write (6,*) 'I_Am, N2, N3, M_s_x :', I_Am, N2, N3, M_s_x

    call PS_STOP ('Topology', ' ','1st Dim : N3-N2+1 < M_s_x')
  endif

  if (M3-M2+1 <  M_s_y) then
    write (6,*) 'I_Am, M2, M3, M_s_y :', I_Am, M2, M3, M_s_y

    call PS_STOP ('Topology', ' ','2st Dim : M3-M2+1 < M_s_y')
  endif

  if (K3-K2+1 <  M_s_z) then
    write (6,*) 'I_Am, N2, N3, M_s_z :', I_Am, N2, N3, M_s_z

    call PS_STOP ('Topology', ' ','3st Dim : K3-K2+1 < M_s_z')
  endif

  END Subroutine Topology_Check_3D 
