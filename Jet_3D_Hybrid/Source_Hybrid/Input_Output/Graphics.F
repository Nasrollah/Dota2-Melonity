  Subroutine TecPlot (Filename)

  Character(LEN=*)   :: Filename
  Character(LEN=256) :: Name

  Character(LEN=16), dimension(NF+NS) :: Variable_Name

  Variable_Name(1:5) = (/ 'Rho', 'RhoU', 'RhoV', 'RhoW', 'E' /)

#if defined (WITH_SPECIES_EQUATION)
  do n = 1,NS-1
    call PS_Make_Filename ('Rhof','',Filename=Name,Number=n)
    Variable_Name(NF+n) = TRIM(Name)
  enddo
#endif

  call PS_Graphics (Graphic_Format, Data_Format,              &
                    N0_g, N5_g, N0, N5, n_0, M_s_x,           &
                    M0_g, M5_g, M0, M5, m_0, M_s_y,           &
                    K0_g, K5_g, K0, K5, k_0, M_s_z,           &
                    NV, x, y, z, Q, Step, Time, dt_Original,  &
                    Filename, Variable_Name, TecPlot_Header, TecPlot_Title)

  End Subroutine TecPlot
