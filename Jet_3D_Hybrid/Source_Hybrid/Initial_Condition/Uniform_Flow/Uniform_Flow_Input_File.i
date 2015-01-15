    if (I_Am == First_Processor) &
      write (6,*) 'Reading .... ', TRIM(Uniform_Input_File)

    lid = 78
    OPEN (UNIT=lid, FILE=TRIM(Uniform_Input_File), STATUS='OLD', IOSTAT=IOS)

    read (lid,'(//)')
    read (lid,*) Restart, Auto_Restart
    read (lid,*)
    read (lid,*) Begin_Time
    read (lid,*) Final_Time
    read (lid,*) Time_Factor
    read (lid,*) CFL_Euler, CFL_NS

    read (lid,'(//)')
    read (lid,*) Mach
    read (lid,*) Gamma
    read (lid,*)
    read (lid,*)  Density_Uniform
    read (lid,*) Pressure_Uniform

    read (lid,'(//)')
    read (lid,*) TecPlot_Title

#if 0
    read (lid,'(//)')
    read (lid,*) Ux, Uy, Uz
    read (lid,*) Source_xc, Source_yc, Source_zc
    read (lid,*) sigmax, sigmay, sigmaz
    read (lid,*) Jet_Strength
#endif

    CLOSE (lid)

    if (I_Am == First_Processor) &
      write (6,*) '....... Done ', TRIM(Uniform_Input_File)
