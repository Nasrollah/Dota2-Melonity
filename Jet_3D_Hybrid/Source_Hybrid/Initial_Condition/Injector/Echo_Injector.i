  if (I_Am == First_Processor) then
    do i = 1,2
      lid=6 ; if (i == 2) lid=60

      write (lid,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
      write (lid,*) 'Injector Specification : '
      write (lid,*) 
      write (lid,*) '  Direction :', Jet_Direction
      write (lid,*) 
      write (lid,*) '  Shape     :', Jet_Shape
      write (lid,*) '  Center    :', Jet_Center
      write (lid,*) '  Circle    :', Jet_Circle
      write (lid,*) '  Square    :', Jet_Square
      write (lid,*) '  Height    :', Jet_Height
      write (lid,*)
      write (lid,*) '  Mach      :', Jet_Mach
      write (lid,*) '  Gamma     :', Jet_Gamma
      write (lid,*)
      write (lid,*) '  Re        :', Jet_Re
      write (lid,*) '  Pr        :', Jet_Pr
      write (lid,*)
      write (lid,*) '  Rho       :', Jet_Rho
      write (lid,*) '  U         :', Jet_U
      write (lid,*) '  P         :', Jet_P
      write (lid,*) '  T         :', Jet_T
      write (lid,*) '  C         :', Jet_C
      write (lid,*) '  f1        :', Jet_f1
      write (lid,*)
      write (lid,*) '  mu        :', Jet_mu
      write (lid,*) '  kappa     :', Jet_kappa
      write (lid,*) '  zeta      :', Jet_zeta
      write (lid,*)
      write (lid,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write (lid,*) 
    enddo
  endif
