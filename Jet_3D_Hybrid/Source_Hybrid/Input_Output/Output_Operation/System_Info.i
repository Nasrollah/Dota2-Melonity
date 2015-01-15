  if (I_Am == First_Processor) then

       write (lid6 ,*) 
       write (lid6 ,*) '--------------------------------'

       write (lid60,*) 
       write (lid60,*) '--------------------------------'

     if (Navier_Stokes) then
       write (lid6 ,*) 'System is the Navier Stokes Equations'
       write (lid6 ,*) 
       if (Viscous_Dissipation) write (lid6 ,*) 'Viscous   Dissipation is ON'
       if (Thermal_Conduction)  write (lid6 ,*) 'Thermal   Conduction  is ON'
       if (Molecular_Diffusion) write (lid6 ,*) 'Molecular Diffusion   is ON'

       write (lid60,*) 'System is the Navier Stokes Equations'
       write (lid60,*) 
       if (Viscous_Dissipation) write (lid60,*) 'Viscous   Dissipation is ON'
       if (Thermal_Conduction)  write (lid60,*) 'Thermal   Conduction  is ON'
       if (Molecular_Diffusion) write (lid60,*) 'Molecular Diffusion   is ON'
     else
       write (lid6 ,*) 'System is the Euler Equations'
       write (lid60,*) 'System is the Euler Equations'
     endif

       write (lid6 ,*) '--------------------------------'
       write (lid6 ,*) 
       write (lid6 ,*) 'Final Time  : ', Final_Time
       write (lid6 ,*) '      Time  : ', Time
       write (lid6 ,*) 'Time Factor : ', Time_Factor
       write (lid6 ,*)
       write (lid6 ,*) '  dt_Euler  : ', Euler_dt
       write (lid6 ,*) '  dt_NS     : ',    NS_dt
       write (lid6 ,*) '  dt        : ', dt_Original
       write (lid6 ,*)

       write (lid60,*) '--------------------------------'
       write (lid60,*) 
       write (lid60,*) 'Final Time  : ', Final_Time
       write (lid60,*) '      Time  : ', Time
       write (lid60,*) 'Time Factor : ', Time_Factor
       write (lid60,*) 
       write (lid60,*) '  dt_Euler  : ', Euler_dt
       write (lid60,*) '  dt_NS     : ',    NS_dt
       write (lid60,*) '  dt        : ', dt_Original
       write (lid60,*) 

  endif
