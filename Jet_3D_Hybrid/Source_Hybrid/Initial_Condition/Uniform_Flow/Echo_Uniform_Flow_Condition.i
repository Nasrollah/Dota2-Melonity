  if (I_Am == First_Processor) then
    do i = 1,2
      lid=6 ; if (i == 2) lid=60

      write (lid,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
      write (lid,*) 'Case Study              : ', Case_Study
      write (lid,*) '     Name               : ', TRIM(Case_Name)
      write (lid,*)
      write (lid,*) 'Gamma Formulation       : ', Gamma_Formulation
      write (lid,*) '   Gamma                : ', Gamma
      write (lid,*)
      write (lid,*) 'Gas Property                  '
      write (lid,*) '   Molecular Weight     : ',    MW_Gas_1, Gas(1)%MW
      write (lid,*) '   Rho                  : ',   Rho_Gas_1, Gas(1)%Rho
      write (lid,*) '   Gamma                : ', Gamma_Gas_1, Gas(1)%Gamma
      write (lid,*) '   Cp                   : ',    Cp_Gas_1, Gas(1)%Cp
      write (lid,*) '   Cv                   : ',    Cv_Gas_1, Gas(1)%Cv
      write (lid,*)
      write (lid,*) '     Mach               :', Mach
      write (lid,*)
      write (lid,*) '     Rho                :', Rho_Gas_1
      write (lid,*) '     U                  :',   U_Gas_1
      write (lid,*) '     P                  :',   P_Gas_1
      write (lid,*) '     C                  :',   C_Gas_1
      write (lid,*)
      write (lid,*) '     mu_0               :',    mu_0
      write (lid,*) '  kappa_0               :', kappa_0
      write (lid,*) '   zeta_0               :',  zeta_0
      write (lid,*)
      write (lid,*) '  R_Universal           :', R_Universal
      write (lid,*) '  R_Specific            :', Gas(1)%R
      write (lid,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write (lid,*) 
    enddo
  endif
