!--------------
#if 1
!--------------




Gamma_Gas = HUGE(ONE)

     Case_Name = 'Uniform Air Flow'

     Gas(1) = Air     ; if (Unit_Conversion) call Gas_Unit_Conversion (Gas(1))
     Gas(2) = Air     ; if (Unit_Conversion) call Gas_Unit_Conversion (Gas(2))

           MW_Gas_1 = Gas(1)%MW           ;    MW_Gas_2 = Gas(2)%MW 
          Rho_Gas_1 = Gas(1)%Rho          ;   Rho_Gas_2 = Gas(2)%Rho
        Gamma_Gas_1 = Gas(1)%Gamma        ; Gamma_Gas_2 = Gas(2)%Gamma
           Cp_Gas_1 = Gas(1)%Cp           ;    Cp_Gas_2 = Gas(2)%Cp

          Rho_Gas_1 =  Density_Uniform    ;   Rho_Gas_2 = Rho_Gas_1
            P_Gas_1 = Pressure_Uniform    ;     P_Gas_2 =   P_Gas_1

!---------
if (Gamma_Formulation == 2) then   ! Constant Gamma
  Gamma_Gas_1 = Gamma ; Gamma_Gas_2 = Gamma 
endif

Cv_Gas_1 = Cp_Gas_1/Gamma_Gas_1 ; Cv_Gas_2 = Cp_Gas_2/Gamma_Gas_2

Gas(1)%f     = ONE                    ; Gas(2)%f     = ONE
Gas(1)%Gamma = Gamma_Gas_1            ; Gas(2)%Gamma = Gamma_Gas_2 
Gas(1)%Cv    = Gas(1)%Cp/Gas(1)%Gamma ; Gas(2)%Cv    = Gas(2)%Cp/Gas(2)%Gamma 

!---------
Gamma_Gas(1,1) = Gas(1)%Gamma ; Gamma_Gas(2,1) = Gas(2)%Gamma
Gamma_Gas(1,2) = Gas(1)%Cp    ; Gamma_Gas(2,2) = Gas(2)%Cp
Gamma_Gas(1,3) = Gas(1)%Cv    ; Gamma_Gas(2,3) = Gas(2)%Cv
Gamma_Gas(1,4) = Gas(1)%MW    ; Gamma_Gas(2,4) = Gas(2)%MW
Gamma_Gas(1,5) = Gas(1)%Rho   ; Gamma_Gas(2,5) = Gas(2)%Rho
!---------

if (Gas_Composition < 0) then  ! Switch the gases
  Gamma_Gas(257,:) = Gamma_Gas(  2,:) ; Gas(9) = Gas(2) 
  Gamma_Gas(  2,:) = Gamma_Gas(  1,:) ; Gas(2) = Gas(1) 
  Gamma_Gas(  1,:) = Gamma_Gas(257,:) ; Gas(1) = Gas(9)
endif








!--------------
#else
!--------------

Gamma_Gas = HUGE(ONE)

     Case_Name = 'Uniform Air Flow'

        Gas(1) = Air              ; Gas(2) = Air

           MW_Gas_1 = Air%MW/1000      ;    MW_Gas_2 =    MW_Gas_1
          Rho_Gas_1 = Density_Uniform  ;   Rho_Gas_2 =   Rho_Gas_1
!         Rho_Gas_1 = Air%Rho          ;   Rho_Gas_2 =   Rho_Gas_1
        Gamma_Gas_1 = Air%Gamma        ; Gamma_Gas_2 = Gamma_Gas_1
            P_Gas_1 = Pressure_Uniform ;     P_Gas_2 =     P_Gas_1
           Cp_Gas_1 = Air%Cp           ;    Cp_Gas_2 =    Cp_Gas_1

if (Gamma_Formulation == 2) then   ! Constant Gamma
!#if defined (NAVIER_STOKES)
! if (Navier_Stokes) Cp_Gas_2 = (MW_Gas_1/MW_Gas_2)*Cp_Gas_1
!#endif

  Gamma_Gas_1 = Gamma ; Gamma_Gas_2 = Gamma 
endif

Cv_Gas_1 = Cp_Gas_1/Gamma_Gas_1 ; Cv_Gas_2 = Cp_Gas_2/Gamma_Gas_2

Gamma_Gas(1,1) = Gamma_Gas_1 ; Gamma_Gas(2,1) = Gamma_Gas_2
Gamma_Gas(1,2) =    Cp_Gas_1 ; Gamma_Gas(2,2) =    Cp_Gas_2
Gamma_Gas(1,3) =    Cv_Gas_1 ; Gamma_Gas(2,3) =    Cv_Gas_2
Gamma_Gas(1,4) =    MW_Gas_1 ; Gamma_Gas(2,4) =    MW_Gas_2
Gamma_Gas(1,5) =   Rho_Gas_1 ; Gamma_Gas(2,5) =   Rho_Gas_2

Gas(1)%f     = ONE            ; Gas(2)%f     = ONE
Gas(1)%Gamma = Gamma_Gas_1    ; Gas(2)%Gamma = Gamma_Gas_2 
Gas(1)%Cp    =    Cp_Gas_1    ; Gas(2)%Cp    =    Cp_Gas_2 
Gas(1)%Cv    =    Cv_Gas_1    ; Gas(2)%Cv    =    Cv_Gas_2
Gas(1)%MW    =    MW_Gas_1    ; Gas(2)%MW    =    MW_Gas_2
Gas(1)%Rho   =   Rho_Gas_1    ; Gas(2)%Rho   =   Rho_Gas_2

if (Gas_Composition < 0) then  ! Switch the gases
  Gamma_Gas(257,:) = Gamma_Gas(  2,:) ; Gas(9) = Gas(2) 
  Gamma_Gas(  2,:) = Gamma_Gas(  1,:) ; Gas(2) = Gas(1) 
  Gamma_Gas(  1,:) = Gamma_Gas(257,:) ; Gas(1) = Gas(9)
endif

!--------------
#endif
!--------------
