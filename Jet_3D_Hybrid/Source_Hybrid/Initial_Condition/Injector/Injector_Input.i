  if (I_Am == First_Processor) &
    write (6,*) 'Reading .... ', TRIM(Injector_Input_File)
  
  lid = 78
  OPEN (UNIT=lid, FILE=TRIM(Injector_Input_File), STATUS='OLD', IOSTAT=IOS)
  
  read (lid,'(//)')
  read (lid,*) Jet_Direction
  read (lid,*)
  read (lid,*) Jet_Shape
  read (lid,*)
  read (lid,*) Jet_Center
  read (lid,*) Jet_Circle
  read (lid,*) Jet_Square
  read (lid,*)
  read (lid,*) Jet_Mach
  read (lid,*) Jet_Re
  read (lid,*)
  read (lid,*) Jet_Rho
  read (lid,*) Jet_U
  read (lid,*) Jet_P
  read (lid,*) Jet_T
  read (lid,*) Jet_f1
  read (lid,*) Jet_Gamma
  read (lid,*) Jet_R

  read (lid,*)
  read (lid,*) Jet_mu
  read (lid,*) Jet_kappa
  read (lid,*) Jet_zeta

  CLOSE (lid)
  
  Jet_C                = SQRT(Jet_Gamma*Jet_R*Jet_T)
  Jet_U(Jet_Direction) = Jet_Mach*Jet_C

if (Jet_Circle <  ZERO) &
  Jet_Circle           = Jet_Re*Jet_mu/(Jet_Rho*Jet_U(Jet_Direction))/2
  
  if (I_Am == First_Processor) &
    write (6,*) '....... Done ', TRIM(Injector_Input_File)

#include "Echo_Injector.i"
