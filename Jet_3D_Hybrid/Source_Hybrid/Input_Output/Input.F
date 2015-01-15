!
!===================================================================
!
  Subroutine Input (Unit)

  IMPLICIT NONE

  integer , OPTIONAL :: Unit
  integer            :: lid, IOS
  integer            :: System

  if (I_Am == First_Processor) &
    write (6,*) 'Reading .... ', TRIM(WENO_Input_File)

  if (PRESENT(Unit)) then ; lid=Unit ; else ; lid = 7  ; endif

  OPEN (UNIT=lid, FILE=TRIM(WENO_Input_File), STATUS='OLD', IOSTAT=IOS)

  read (lid,'(//)')
  read (lid,*) System
  read (lid,*) 
  read (lid,*) Case_Study

                           Euler = .TRUE. 
                   Navier_Stokes = .FALSE.
  if (System == 1) Navier_Stokes = .TRUE.

#include "WENO_Input.i"

  read (lid,'(//)')
  read (lid,*) Distributed_x, Periodicity_x, Topology_Grid_x
  read (lid,*) Distributed_y, Periodicity_y, Topology_Grid_y
  read (lid,*) Distributed_z, Periodicity_z, Topology_Grid_z

  read (lid,'(///)')
  read (lid,*) BC_On_Off_x, BC_Symmetry_x,  BC_Type_x
  read (lid,*) BC_On_Off_y, BC_Symmetry_y,  BC_Type_y
  read (lid,*) BC_On_Off_z, BC_Symmetry_z,  BC_Type_z

  read (lid,'(//)')
  read (lid,*) D_F_Zero, Adjust_D_Flux
  read (lid,*)   F_Zero, Adjust_Mass_Fraction, Check_Mass_Fraction

  CLOSE (lid)

  Filename = TRIM(Input_Dir)//'Gas_Properties_Database.input'
  call Read_Gas_Properties (Filename)

#include "Output_Control.i"
#include "Open_Output_File.i"

#include "Uniform_Flow_Input_File.i"

#if 1
#include "Injector_Input.i"
#else
   call Injector_Input (Injector_Input_File)
#endif


#if defined (MULTI_RESOLUTION)
#include "MR_Input_File.i"
#endif

  if (CFL_Euler >  0.450d0) then
    if (I_Am == First_Processor) write (lid6,*) 'WARNING :  CFL_Euler > 0.45'
  endif

#if defined (NAVIER_STOKES)
  if ((Navier_Stokes) .AND. (CFL_NS    >  1.0d0)) then
    if (I_Am == First_Processor) write (lid6,*) 'WARNING :  CFL_NS    > 0.50'
  endif
#endif

#include "Read_Auto_Restart_File.i"

  if (I_Am == First_Processor) &
    write (6,*) '....... Done ', TRIM(WENO_Input_File)

  End Subroutine Input
