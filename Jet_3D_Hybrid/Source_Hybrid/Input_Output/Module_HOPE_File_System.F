MODULE HOPE_File_System
  
  USE PseudoPack, ONLY : I_Am, First_Processor, lid6
  
IMPLICIT NONE

Character(LEN=256) ::       Root_Dir
Character(LEN=256) ::      Input_Dir
  
Character(LEN=256) :: Input_Root_Dir
Character(LEN=256) ::       Case_Dir
Character(LEN=256) ::      WENO_Dir,      WENO_Input,      WENO_Input_File
Character(LEN=256) ::    Output_Dir,    Output_Input,    Output_Input_File
Character(LEN=256) ::     Shock_Dir,     Shock_Input,     Shock_Input_File
Character(LEN=256) ::        RT_Dir,        RT_Input,        RT_Input_File
Character(LEN=256) ::        RM_Dir,        RM_Input,        RM_Input_File
Character(LEN=256) ::  Particle_Dir,  Particle_Input,  Particle_Input_File
Character(LEN=256) :: Advection_Dir, Advection_Input, Advection_Input_File
Character(LEN=256) ::   Uniform_Dir,   Uniform_Input,   Uniform_Input_File
Character(LEN=256) ::  Injector_Dir,  Injector_Input,  Injector_Input_File
Character(LEN=256) :: Membrance_Dir, Membrance_Input, Membrance_Input_File
Character(LEN=256) ::     Notch_Dir,     Notch_Input,     Notch_Input_File
Character(LEN=256) ::      Mode_Dir,      Mode_Input,      Mode_Input_File
Character(LEN=256) ::    Random_Dir,    Random_Input,    Random_Input_File
Character(LEN=256) ::   Bubbles_Dir,   Bubbles_Input,   Bubbles_Input_File
Character(LEN=256) ::        MR_Dir,        MR_Input,        MR_Input_File
Character(LEN=256) :: xAnalysis_Dir,  Analysis_Input,  Analysis_Input_File
  
Character(LEN=256) :: Output_Root_Dir
Character(LEN=256) :: Prefix
Character(LEN=256) ::      Plot_Dir,      Plot_File ,   Plot_BlowUp_File
Character(LEN=256) ::                  Plot_IC_File ,    Plot_Final_File
Character(LEN=256) ::  Topology_Dir,                       Topology_File
Character(LEN=256) ::                            Processor_Topology_File
Character(LEN=256) ::      Misc_Dir,      Misc_File ,        Output_File
Character(LEN=256) ::     Movie_Dir,     Movie_File
Character(LEN=256) ::   MR_Flag_Dir,   MR_Flag_File ,    MR_Flag_IC_File
Character(LEN=256) ::                                 MR_Flag_Final_File
Character(LEN=256) ::  Analysis_Dir,  Analysis_File 
Character(LEN=256) ::             Mixing_Layer_File
Character(LEN=256) :: Flowfield_Dir, Flowfield_File 
Character(LEN=256) :: Statistic_Dir, Statistic_File 
Character(LEN=256) ::  Spectrum_Dir,  Spectrum_File 
Character(LEN=256) ::   Profile_Dir,   Profile_File 
Character(LEN=256) ::       PDF_Dir,       PDF_File 
Character(LEN=256) ::   Restart_Dir,   Restart_File ,  Auto_Restart_File
Character(LEN=256) ::                                  Restart_List_File
Character(LEN=256) ::     Restart_Data_Dir,            Restart_Data_File
Character(LEN=256) ::  Adaptive_Domain_Dir,         Adaptive_Domain_File
Character(LEN=256) ::    Analysis_Data_Dir,           Analysis_Data_File
  
PUBLIC :: File_System_Setup
  
CONTAINS
  
Subroutine File_System_Setup (Dimensions)
  
  integer :: Dimensions
  
  call  Input_Filenames (Dimensions)
  
  call Output_Filenames 
  
END Subroutine File_System_Setup
!
! ==========================================================================
!
Subroutine Input_Filenames (Dimensions)
  
  integer            :: Dimensions
  integer            :: lid, IOS
  Character(LEN=256) :: Mode_2D, Mode_3D, Notch_2D, Notch_3D, Filename
  
  Filename = 'Input/Input_Filenames.input'

  if (I_Am == First_Processor) &
    write (6,*) 'Reading .... ', TRIM(Filename)

  lid = 101 
  OPEN (Unit=lid, FILE=TRIM(Filename), STATUS='OLD', IOSTAT=IOS)

  read (lid,'(/////)') 
  read (lid,*) Input_Root_Dir
  read (lid,*)       Case_Dir
  read (lid,'(//)')    
  read (lid,*)      WENO_Dir ; read (lid,*)      WENO_Input  
  read (lid,'(/)')     
  read (lid,*)    Output_Dir ; read (lid,*)    Output_Input  
  read (lid,'(////)')  
  read (lid,*)     Shock_Dir ; read (lid,*)     Shock_Input  
  read (lid,'(/)')     
  read (lid,*)        RT_Dir ; read (lid,*)        RT_Input  
  read (lid,'(/)')     
  read (lid,*)        RM_Dir ; read (lid,*)        RM_Input  
  read (lid,'(/)')     
  read (lid,*)  Particle_Dir ; read (lid,*)  Particle_Input  
  read (lid,'(/)')     
  read (lid,*) Advection_Dir ; read (lid,*) Advection_Input  
  read (lid,'(/)')
  read (lid,*)   Uniform_Dir ; read (lid,*)   Uniform_Input
  read (lid,'(/)')
  read (lid,*)  Injector_Dir ; read (lid,*)  Injector_Input
  read (lid,'(/////)') 
  read (lid,*) Membrance_Dir ; read (lid,*) Membrance_Input  
  read (lid,'(/)')     
  read (lid,*)     Notch_Dir ; read (lid,*) Notch_2D ; read (lid,*) Notch_3D 
  read (lid,'(/)')
  read (lid,*)      Mode_Dir ; read (lid,*)  Mode_2D ; read (lid,*)  Mode_3D 
  read (lid,'(/)')
  read (lid,*)    Random_Dir ; read (lid,*)    Random_Input
  read (lid,'(/)')
  read (lid,*)   Bubbles_Dir ; read (lid,*)   Bubbles_Input
  read (lid,'(/)')
  read (lid,*)        MR_Dir ; read (lid,*)        MR_Input 
  read (lid,'(/////)') 
  read (lid,*)  Analysis_Dir ; read (lid,*)  Analysis_Input  
  read (lid,'(/)')     
  
  Notch_Input = Notch_2D ; if (Dimensions == 3) Notch_Input = Notch_3D
   Mode_Input =  Mode_2D ; if (Dimensions == 3)  Mode_Input =  Mode_3D
  
                                   Root_Dir = '' 
  if (TRIM(Input_Root_Dir) /= ' ') Root_Dir = TRIM(Input_Root_Dir)//'/'  

  SELECT CASE (Dimensions)
    CASE DEFAULT
      Root_Dir = TRIM(Root_Dir)//'2D/'

    CASE (1)
      Root_Dir = TRIM(Root_Dir)//'1D/'

    CASE (3)
      Root_Dir = TRIM(Root_Dir)//'3D/'

  END SELECT

  if (TRIM(      Case_Dir) /= ' ') Root_Dir = TRIM(Root_Dir)//TRIM(Case_Dir)//'/'

  Input_Dir = TRIM(Root_Dir)

  call Form_Path (Root_Dir,      WENO_Dir,      WENO_Input,      WENO_Input_File)
  call Form_Path (Root_Dir,    Output_Dir,    Output_Input,    Output_Input_File)
  call Form_Path (Root_Dir,     Shock_Dir,     Shock_Input,     Shock_Input_File)
  call Form_Path (Root_Dir,        RT_Dir,        RT_Input,        RT_Input_File)
  call Form_Path (Root_Dir,        RM_Dir,        RM_Input,        RM_Input_File)
  call Form_Path (Root_Dir,  Particle_Dir,  Particle_Input,  Particle_Input_File)
  call Form_Path (Root_Dir, Advection_Dir, Advection_Input, Advection_Input_File)
  call Form_Path (Root_Dir,   Uniform_Dir,   Uniform_Input,   Uniform_Input_File)
  call Form_Path (Root_Dir,  Injector_Dir,  Injector_Input,  Injector_Input_File)
  call Form_Path (Root_Dir, Membrance_Dir, Membrance_Input, Membrance_Input_File)
  call Form_Path (Root_Dir,     Notch_Dir,     Notch_Input,     Notch_Input_File)
  call Form_Path (Root_Dir,      Mode_Dir,      Mode_Input,      Mode_Input_File)
  call Form_Path (Root_Dir,    Random_Dir,    Random_Input,    Random_Input_File)
  call Form_Path (Root_Dir,   Bubbles_Dir,   Bubbles_Input,   Bubbles_Input_File)
  call Form_Path (Root_Dir,        MR_Dir,        MR_Input,        MR_Input_File)
  call Form_Path (Root_Dir,  Analysis_Dir,  Analysis_Input,  Analysis_Input_File)
  
  CLOSE (lid)
  
  if (I_Am == First_Processor) &
    write (6,*) '....... Done ', TRIM(Filename)

  call Input_File_List (lid6)
  
  CONTAINS
  
  Subroutine Form_Path (Root, Dir, Input, Filename)
  
  Character(LEN=*)   :: Root, Dir, Input, Filename
  Character(LEN=256) :: Path
  
                            Path = TRIM(Root)
    if (TRIM(Dir  ) /= ' ') Path = TRIM(Path)//TRIM(Dir  )//'/' 
    if (TRIM(Input) /= ' ') Path = TRIM(Path)//TRIM(Input)
  
    Filename = TRIM(Path)
  
  END Subroutine Form_Path 
!
! -----------------------------------------------------------
!
  Subroutine Input_File_List (lid)
   
    integer  :: lid
  
    if (I_Am /= First_Processor) RETURN
  
    write (lid,*) '---------------------------------------------------'
    write (lid,*) 'Input File Directory and Name'
    write (lid,*)
    write (lid,*) TRIM(Root_Dir)
    write (lid,*) TRIM(     WENO_Input_File)
    write (lid,*) TRIM(   Output_Input_File)
    write (lid,*) TRIM(    Shock_Input_File)
    write (lid,*) TRIM(       RT_Input_File)
    write (lid,*) TRIM(       RM_Input_File)
    write (lid,*) TRIM( Particle_Input_File)
    write (lid,*) TRIM(Advection_Input_File)
    write (lid,*) TRIM(  Uniform_Input_File)
    write (lid,*) TRIM( Injector_Input_File)
    write (lid,*) TRIM(Membrance_Input_File)
    write (lid,*) TRIM(    Notch_Input_File)
    write (lid,*) TRIM(     Mode_Input_File)
    write (lid,*) TRIM(   Random_Input_File)
    write (lid,*) TRIM(  Bubbles_Input_File)
    write (lid,*) TRIM(       MR_Input_File)
    write (lid,*) TRIM( Analysis_Input_File)
    write (lid,*)
    write (lid,*) '---------------------------------------------------'
  
  END Subroutine Input_File_List 
  
END Subroutine Input_Filenames 
!
!--------------------------------------------------------------------------
!
Subroutine Output_Filenames 
  
  integer            :: lid, IOS
  Character(LEN=256) :: Tmp_File, Filename
  
  Filename = TRIM(Input_Dir)//'Output_Filenames.input'

  if (I_Am == First_Processor) &
    write (6,*) 'Reading .... ', TRIM(Filename)

  lid = 102 
  OPEN (Unit=lid, FILE=TRIM(Filename), STATUS='OLD', IOSTAT=IOS)
  
  read (lid,'(/////)') ; read (lid,*) Output_Root_Dir
  read (lid,'(//)')    ; read (lid,*) Prefix 
  read (lid,'(/)')     ; read (lid,*)            Plot_Dir 
                         read (lid,*)           Movie_Dir 
                         read (lid,*)         MR_Flag_Dir 
  read (lid,'(/)')     ; read (lid,*)        Analysis_Dir  
                         read (lid,*)   Analysis_Data_Dir 
                         read (lid,*)       Flowfield_Dir 
                         read (lid,*)       Statistic_Dir 
                         read (lid,*)             PDF_Dir 
                         read (lid,*)         Profile_Dir 
                         read (lid,*)        Spectrum_Dir 
  read (lid,'(/)')     ; read (lid,*)         Restart_Dir 
                         read (lid,*)    Restart_Data_Dir
  read (lid,'(/)')     ; read (lid,*)        Topology_Dir
  read (lid,'(/)')     ; read (lid,*) Adaptive_Domain_Dir
  read (lid,'(/)')     ; read (lid,*)            Misc_Dir 
  
                                    Root_Dir = '' 
  if (TRIM(Output_Root_Dir) /= ' ') Root_Dir = TRIM(Output_Root_Dir)//'/'  
  
  Output_Dir = Root_Dir

  call Form_PATH (Root_Dir,            Plot_Dir, Prefix,            Plot_File)
  call Form_PATH (Root_Dir,           Movie_Dir, Prefix,           Movie_File)
  call Form_PATH (Root_Dir,         MR_Flag_Dir, Prefix,         MR_Flag_File)
  call Form_PATH (Root_Dir,        Analysis_Dir, Prefix,        Analysis_File)
  call Form_PATH (Root_Dir,   Analysis_Data_Dir, Prefix,   Analysis_Data_File)
  call Form_PATH (Root_Dir,       Flowfield_Dir, Prefix,       Flowfield_File)
  call Form_PATH (Root_Dir,       Statistic_Dir, Prefix,       Statistic_File)
  call Form_PATH (Root_Dir,             PDF_Dir, Prefix,             PDF_File)
  call Form_PATH (Root_Dir,         Profile_Dir, Prefix,         Profile_File)
  call Form_PATH (Root_Dir,        Spectrum_Dir, Prefix,        Spectrum_File)
  call Form_PATH (Root_Dir,         Restart_Dir, Prefix,         Restart_File)
  call Form_PATH (Root_Dir,    Restart_Data_Dir, Prefix,    Restart_Data_File)
  call Form_PATH (Root_Dir, Adaptive_Domain_Dir, Prefix, Adaptive_Domain_File)
  call Form_PATH (Root_Dir,        Topology_Dir, Prefix,        Topology_File)
  call Form_PATH (Root_Dir,            Misc_Dir, Prefix,            Misc_File)
  
             Plot_IC_File = TRIM(           Plot_File)//'IC'
          Plot_Final_File = TRIM(           Plot_File)//'Final'
         Plot_BlowUp_File = TRIM(           Plot_File)//'BlowUp'
          MR_Flag_IC_File = TRIM(        MR_Flag_File)//'IC'
       MR_Flag_Final_File = TRIM(        MR_Flag_File)//'Final'
            Analysis_File = TRIM(       Analysis_File)//'Analysis'
       Analysis_Data_File = TRIM(  Analysis_Data_File)//'Analysis'
        Mixing_Layer_File = TRIM(       Analysis_File)//'_Mixing_Layer.dat'
           Flowfield_File = TRIM(      Flowfield_File)//'Flowfield.dat'
                 PDF_File = TRIM(            PDF_File)//'PDF'
           Statistic_File = TRIM(      Statistic_File)//'Statistic'
             Profile_File = TRIM(        Profile_File)//'Profile'
            Spectrum_File = TRIM(       Spectrum_File)//'Spectrum'
        Auto_Restart_File = TRIM(        Restart_File)//'Auto_Restart_File.dat'
        Restart_List_File = TRIM(        Restart_File)//'Restart_List_File.dat'
        Restart_Data_File = TRIM(   Restart_Data_File)//'Restart'
  Processor_Topology_File = TRIM(       Topology_File)//'Processor_Topology.dat'
     Adaptive_Domain_File = TRIM(Adaptive_Domain_File)//'Adaptive_Domain'
              Output_File = TRIM(           Misc_File)//'Output'
                Misc_File = TRIM(           Misc_File)//' '
  
  CLOSE (lid)
  
  if (I_Am == First_Processor) &
    write (6,*) '....... Done ', TRIM(Filename)

  call Output_File_List (lid6)
  
  CONTAINS
  
  Subroutine Form_Path (Root, Dir, Prefix, Filename)
  
    Character(LEN=*)   :: Root, Dir, Prefix, Filename
    Character(LEN=256) :: Path
  
                             Path = TRIM(Root)
    if (TRIM(Dir   ) /= ' ') Path = TRIM(Path)//TRIM(Dir)//'/' 
    if (TRIM(Prefix) /= ' ') Path = TRIM(Path)//TRIM(Prefix)//'_'
  
    Filename = TRIM(Path)
  
  END Subroutine Form_Path
!
! -----------------------------------------------------------
!
  Subroutine Output_File_List (lid)
   
    integer  :: lid
  
    if (I_Am /= First_Processor) RETURN
  
    write (lid,*) '---------------------------------------------------'
    write (lid,*) 'Output File Directory and Prefix'
    write (lid,*)
    write (lid,*) TRIM(Root_Dir)
    write (lid,*) TRIM(Plot_File)
    write (lid,*) TRIM(Plot_IC_File)
    write (lid,*) TRIM(Plot_Final_File)
    write (lid,*) TRIM(Plot_BlowUp_File)
    write (lid,*) TRIM(Movie_File)
    write (lid,*) TRIM(MR_Flag_File)
    write (lid,*) TRIM(MR_Flag_IC_File)
    write (lid,*) TRIM(MR_Flag_Final_File)
    write (lid,*) TRIM(Analysis_File)
    write (lid,*) TRIM(Analysis_Data_File)
    write (lid,*) TRIM(Mixing_Layer_File)
    write (lid,*) TRIM(Flowfield_File)
    write (lid,*) TRIM(PDF_File)
    write (lid,*) TRIM(Statistic_File)
    write (lid,*) TRIM(Profile_File)
    write (lid,*) TRIM(Spectrum_File)
    write (lid,*) TRIM(Auto_Restart_File)
    write (lid,*) TRIM(Restart_List_File)
    write (lid,*) TRIM(Restart_Data_File)
    write (lid,*) TRIM(Processor_Topology_File)
    write (lid,*) TRIM(Adaptive_Domain_File)
    write (lid,*) TRIM(Output_File)
    write (lid,*) TRIM(Misc_File)
    write (lid,*)
    write (lid,*) '---------------------------------------------------'
  
  END Subroutine Output_File_List 
  
END Subroutine Output_Filenames 
  
END MODULE HOPE_File_System
