#include "pseudopack.h"


MODULE Graphic_Parameter

implicit NONE

logical  :: Print_IC    = .TRUE.
logical  :: Print_Final = .TRUE.

logical  :: Movie       = .FALSE.
logical  :: Movie_Q     = .FALSE. , Start_Movie_Q = .TRUE.
logical  :: Movie_V     = .FALSE. , Start_Movie_V = .TRUE.

integer  :: Movie_Interval = 5
integer  :: Movie_Step
REALTYPE :: Movie_Time, Time_Movie

integer,  parameter, PRIVATE :: N_Fluid = 5
integer,             PRIVATE :: i

REALTYPE :: x0_Plot , x1_Plot , y0_Plot , y1_Plot , z0_Plot , z1_Plot
REALTYPE :: x0_Movie, x1_Movie, y0_Movie, y1_Movie, z0_Movie, z1_Movie

integer, dimension(40) :: Variable_Map

Character(LEN=40)      ::      Plot_Kind     = 'Primitive'
Character(LEN=256)     ::      Plot_Title    = 'Combustion'
Logical, dimension(40) ::      Plot_Variable = (/ (.TRUE., i=1,40) /)

Character(LEN=40)      ::     Movie_Kind     = 'Primitive'
Character(LEN=256)     ::     Movie_Title    = 'Combustion'
Logical, dimension(40) ::     Movie_Variable = (/ (.TRUE., i=1,40) /)

Character(LEN=40)      :: Vorticity_Kind     = ' '
Character(LEN=256)     :: Vorticity_Title    = 'Vorticity - The Movie'
Logical, dimension(3)  :: Vorticity_Variable = (/ (.TRUE., i=1,3)/)

Character(LEN=40)      :: Viscosity_Kind     = ' '
Character(LEN=256)     :: Viscosity_Title    = 'Mu, Kappa, T'
Logical, dimension(3)  :: Viscosity_Variable = (/ (.TRUE., i=1,3)/)

Character(LEN=40)      ::     Debug_Kind     = ' '
Character(LEN=256)     ::     Debug_Title    = 'Debug - XX'
Logical, dimension(1)  ::     Debug_Variable = (/ (.TRUE., i=1,1)/)

Character(LEN=7), parameter :: Ch_Rho  = 'Rho'    
Character(LEN=7), parameter :: Ch_RhoU = 'Rho_U' , Ch_U = 'U' 
Character(LEN=7), parameter :: Ch_RhoV = 'Rho_V' , Ch_V = 'V' 
Character(LEN=7), parameter :: Ch_RhoW = 'Rho_W' , Ch_W = 'W' 
Character(LEN=7), parameter :: Ch_E    = 'E'     , Ch_T = 'T'
Character(LEN=7), parameter :: Ch_Rhof = 'Rho_f_', Ch_f = 'f_'

Character(LEN=7), parameter :: Ch_Mu   = 'Mu'
Character(LEN=7), parameter :: Ch_Kappa= 'Kappa'

Character(LEN=7), parameter :: Ch_X    = 'XX'

Character(LEN=7), parameter :: Ch_Omega_x = 'Omega_x'
Character(LEN=7), parameter :: Ch_Omega_y = 'Omega_y'
Character(LEN=7), parameter :: Ch_Omega_z = 'Omega_z'

Character(LEN=7), dimension(N_Fluid+1), parameter ::   &
  Variable_Name_C = (/ Ch_Rho, Ch_RhoU, Ch_RhoV, Ch_RhoW, Ch_E   ,   &
                       Ch_Rhof /)

Character(LEN=7), dimension(N_Fluid+1), parameter ::   &
  Variable_Name_P = (/ Ch_Rho, Ch_U   , Ch_V   , Ch_W   , Ch_T   ,   &
                       Ch_f    /)

Character(LEN=7), dimension(3), parameter ::   &
  Variable_Name_M = (/ Ch_Mu, Ch_Kappa, Ch_T /)

Character(LEN=7), dimension(1), parameter ::   &
  Variable_Name_D = (/ Ch_X /)

#if defined (THREE_DIMENSION) 
Character(LEN=7), dimension(3), parameter ::   &
  Variable_Name_V = (/ Ch_Omega_z, Ch_Omega_y, Ch_Omega_x /)
#endif
#if defined (TWO_DIMENSION) 
Character(LEN=7), dimension(1), parameter ::   &
  Variable_Name_V = (/ Ch_Omega_z /)
#endif
#if defined (ONE_DIMENSION) 
Character(LEN=7), dimension(3), parameter ::   &
  Variable_Name_V = (/ Ch_Omega_z /)
#endif

Character(LEN=7), dimension(N_Fluid+1)            ::   &
  Variable_Name_G = (/ Ch_Rho, Ch_U   , Ch_V   , Ch_W   , Ch_T   ,  &
                       Ch_f    /)

CONTAINS
!
!=======================================================================
!
Subroutine PS_TecPlot_Title (Dimension, NCS, NV,              &
                             V_Title, V_Kind, V_Plot, V_Name, &
                             Step, Time, lid, New_Graphic_File)

REALTYPE                          :: Time
integer                           :: Dimension, NCS, NV, Step, lid
character(LEN=40)                 :: V_Kind
character(LEN=256)                :: V_Title
character(LEN=7)  , dimension(*)  :: V_Name
logical           , dimension(*)  :: V_Plot
character(LEN=256)                :: Title, V_List, Var_Name
logical , OPTIONAL                :: New_Graphic_File

integer  :: V, L

if (.NOT. PRESENT(New_Graphic_File)) RETURN
if (.NOT.         New_Graphic_File ) RETURN 

New_Graphic_File = .FALSE. 

SELECT CASE (Dimension)
  CASE (1)
    Title='Title="1-D' ; V_List='Variables=x,'

  CASE (2)
    Title='Title="2-D' ; V_List='Variables=x, y,'

  CASE (3)
    Title='Title="3-D' ; V_List='Variables=x, y, z,'
END SELECT

Title = TRIM(Title)//' '//TRIM(V_Title)//' - '//TRIM(V_Kind)//' : '

if ((V_Kind == 'Conservative') .OR. (V_Kind == 'Primitive')) then
  do V = 1,N_Fluid+NCS
    if (.NOT. V_Plot(V)) CYCLE
 
    if (V <= 5) then ; L = V ; else ; L = 6 ; endif

    if (V_Kind == 'Conservative') Var_Name = TRIM(Variable_Name_C(L))
    if (V_Kind == 'Primitive'   ) Var_Name = TRIM(Variable_Name_P(L))

    if (V >  5) Var_Name = TRIM(Var_Name)//CHAR_N(V-5)

    V_List = TRIM(V_List)//' '//TRIM(Var_Name)//','
  enddo
else
  do V = 1,NV
    if (.NOT. V_Plot(V)) CYCLE
    
    Var_Name = TRIM(V_Name(V))

    V_List = TRIM(V_List)//' '//TRIM(Var_Name)//','
  enddo
endif

L = LEN_TRIM(V_List) ; V_List(L:L) = ' '

write (lid,1000) TRIM(Title), Step, Time, TRIM(V_List)

 1000 format (1x,(a),' Step=',i5,' Time=',f10.7,'"'/1x,(a))

END Subroutine PS_TecPlot_Title 
!
!=======================================================================
!
Subroutine PS_TecPlot_Zone_Header (Dimension, Time, x, y, z, lid,  &
                                   Domain_Number, New_Graphic_File)

integer  :: Dimension, lid
REALTYPE :: Time

REALTYPE, dimension(:,:,:) :: x, y, z

logical , OPTIONAL :: New_Graphic_File
integer , OPTIONAL :: Domain_Number

integer          :: i, j, l, k, i0, i1, j0, j1, k0, k1
logical          :: L_New_Graphic_File
Character(LEN=9) :: Repeat_Factor 

i0 = LBOUND(x,DIM=1) ; i1 = UBOUND(x,DIM=1)
j0 = LBOUND(x,DIM=2) ; j1 = UBOUND(x,DIM=2)
k0 = LBOUND(x,DIM=3) ; k1 = UBOUND(x,DIM=3)

k = 0 ; if (PRESENT(Domain_Number)) k = Domain_Number 

L_New_Graphic_File = PRESENT(New_Graphic_File)

if ( (.NOT. L_New_Graphic_File) .OR.                    &
     (L_New_Graphic_File .AND. New_Graphic_File) ) then

  if (L_New_Graphic_File) New_Graphic_File = .FALSE.

    write (lid,1010) k, Time, ABS(i1-i0+1), ABS(j1-j0+1), ABS(k1-k0+1)

                      write (lid,1100) (((x(i,j,l) ,i=i0,i1),j=j0,j1),l=k0,k1)
  if (Dimension >= 2) write (lid,1100) (((y(i,j,l) ,i=i0,i1),j=j0,j1),l=k0,k1)
  if (Dimension == 3) write (lid,1100) (((z(i,j,l) ,i=i0,i1),j=j0,j1),l=k0,k1)

else

  if (Dimension == 1) Repeat_Factor = 'D=(1)'
  if (Dimension == 2) Repeat_Factor = 'D=(1,2)'
  if (Dimension == 3) Repeat_Factor = 'D=(1,2,3)'

  write (lid,1010) k, Time, ABS(i1-i0+1), ABS(j1-j0+1), ABS(k1-k0+1), &
                      TRIM(Repeat_Factor)

endif

 1010 format (1x,'Zone T="',i4,1x,G14.7,'", I=',i4,' J=',i4,' K=',i4,  &
           ' F=Block ',:,(a))
 1100 format (  8(1x,G14.7))

END Subroutine PS_TecPlot_Zone_Header 

END MODULE Graphic_Parameter
