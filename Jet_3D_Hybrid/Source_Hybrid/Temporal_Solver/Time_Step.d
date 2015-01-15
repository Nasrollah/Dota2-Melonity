REALTYPE :: Begin_Time    = ZERO                   ! s
REALTYPE :: Final_Time    = 8.0d0                  ! s

REALTYPE :: Hybrid_Start_Time = ZERO               ! s

integer  :: Step          = 0
REALTYPE :: Time          = ZERO                   ! s
integer  :: Restart       = 0
integer  :: Auto_Restart  = 0

REALTYPE :: CFL_Euler, CFL_NS
REALTYPE :: dt, dt_Original, dt_Restart, dt_Movie, dt_Analysis, dt_Output

REALTYPE :: Euler_dt, NS_dt

logical  :: Begin         = .TRUE. 
logical  :: Final         = .FALSE.

logical  :: Euler         = .TRUE.
logical  :: Scalar        = .FALSE.

REALTYPE :: Wave_Speed_x = ONE
REALTYPE :: Wave_Speed_y = ONE
REALTYPE :: Wave_Speed_z = ONE
