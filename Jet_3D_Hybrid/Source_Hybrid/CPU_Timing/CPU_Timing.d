REALTYPE :: Time_Factor   = 1.0d6                  ! micro-second

REALTYPE :: CPU_Start, CPU_RK_Begin, CPU_RK_End, CPU_Usage

integer , dimension(N_Clock) :: CPU_Used

integer  :: Max_Output_Step = 10000
integer  :: CPU_Step_Factor = 10
