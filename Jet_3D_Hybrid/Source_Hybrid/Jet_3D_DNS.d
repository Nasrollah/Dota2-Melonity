integer            :: Spatial_Order = 5
integer            :: Case_Study    = 1           ! 1=RM  2=RT
Character(LEN=60)  :: Case_Name     = 'RM'
integer            :: Flux_Type     = 1
integer            :: Source        = 0

integer            :: Write_IC              = 1
integer            :: Write_Final           = 1
integer            :: Write_Initial_Movie   = 0
integer            :: Write_Initial_MR_Flag = 0

integer            ::      MR_Output_Screen     = 1
integer            :: Restart_Distributed_Read  = 0
integer            :: Restart_Distributed_Write = 0

REALTYPE           ::    Movie_Resolution = 1.0d-5
REALTYPE           ::  Restart_Resolution = 1.0d-5

REALTYPE           :: D_F_Zero 
integer            :: Adjust_D_Flux

REALTYPE           ::   F_Zero 
integer            :: Adjust_Mass_Fraction, Check_Mass_Fraction

REALTYPE           :: Atwood_Number
Character(LEN=256) :: Filename

integer  :: Number_Of_Species = 1
integer  :: i, j, k, N, M, ND, NF, NS, NV, N_Inc

REALTYPE :: CPU_Euler, CPU_Navier_Stokes, CPU_t0, CPU_t1

logical  :: Error

#include "Declaration.d"

integer  :: MPI_Comm_1D

TYPE (Processor_Topology)   :: Topology_1D

