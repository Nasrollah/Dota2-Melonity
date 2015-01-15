MODULE MPIF

implicit NONE

#if defined (PARALLEL_MPI)
! include "mpif.h"
#else
integer  :: MPI_Comm_World = -1
integer  :: MPI_Comm_Null  = -1
#endif

END MODULE MPIF
