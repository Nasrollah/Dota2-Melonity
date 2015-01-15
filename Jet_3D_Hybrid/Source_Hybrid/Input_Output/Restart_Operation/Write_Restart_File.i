#include "Restart_File_Write.i"

#include "Write_Auto_Restart_File.i"
#include "Write_Restart_List_File.i"

#if defined (PARALLEL_MPI)
call PS_MPI_Synchronize ()
#endif
