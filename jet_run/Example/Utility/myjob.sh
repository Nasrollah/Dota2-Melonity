#!/bin/tcsh
limit stacksize unlimited
mpirun -r ssh -f mpd.hosts -machinefile machfile -np 32 ./Executable/sp*
#~/MPICH/MPICH2/bin/mpirun -np 8 ./Executable/sp*
