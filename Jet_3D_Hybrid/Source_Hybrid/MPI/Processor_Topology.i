                   Distributed(Index_x) =   Distributed_x 
                   Periodicity(Index_x) =   Periodicity_x
                          Grid(Index_x) = Topology_Grid_x

                   Distributed(Index_y) =   Distributed_y 
                   Periodicity(Index_y) =   Periodicity_y
                          Grid(Index_y) = Topology_Grid_y

                   Distributed(Index_z) =   Distributed_z 
                   Periodicity(Index_z) =   Periodicity_z
                          Grid(Index_z) = Topology_Grid_z

#if defined (PARALLEL_MPI)

#if 1

if (.NOT. Distributed_x) then
  Grid(Index_x) = 1 ; Grid(Index_y) = 0 ; Grid(Index_z) = 0

  do i = 1,8
    call MPI_DIMS_CREATE (N_Processor, Dimensions, Grid, MPI_Error_Status)

    j = Grid(Index_y) ; if (j  <=  N_y-1) EXIT

    Grid(Index_x) = 0 ; Grid(Index_y) = j/2 ; Grid(Index_z) = 0

    Distributed(Index_x) = .TRUE. 
  enddo
    Grid(Index_y) = 0 ; Grid(Index_z) = 0
else
  if ((Dimensions > 1) .AND. (MAXVAL(Grid(1:Dimensions)) == 0)) then
    call MPI_DIMS_CREATE (N_Processor, Dimensions, Grid, MPI_Error_Status)

    Grid(Index_x) = MINVAL(Grid) ; Grid(Index_y) = 0 ; Grid(Index_z) = 0
  endif
endif

#else

if (Distributed_x) then
  if ((Dimensions > 1) .AND. (MAXVAL(Grid(1:Dimensions)) == 0)) then
    call MPI_DIMS_CREATE (N_Processor, Dimensions, Grid, MPI_Error_Status)

    Grid(Index_x) = MINVAL(Grid) ; Grid(Index_y) = 0 ; Grid(Index_z) = 0
  endif
endif

#endif

#endif

#if 0
  if (I_Am == First_Processor) then
    write (6,*) '================================================='
    write (6,*) 'Info : Processor Cartesian Topology              '
    write (6,*) '          Grid(1)       : ', Grid(Index_x)
    write (6,*) '          Grid(2)       : ', Grid(Index_y)
    write (6,*) '================================================='
  endif
#endif

! ------------------------------------------------------------------
if ((Restart /= 0) .AND. (Restart_Distributed_Read == 1)) then
  call PS_MPI_Read_Processor_Topology (Dimensions, Grid, Periodicity,     &
                                       Processor_Topology_File, Error    )
endif

call PS_MPI_Processor_Topology_Setup (Topology, Dimensions, &
                                      Distributed, Periodicity, Grid)

Distributed_x = Distributed(Index_x)
Distributed_y = Distributed(Index_y)
Distributed_z = Distributed(Index_z)

MPI_CART_Communicator = Topology%MPI_Communicator

#if defined (PARALLEL_MPI)
call MPI_CART_GET    (Topology%MPI_Communicator, Topology%Dimension,        &
                      Topology%Grid            , Topology%Periodicity,      &
                      Topology%Coordination    , MPI_Error_Status)
#endif

call PS_MPI_Write_Processor_Topology (Topology%Dimension, Topology%Grid, &
                                      Topology%Periodicity,              &
                                      Topology%Coordination,             &
                                      Processor_Topology_File,           &
                                      Topology%MPI_Communicator          )

! ------------------------------------------------------------------
if (I_Am == First_Processor) then
  write (6    ,298) Topology%Grid
  write (6    ,299) Topology%Periodicity
  write (lid60,298) Topology%Grid
  write (lid60,299) Topology%Periodicity

 298 format (/1x,40('=')/3x,'Processor Cartesian Topology : ',3(i5)/)
 299 format (            3x,'Processor Periodicity        : ',3(i5)/1x,40('=')/)
endif
