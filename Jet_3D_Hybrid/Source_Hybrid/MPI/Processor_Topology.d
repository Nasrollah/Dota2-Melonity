#if 0
logical  :: Distributed_x
logical  :: Distributed_y
logical  :: Distributed_z
#endif

logical  :: Periodicity_x
logical  :: Periodicity_y
logical  :: Periodicity_z

integer  :: Topology_Grid_x
integer  :: Topology_Grid_y
integer  :: Topology_Grid_z

logical  , dimension(Dimensions) :: Distributed, Periodicity
integer  , dimension(Dimensions) :: Grid

TYPE (Processor_Topology) :: Topology

integer  :: MPI_CART_Communicator
integer  :: K_Processor
