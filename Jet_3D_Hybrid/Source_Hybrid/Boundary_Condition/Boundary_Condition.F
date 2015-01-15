  Subroutine BC (Q,time,dt3)

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)          :: Q
  REALTYPE, dimension(3,NV)                          :: BC_IEO
  real(kind=8) :: time, dt3

  integer  :: M_s, M_s_0, n, Index
  integer  :: MPI_Communicator

  MPI_Communicator = Topology%MPI_Communicator



#if defined (PARALLEL_MPI)
#include "Ghost_Cell_MPI_Exchange.i"
#endif

  BC_IEO =  ONE 

  if (Euler .OR. Navier_Stokes) then
    BC_IEO(1,2) = -ONE ; BC_IEO(2,3) = -ONE ; BC_IEO(3,4) = -ONE
  endif

  call PS_WENO_Boundary_Condition (N0, N5, N2, N3,                             &
                                   M0, M5, M2, M3,                             &
                                   K0, K5, K2, K3,  NV, Q,                     &
                                   BC_On_Off, BC_Distributed, BC_Order,        &
                                   BC_Method, BC_Symmetry   , BC_Type, BC_IEO, &
                                   MPI_Communicator)

! Zero Velocity at the lower z domain 
  call PS_Boundary_Condition_z    (-1,                                         &
                                   N0, N5, N2, N3,                             &
                                   M0, M5, M2, M3,                             &
                                   K0, K5, K2, K3,  NV, Q,                     &
                                   MPI_Communicator)

! Injector at the lower z domain
  if ((Jet_Direction >= 2) .AND. (Jet_Direction <= 3)) then
    call Injector                 (N0, N5, N2, N3,                             &
                                   M0, M5, M2, M3,                             &
                                   K0, K5, K2, K3,  NV, x, y, z, Q,            &
                                   MPI_Communicator,time,dt3)
  endif

  END Subroutine BC
