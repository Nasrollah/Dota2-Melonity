#if defined (VARIABLE_GAMMA)
#define _EULER_GAMMA	Gamma_Gas, Gamma
#else
#define _EULER_GAMMA	Gamma
#endif

  Subroutine WENO_Flux (Q, D_Flux)

  REALTYPE            :: D_F_Epsilon

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, D_Flux
#if defined (NAVIER_STOKES)
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, P, T
#endif
#if defined (VARIABLE_GAMMA)
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: Gamma

  call PS_Specific_Heat_Gamma  (Gamma_Formulation, Q, _EULER_GAMMA)
#endif

#include "Adjust_D_Flux_Epsilon.i"

    call PS_WENO_Euler_Flux    (Index_x, Distributed_x, WENO_Order_x, dx, &
                                         N0, N5, N2, N3,                  &
                                Index_y, Distributed_y, WENO_Order_y, dy, &
                                         M0, M5, M2, M3,                  &
                                Index_z, Distributed_z, WENO_Order_z, dz, &
                                         K0, K5, K2, K3, NV,              &
                                _EULER_GAMMA,                             &
                                Q, D_Flux,                                &
                                  D_F_Epsilon=D_F_Epsilon,                &
                                MPI_Comm_Type=Topology%MPI_Communicator)

#if defined (NAVIER_STOKES)
  if (Navier_Stokes) then
    call PS_Primitive_Variable (N0, N5, M0, M5, K0, K5, NV,                    &
                                Gamma, Q, U, V, W, P)

    call PS_Temperature        (N0, N5, M0, M5, K0, K5, NV,                    &
                                       Q, P, T)

    call PS_CD_NS_Flux         (        N0, N5, N2, N3,                        &
                                        M0, M5, M2, M3,                        &
                                        K0, K5, K2, K3,                        &
                                        NV, NF, NS,                            &
                                Q, U, V, W, T, D_F, D_F_Epsilon,               &
                                MPI_Comm_Type=Topology%MPI_Communicator)

    D_Flux = D_Flux + D_F
  endif
#endif

  END Subroutine WENO_Flux

#if defined (_EULER_GAMMA)
#undef _EULER_GAMMA
#endif

