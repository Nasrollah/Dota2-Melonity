#if defined (VARIABLE_GAMMA)
#define _EULER_GAMMA	Gamma_Formulation, Gamma_Gas 
#else
#define _EULER_GAMMA	Gamma
#endif

  Subroutine CD_Euler_Flux (Q, D_Flux)

  REALTYPE            :: D_F_Epsilon 

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, D_Flux

#include "Adjust_D_Flux_Epsilon.i"

  call PS_CD_Euler_Flux    (N0, N5, N2, N3,                        &
                            M0, M5, M2, M3,                        &
                            K0, K5, K2, K3, NV,                    &
                            _EULER_GAMMA ,                         &
                            Q,                                     &
                            CD_x, CD_y, CD_z, D_Flux,              &
                            -ONE, D_F_Epsilon)

  END Subroutine CD_Euler_Flux

#if defined (_EULER_GAMMA)
#undef _EULER_GAMMA
#endif
