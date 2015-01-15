#include "pseudopack.h"

#if defined (MULTI_RESOLUTION) 
#define _RK_MR_FLAG(n)  , MR_Flag(N0:N5,M0:M5,K0:K5,n)                     
#else                          
#define _RK_MR_FLAG(n)          
#endif   

  Subroutine Central_Smoothing (Q)

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q

#if defined (MULTI_RESOLUTION)
  if (MR_Activation == 1) then               
#endif     
    call PS_Central_Smooth (1, N0, N5, N2, N3,                             &
                               M0, M5, M2, M3,                             &
                               K0, K5, K2, K3, NV, Q,                      &
                            Filter_Choice_x, Filter_Order_x,               &
                                             Filter_Viscosity_x            &
                            _RK_MR_FLAG(1) )


    call PS_Central_Smooth (2, N0, N5, N2, N3,                             &
                               M0, M5, M2, M3,                             &
                               K0, K5, K2, K3, NV, Q,                      &
                            Filter_Choice_y, Filter_Order_y,               &
                                             Filter_Viscosity_y            &
                            _RK_MR_FLAG(2) )

    call PS_Central_Smooth (3, N0, N5, N2, N3,                             &
                               M0, M5, M2, M3,                             &
                               K0, K5, K2, K3, NV, Q,                      &
                            Filter_Choice_z, Filter_Order_z,               &
                                             Filter_Viscosity_z            &
                            _RK_MR_FLAG(3) )

#if defined (MULTI_RESOLUTION)
  endif
#endif

  END Subroutine Central_Smoothing 

#if defined (MULTI_RESOLUTION)
#undef _RK_MR_FLAG
#endif
