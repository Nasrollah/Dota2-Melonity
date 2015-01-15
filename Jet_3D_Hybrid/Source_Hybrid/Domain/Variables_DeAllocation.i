  if (ALLOCATED(Q  )) DEALLOCATE (Q )

  if (ALLOCATED(Specific_Heat)) DEALLOCATE (Specific_Heat)

  if (ALLOCATED(Q_0)) DEALLOCATE (Q_0)

  if (ALLOCATED(xi  )) DEALLOCATE (xi  )
  if (ALLOCATED(x   )) DEALLOCATE (x   )

  if (ALLOCATED(x_1D)) DEALLOCATE (x_1D)

  if (ALLOCATED(Q0_N0)) DEALLOCATE (Q0_N0)
  if (ALLOCATED(Q0_N5)) DEALLOCATE (Q0_N5)

  if (ALLOCATED(Normal_N0)) DEALLOCATE (Normal_N0)
  if (ALLOCATED(Normal_N5)) DEALLOCATE (Normal_N5)

  if (ALLOCATED(eta )) DEALLOCATE (eta )
  if (ALLOCATED(y   )) DEALLOCATE (y   )

  if (ALLOCATED(y_1D)) DEALLOCATE (y_1D)

  if (ALLOCATED(Q0_M0)) DEALLOCATE (Q0_M0)
  if (ALLOCATED(Q0_M5)) DEALLOCATE (Q0_M5)

  if (ALLOCATED(Normal_M0)) DEALLOCATE (Normal_M0)
  if (ALLOCATED(Normal_M5)) DEALLOCATE (Normal_M5)

  if (ALLOCATED(fi  )) DEALLOCATE (fi  )
  if (ALLOCATED(z   )) DEALLOCATE (z   )
 
  if (ALLOCATED(z_1D)) DEALLOCATE (z_1D)

  if (ALLOCATED(Q0_K0)) DEALLOCATE (Q0_K0)
  if (ALLOCATED(Q0_K5)) DEALLOCATE (Q0_K5)

  if (ALLOCATED(Normal_K0)) DEALLOCATE (Normal_K0)
  if (ALLOCATED(Normal_K5)) DEALLOCATE (Normal_K5)

#if defined (NAVIER_STOKES)
#include "DeAllocate_Navier_Stokes_Coefficients.i"
#endif

#if defined (MULTI_RESOLUTION)
if (ALLOCATED(Flag_WENO)) DEALLOCATE (Flag_WENO)
#endif

