  Allocate ( xi(N0_g:N5_g), x_1D(N0_g:N5_g))
  Allocate (eta(M0_g:M5_g), y_1D(M0_g:M5_g))
  Allocate ( fi(K0_g:K5_g), z_1D(K0_g:K5_g))

  Allocate (x(N0:N5,M0:M5,K0:K5), y(N0:N5,M0:M5,K0:K5), z(N0:N5,M0:M5,K0:K5))

  Allocate (Q(N0:N5,M0:M5,K0:K5,NV))

  Allocate (S1(N0:N5,M0:M5,K0:K5,NV))

  Allocate (Specific_Heat(N0:N5,M0:M5,K0:K5))

  Allocate (Normal_N0(N0:N2,M0:M5,K0:K5,ND), Normal_N5(N3:N5,M0:M5,K0:K5,ND))
  Allocate (Normal_M0(N0:N5,M0:M2,K0:K5,ND), Normal_M5(N0:N5,M3:M5,K0:K5,ND))
  Allocate (Normal_K0(N0:N5,M0:M5,K0:K2,ND), Normal_K5(N0:N5,M0:M5,K3:K5,ND))
  Allocate (    Q0_N0(N0:N2,M0:M5,K0:K5,NV),     Q0_N5(N3:N5,M0:M5,K0:K5,NV))
  Allocate (    Q0_M0(N0:N5,M0:M2,K0:K5,NV),     Q0_M5(N0:N5,M3:M5,K0:K5,NV))
  Allocate (    Q0_K0(N0:N5,M0:M5,K0:K2,NV),     Q0_K5(N0:N5,M0:M5,K3:K5,NV))

  x = HUGE(ONE) ; y = HUGE(ONE) ; z = HUGE(ONE)

#if defined (NAVIER_STOKES)
#include "Allocate_Navier_Stokes_Coefficients.i"
#endif

  Allocate (Q_0(NV,0:1))

  Q = HUGE(ONE) ; Specific_Heat = HUGE(ONE)

#if defined (MULTI_RESOLUTION)
  MR_N0 = N0 ; MR_N2 = N2 ; MR_N3 = N3 ; MR_N5 = N5
  MR_M0 = M0 ; MR_M2 = M2 ; MR_M3 = M3 ; MR_M5 = M5
  MR_K0 = K0 ; MR_K2 = K2 ; MR_K3 = K3 ; MR_K5 = K5
  
  call PS_MR_Adjust_Ghost_Index (MR_Order, MR_N2, MR_N3, MR_N0, MR_N5)
  call PS_MR_Adjust_Ghost_Index (MR_Order, MR_M2, MR_M3, MR_M0, MR_M5)
  call PS_MR_Adjust_Ghost_Index (MR_Order, MR_K2, MR_K3, MR_K0, MR_K5)
  
  call PS_MR_Allocate_Flag_And_Scale (Dimensions, .TRUE., .TRUE., .TRUE.)
  
  Allocate (Flag_WENO(N0:N5,M0:M5,K0:K5,3))
#endif
