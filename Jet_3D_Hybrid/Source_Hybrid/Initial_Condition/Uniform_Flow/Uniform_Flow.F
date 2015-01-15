  Subroutine Uniform_Flow(time,dt3)

  REALTYPE :: U_Gas_1, C_Gas_1, f
  real(kind=8) :: time, dt3 

  C_Gas_1 = SQRT(Gamma_Gas_1*P_Gas_1/Rho_Gas_1)
  U_Gas_1 = Mach*C_Gas_1

#if defined (ECHO_FLOW_CONDITION)
#include "Echo_Uniform_Flow_Condition.i"
#endif

  f = ONE

  Q(:,:,:,1) = Rho_Gas_1*f
  Q(:,:,:,2) = Rho_Gas_1*U_Gas_1
  Q(:,:,:,3) =      ZERO
  Q(:,:,:,4) =      ZERO
  Q(:,:,:,5) =   P_Gas_1/(Gamma_Gas_1-ONE) + HALF*Rho_Gas_1*U_Gas_1**2
  Q(:,:,:,6) = Rho_Gas_1*f

  call BC (Q,time,dt3)

  Specific_Heat = Gamma

  Q0_N0 = Q(N0:N2,:,:,:) ; Q0_N5 = Q(N3:N5,:,:,:)
  Q0_M0 = Q(:,M0:M2,:,:) ; Q0_M5 = Q(:,M3:M5,:,:)
  Q0_K0 = Q(:,:,K0:K2,:) ; Q0_K5 = Q(:,:,K3:K5,:)

  Q_0(:,0) = Q(N0,M0,K0,:) 
  Q_0(:,1) = Q(N5,M0,K0,:) 

  END Subroutine Uniform_Flow
