#include "WENO_Zico_Tau_1.i"

#ifdef NO_ABS_BETA_ZICO
        a(0) = d(0)*(ONE+ Tau/(    b(0)  + Epsilon_WENO))**Power
#else
#ifdef ZICO_ONE_PLUS_TAU_OVER_BETA_POWER
        a(0) = d(0)*(ONE+ Tau/(ABS(b(0)) + Epsilon_WENO))**Power
#else
        a(0) = d(0)*(ONE+(Tau/(ABS(b(0)) + Epsilon_WENO))**Power)
#endif
#endif
