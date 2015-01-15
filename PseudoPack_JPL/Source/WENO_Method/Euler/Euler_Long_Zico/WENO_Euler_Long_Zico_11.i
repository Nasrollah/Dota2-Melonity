#include "WENO_Zico_Tau_11.i"

#ifdef NO_ABS_BETA_ZICO
        a(0) = d(0)*(ONE+ Tau/(    b(0)  + Epsilon_WENO))**Power
        a(1) = d(1)*(ONE+ Tau/(    b(1)  + Epsilon_WENO))**Power
        a(2) = d(2)*(ONE+ Tau/(    b(2)  + Epsilon_WENO))**Power
        a(3) = d(3)*(ONE+ Tau/(    b(3)  + Epsilon_WENO))**Power
        a(4) = d(4)*(ONE+ Tau/(    b(4)  + Epsilon_WENO))**Power
        a(5) = d(5)*(ONE+ Tau/(    b(5)  + Epsilon_WENO))**Power
#else
#ifdef ZICO_ONE_PLUS_TAU_OVER_BETA_POWER
        a(0) = d(0)*(ONE+ Tau/(ABS(b(0)) + Epsilon_WENO))**Power
        a(1) = d(1)*(ONE+ Tau/(ABS(b(1)) + Epsilon_WENO))**Power
        a(2) = d(2)*(ONE+ Tau/(ABS(b(2)) + Epsilon_WENO))**Power
        a(3) = d(3)*(ONE+ Tau/(ABS(b(3)) + Epsilon_WENO))**Power
        a(4) = d(4)*(ONE+ Tau/(ABS(b(4)) + Epsilon_WENO))**Power
        a(5) = d(5)*(ONE+ Tau/(ABS(b(5)) + Epsilon_WENO))**Power
#else
        a(0) = d(0)*(ONE+(Tau/(ABS(b(0)) + Epsilon_WENO))**Power)
        a(1) = d(1)*(ONE+(Tau/(ABS(b(1)) + Epsilon_WENO))**Power)
        a(2) = d(2)*(ONE+(Tau/(ABS(b(2)) + Epsilon_WENO))**Power)
        a(3) = d(3)*(ONE+(Tau/(ABS(b(3)) + Epsilon_WENO))**Power)
        a(4) = d(4)*(ONE+(Tau/(ABS(b(4)) + Epsilon_WENO))**Power)
        a(5) = d(5)*(ONE+(Tau/(ABS(b(5)) + Epsilon_WENO))**Power)
#endif
#endif
