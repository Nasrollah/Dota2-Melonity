#ifdef NO_ABS_BETA_ZICO
        a = d*(ONE + Tau/(    b +Epsilon_WENO))**Power
#else
#ifdef ZICO_ONE_PLUS_TAU_OVER_BETA_POWER
        a = d*(ONE + Tau/(ABS(b)+Epsilon_WENO))**Power
#else
        a = d*(ONE + (Tau/(ABS(b)+Epsilon_WENO))**Power)
#endif
#endif
