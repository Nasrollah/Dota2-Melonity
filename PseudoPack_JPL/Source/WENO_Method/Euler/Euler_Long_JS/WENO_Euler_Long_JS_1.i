#ifdef NO_ABS_BETA_JS
        a(0) = d(0)/(    b(0)  + Epsilon_WENO)**Power
#else
        a(0) = d(0)/(ABS(b(0)) + Epsilon_WENO)**Power
#endif
