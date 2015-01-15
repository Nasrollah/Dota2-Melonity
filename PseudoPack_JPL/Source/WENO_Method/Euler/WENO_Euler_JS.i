#ifdef NO_ABS_BETA_JS
        a = d/(    b  + Epsilon_WENO)**Power
#else
        a = d/(ABS(b) + Epsilon_WENO)**Power
#endif
