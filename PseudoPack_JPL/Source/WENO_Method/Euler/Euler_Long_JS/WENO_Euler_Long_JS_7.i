#ifdef NO_ABS_BETA_JS
        a(0) = d(0)/(    b(0)  + Epsilon_WENO)**Power
        a(1) = d(1)/(    b(1)  + Epsilon_WENO)**Power
        a(2) = d(2)/(    b(2)  + Epsilon_WENO)**Power
        a(3) = d(3)/(    b(3)  + Epsilon_WENO)**Power
#else
        a(0) = d(0)/(ABS(b(0)) + Epsilon_WENO)**Power
        a(1) = d(1)/(ABS(b(1)) + Epsilon_WENO)**Power
        a(2) = d(2)/(ABS(b(2)) + Epsilon_WENO)**Power
        a(3) = d(3)/(ABS(b(3)) + Epsilon_WENO)**Power
#endif
