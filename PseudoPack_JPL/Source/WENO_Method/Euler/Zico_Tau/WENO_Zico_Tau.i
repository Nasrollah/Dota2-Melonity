      SELECT CASE (Order)
        CASE (1)
#include "WENO_Zico_Tau_1.i"

        CASE (3)
#include "WENO_Zico_Tau_3.i"

        CASE (5)
#include "WENO_Zico_Tau_5.i"

        CASE (7)
#include "WENO_Zico_Tau_7.i"

        CASE (9)
#include "WENO_Zico_Tau_9.i"

        CASE (11)
#include "WENO_Zico_Tau_11.i"

      END SELECT
