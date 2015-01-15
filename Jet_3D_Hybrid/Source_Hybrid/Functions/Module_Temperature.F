#include "pseudopack.h"


MODULE Temperature

  USE Gas_Properties
  USE Gas_Mixture_Formulation

implicit NONE

INTERFACE PS_Temperature
  MODULE PROCEDURE Temperature_0D
  MODULE PROCEDURE Temperature_1D
  MODULE PROCEDURE Temperature_2D
  MODULE PROCEDURE Temperature_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Temperature

CONTAINS
!
!===================================================================
!
  Subroutine Temperature_0D (NV, Q, P, T, EOS)

  integer                       :: NV 

  REALTYPE, dimension(NV)                   :: Q
  REALTYPE                                  :: P, T
  REALTYPE                                  :: f_1, f_2, MW

  integer                        , OPTIONAL :: EOS
  integer                                   :: Equation_Of_state

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
      f_1 = Q(4)/Q(1) ; f_2 = ONE-f_1

      call PS_Mixture_Molecular_Weight (f_1, Gas(1)%MW, f_2, Gas(2)%MW, MW)

      T = P*MW/(R_Universal*Q(1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Temperature_0D
!
!===================================================================
!
  Subroutine Temperature_1D (N0,N5, NV, Q, P, T, EOS)

  integer                       :: N0,N5, NV 

  REALTYPE, dimension(N0:N5,NV)             :: Q
  REALTYPE, dimension(N0:N5)                :: P, T
  REALTYPE, dimension(N0:N5)                :: f_1, f_2, MW

  integer                        , OPTIONAL :: EOS
  integer                                   :: Equation_Of_state

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
      f_1 = Q(:,4)/Q(:,1) ; f_2 = ONE-f_1

      call PS_Mixture_Molecular_Weight (f_1, Gas(1)%MW, f_2, Gas(2)%MW, MW)

      T = P*MW/(R_Universal*Q(:,1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Temperature_1D
!
!===================================================================
!
  Subroutine Temperature_2D (N0,N5, M0,M5, NV, Q, P, T, EOS)

  integer                       :: N0,N5, M0,M5, NV 

  REALTYPE, dimension(N0:N5,M0:M5,NV)       :: Q
  REALTYPE, dimension(N0:N5,M0:M5)          :: P, T
  REALTYPE, dimension(N0:N5,M0:M5)          :: f_1, f_2, MW

  integer                        , OPTIONAL :: EOS
  integer                                   :: Equation_Of_state

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
      f_1 = Q(:,:,5)/Q(:,:,1) ; f_2 = ONE-f_1

      call PS_Mixture_Molecular_Weight (f_1, Gas(1)%MW, f_2, Gas(2)%MW, MW)

      T = P*MW/(R_Universal*Q(:,:,1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Temperature_2D
!
!===================================================================
!
  Subroutine Temperature_3D (N0,N5, M0,M5, K0,K5, NV, Q, P, T, EOS)

  integer                       :: N0,N5, M0,M5, K0,K5, NV 

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: P, T
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: f_1, f_2, MW

  integer                        , OPTIONAL :: EOS
  integer                                   :: Equation_Of_state

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
        f_1 = Q(:,:,:,6)/Q(:,:,:,1) ; f_2 = ONE-f_1

        call PS_Mixture_Molecular_Weight (f_1, Gas(1)%MW, f_2, Gas(2)%MW, MW)

        T = P*MW/(R_Universal*Q(:,:,:,1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Temperature_3D

END MODULE Temperature
