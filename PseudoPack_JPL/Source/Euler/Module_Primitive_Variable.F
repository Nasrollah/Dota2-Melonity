#include "pseudopack.h"


MODULE Primitive_Variable

implicit NONE

INTERFACE PS_Primitive_Variable
  MODULE PROCEDURE Primitive_Variable_11
  MODULE PROCEDURE Primitive_Variable_21
  MODULE PROCEDURE Primitive_Variable_31

  MODULE PROCEDURE Primitive_Variable_12
  MODULE PROCEDURE Primitive_Variable_22
  MODULE PROCEDURE Primitive_Variable_32
END INTERFACE

REALTYPE, parameter :: R_Universal = 8.314410d7        ! erg/K-mol
REALTYPE, parameter :: R_Specific  = 287.0d0           !   J/K-kg

PRIVATE
PUBLIC :: PS_Primitive_Variable

CONTAINS
!
!===================================================================
!
  Subroutine Primitive_Variable_11 (N0,N5, NV, Gamma,   &
                                    Q, U, P, h, C, S, T, EOS)

  integer                       :: N0,N5, NV 
  logical                       :: PhCST
  REALTYPE                      :: Gamma

  REALTYPE, dimension(N0:N5,NV)           :: Q
  REALTYPE, dimension(N0:N5)              :: U, PP
  REALTYPE, dimension(N0:N5)   , OPTIONAL :: P, h, C, S, T

  integer                      , OPTIONAL :: EOS
  integer                                 :: Equation_Of_state

  PhCST = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C)  &
                     .OR. PRESENT(S) .OR. PRESENT(T)

                  U = Q(:,2)/Q(:,1)
  if (PhCST)     PP = (Gamma-ONE)*(Q(:,3) - HALF*Q(:,1)*(U**2))

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
      if (PRESENT(P)) P = PP
      if (PRESENT(h)) h = (Q(:,3)+PP)/Q(:,1)
      if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,1)))
      if (PRESENT(S)) S = LOG(PP/Q(:,1)**Gamma)
      if (PRESENT(T)) T = PP/(R_Universal*Q(:,1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Primitive_Variable_11
!
!===================================================================
!
  Subroutine Primitive_Variable_21 (N0,N5, M0,M5, NV, Gamma,  &
                                    Q, U, V, P, h, C, S, T, EOS)

  integer                       :: N0,N5, M0,M5, NV 
  logical                       :: PhCST
  REALTYPE                      :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Q
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V, PP
  REALTYPE, dimension(N0:N5,M0:M5)   , OPTIONAL :: P, h, C, S, T

  integer                            , OPTIONAL :: EOS
  integer                                       :: Equation_Of_state

  PhCST = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C)  &
                     .OR. PRESENT(S) .OR. PRESENT(T)

                  U = Q(:,:,2)/Q(:,:,1)
                  V = Q(:,:,3)/Q(:,:,1)

  if (PhCST)     PP = (Gamma-ONE)*(Q(:,:,4) - HALF*Q(:,:,1)*(U*U+V*V))

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
      if (PRESENT(P)) P = PP
      if (PRESENT(h)) h = (Q(:,:,4)+PP)/Q(:,:,1)
      if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,:,1)))
      if (PRESENT(S)) S = LOG(PP/Q(:,:,1)**Gamma)
      if (PRESENT(T)) T = PP/(R_Universal*Q(:,:,1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Primitive_Variable_21
!
!===================================================================
!
  Subroutine Primitive_Variable_31 (N0,N5, M0,M5, K0,K5, NV,          &
                                    Gamma, Q, U, V, W, P, h, C, S, T, EOS)

  integer                       :: N0,N5, M0,M5, K0,K5, NV 
  logical                       :: PhCST
  REALTYPE                      :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, PP
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , OPTIONAL :: P, h, C, S, T

  integer                                  , OPTIONAL :: EOS
  integer                                             :: Equation_Of_state

  PhCST = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C)  &
                     .OR. PRESENT(S) .OR. PRESENT(T)

                  U = Q(:,:,:,2)/Q(:,:,:,1)
                  V = Q(:,:,:,3)/Q(:,:,:,1)
                  W = Q(:,:,:,4)/Q(:,:,:,1)
  if (PhCST)     PP = (Gamma-ONE)*(Q(:,:,:,5) - HALF*Q(:,:,:,1)*(U*U+V*V+W*W))

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
      if (PRESENT(P)) P = PP
      if (PRESENT(h)) h = (Q(:,:,:,5)+PP)/Q(:,:,:,1)
      if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,:,:,1)))
      if (PRESENT(S)) S = LOG(PP/Q(:,:,:,1)**Gamma)
      if (PRESENT(T)) T = PP/(R_Universal*Q(:,:,:,1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Primitive_Variable_31
!
!===================================================================
!
  Subroutine Primitive_Variable_12 (N0,N5, NV,                 &
                                    Gamma, Q, U, P, h, C, S, T, EOS)

  integer                       :: N0,N5, NV 
  logical                       :: PhCST

  REALTYPE, dimension(N0:N5,NV)           :: Q
  REALTYPE, dimension(N0:N5)              :: U, PP, Gamma
  REALTYPE, dimension(N0:N5)   , OPTIONAL :: P, h, C, S, T

  integer                      , OPTIONAL :: EOS
  integer                                 :: Equation_Of_state

  PhCST = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C)  &
                     .OR. PRESENT(S) .OR. PRESENT(T)

                  U = Q(:,2)/Q(:,1)
  if (PhCST)     PP = (Gamma-ONE)*(Q(:,3) - HALF*Q(:,1)*(U**2))

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
      if (PRESENT(P)) P = PP
      if (PRESENT(h)) h = (Q(:,3)+PP)/Q(:,1)
      if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,1)))
      if (PRESENT(S)) S = LOG(PP/Q(:,1)**Gamma)
      if (PRESENT(T)) T = PP/(R_Universal*Q(:,1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Primitive_Variable_12
!
!===================================================================
!
  Subroutine Primitive_Variable_22 (N0,N5, M0,M5, NV, &
                                    Gamma, Q, U, V, P, h, C, S, T, EOS)

  integer                       :: N0,N5, M0,M5, NV 
  logical                       :: PhCST

  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Q
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V, PP, Gamma
  REALTYPE, dimension(N0:N5,M0:M5)   , OPTIONAL :: P, h, C, S, T

  integer                            , OPTIONAL :: EOS
  integer                                       :: Equation_Of_state

  PhCST = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C)  &
                     .OR. PRESENT(S) .OR. PRESENT(T)

                  U = Q(:,:,2)/Q(:,:,1)
                  V = Q(:,:,3)/Q(:,:,1)

  if (PhCST)     PP = (Gamma-ONE)*(Q(:,:,4) - HALF*Q(:,:,1)*(U*U+V*V))

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
      if (PRESENT(P)) P = PP
      if (PRESENT(h)) h = (Q(:,:,4)+PP)/Q(:,:,1)
      if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,:,1)))
      if (PRESENT(S)) S = LOG(PP/Q(:,:,1)**Gamma)
      if (PRESENT(T)) T = PP/(R_Universal*Q(:,:,1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Primitive_Variable_22
!
!===================================================================
!
  Subroutine Primitive_Variable_32 (N0,N5, M0,M5, K0,K5, NV,          &
                                    Gamma, Q, U, V, W, P, h, C, S, T, EOS)

  integer                       :: N0,N5, M0,M5, K0,K5, NV 
  logical                       :: PhCST

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, PP, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , OPTIONAL :: P, h, C, S, T

  integer                                  , OPTIONAL :: EOS
  integer                                             :: Equation_Of_state

  PhCST = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C)  &
                     .OR. PRESENT(S) .OR. PRESENT(T)

                  U = Q(:,:,:,2)/Q(:,:,:,1)
                  V = Q(:,:,:,3)/Q(:,:,:,1)
                  W = Q(:,:,:,4)/Q(:,:,:,1)

  if (PhCST)     PP = (Gamma-ONE)*(Q(:,:,:,5) - HALF*Q(:,:,:,1)*(U*U+V*V+W*W))

                    Equation_Of_State = 0
  if (PRESENT(EOS)) Equation_Of_State = EOS
 
  SELECT CASE (Equation_Of_State)
    CASE DEFAULT
      if (PRESENT(P)) P = PP
      if (PRESENT(h)) h = (Q(:,:,:,5)+PP)/Q(:,:,:,1)
      if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,:,:,1)))
      if (PRESENT(S)) S = LOG(PP/Q(:,:,:,1)**Gamma)
      if (PRESENT(T)) T = PP/(R_Universal*Q(:,:,:,1))

    CASE (1)
      RETURN

  END SELECT

  END Subroutine Primitive_Variable_32

END MODULE Primitive_Variable
