#include "pseudopack.h"


MODULE Euler_Function

  USE MPI_AllReduce_Operation

implicit NONE

INTERFACE PS_Euler_Primitive_Variable
  MODULE PROCEDURE Primitive_Variable_11
  MODULE PROCEDURE Primitive_Variable_21
  MODULE PROCEDURE Primitive_Variable_31

  MODULE PROCEDURE Primitive_Variable_12
  MODULE PROCEDURE Primitive_Variable_22
  MODULE PROCEDURE Primitive_Variable_32
END INTERFACE

INTERFACE PS_Euler_Flux
  MODULE PROCEDURE Flux_11
  MODULE PROCEDURE Flux_21
  MODULE PROCEDURE Flux_31
END INTERFACE

INTERFACE PS_Euler_Roe_Eigenvector
  MODULE PROCEDURE Roe_Eigenvector_10
  MODULE PROCEDURE Roe_Eigenvector_20
  MODULE PROCEDURE Roe_Eigenvector_30

  MODULE PROCEDURE Roe_Eigenvector_12
  MODULE PROCEDURE Roe_Eigenvector_22
  MODULE PROCEDURE Roe_Eigenvector_32
END INTERFACE

INTERFACE PS_Euler_Eigenvalue
  MODULE PROCEDURE Euler_Eigenvalue_11
  MODULE PROCEDURE Euler_Eigenvalue_21
  MODULE PROCEDURE Euler_Eigenvalue_31
END INTERFACE

INTERFACE PS_Maximum_Eigenvalue
  MODULE PROCEDURE Maximum_Eigenvalue_1D
  MODULE PROCEDURE Maximum_Eigenvalue_2D
  MODULE PROCEDURE Maximum_Eigenvalue_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Euler_Primitive_Variable
PUBLIC :: PS_Euler_Flux
PUBLIC :: PS_Euler_Roe_Eigenvector
PUBLIC :: PS_Euler_Eigenvalue
PUBLIC :: PS_Maximum_Eigenvalue

CONTAINS
!
!===================================================================
!
  Subroutine Roe_Eigenvector_10 (N0, N5, i0, i1, NV,      &
                                 Normal, Gamma,           &
                                 Rho_m, U_m, h_m, f1_m,   &
                                 E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: Gamma

  REALTYPE, dimension(1)                     :: Normal
  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, h_m, f1_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R

  SELECT CASE (NV)
    CASE (3)
      call Roe_Eigenvector_14 (N0, N5, i0, i1, NV,                   &
                               Normal, Gamma, Rho_m, U_m, h_m,       &
                               E_L, E_R)

    CASE (4)
      call Roe_Eigenvector_15 (N0, N5, i0, i1, NV,                   &
                               Normal, Gamma, Rho_m, U_m, h_m, f1_m, &
                               E_L, E_R)

  END SELECT

  END Subroutine Roe_Eigenvector_10 
!
!===================================================================
!
  Subroutine Roe_Eigenvector_20 (N0, N5, i0, i1, NV,          &
                                 Normal, Gamma,               &
                                 Rho_m, U_m, V_m, h_m, f1_m,  &
                                 E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: Gamma

  REALTYPE, dimension(2)                     :: Normal
  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, V_m, h_m, f1_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R

  SELECT CASE (NV)
    CASE (4)
      call Roe_Eigenvector_24 (N0, N5, i0, i1, NV,                        &
                               Normal, Gamma, Rho_m, U_m, V_m, h_m,       &
                               E_L, E_R)

    CASE (5)
      call Roe_Eigenvector_25 (N0, N5, i0, i1, NV,                        &
                               Normal, Gamma, Rho_m, U_m, V_m, h_m, f1_m, &
                               E_L, E_R)

  END SELECT

  END Subroutine Roe_Eigenvector_20 
!
!===================================================================
!
  Subroutine Roe_Eigenvector_30 (N0, N5, i0, i1, NV,              &
                                 Normal, Gamma,                   &
                                 Rho_m, U_m, V_m, W_m, h_m, f1_m, &
                                 E_L, E_R)

  integer  :: N0,N5, NV, i0,i1
  REALTYPE :: Gamma

  REALTYPE, dimension(3)           :: Normal
  REALTYPE, dimension(N0:N5)       :: Rho_m, U_m, V_m, W_m, h_m, f1_m
  REALTYPE, dimension(N0:N5,NV,NV) :: E_L, E_R

  SELECT CASE (NV)
    CASE (5)
      call Roe_Eigenvector_34 (N0, N5, i0, i1, NV,                             &
                               Normal, Gamma, Rho_m, U_m, V_m, W_m, h_m,       &
                               E_L, E_R)

    CASE (6)
      call Roe_Eigenvector_35 (N0, N5, i0, i1, NV,                             &
                               Normal, Gamma, Rho_m, U_m, V_m, W_m, h_m, f1_m, &
                               E_L, E_R)

  END SELECT

  END Subroutine Roe_Eigenvector_30 
!
!===================================================================
!
  Subroutine Primitive_Variable_11 (N0,N5, NV, Gamma,   &
                                    Q, U, P, h, C, S)

  integer                       :: N0,N5, NV 
  logical                       :: PhCS
  REALTYPE                      :: Gamma

  REALTYPE, dimension(N0:N5,NV)           :: Q
  REALTYPE, dimension(N0:N5)              :: U, PP
  REALTYPE, dimension(N0:N5)   , OPTIONAL :: P, h, C, S

  PhCS = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C) .OR. PRESENT(S)

                  U = Q(:,2)/Q(:,1)
  if (PhCS)      PP = (Gamma-ONE)*(Q(:,3) - HALF*Q(:,1)*(U**2))

  if (PRESENT(P)) P = PP
  if (PRESENT(h)) h = (Q(:,3)+PP)/Q(:,1)
  if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,1)))
  if (PRESENT(S)) S = LOG(PP/Q(:,1)**Gamma)

  END Subroutine Primitive_Variable_11
!
!===================================================================
!
  Subroutine Flux_11 (Index, N0,N5, NV, Q, U, P, F, S)

  integer                       :: Index, N0,N5, NV 

  REALTYPE, dimension(N0:N5,NV)           :: Q, F
  REALTYPE, dimension(N0:N5)              :: U, P
  REALTYPE, dimension(1)       , OPTIONAL :: S

  integer  :: n

  Index =1

        F(:,1) =   Q(:,1)      * U    
        F(:,2) =   Q(:,2)      * U + P
        F(:,3) =  (Q(:,3) + P) * U

      do n = 4,NV
        F(:,n) =   Q(:,n)      * U
      enddo

  if (.NOT. PRESENT(S)) RETURN
  if (S(Index) == ZERO) RETURN

  F = F - S(Index)*Q

  END Subroutine Flux_11
!
!===================================================================
!
  Subroutine Roe_Eigenvector_14 (N0, N5, i0, i1, NV,                &
                                 Normal, Gamma, Rho_m, U_m, h_m,    &
                                 E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: N_x, Gamma, gm1

  integer  :: i
  REALTYPE :: s0, s1, U, h, q, c, UdotN, M_x
  REALTYPE :: b1, b2, c1, u1

  REALTYPE, dimension(1)                     :: Normal
  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, h_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R

  N_x = Normal(1)

#if defined (DEBUG)
  if (NV >  3) then
    write (6 ,*) 'STOP!  Roe_Eigenvector_14, NV > 3'
    write (99,*) 'STOP!  Roe_Eigenvector_14, NV > 3'
    STOP
  endif
  
  if (N_x*N_x /= ONE) then
    write (6 ,*) 'STOP!  The normal N_x was not normalized to 1.'
    write (99,*) 'STOP!  The normal N_x was not normalized to 1.'
    STOP
  endif
#endif

  gm1 = Gamma-ONE

  do i = i0,i1
    s0 = Rho_m(i)/(Rho_m(i) + Rho_m(i+1)) ; s1 = ONE-s0

    U = s0*U_m(i) + s1*U_m(i+1)
    h = s0*h_m(i) + s1*h_m(i+1)

    q = HALF*(U*U)
    c = SQRT(ABS(gm1*(h - q)))

      UdotN = U*N_x

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1
      c1 = HALF/c        ; u1 = U*b1

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x

      E_R(i,2,1) = U    + c*N_x
      E_R(i,2,2) = U    - c*N_x
      E_R(i,2,3) = ONE

      E_R(i,3,1) =  h     + c*UdotN
      E_R(i,3,2) =  h     - c*UdotN
      E_R(i,3,3) = -q*M_x + U

      E_L(i,1,1) =     - c1*UdotN + b2
      E_L(i,1,2) =       c1*N_x   - u1
      E_L(i,1,3) =                  b1

      E_L(i,2,1) =       c1*UdotN + b2
      E_L(i,2,2) =     - c1*N_x   - u1
      E_L(i,2,3) =                  b1

      E_L(i,3,1) =  UdotN*N_x-TWO*U*b2
      E_L(i,3,2) =            TWO*U*u1
      E_L(i,3,3) =           -TWO*U*b1
  enddo

  END Subroutine Roe_Eigenvector_14 
!
!===================================================================
!
  Subroutine Euler_Eigenvalue_11 (N0,N5, NV, Normal, U, C, Lambda, S)

  integer  :: N0,N5, NV
  REALTYPE :: N_x

  REALTYPE, dimension(N0:N5)                 :: U, C, UdotN
  REALTYPE, dimension(N0:N5,NV)              :: Lambda
  REALTYPE, dimension(1)                     :: Normal
  REALTYPE, dimension(1)          , OPTIONAL :: S

  integer  :: n

  N_x = Normal(1)

#if defined (DEBUG)
  if (N_x*N_x /= ONE) then
    write (6 ,*) 'STOP!  The normal N_x was not normalized to 1.'
    write (99,*) 'STOP!  The normal N_x was not normalized to 1.'
    STOP
  endif
#endif

  if (PRESENT(S)) then
    UdotN = (U-S(1))*N_x 
  else
    UdotN = U*N_x 
  endif

    Lambda(:,1) = UdotN + C
    Lambda(:,2) = UdotN - C

  do n = 3,NV
    Lambda(:,n) = UdotN
  enddo

  END Subroutine Euler_Eigenvalue_11 
!
!===================================================================
!
  Subroutine Primitive_Variable_21 (N0,N5, M0,M5, NV, Gamma,  &
                                    Q, U, V, P, h, C, S)

  integer                       :: N0,N5, M0,M5, NV 
  logical                       :: PhCS
  REALTYPE                      :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Q
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V, PP
  REALTYPE, dimension(N0:N5,M0:M5)   , OPTIONAL :: P, h, C, S

  PhCS = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C) .OR. PRESENT(S)

                  U = Q(:,:,2)/Q(:,:,1)
                  V = Q(:,:,3)/Q(:,:,1)

  if (PhCS)      PP = (Gamma-ONE)*(Q(:,:,4) - HALF*Q(:,:,1)*(U*U+V*V))

  if (PRESENT(P)) P = PP
  if (PRESENT(h)) h = (Q(:,:,4)+PP)/Q(:,:,1)
  if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,:,1)))
  if (PRESENT(S)) S = LOG(PP/Q(:,:,1)**Gamma)

  END Subroutine Primitive_Variable_21
!
!===================================================================
!
  Subroutine Flux_21 (Index, N0,N5, M0,M5, NV, Q, U, V, P, F, S)

  integer                       :: Index, N0,N5, M0,M5, NV 

  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Q, F
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V, P
  REALTYPE, dimension(2)             , OPTIONAL :: S

  integer  :: n

  SELECT CASE (Index)
    CASE (1)
        F(:,:,1) =  Q(:,:,1)      * U    
        F(:,:,2) =  Q(:,:,2)      * U + P
        F(:,:,3) =  Q(:,:,3)      * U 
        F(:,:,4) = (Q(:,:,4) + P) * U

      do n = 5,NV
        F(:,:,n) =  Q(:,:,n)      * U 
      enddo

    CASE (2)
        F(:,:,1) =  Q(:,:,1)      * V    
        F(:,:,2) =  Q(:,:,2)      * V 
        F(:,:,3) =  Q(:,:,3)      * V + P
        F(:,:,4) = (Q(:,:,4) + P) * V

      do n = 5,NV
        F(:,:,n) =  Q(:,:,n)      * V 
      enddo

  END SELECT

  if (.NOT. PRESENT(S)) RETURN
  if (S(Index) == ZERO) RETURN

  F = F - S(Index)*Q

  END Subroutine Flux_21
!
!===================================================================
!
  Subroutine Roe_Eigenvector_24 (N0, N5, i0, i1, NV,                     &
                                 Normal, Gamma, Rho_m, U_m, V_m, h_m,    &
                                 E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: N_x, N_y, Gamma, gm1

  integer  :: i
  REALTYPE :: s0, s1, U, V, h, q, c, UdotN, M_x, M_y
  REALTYPE :: b1, b2, c1, u1, v1

  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, V_m, h_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R
  REALTYPE, dimension(2)                     :: Normal

  N_x = Normal(1) ; N_y = Normal(2)

#if defined (DEBUG)
  if (NV >  4) then
    write (6 ,*) 'STOP!  Roe_Eigenvector_24, NV > 4'
    write (99,*) 'STOP!  Roe_Eigenvector_24, NV > 4'
    STOP
  endif

  if (N_x*N_x + N_y*N_y /= ONE) then
    write (6 ,*) 'STOP!  The normal (N_x,N_y) was not normalized to 1.'
    write (99,*) 'STOP!  The normal (N_x,N_y) was not normalized to 1.'
    STOP
  endif
#endif

  gm1 = Gamma-ONE

  do i = i0,i1
    s0 = Rho_m(i)/(Rho_m(i) + Rho_m(i+1)) ; s1 = ONE-s0

    U = s0*U_m(i) + s1*U_m(i+1)
    V = s0*V_m(i) + s1*V_m(i+1)
    h = s0*h_m(i) + s1*h_m(i+1)

    q = HALF*(U*U + V*V)
    c = SQRT(ABS(gm1*(h - q)))

      UdotN = U*N_x + V*N_y

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = M_y

      E_R(i,2,1) =  U     + c*N_x
      E_R(i,2,2) =  U     - c*N_x
      E_R(i,2,3) = ONE
      E_R(i,2,4) = ZERO

      E_R(i,3,1) =  V     + c*N_y
      E_R(i,3,2) =  V     - c*N_y
      E_R(i,3,3) = ZERO
      E_R(i,3,4) = ONE

      E_R(i,4,1) =  h     + c*UdotN
      E_R(i,4,2) =  h     - c*UdotN
      E_R(i,4,3) = -q*M_x + U
      E_R(i,4,4) = -q*M_y + V

      E_L(i,1,1) =     - c1*UdotN + b2
      E_L(i,1,2) =       c1*N_x   - u1
      E_L(i,1,3) =       c1*N_y   - v1
      E_L(i,1,4) =                  b1

      E_L(i,2,1) =       c1*UdotN + b2
      E_L(i,2,2) =     - c1*N_x   - u1
      E_L(i,2,3) =     - c1*N_y   - v1
      E_L(i,2,4) =                  b1

      E_L(i,3,1) =  UdotN*N_x-TWO*U*b2
      E_L(i,3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(i,3,3) =            TWO*U*v1 - N_x*N_y
      E_L(i,3,4) =           -TWO*U*b1

      E_L(i,4,1) =  UdotN*N_y-TWO*V*b2
      E_L(i,4,2) =            TWO*V*u1 - N_y*N_x
      E_L(i,4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(i,4,4) =           -TWO*V*b1

  enddo

  END Subroutine Roe_Eigenvector_24 
!
!===================================================================
!
  Subroutine Euler_Eigenvalue_21 (N0,N5, M0,M5, NV, Normal, U, V, C,  &
                                         Lambda, S)

  integer  :: N0,N5, M0,M5, NV
  REALTYPE :: N_x, N_y

  REALTYPE, dimension(N0:N5,M0:M5)                 :: U, V, C, UdotN
  REALTYPE, dimension(N0:N5,M0:M5,NV)              :: Lambda
  REALTYPE, dimension(2)                           :: Normal
  REALTYPE, dimension(2)                , OPTIONAL :: S

  integer  :: n

  N_x = Normal(1); N_y = Normal(2)

#if defined (DEBUG)
  if (N_x*N_x + N_y*N_y /= ONE) then
    write (6 ,*) 'STOP!  The normal (N_x,N_y) was not normalized to 1.'
    write (99,*) 'STOP!  The normal (N_x,N_y) was not normalized to 1.'
    STOP
  endif
#endif

  if (PRESENT(S)) then
    UdotN = (U-S(1))*N_x + (V-S(2))*N_y 
  else
    UdotN = U*N_x + V*N_y 
  endif

    Lambda(:,:,1) = UdotN + C
    Lambda(:,:,2) = UdotN - C

  do n = 3,NV
    Lambda(:,:,n) = UdotN
  enddo

  END Subroutine Euler_Eigenvalue_21 
!
!===================================================================
!
  Subroutine Primitive_Variable_31 (N0,N5, M0,M5, K0,K5, NV,          &
                                    Gamma, Q, U, V, W, P, h, C, S)

  integer                       :: N0,N5, M0,M5, K0,K5, NV 
  logical                       :: PhCS
  REALTYPE                      :: Gamma

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, PP
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , OPTIONAL :: P, h, C, S

  PhCS = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C) .OR. PRESENT(S)

                  U = Q(:,:,:,2)/Q(:,:,:,1)
                  V = Q(:,:,:,3)/Q(:,:,:,1)
                  W = Q(:,:,:,4)/Q(:,:,:,1)
  if (PhCS)      PP = (Gamma-ONE)*(Q(:,:,:,5) - HALF*Q(:,:,:,1)*(U*U+V*V+W*W))

  if (PRESENT(P)) P = PP
  if (PRESENT(h)) h = (Q(:,:,:,5)+PP)/Q(:,:,:,1)
  if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,:,:,1)))
  if (PRESENT(S)) S = LOG(PP/Q(:,:,:,1)**Gamma)

  END Subroutine Primitive_Variable_31
!
!===================================================================
!
  Subroutine Flux_31 (Index, N0,N5, M0,M5, K0,K5, NV, Q, U, V, W, P, F, S)

  integer                       :: Index, N0,N5, M0,M5, K0,K5, NV 

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, P
  REALTYPE, dimension(3)                   , OPTIONAL :: S

  integer  :: n

  SELECT CASE (Index)
    CASE (1)
        F(:,:,:,1) =  Q(:,:,:,1)      * U    
        F(:,:,:,2) =  Q(:,:,:,2)      * U + P
        F(:,:,:,3) =  Q(:,:,:,3)      * U 
        F(:,:,:,4) =  Q(:,:,:,4)      * U 
        F(:,:,:,5) = (Q(:,:,:,5) + P) * U

      do n = 6,NV
        F(:,:,:,n) =  Q(:,:,:,n)      * U 
      enddo

    CASE (2)
        F(:,:,:,1) =  Q(:,:,:,1)      * V    
        F(:,:,:,2) =  Q(:,:,:,2)      * V 
        F(:,:,:,3) =  Q(:,:,:,3)      * V + P
        F(:,:,:,4) =  Q(:,:,:,4)      * V 
        F(:,:,:,5) = (Q(:,:,:,5) + P) * V

      do n = 6,NV
        F(:,:,:,n) =  Q(:,:,:,n)      * V 
      enddo

    CASE (3)
        F(:,:,:,1) =  Q(:,:,:,1)      * W
        F(:,:,:,2) =  Q(:,:,:,2)      * W
        F(:,:,:,3) =  Q(:,:,:,3)      * W 
        F(:,:,:,4) =  Q(:,:,:,4)      * W + P
        F(:,:,:,5) = (Q(:,:,:,5) + P) * W

      do n = 6,NV
        F(:,:,:,n) =  Q(:,:,:,n)      * W
      enddo

  END SELECT

  if (.NOT. PRESENT(S)) RETURN
  if (S(Index) == ZERO) RETURN

  F = F - S(Index)*Q

  END Subroutine Flux_31
!
!===================================================================
!
  Subroutine Roe_Eigenvector_34 (N0, N5, i0, i1, NV,                       &
                                 Normal, Gamma, Rho_m, U_m, V_m, W_m, h_m, &
                                 E_L, E_R)

  integer  :: N0,N5, NV, i0,i1
  REALTYPE :: N_x, N_y, N_z, Gamma, gm1

  integer  :: i
  REALTYPE :: s0, s1, U, V, W, h, q, c, UdotN, M_x, M_y, M_z
  REALTYPE :: b1, b2, c1, u1, v1, w1

  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, V_m, W_m, h_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R
  REALTYPE, dimension(3)                     :: Normal

  N_x = Normal(1); N_y = Normal(2) ; N_z = Normal(3)

#if defined (DEBUG)
  if (NV >  5) then
    write (6 ,*) 'STOP!  Roe_Eigenvector_34, NV > 5'
    write (99,*) 'STOP!  Roe_Eigenvector_34, NV > 5'
    STOP
  endif

  if (N_x*N_x + N_y*N_y + N_z*N_z/= ONE) then
    write (6 ,*) 'STOP!  The normal (N_x,N_y,N_z) was not normalized to 1.'
    write (99,*) 'STOP!  The normal (N_x,N_y,N_z) was not normalized to 1.'
    STOP
  endif
#endif

  gm1 = Gamma-ONE

  do i = i0,i1
    s0 = Rho_m(i)/(Rho_m(i) + Rho_m(i+1)) ; s1 = ONE-s0

    U = s0*U_m(i) + s1*U_m(i+1)
    V = s0*V_m(i) + s1*V_m(i+1)
    W = s0*W_m(i) + s1*W_m(i+1)
    h = s0*h_m(i) + s1*h_m(i+1)

    q = HALF*(U*U + V*V + W*W)
    c = SQRT(ABS(gm1*(h - q)))

      UdotN = U*N_x + V*N_y + W*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 ; w1 = W*b1

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = M_y
      E_R(i,1,5) = M_z

      E_R(i,2,1) =  U     + c*N_x
      E_R(i,2,2) =  U     - c*N_x
      E_R(i,2,3) = ONE
      E_R(i,2,4) = ZERO
      E_R(i,2,5) = ZERO

      E_R(i,3,1) =  V     + c*N_y
      E_R(i,3,2) =  V     - c*N_y
      E_R(i,3,3) = ZERO
      E_R(i,3,4) = ONE
      E_R(i,3,5) = ZERO

      E_R(i,4,1) =  W     + c*N_z
      E_R(i,4,2) =  W     - c*N_z
      E_R(i,4,3) = ZERO
      E_R(i,4,4) = ZERO
      E_R(i,4,5) = ONE

      E_R(i,5,1) =  h     + c*UdotN
      E_R(i,5,2) =  h     - c*UdotN
      E_R(i,5,3) = -q*M_x + U
      E_R(i,5,4) = -q*M_y + V
      E_R(i,5,5) = -q*M_z + W

      E_L(i,1,1) =     - c1*UdotN + b2
      E_L(i,1,2) =       c1*N_x   - u1
      E_L(i,1,3) =       c1*N_y   - v1
      E_L(i,1,4) =       c1*N_z   - w1
      E_L(i,1,5) =                  b1

      E_L(i,2,1) =       c1*UdotN + b2
      E_L(i,2,2) =     - c1*N_x   - u1
      E_L(i,2,3) =     - c1*N_y   - v1
      E_L(i,2,4) =     - c1*N_z   - w1
      E_L(i,2,5) =                  b1

      E_L(i,3,1) =  UdotN*N_x-TWO*U*b2
      E_L(i,3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(i,3,3) =            TWO*U*v1 - N_x*N_y
      E_L(i,3,4) =            TWO*U*w1 - N_x*N_z
      E_L(i,3,5) =           -TWO*U*b1

      E_L(i,4,1) =  UdotN*N_y-TWO*V*b2
      E_L(i,4,2) =            TWO*V*u1 - N_y*N_x
      E_L(i,4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(i,4,4) =            TWO*V*w1 - N_y*N_z
      E_L(i,4,5) =           -TWO*V*b1

      E_L(i,5,1) =  UdotN*N_z-TWO*W*b2
      E_L(i,5,2) =            TWO*W*u1 - N_z*N_x
      E_L(i,5,3) =            TWO*W*v1 - N_z*N_y
      E_L(i,5,4) =            TWO*W*w1 - N_z*N_z + ONE
      E_L(i,5,5) =           -TWO*W*b1

  enddo

  END Subroutine Roe_Eigenvector_34 
!
!===================================================================
!
  Subroutine Euler_Eigenvalue_31 (N0,N5, M0,M5, K0,K5, NV,      &
                                  Normal, U, V, W, C, Lambda, S)

  integer  :: N0,N5, M0,M5, K0,K5, NV
  REALTYPE :: N_x, N_y, N_z

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)                 :: U, V, W, C, UdotN
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)              :: Lambda
  REALTYPE, dimension(3)                                 :: Normal
  REALTYPE, dimension(3)                      , OPTIONAL :: S

  integer  :: n

  N_x = Normal(1); N_y = Normal(2) ; N_z = Normal(3)

#if defined (DEBUG)
  if (N_x*N_x + N_y*N_y + N_z*N_z/= ONE) then
    write (6 ,*) 'STOP!  The normal (N_x,N_y,N_z) was not normalized to 1.'
    write (99,*) 'STOP!  The normal (N_x,N_y,N_z) was not normalized to 1.'
    STOP
  endif
#endif

  if (PRESENT(S)) then
    UdotN = (U-S(1))*N_x + (V-S(2))*N_y + (W-S(3))*N_z
  else
    UdotN = U*N_x + V*N_y + W*N_z
  endif

    Lambda(:,:,:,1) = UdotN + C
    Lambda(:,:,:,2) = UdotN - C

  do n = 3,NV
    Lambda(:,:,:,n) = UdotN
  enddo

  END Subroutine Euler_Eigenvalue_31 
!
!===================================================================
!
  Subroutine Primitive_Variable_12 (N0,N5, NV,                 &
                                    Gamma, Q, U, P, h, C, S)

  integer                       :: N0,N5, NV 
  logical                       :: PhCS

  REALTYPE, dimension(N0:N5,NV)           :: Q
  REALTYPE, dimension(N0:N5)              :: U, PP, Gamma
  REALTYPE, dimension(N0:N5)   , OPTIONAL :: P, h, C, S

  PhCS = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C) .OR. PRESENT(S)

                  U = Q(:,2)/Q(:,1)
  if (PhCS)      PP = (Gamma-ONE)*(Q(:,3) - HALF*Q(:,1)*(U**2))

  if (PRESENT(P)) P = PP
  if (PRESENT(h)) h = (Q(:,3)+PP)/Q(:,1)
  if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,1)))
  if (PRESENT(S)) S = LOG(PP/Q(:,1)**Gamma)

  END Subroutine Primitive_Variable_12
!
!===================================================================
!
  Subroutine Primitive_Variable_22 (N0,N5, M0,M5, NV, &
                                    Gamma, Q, U, V, P, h, C, S)

  integer                       :: N0,N5, M0,M5, NV 
  logical                       :: PhCS

  REALTYPE, dimension(N0:N5,M0:M5,NV)           :: Q
  REALTYPE, dimension(N0:N5,M0:M5)              :: U, V, PP, Gamma
  REALTYPE, dimension(N0:N5,M0:M5)   , OPTIONAL :: P, h, C, S

  PhCS = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C) .OR. PRESENT(S)

                  U = Q(:,:,2)/Q(:,:,1)
                  V = Q(:,:,3)/Q(:,:,1)

  if (PhCS)      PP = (Gamma-ONE)*(Q(:,:,4) - HALF*Q(:,:,1)*(U*U+V*V))

  if (PRESENT(P)) P = PP
  if (PRESENT(h)) h = (Q(:,:,4)+PP)/Q(:,:,1)
  if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,:,1)))
  if (PRESENT(S)) S = LOG(PP/Q(:,:,1)**Gamma)

  END Subroutine Primitive_Variable_22
!
!===================================================================
!
  Subroutine Primitive_Variable_32 (N0,N5, M0,M5, K0,K5, NV,          &
                                    Gamma, Q, U, V, W, P, h, C, S)

  integer                       :: N0,N5, M0,M5, K0,K5, NV 
  logical                       :: PhCS

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)              :: U, V, W, PP, Gamma
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   , OPTIONAL :: P, h, C, S

  PhCS = PRESENT(P) .OR. PRESENT(h) .OR. PRESENT(C) .OR. PRESENT(S)

                  U = Q(:,:,:,2)/Q(:,:,:,1)
                  V = Q(:,:,:,3)/Q(:,:,:,1)
                  W = Q(:,:,:,4)/Q(:,:,:,1)

  if (PhCS)      PP = (Gamma-ONE)*(Q(:,:,:,5) - HALF*Q(:,:,:,1)*(U*U+V*V+W*W))

  if (PRESENT(P)) P = PP
  if (PRESENT(h)) h = (Q(:,:,:,5)+PP)/Q(:,:,:,1)
  if (PRESENT(C)) C = SQRT(ABS(Gamma*PP/Q(:,:,:,1)))
  if (PRESENT(S)) S = LOG(PP/Q(:,:,:,1)**Gamma)

  END Subroutine Primitive_Variable_32
!
!===================================================================
!
  Subroutine Roe_Eigenvector_12 (N0, N5, i0, i1, NV,                &
                                 Normal, Gamma_Gas,                 &
                                 Rho_m, U_m, h_m, f1_m, E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: N_x
  REALTYPE :: M_x

  integer  :: i
  REALTYPE :: s0, s1, U, h, q, e, c, f1, f2, X, Y
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1, UdotN
  REALTYPE :: Cp_1, Cv_1
  REALTYPE :: Cp_2, Cv_2
  REALTYPE :: Gamma, Gamma_Prime, gm1

  REALTYPE, dimension(:,:)                   :: Gamma_Gas
  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, h_m, f1_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R
  REALTYPE, dimension(1)                     :: Normal

  N_x = Normal(1)

#if defined (DEBUG)
  if (NV >  4) then
    write (6 ,*) 'STOP!  Roe_Eigenvector_12, NV > 4'
    write (99,*) 'STOP!  Roe_Eigenvector_12, NV > 4'
    STOP
  endif

  if (N_x*N_x /= ONE) then
    write (6 ,*) 'STOP!  The normal N_x was not normalized to 1.'
    write (99,*) 'STOP!  The normal N_x was not normalized to 1.'
    STOP
  endif
#endif

    Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
    Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

    do i = i0-1, i1
      s0 = Rho_m(i)/(Rho_m(i) + Rho_m(i+1)) ; s1 = ONE-s0

       U = s0* U_m(i) + s1* U_m(i+1)
       h = s0* h_m(i) + s1* h_m(i+1)
      f1 = s0*f1_m(i) + s1*f1_m(i+1)

      f1 = MIN(MAX(ZERO,f1),ONE) ; f2 = ONE - f1

      gm1         =  Cv_1*f1   + Cv_2*f2
      Gamma       = (Cp_1*f1   + Cp_2*f2  )/gm1
      Gamma_Prime = (Cp_1*Cv_2 - Cp_2*Cv_1)/gm1**2

      gm1         = Gamma-ONE

      q = HALF*(U*U)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = Gamma_Prime*e ; Y = X/gm1

      UdotN = U*N_x

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN 

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = ZERO

      E_R(i,2,1) =  U     + c*N_x
      E_R(i,2,2) =  U     - c*N_x
      E_R(i,2,3) = ONE
      E_R(i,2,4) = ZERO

      E_R(i,3,1) =  h     + c*UdotN
      E_R(i,3,2) =  h     - c*UdotN
      E_R(i,3,3) = -q*M_x + U         
      E_R(i,3,4) =        - Y

      E_R(i,4,1) = f1
      E_R(i,4,2) = f1
      E_R(i,4,3) = f1*M_x
      E_R(i,4,4) = ONE

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 

      E_L(i,1,1) =     - c1*UdotN  + b4
      E_L(i,1,2) =       c1*N_x    - u1
      E_L(i,1,3) =                   b1
      E_L(i,1,4) =                   b3

      E_L(i,2,1) =       c1*UdotN  + b4
      E_L(i,2,2) =     - c1*N_x    - u1
      E_L(i,2,3) =                   b1
      E_L(i,2,4) =                   b3

      E_L(i,3,1) =  UdotN*N_x-TWO*U*b4
      E_L(i,3,2) =            TWO*U*u1
      E_L(i,3,3) =           -TWO*U*b1
      E_L(i,3,4) =           -TWO*U*b3

      E_L(i,4,1) = -f1
      E_L(i,4,2) =  ZERO
      E_L(i,4,3) =  ZERO
      E_L(i,4,4) =  ONE

    enddo

  END Subroutine Roe_Eigenvector_12 
!
!===================================================================
!
  Subroutine Roe_Eigenvector_22 (N0, N5, i0, i1, NV,                &
                                 Normal, Gamma_Gas,                 &
                                 Rho_m, U_m, V_m, h_m, f1_m, E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: N_x, N_y
  REALTYPE :: M_x, M_y

  integer  :: i
  REALTYPE :: s0, s1, U, V, h, q, e, c, f1, f2, X, Y
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1, UdotN
  REALTYPE :: Cp_1, Cv_1
  REALTYPE :: Cp_2, Cv_2
  REALTYPE :: Gamma, Gamma_Prime, gm1

  REALTYPE, dimension(:,:)                   :: Gamma_Gas
  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, V_m, h_m, f1_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R
  REALTYPE, dimension(2)                     :: Normal

  N_x = Normal(1) ; N_y = Normal(2)

#if defined (DEBUG)
  if (NV >  5) then
    write (6 ,*) 'STOP!  Roe_Eigenvector_22, NV > 5'
    write (99,*) 'STOP!  Roe_Eigenvector_22, NV > 5'
    STOP
  endif

  if (N_x*N_x + N_y*N_y /= ONE) then
    write (6 ,*) 'STOP!  The normal (N_x,N_y) was not normalized to 1.'
    write (99,*) 'STOP!  The normal (N_x,N_y) was not normalized to 1.'
    STOP
  endif
#endif

    Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
    Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

    do i = i0-1, i1
      s0 = Rho_m(i)/(Rho_m(i) + Rho_m(i+1)) ; s1 = ONE-s0

       U = s0* U_m(i) + s1* U_m(i+1)
       V = s0* V_m(i) + s1* V_m(i+1)
       h = s0* h_m(i) + s1* h_m(i+1)
      f1 = s0*f1_m(i) + s1*f1_m(i+1)

      f1 = MIN(MAX(ZERO,f1),ONE) ; f2 = ONE - f1

      gm1         =  Cv_1*f1   + Cv_2*f2
      Gamma       = (Cp_1*f1   + Cp_2*f2  )/gm1
      Gamma_Prime = (Cp_1*Cv_2 - Cp_2*Cv_1)/gm1**2

      gm1         = Gamma-ONE

      q = HALF*(U*U + V*V)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = Gamma_Prime*e ; Y = X/gm1

      UdotN = U*N_x + V*N_y

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = M_y
      E_R(i,1,5) = ZERO

      E_R(i,2,1) =  U     + c*N_x
      E_R(i,2,2) =  U     - c*N_x
      E_R(i,2,3) = ONE                  ! UdotN
      E_R(i,2,4) = ZERO
      E_R(i,2,5) = ZERO

      E_R(i,3,1) =  V     + c*N_y
      E_R(i,3,2) =  V     - c*N_y
      E_R(i,3,3) = ZERO
      E_R(i,3,4) = ONE                  ! UdotN
      E_R(i,3,5) = ZERO

      E_R(i,4,1) =  h     + c*UdotN
      E_R(i,4,2) =  h     - c*UdotN
      E_R(i,4,3) = -q*M_x + U         ! -q*N_x + U*UdotN
      E_R(i,4,4) = -q*M_y + V         ! -q*N_y + V*UdotN
      E_R(i,4,5) =         - Y

      E_R(i,5,1) = f1
      E_R(i,5,2) = f1
      E_R(i,5,3) = f1*M_x              ! f1*N_x
      E_R(i,5,4) = f1*M_y              ! f1*N_y
      E_R(i,5,5) = ONE

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1

      E_L(i,1,1) =     - c1*UdotN  + b4
      E_L(i,1,2) =       c1*N_x    - u1
      E_L(i,1,3) =       c1*N_y    - v1
      E_L(i,1,4) =                   b1
      E_L(i,1,5) =                   b3

      E_L(i,2,1) =       c1*UdotN  + b4
      E_L(i,2,2) =     - c1*N_x    - u1
      E_L(i,2,3) =     - c1*N_y    - v1
      E_L(i,2,4) =                   b1
      E_L(i,2,5) =                   b3

      E_L(i,3,1) =  UdotN*N_x-TWO*U*b4
      E_L(i,3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(i,3,3) =            TWO*U*v1 - N_x*N_y
      E_L(i,3,4) =           -TWO*U*b1
      E_L(i,3,5) =           -TWO*U*b3

      E_L(i,4,1) =  UdotN*N_y-TWO*V*b4
      E_L(i,4,2) =            TWO*V*u1 - N_y*N_x
      E_L(i,4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(i,4,4) =           -TWO*V*b1
      E_L(i,4,5) =           -TWO*V*b3

      E_L(i,5,1) = -f1
      E_L(i,5,2) =  ZERO
      E_L(i,5,3) =  ZERO
      E_L(i,5,4) =  ZERO
      E_L(i,5,5) =  ONE

    enddo

  END Subroutine Roe_Eigenvector_22 
!
!===================================================================
!
  Subroutine Roe_Eigenvector_32 (N0, N5, i0, i1, NV,                &
                                 Normal, Gamma_Gas,                 &
                                 Rho_m, U_m, V_m, W_m, h_m, f1_m, E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: N_x, N_y, N_z
  REALTYPE :: M_x, M_y, M_z

  integer  :: i
  REALTYPE :: s0, s1, U, V, W, h, q, e, c, f1, f2, X, Y
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1, w1, UdotN
  REALTYPE :: Cp_1, Cv_1
  REALTYPE :: Cp_2, Cv_2
  REALTYPE :: Gamma, Gamma_Prime, gm1

  REALTYPE, dimension(:,:)                   :: Gamma_Gas
  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, V_m, W_m, h_m, f1_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R
  REALTYPE, dimension(3)                     :: Normal

  N_x = Normal(1) ; N_y = Normal(2) ; N_z = Normal(3)

#if defined (DEBUG)
  if (NV >  6) then
    write (6 ,*) 'STOP!  Roe_Eigenvector_32, NV > 6'
    write (99,*) 'STOP!  Roe_Eigenvector_32, NV > 6'
    STOP
  endif

  if (N_x*N_x + N_y*N_y + N_z*N_z /= ONE) then
    write (6 ,*) 'STOP!  The normal (N_x,N_y,N_z) was not normalized to 1.'
    write (99,*) 'STOP!  The normal (N_x,N_y,N_z) was not normalized to 1.'
    STOP
  endif
#endif

    Cp_1 = Gamma_Gas(1,2) ; Cv_1 = Gamma_Gas(1,3)
    Cp_2 = Gamma_Gas(2,2) ; Cv_2 = Gamma_Gas(2,3)

    do i = i0-1, i1
      s0 = Rho_m(i)/(Rho_m(i) + Rho_m(i+1)) ; s1 = ONE-s0

       U = s0* U_m(i) + s1* U_m(i+1)
       V = s0* V_m(i) + s1* V_m(i+1)
       W = s0* W_m(i) + s1* W_m(i+1)
       h = s0* h_m(i) + s1* h_m(i+1)
      f1 = s0*f1_m(i) + s1*f1_m(i+1)

      f1 = MIN(MAX(ZERO,f1),ONE) ; f2 = ONE - f1

      gm1         =  Cv_1*f1   + Cv_2*f2
      Gamma       = (Cp_1*f1   + Cp_2*f2  )/gm1
      Gamma_Prime = (Cp_1*Cv_2 - Cp_2*Cv_1)/gm1**2

      gm1         = Gamma-ONE

      q = HALF*(U*U + V*V + W*W)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = Gamma_Prime*e ; Y = X/gm1

      UdotN = U*N_x + V*N_y + W*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN 

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = M_y
      E_R(i,1,5) = M_z
      E_R(i,1,6) = ZERO

      E_R(i,2,1) =  U     + c*N_x
      E_R(i,2,2) =  U     - c*N_x
      E_R(i,2,3) = ONE
      E_R(i,2,4) = ZERO
      E_R(i,2,5) = ZERO
      E_R(i,2,6) = ZERO

      E_R(i,3,1) =  V     + c*N_y
      E_R(i,3,2) =  V     - c*N_y
      E_R(i,3,3) = ZERO
      E_R(i,3,4) = ONE
      E_R(i,3,5) = ZERO
      E_R(i,3,6) = ZERO

      E_R(i,4,1) =  W     + c*N_z
      E_R(i,4,2) =  W     - c*N_z
      E_R(i,4,3) = ZERO
      E_R(i,4,4) = ZERO
      E_R(i,4,5) = ONE
      E_R(i,4,6) = ZERO

      E_R(i,5,1) =  h     + c*UdotN
      E_R(i,5,2) =  h     - c*UdotN
      E_R(i,5,3) = -q*M_x + U
      E_R(i,5,4) = -q*M_y + V
      E_R(i,5,5) = -q*M_z + W
      E_R(i,5,6) =        - Y

      E_R(i,6,1) = f1
      E_R(i,6,2) = f1
      E_R(i,6,3) = f1*M_x
      E_R(i,6,4) = f1*M_y
      E_R(i,6,5) = f1*M_z
      E_R(i,6,6) = ONE

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 ; w1 = W*b1

      E_L(i,1,1) =     - c1*UdotN  + b4
      E_L(i,1,2) =       c1*N_x    - u1
      E_L(i,1,3) =       c1*N_y    - v1
      E_L(i,1,4) =       c1*N_z    - w1
      E_L(i,1,5) =                   b1
      E_L(i,1,6) =                   b3

      E_L(i,2,1) =       c1*UdotN  + b4
      E_L(i,2,2) =     - c1*N_x    - u1
      E_L(i,2,3) =     - c1*N_y    - v1
      E_L(i,2,4) =     - c1*N_z    - w1
      E_L(i,2,5) =                   b1
      E_L(i,2,6) =                   b3

      E_L(i,3,1) =  UdotN*N_x-TWO*U*b4
      E_L(i,3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(i,3,3) =            TWO*U*v1 - N_x*N_y
      E_L(i,3,4) =            TWO*U*w1 - N_x*N_z
      E_L(i,3,5) =           -TWO*U*b1
      E_L(i,3,6) =           -TWO*U*b3

      E_L(i,4,1) =  UdotN*N_y-TWO*V*b4
      E_L(i,4,2) =            TWO*V*u1 - N_y*N_x
      E_L(i,4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(i,4,4) =            TWO*V*w1 - N_y*N_z
      E_L(i,4,5) =           -TWO*V*b1
      E_L(i,4,6) =           -TWO*V*b3

      E_L(i,5,1) =  UdotN*N_z-TWO*W*b4
      E_L(i,5,2) =            TWO*W*u1 - N_z*N_x
      E_L(i,5,3) =            TWO*W*v1 - N_z*N_y
      E_L(i,5,4) =            TWO*W*w1 - N_z*N_z + ONE
      E_L(i,5,5) =           -TWO*W*b1
      E_L(i,5,6) =           -TWO*W*b3

      E_L(i,6,1) = -f1
      E_L(i,6,2) =  ZERO
      E_L(i,6,3) =  ZERO
      E_L(i,6,4) =  ZERO
      E_L(i,6,5) =  ZERO
      E_L(i,6,6) =  ONE

    enddo

  END Subroutine Roe_Eigenvector_32 
#if 1
!
!===================================================================
!
  Subroutine Roe_Eigenvector_15 (N0, N5, i0, i1, NV,                    &
                                 Normal, Gamma,                         &
                                 Rho_m, U_m, h_m, f1_m, E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: N_x
  REALTYPE :: M_x

  integer  :: i
  REALTYPE :: s0, s1, U, h, q, e, c, f1, f2, X, Y
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1, UdotN
  REALTYPE :: Gamma, Gamma_Prime, gm1

  REALTYPE, dimension(1)                     :: Normal
  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, h_m, f1_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R

  N_x = Normal(1)

#if defined (DEBUG)
  if (NV >  4) then
    write (6 ,*) 'STOP!  Roe_Eigenvector_15, NV > 4'
    write (99,*) 'STOP!  Roe_Eigenvector_15, NV > 4'
    STOP
  endif

  if (N_x*N_x /= ONE) then
    write (6 ,*) 'STOP!  The normal N_x was not normalized to 1.'
    write (99,*) 'STOP!  The normal N_x was not normalized to 1.'
    STOP
  endif
#endif

    do i = i0-1, i1
      s0 = Rho_m(i)/(Rho_m(i) + Rho_m(i+1)) ; s1 = ONE-s0

       U = s0* U_m(i) + s1* U_m(i+1)
       h = s0* h_m(i) + s1* h_m(i+1)
      f1 = s0*f1_m(i) + s1*f1_m(i+1)

      f1 = MIN(MAX(ZERO,f1),ONE) ; f2 = ONE - f1

      Gamma_Prime = ZERO
      gm1         = Gamma-ONE

      q = HALF*(U*U)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = Gamma_Prime*e ; Y = X/gm1

      UdotN = U*N_x

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN 

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = ZERO

      E_R(i,2,1) =  U     + c*N_x
      E_R(i,2,2) =  U     - c*N_x
      E_R(i,2,3) = ONE
      E_R(i,2,4) = ZERO

      E_R(i,3,1) =  h     + c*UdotN
      E_R(i,3,2) =  h     - c*UdotN
      E_R(i,3,3) = -q*M_x + U         
      E_R(i,3,4) =        - Y

      E_R(i,4,1) = f1
      E_R(i,4,2) = f1
      E_R(i,4,3) = f1*M_x
      E_R(i,4,4) = ONE

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 

      E_L(i,1,1) =     - c1*UdotN  + b4
      E_L(i,1,2) =       c1*N_x    - u1
      E_L(i,1,3) =                   b1
      E_L(i,1,4) =                   b3

      E_L(i,2,1) =       c1*UdotN  + b4
      E_L(i,2,2) =     - c1*N_x    - u1
      E_L(i,2,3) =                   b1
      E_L(i,2,4) =                   b3

      E_L(i,3,1) =  UdotN*N_x-TWO*U*b4
      E_L(i,3,2) =            TWO*U*u1
      E_L(i,3,3) =           -TWO*U*b1
      E_L(i,3,4) =           -TWO*U*b3

      E_L(i,4,1) = -f1
      E_L(i,4,2) =  ZERO
      E_L(i,4,3) =  ZERO
      E_L(i,4,4) =  ONE

    enddo

  END Subroutine Roe_Eigenvector_15 
!
!===================================================================
!
  Subroutine Roe_Eigenvector_25 (N0, N5, i0, i1, NV,                    &
                                 Normal, Gamma,                         &
                                 Rho_m, U_m, V_m, h_m, f1_m, E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: N_x, N_y
  REALTYPE :: M_x, M_y

  integer  :: i
  REALTYPE :: s0, s1, U, V, h, q, e, c, f1, f2, X, Y
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1, UdotN
  REALTYPE :: Gamma, Gamma_Prime, gm1

  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, V_m, h_m, f1_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R
  REALTYPE, dimension(2)                     :: Normal

  N_x = Normal(1) ; N_y = Normal(2)

#if defined (DEBUG)
  if (NV >  5) then
    write (6 ,*) 'STOP!  Roe_Eigenvector_25, NV > 5'
    write (99,*) 'STOP!  Roe_Eigenvector_25, NV > 5'
    STOP
  endif

  if (N_x*N_x + N_y*N_y /= ONE) then
    write (6 ,*) 'STOP!  The normal (N_x,N_y) was not normalized to 1.'
    write (99,*) 'STOP!  The normal (N_x,N_y) was not normalized to 1.'
    STOP
  endif
#endif

    do i = i0-1, i1
      s0 = Rho_m(i)/(Rho_m(i) + Rho_m(i+1)) ; s1 = ONE-s0

       U = s0* U_m(i) + s1* U_m(i+1)
       V = s0* V_m(i) + s1* V_m(i+1)
       h = s0* h_m(i) + s1* h_m(i+1)
      f1 = s0*f1_m(i) + s1*f1_m(i+1)

      f1 = MIN(MAX(ZERO,f1),ONE) ; f2 = ONE - f1

      Gamma_Prime = ZERO

      gm1         = Gamma-ONE

      q = HALF*(U*U + V*V)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = Gamma_Prime*e ; Y = X/gm1

      UdotN = U*N_x + V*N_y

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = M_y
      E_R(i,1,5) = ZERO

      E_R(i,2,1) =  U     + c*N_x
      E_R(i,2,2) =  U     - c*N_x
      E_R(i,2,3) = ONE                  ! UdotN
      E_R(i,2,4) = ZERO
      E_R(i,2,5) = ZERO

      E_R(i,3,1) =  V     + c*N_y
      E_R(i,3,2) =  V     - c*N_y
      E_R(i,3,3) = ZERO
      E_R(i,3,4) = ONE                  ! UdotN
      E_R(i,3,5) = ZERO

      E_R(i,4,1) =  h     + c*UdotN
      E_R(i,4,2) =  h     - c*UdotN
      E_R(i,4,3) = -q*M_x + U         ! -q*N_x + U*UdotN
      E_R(i,4,4) = -q*M_y + V         ! -q*N_y + V*UdotN
      E_R(i,4,5) =         - Y

      E_R(i,5,1) = f1
      E_R(i,5,2) = f1
      E_R(i,5,3) = f1*M_x              ! f1*N_x
      E_R(i,5,4) = f1*M_y              ! f1*N_y
      E_R(i,5,5) = ONE

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1

      E_L(i,1,1) =     - c1*UdotN  + b4
      E_L(i,1,2) =       c1*N_x    - u1
      E_L(i,1,3) =       c1*N_y    - v1
      E_L(i,1,4) =                   b1
      E_L(i,1,5) =                   b3

      E_L(i,2,1) =       c1*UdotN  + b4
      E_L(i,2,2) =     - c1*N_x    - u1
      E_L(i,2,3) =     - c1*N_y    - v1
      E_L(i,2,4) =                   b1
      E_L(i,2,5) =                   b3

      E_L(i,3,1) =  UdotN*N_x-TWO*U*b4
      E_L(i,3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(i,3,3) =            TWO*U*v1 - N_x*N_y
      E_L(i,3,4) =           -TWO*U*b1
      E_L(i,3,5) =           -TWO*U*b3

      E_L(i,4,1) =  UdotN*N_y-TWO*V*b4
      E_L(i,4,2) =            TWO*V*u1 - N_y*N_x
      E_L(i,4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(i,4,4) =           -TWO*V*b1
      E_L(i,4,5) =           -TWO*V*b3

      E_L(i,5,1) = -f1
      E_L(i,5,2) =  ZERO
      E_L(i,5,3) =  ZERO
      E_L(i,5,4) =  ZERO
      E_L(i,5,5) =  ONE

    enddo

  END Subroutine Roe_Eigenvector_25 
!
!===================================================================
!
  Subroutine Roe_Eigenvector_35 (N0, N5, i0, i1, NV,                        &
                                 Normal, Gamma,                             &
                                 Rho_m, U_m, V_m, W_m, h_m, f1_m, E_L, E_R)

  integer  :: N0, N5, i0, i1, NV
  REALTYPE :: N_x, N_y, N_z
  REALTYPE :: M_x, M_y, M_z

  integer  :: i
  REALTYPE :: s0, s1, U, V, W, h, q, e, c, f1, f2, X, Y
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1, w1, UdotN
  REALTYPE :: Gamma, Gamma_Prime, gm1

  REALTYPE, dimension(N0:N5)                 :: Rho_m, U_m, V_m, W_m, h_m, f1_m
  REALTYPE, dimension(N0:N5,NV,NV)           :: E_L, E_R
  REALTYPE, dimension(3)                     :: Normal

  N_x = Normal(1) ; N_y = Normal(2) ; N_z = Normal(3)

#if defined (DEBUG)
  if (NV >  6) then
    write (6 ,*) 'STOP!  Roe_Eigenvector_35, NV > 6'
    write (99,*) 'STOP!  Roe_Eigenvector_35, NV > 6'
    STOP
  endif

  if (N_x*N_x + N_y*N_y + N_z*N_z /= ONE) then
    write (6 ,*) 'STOP!  The normal (N_x,N_y,N_z) was not normalized to 1.'
    write (99,*) 'STOP!  The normal (N_x,N_y,N_z) was not normalized to 1.'
    STOP
  endif
#endif

    do i = i0-1, i1
      s0 = Rho_m(i)/(Rho_m(i) + Rho_m(i+1)) ; s1 = ONE-s0

       U = s0* U_m(i) + s1* U_m(i+1)
       V = s0* V_m(i) + s1* V_m(i+1)
       W = s0* W_m(i) + s1* W_m(i+1)
       h = s0* h_m(i) + s1* h_m(i+1)
      f1 = s0*f1_m(i) + s1*f1_m(i+1)

      f1 = MIN(MAX(ZERO,f1),ONE) ; f2 = ONE - f1

      Gamma_Prime = ZERO
      gm1         = Gamma-ONE

      q = HALF*(U*U + V*V + W*W)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = Gamma_Prime*e ; Y = X/gm1

      UdotN = U*N_x + V*N_y + W*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)

      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN 

      E_R(i,1,1) = ONE
      E_R(i,1,2) = ONE
      E_R(i,1,3) = M_x
      E_R(i,1,4) = M_y
      E_R(i,1,5) = M_z
      E_R(i,1,6) = ZERO

      E_R(i,2,1) =  U     + c*N_x
      E_R(i,2,2) =  U     - c*N_x
      E_R(i,2,3) = ONE
      E_R(i,2,4) = ZERO
      E_R(i,2,5) = ZERO
      E_R(i,2,6) = ZERO

      E_R(i,3,1) =  V     + c*N_y
      E_R(i,3,2) =  V     - c*N_y
      E_R(i,3,3) = ZERO
      E_R(i,3,4) = ONE
      E_R(i,3,5) = ZERO
      E_R(i,3,6) = ZERO

      E_R(i,4,1) =  W     + c*N_z
      E_R(i,4,2) =  W     - c*N_z
      E_R(i,4,3) = ZERO
      E_R(i,4,4) = ZERO
      E_R(i,4,5) = ONE
      E_R(i,4,6) = ZERO

      E_R(i,5,1) =  h     + c*UdotN
      E_R(i,5,2) =  h     - c*UdotN
      E_R(i,5,3) = -q*M_x + U
      E_R(i,5,4) = -q*M_y + V
      E_R(i,5,5) = -q*M_z + W
      E_R(i,5,6) =        - Y

      E_R(i,6,1) = f1
      E_R(i,6,2) = f1
      E_R(i,6,3) = f1*M_x
      E_R(i,6,4) = f1*M_y
      E_R(i,6,5) = f1*M_z
      E_R(i,6,6) = ONE

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 ; w1 = W*b1

      E_L(i,1,1) =     - c1*UdotN  + b4
      E_L(i,1,2) =       c1*N_x    - u1
      E_L(i,1,3) =       c1*N_y    - v1
      E_L(i,1,4) =       c1*N_z    - w1
      E_L(i,1,5) =                   b1
      E_L(i,1,6) =                   b3

      E_L(i,2,1) =       c1*UdotN  + b4
      E_L(i,2,2) =     - c1*N_x    - u1
      E_L(i,2,3) =     - c1*N_y    - v1
      E_L(i,2,4) =     - c1*N_z    - w1
      E_L(i,2,5) =                   b1
      E_L(i,2,6) =                   b3

      E_L(i,3,1) =  UdotN*N_x-TWO*U*b4
      E_L(i,3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(i,3,3) =            TWO*U*v1 - N_x*N_y
      E_L(i,3,4) =            TWO*U*w1 - N_x*N_z
      E_L(i,3,5) =           -TWO*U*b1
      E_L(i,3,6) =           -TWO*U*b3

      E_L(i,4,1) =  UdotN*N_y-TWO*V*b4
      E_L(i,4,2) =            TWO*V*u1 - N_y*N_x
      E_L(i,4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(i,4,4) =            TWO*V*w1 - N_y*N_z
      E_L(i,4,5) =           -TWO*V*b1
      E_L(i,4,6) =           -TWO*V*b3

      E_L(i,5,1) =  UdotN*N_z-TWO*W*b4
      E_L(i,5,2) =            TWO*W*u1 - N_z*N_x
      E_L(i,5,3) =            TWO*W*v1 - N_z*N_y
      E_L(i,5,4) =            TWO*W*w1 - N_z*N_z + ONE
      E_L(i,5,5) =           -TWO*W*b1
      E_L(i,5,6) =           -TWO*W*b3

      E_L(i,6,1) = -f1
      E_L(i,6,2) =  ZERO
      E_L(i,6,3) =  ZERO
      E_L(i,6,4) =  ZERO
      E_L(i,6,5) =  ZERO
      E_L(i,6,6) =  ONE

    enddo

  END Subroutine Roe_Eigenvector_35 
#endif
!
!==========================================================================
!
  Subroutine Maximum_Eigenvalue_1D (Index, N0, N5, NV,                   &
                                           U, C,                         &
                                    Lambda_Max, V_Moving_Frame,          &
                                    MPI_Comm_Type)

  integer  :: Index
  integer  :: N0, N5

  integer  :: NV

  REALTYPE, dimension(N0:N5)          :: U, C
  REALTYPE, dimension(N0:N5,NV)       :: Lambda
  REALTYPE, dimension(:)              :: Lambda_Max
  REALTYPE, dimension(:), ALLOCATABLE :: Lambda_p

  REALTYPE, dimension(1), OPTIONAL    :: V_Moving_Frame
  REALTYPE, dimension(1)              :: Normal

  integer ,               OPTIONAL    :: MPI_Comm_Type

  integer  :: i, n, NN

  character(LEN=32) :: Operation

      Normal = (/ ONE /)

      call PS_Euler_Eigenvalue (N0, N5, NV, Normal, U, C, Lambda,  &
                                V_Moving_Frame)

      ALLOCATE (Lambda_p(1:NV))
 
      do n = 1,NV
        Lambda_p(n) = MAXVAL(ABS(Lambda(:,n)))
      enddo

  Operation = 'MAX' ; NN = SIZE(Lambda_p)
#if 0
  call PS_MPI_AllReduce (Operation, Index, NN, Lambda_p, Lambda_Max, MPI_Comm_Type)
#else
  call PS_MPI_AllReduce (Operation, Index,     Lambda_p, Lambda_Max, MPI_Comm_Type)
#endif

  DEALLOCATE (Lambda_p)

  END Subroutine Maximum_Eigenvalue_1D 
!
!==========================================================================
!
  Subroutine Maximum_Eigenvalue_2D (Index, N0, N5, M0, M5, NV,           &
                                           U, V, C,                      &
                                    Lambda_Max, V_Moving_Frame,          &
                                    MPI_Comm_Type)

  integer  :: Index
  integer  :: N0, N5
  integer  :: M0, M5

  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5)      :: U, V, C
  REALTYPE, dimension(N0:N5,M0:M5,NV)   :: Lambda
  REALTYPE, dimension(:,:)              :: Lambda_Max
  REALTYPE, dimension(:,:), ALLOCATABLE :: Lambda_p

  REALTYPE, dimension(2),  OPTIONAL     :: V_Moving_Frame
  REALTYPE, dimension(2)                :: Normal

  integer , OPTIONAL :: MPI_Comm_Type

  integer  :: i, j, n, NN

  character(LEN=32) :: Operation

  SELECT CASE (Index)
    CASE (1)
      Normal = (/ ONE, ZERO /)

      call PS_Euler_Eigenvalue (N0, N5, M0, M5, NV, Normal, U, V, C, Lambda,  &
                                V_Moving_Frame)

      ALLOCATE (Lambda_p(M0:M5,1:NV))
 
      do n = 1,NV
        do j = M0,M5
          Lambda_p(j,n) = MAXVAL(ABS(Lambda(:,j,n)))
        enddo
      enddo

    CASE (2)
      Normal = (/ ZERO, ONE /)

      call PS_Euler_Eigenvalue (N0, N5, M0, M5, NV, Normal, U, V, C, Lambda,  &
                                V_Moving_Frame)

      ALLOCATE (Lambda_p(N0:N5,1:NV))

      do n = 1,NV
        do i = N0,N5
          Lambda_p(i,n) = MAXVAL(ABS(Lambda(i,:,n)))
        enddo
      enddo

  END SELECT

  Operation = 'MAX' ; NN = SIZE(Lambda_p)
#if 0
  call PS_MPI_AllReduce (Operation, Index, NN, Lambda_p, Lambda_Max, MPI_Comm_Type)
#else
  call PS_MPI_AllReduce (Operation, Index,     Lambda_p, Lambda_Max, MPI_Comm_Type)
#endif

  DEALLOCATE (Lambda_p)

  END Subroutine Maximum_Eigenvalue_2D 
!
!==========================================================================
!
  Subroutine Maximum_Eigenvalue_3D (Index, N0, N5, M0, M5, K0, K5, NV,   &
                                           U, V, W, C,                   &
                                    Lambda_Max, V_Moving_Frame,          &
                                    MPI_Comm_Type)

  integer  :: Index
  integer  :: N0, N5
  integer  :: M0, M5
  integer  :: K0, K5

  integer  :: NV

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)      :: U, V, W, C
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)   :: Lambda
  REALTYPE, dimension(:,:,:)                  :: Lambda_Max
  REALTYPE, dimension(:,:,:), ALLOCATABLE     :: Lambda_p

  REALTYPE, dimension(3), OPTIONAL            :: V_Moving_Frame
  REALTYPE, dimension(3)                      :: Normal

  integer , OPTIONAL :: MPI_Comm_Type

  integer  :: i, j, k, n, NN

  character(LEN=32) :: Operation

  SELECT CASE (Index)
    CASE (1)
      Normal = (/ ONE, ZERO, ZERO /)

      call PS_Euler_Eigenvalue (N0, N5, M0, M5, K0, K5, NV,               &
                                Normal, U, V, W, C, Lambda, V_Moving_Frame)

      ALLOCATE (Lambda_p(M0:M5,K0:K5,1:NV))
 
      do n = 1,NV
        do k = K0,K5
          do j = M0,M5
            Lambda_p(j,k,n) = MAXVAL(ABS(Lambda(:,j,k,n)))
          enddo
        enddo
      enddo

    CASE (2)
      Normal = (/ ZERO, ONE, ZERO /)

      call PS_Euler_Eigenvalue (N0, N5, M0, M5, K0, K5, NV,               &
                                Normal, U, V, W, C, Lambda, V_Moving_Frame)

      ALLOCATE (Lambda_p(N0:N5,K0:K5,1:NV))

      do n = 1,NV
        do k = K0,K5
          do i = N0,N5
            Lambda_p(i,k,n) = MAXVAL(ABS(Lambda(i,:,k,n)))
          enddo
        enddo
      enddo

    CASE (3)
      Normal = (/ ZERO, ZERO, ONE /)

      call PS_Euler_Eigenvalue (N0, N5, M0, M5, K0, K5, NV,               &
                                Normal, U, V, W, C, Lambda, V_Moving_Frame)

      ALLOCATE (Lambda_p(N0:N5,M0:M5,1:NV))

      do n = 1,NV
        do j = M0,M5
          do i = N0,N5
            Lambda_p(i,j,n) = MAXVAL(ABS(Lambda(i,j,:,n)))
          enddo
        enddo
      enddo

  END SELECT

  Operation = 'MAX' ; NN = SIZE(Lambda_p)
#if 0
  call PS_MPI_AllReduce (Operation, Index, NN, Lambda_p, Lambda_Max, MPI_Comm_Type)
#else
  call PS_MPI_AllReduce (Operation, Index,     Lambda_p, Lambda_Max, MPI_Comm_Type)
#endif

  DEALLOCATE (Lambda_p)

  END Subroutine Maximum_Eigenvalue_3D 

END MODULE Euler_Function
