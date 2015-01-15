#include "pseudopack.h"


Module Euler_Eigensystem 

#if defined (PARALLEL_MPI)
  USE Processor
#endif

implicit NONE

INTERFACE PS_Euler_Eigensystem_1_Species
  MODULE PROCEDURE Eigensystem_1D_1S
  MODULE PROCEDURE Eigensystem_2D_1S
  MODULE PROCEDURE Eigensystem_3D_1S
END INTERFACE

INTERFACE PS_Euler_Eigensystem_2_Species
  MODULE PROCEDURE Eigensystem_1D_2S
  MODULE PROCEDURE Eigensystem_2D_2S
  MODULE PROCEDURE Eigensystem_3D_2S
END INTERFACE

PRIVATE
PUBLIC  :: PS_Euler_Eigensystem_1_Species
PUBLIC  :: PS_Euler_Eigensystem_2_Species
#if 0
PUBLIC  :: PS_MPI_Euler_Eigenvalue
#endif

CONTAINS

#if 0

  Subroutine Eigensystem_1D_1S (NV, Normal,                            &
                                r, U, h, Gamma,                        &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE, dimension(1) :: Normal

  REALTYPE :: N_x, r, U, h, Gamma
  REALTYPE :: gm1, q, c, e, UdotN, M_x
  REALTYPE :: b1, b2, c1, u1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

  N_x = Normal(1)

#if defined (DEBUG)
  if (N_x /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      UdotN = U*N_x 

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1 
      c1 = HALF/c        ; u1 = U*b1 

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE

      E_R(3,1) =  h     + c*UdotN
      E_R(3,2) =  h     - c*UdotN
      E_R(3,3) = -q*M_x + U

      E_L(1,1) =     - c1*UdotN + b2
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =                  b1

      E_L(2,1) =       c1*UdotN + b2
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =                  b1

      E_L(3,1) =  UdotN*N_x-TWO*U*b2
      E_L(3,2) =            TWO*U*u1 
      E_L(3,3) =           -TWO*U*b1

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
      endif

  END Subroutine Eigensystem_1D_1S
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_2D_1S (NV, Normal,                          &
                                r, U, V, h, Gamma,                   &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE, dimension(2) :: Normal

  REALTYPE :: N_x, N_y, r, U, V, h, Gamma
  REALTYPE :: gm1, q, c, e, UdotN, M_x, M_y
  REALTYPE :: b1, b2, c1, u1, v1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

  N_x = Normal(1) ; N_y = Normal(2) 

#if defined (DEBUG)
  if (N_x**2 + N_y**2 /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U + V*V)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      UdotN = U*N_x + V*N_y 

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN ; M_y = N_y/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1 
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = M_y

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO

      E_R(3,1) = V    + c*N_y
      E_R(3,2) = V    - c*N_y
      E_R(3,3) = ZERO
      E_R(3,4) = ONE

      E_R(4,1) =  h     + c*UdotN
      E_R(4,2) =  h     - c*UdotN
      E_R(4,3) = -q*M_x + U
      E_R(4,4) = -q*M_y + V

      E_L(1,1) =     - c1*UdotN + b2
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =       c1*N_y   - v1
      E_L(1,4) =                  b1

      E_L(2,1) =       c1*UdotN + b2
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =     - c1*N_y   - v1
      E_L(2,4) =                  b1

      E_L(3,1) =  UdotN*N_x-TWO*U*b2
      E_L(3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(3,3) =            TWO*U*v1 - N_x*N_y
      E_L(3,4) =           -TWO*U*b1

      E_L(4,1) =  UdotN*N_y-TWO*V*b2
      E_L(4,2) =            TWO*V*u1 - N_y*N_x
      E_L(4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(4,4) =           -TWO*V*b1

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
      endif

  END Subroutine Eigensystem_2D_1S
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_3D_1S (NV, Normal,                          &
                                r, U, V, W, h, Gamma,                &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE, dimension(3) :: Normal

  REALTYPE :: N_x, N_y, N_z, r, U, V, W, h, Gamma
  REALTYPE :: gm1, q, c, e, UdotN, M_x, M_y, M_z
  REALTYPE :: b1, b2, c1, u1, v1, w1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

  N_x = Normal(1) ; N_y = Normal(2) ; N_z = Normal(3)

#if defined (DEBUG)
  if (N_x**2 + N_y**2 + N_z**2 /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U + V*V + W*W)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      UdotN = U*N_x + V*N_y + W*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1 
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 ; w1 = W*b1

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = M_y
      E_R(1,5) = M_z

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO
      E_R(2,5) = ZERO

      E_R(3,1) = V    + c*N_y
      E_R(3,2) = V    - c*N_y
      E_R(3,3) = ZERO
      E_R(3,4) = ONE
      E_R(3,5) = ZERO

      E_R(4,1) = W    + c*N_z
      E_R(4,2) = W    - c*N_z
      E_R(4,3) = ZERO
      E_R(4,4) = ZERO
      E_R(4,5) = ONE

      E_R(5,1) =  h     + c*UdotN
      E_R(5,2) =  h     - c*UdotN
      E_R(5,3) = -q*M_x + U
      E_R(5,4) = -q*M_y + V
      E_R(5,5) = -q*M_z + W

      E_L(1,1) =     - c1*UdotN + b2
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =       c1*N_y   - v1
      E_L(1,4) =       c1*N_z   - w1
      E_L(1,5) =                  b1

      E_L(2,1) =       c1*UdotN + b2
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =     - c1*N_y   - v1
      E_L(2,4) =     - c1*N_z   - w1
      E_L(2,5) =                  b1

      E_L(3,1) =  UdotN*N_x-TWO*U*b2
      E_L(3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(3,3) =            TWO*U*v1 - N_x*N_y
      E_L(3,4) =            TWO*U*w1 - N_x*N_z
      E_L(3,5) =           -TWO*U*b1

      E_L(4,1) =  UdotN*N_y-TWO*V*b2
      E_L(4,2) =            TWO*V*u1 - N_y*N_x
      E_L(4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(4,4) =            TWO*V*w1 - N_y*N_z
      E_L(4,5) =           -TWO*V*b1

      E_L(5,1) =  UdotN*N_z-TWO*W*b2
      E_L(5,2) =            TWO*W*u1 - N_z*N_x
      E_L(5,3) =            TWO*W*v1 - N_z*N_y 
      E_L(5,4) =            TWO*W*w1 - N_z*N_z + ONE
      E_L(5,5) =           -TWO*W*b1

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
        Lambda(5) = UdotN
      endif

  END Subroutine Eigensystem_3D_1S 
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_1D_2S (NV, Normal,                            &
                                r, U, h, f1, Gamma, D_Gamma,           &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE, dimension(1) :: Normal

  REALTYPE :: N_x, r, U, h, f1, Gamma, D_Gamma
  REALTYPE :: gm1, q, c, e, X, Y, UdotN, M_x
  REALTYPE :: b1, b2, b3, b4, c1, u1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

  N_x = Normal(1) 

#if defined (DEBUG)
  if (N_x /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = D_Gamma*e ; Y = X/gm1

      UdotN = U*N_x 

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN 

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = ZERO

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO

      E_R(3,1) =  h     + c*UdotN
      E_R(3,2) =  h     - c*UdotN
      E_R(3,3) = -q*M_x + U
      E_R(3,4) =        - Y

      E_R(4,1) = f1
      E_R(4,2) = f1
      E_R(4,3) = f1*M_x
      E_R(4,4) = ONE

      E_L(1,1) =     - c1*UdotN + b4 
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =                  b1
      E_L(1,4) =                  b3

      E_L(2,1) =       c1*UdotN + b4
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =                  b1
      E_L(2,4) =                  b3

      E_L(3,1) =  UdotN*N_x-TWO*U*b4
      E_L(3,2) =            TWO*U*u1 
      E_L(3,3) =           -TWO*U*b1
      E_L(3,4) =           -TWO*U*b3

      E_L(4,1) = -f1
      E_L(4,2) =  ZERO
      E_L(4,3) =  ZERO
      E_L(4,4) =  ONE

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
      endif

  END Subroutine Eigensystem_1D_2S
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_2D_2S (NV, Normal,                          &
                                r, U, V, h, f1, Gamma, D_Gamma,      &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE, dimension(2) :: Normal

  REALTYPE :: N_x, N_y, r, U, V, h, f1, Gamma, D_Gamma
  REALTYPE :: gm1, q, c, e, X, Y, UdotN, M_x, M_y
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

  N_x = Normal(1) ; N_y = Normal(2) 

#if defined (DEBUG)
  if (N_x**2 + N_y**2 /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U + V*V)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = D_Gamma*e ; Y = X/gm1

      UdotN = U*N_x + V*N_y 

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN ; M_y = N_y/UdotN 

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = M_y
      E_R(1,5) = ZERO

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO
      E_R(2,5) = ZERO

      E_R(3,1) = V    + c*N_y
      E_R(3,2) = V    - c*N_y
      E_R(3,3) = ZERO
      E_R(3,4) = ONE
      E_R(3,5) = ZERO

      E_R(4,1) =  h     + c*UdotN
      E_R(4,2) =  h     - c*UdotN
      E_R(4,3) = -q*M_x + U
      E_R(4,4) = -q*M_y + V
      E_R(4,5) =        - Y

      E_R(5,1) = f1
      E_R(5,2) = f1
      E_R(5,3) = f1*M_x
      E_R(5,4) = f1*M_y
      E_R(5,5) = ONE

      E_L(1,1) =     - c1*UdotN + b4
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =       c1*N_y   - v1
      E_L(1,4) =                  b1
      E_L(1,5) =                  b3

      E_L(2,1) =       c1*UdotN + b4
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =     - c1*N_y   - v1
      E_L(2,4) =                  b1
      E_L(2,5) =                  b3

      E_L(3,1) =  UdotN*N_x-TWO*U*b4
      E_L(3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(3,3) =            TWO*U*v1 - N_x*N_y
      E_L(3,4) =           -TWO*U*b1
      E_L(3,5) =           -TWO*U*b3

      E_L(4,1) =  UdotN*N_y-TWO*V*b4
      E_L(4,2) =            TWO*V*u1 - N_y*N_x
      E_L(4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(4,4) =           -TWO*V*b1
      E_L(4,5) =           -TWO*V*b3

      E_L(5,1) = -f1
      E_L(5,2) =  ZERO
      E_L(5,3) =  ZERO
      E_L(5,4) =  ZERO
      E_L(5,5) =  ONE

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
        Lambda(5) = UdotN
      endif

  END Subroutine Eigensystem_2D_2S
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_3D_2S (NV, Normal,                          &
                                r, U, V, W, h, f1, Gamma, D_Gamma,   &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE, dimension(3) :: Normal

  REALTYPE :: N_x, N_y, N_z, r, U, V, W, h, f1, Gamma, D_Gamma
  REALTYPE :: gm1, q, c, e, X, Y, UdotN, M_x, M_y, M_z
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1, w1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

  N_x = Normal(1) ; N_y = Normal(2) ; N_z = Normal(3)

#if defined (DEBUG)
  if (N_x**2 + N_y**2 + N_z**2 /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U + V*V + W*W)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = D_Gamma*e ; Y = X/gm1

      UdotN = U*N_x + V*N_y + W*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 ; w1 = W*b1

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = M_y
      E_R(1,5) = M_z
      E_R(1,6) = ZERO

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO
      E_R(2,5) = ZERO
      E_R(2,6) = ZERO

      E_R(3,1) = V    + c*N_y
      E_R(3,2) = V    - c*N_y
      E_R(3,3) = ZERO
      E_R(3,4) = ONE
      E_R(3,5) = ZERO
      E_R(3,6) = ZERO

      E_R(4,1) = W    + c*N_z
      E_R(4,2) = W    - c*N_z
      E_R(4,3) = ZERO
      E_R(4,4) = ZERO
      E_R(4,5) = ONE
      E_R(4,6) = ZERO

      E_R(5,1) =  h     + c*UdotN
      E_R(5,2) =  h     - c*UdotN
      E_R(5,3) = -q*M_x + U
      E_R(5,4) = -q*M_y + V
      E_R(5,5) = -q*M_z + W
      E_R(5,6) =        - Y

      E_R(6,1) = f1
      E_R(6,2) = f1
      E_R(6,3) = f1*M_x
      E_R(6,4) = f1*M_y
      E_R(6,5) = f1*M_z
      E_R(6,6) = ONE

      E_L(1,1) =     - c1*UdotN + b4
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =       c1*N_y   - v1
      E_L(1,4) =       c1*N_z   - w1
      E_L(1,5) =                  b1
      E_L(1,6) =                  b3

      E_L(2,1) =       c1*UdotN + b4
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =     - c1*N_y   - v1
      E_L(2,4) =     - c1*N_z   - w1
      E_L(2,5) =                  b1
      E_L(2,6) =                  b3

      E_L(3,1) =  UdotN*N_x-TWO*U*b4
      E_L(3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(3,3) =            TWO*U*v1 - N_x*N_y
      E_L(3,4) =            TWO*U*w1 - N_x*N_z
      E_L(3,5) =           -TWO*U*b1
      E_L(3,6) =           -TWO*U*b3

      E_L(4,1) =  UdotN*N_y-TWO*V*b4
      E_L(4,2) =            TWO*V*u1 - N_y*N_x
      E_L(4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(4,4) =            TWO*V*w1 - N_y*N_z
      E_L(4,5) =           -TWO*V*b1
      E_L(4,6) =           -TWO*V*b3

      E_L(5,1) =  UdotN*N_z-TWO*W*b4
      E_L(5,2) =            TWO*W*u1 - N_z*N_x
      E_L(5,3) =            TWO*W*v1 - N_z*N_y 
      E_L(5,4) =            TWO*W*w1 - N_z*N_z + ONE
      E_L(5,5) =           -TWO*W*b1
      E_L(5,6) =           -TWO*W*b3

      E_L(6,1) = -f1
      E_L(6,2) =  ZERO
      E_L(6,3) =  ZERO
      E_L(6,4) =  ZERO
      E_L(6,5) =  ZERO
      E_L(6,6) =  ONE

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
        Lambda(5) = UdotN
        Lambda(6) = UdotN
      endif

  END Subroutine Eigensystem_3D_2S 

END Module Euler_Eigensystem 

!********************************************************************
#else
!********************************************************************

  Subroutine Eigensystem_1D_1S (NV, N_x,                               &
                                r, U, h, Gamma,                        &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE :: N_x, r, U, h, Gamma
  REALTYPE :: gm1, q, c, e, UdotN, M_x
  REALTYPE :: b1, b2, c1, u1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

#if defined (DEBUG)
  if (N_x /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      UdotN = U*N_x 

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1 
      c1 = HALF/c        ; u1 = U*b1 

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE

      E_R(3,1) =  h     + c*UdotN
      E_R(3,2) =  h     - c*UdotN
      E_R(3,3) = -q*M_x + U

      E_L(1,1) =     - c1*UdotN + b2
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =                  b1

      E_L(2,1) =       c1*UdotN + b2
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =                  b1

      E_L(3,1) =  UdotN*N_x-TWO*U*b2
      E_L(3,2) =            TWO*U*u1 
      E_L(3,3) =           -TWO*U*b1

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
      endif

  END Subroutine Eigensystem_1D_1S
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_2D_1S (NV, N_x, N_y,                        &
                                r, U, V, h, Gamma,                   &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE :: N_x, N_y, r, U, V, h, Gamma
  REALTYPE :: gm1, q, c, e, UdotN, M_x, M_y
  REALTYPE :: b1, b2, c1, u1, v1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

#if defined (DEBUG)
  if (N_x**2 + N_y**2 /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U + V*V)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      UdotN = U*N_x + V*N_y 

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN ; M_y = N_y/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1 
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = M_y

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO

      E_R(3,1) = V    + c*N_y
      E_R(3,2) = V    - c*N_y
      E_R(3,3) = ZERO
      E_R(3,4) = ONE

      E_R(4,1) =  h     + c*UdotN
      E_R(4,2) =  h     - c*UdotN
      E_R(4,3) = -q*M_x + U
      E_R(4,4) = -q*M_y + V

      E_L(1,1) =     - c1*UdotN + b2
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =       c1*N_y   - v1
      E_L(1,4) =                  b1

      E_L(2,1) =       c1*UdotN + b2
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =     - c1*N_y   - v1
      E_L(2,4) =                  b1

      E_L(3,1) =  UdotN*N_x-TWO*U*b2
      E_L(3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(3,3) =            TWO*U*v1 - N_x*N_y
      E_L(3,4) =           -TWO*U*b1

      E_L(4,1) =  UdotN*N_y-TWO*V*b2
      E_L(4,2) =            TWO*V*u1 - N_y*N_x
      E_L(4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(4,4) =           -TWO*V*b1

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
      endif

  END Subroutine Eigensystem_2D_1S
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_3D_1S (NV, N_x, N_y, N_z,                   &
                                r, U, V, W, h, Gamma,                &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE :: N_x, N_y, N_z, r, U, V, W, h, Gamma
  REALTYPE :: gm1, q, c, e, UdotN, M_x, M_y, M_z
  REALTYPE :: b1, b2, c1, u1, v1, w1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

#if defined (DEBUG)
  if (N_x**2 + N_y**2 + N_z**2 /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U + V*V + W*W)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      UdotN = U*N_x + V*N_y + W*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1 
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 ; w1 = W*b1

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = M_y
      E_R(1,5) = M_z

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO
      E_R(2,5) = ZERO

      E_R(3,1) = V    + c*N_y
      E_R(3,2) = V    - c*N_y
      E_R(3,3) = ZERO
      E_R(3,4) = ONE
      E_R(3,5) = ZERO

      E_R(4,1) = W    + c*N_z
      E_R(4,2) = W    - c*N_z
      E_R(4,3) = ZERO
      E_R(4,4) = ZERO
      E_R(4,5) = ONE

      E_R(5,1) =  h     + c*UdotN
      E_R(5,2) =  h     - c*UdotN
      E_R(5,3) = -q*M_x + U
      E_R(5,4) = -q*M_y + V
      E_R(5,5) = -q*M_z + W

      E_L(1,1) =     - c1*UdotN + b2
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =       c1*N_y   - v1
      E_L(1,4) =       c1*N_z   - w1
      E_L(1,5) =                  b1

      E_L(2,1) =       c1*UdotN + b2
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =     - c1*N_y   - v1
      E_L(2,4) =     - c1*N_z   - w1
      E_L(2,5) =                  b1

      E_L(3,1) =  UdotN*N_x-TWO*U*b2
      E_L(3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(3,3) =            TWO*U*v1 - N_x*N_y
      E_L(3,4) =            TWO*U*w1 - N_x*N_z
      E_L(3,5) =           -TWO*U*b1

      E_L(4,1) =  UdotN*N_y-TWO*V*b2
      E_L(4,2) =            TWO*V*u1 - N_y*N_x
      E_L(4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(4,4) =            TWO*V*w1 - N_y*N_z
      E_L(4,5) =           -TWO*V*b1

      E_L(5,1) =  UdotN*N_z-TWO*W*b2
      E_L(5,2) =            TWO*W*u1 - N_z*N_x
      E_L(5,3) =            TWO*W*v1 - N_z*N_y 
      E_L(5,4) =            TWO*W*w1 - N_z*N_z + ONE
      E_L(5,5) =           -TWO*W*b1

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
        Lambda(5) = UdotN
      endif

  END Subroutine Eigensystem_3D_1S 
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_1D_2S (NV, N_x,                               &
                                r, U, h, f1, Gamma, D_Gamma,           &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE :: N_x, r, U, h, f1, Gamma, D_Gamma
  REALTYPE :: gm1, q, c, e, X, Y, UdotN, M_x
  REALTYPE :: b1, b2, b3, b4, c1, u1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

#if defined (DEBUG)
  if (N_x /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = D_Gamma*e ; Y = X/gm1

      UdotN = U*N_x 

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN 

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = ZERO

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO

      E_R(3,1) =  h     + c*UdotN
      E_R(3,2) =  h     - c*UdotN
      E_R(3,3) = -q*M_x + U
      E_R(3,4) =        - Y

      E_R(4,1) = f1
      E_R(4,2) = f1
      E_R(4,3) = f1*M_x
      E_R(4,4) = ONE

      E_L(1,1) =     - c1*UdotN + b4 
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =                  b1
      E_L(1,4) =                  b3

      E_L(2,1) =       c1*UdotN + b4
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =                  b1
      E_L(2,4) =                  b3

      E_L(3,1) =  UdotN*N_x-TWO*U*b4
      E_L(3,2) =            TWO*U*u1 
      E_L(3,3) =           -TWO*U*b1
      E_L(3,4) =           -TWO*U*b3

      E_L(4,1) = -f1
      E_L(4,2) =  ZERO
      E_L(4,3) =  ZERO
      E_L(4,4) =  ONE

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
      endif

  END Subroutine Eigensystem_1D_2S
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_2D_2S (NV, N_x, N_y,                        &
                                r, U, V, h, f1, Gamma, D_Gamma,      &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE :: N_x, N_y, r, U, V, h, f1, Gamma, D_Gamma
  REALTYPE :: gm1, q, c, e, X, Y, UdotN, M_x, M_y
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

#if defined (DEBUG)
  if (N_x**2 + N_y**2 /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U + V*V)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = D_Gamma*e ; Y = X/gm1

      UdotN = U*N_x + V*N_y 

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN ; M_y = N_y/UdotN 

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = M_y
      E_R(1,5) = ZERO

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO
      E_R(2,5) = ZERO

      E_R(3,1) = V    + c*N_y
      E_R(3,2) = V    - c*N_y
      E_R(3,3) = ZERO
      E_R(3,4) = ONE
      E_R(3,5) = ZERO

      E_R(4,1) =  h     + c*UdotN
      E_R(4,2) =  h     - c*UdotN
      E_R(4,3) = -q*M_x + U
      E_R(4,4) = -q*M_y + V
      E_R(4,5) =        - Y

      E_R(5,1) = f1
      E_R(5,2) = f1
      E_R(5,3) = f1*M_x
      E_R(5,4) = f1*M_y
      E_R(5,5) = ONE

      E_L(1,1) =     - c1*UdotN + b4
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =       c1*N_y   - v1
      E_L(1,4) =                  b1
      E_L(1,5) =                  b3

      E_L(2,1) =       c1*UdotN + b4
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =     - c1*N_y   - v1
      E_L(2,4) =                  b1
      E_L(2,5) =                  b3

      E_L(3,1) =  UdotN*N_x-TWO*U*b4
      E_L(3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(3,3) =            TWO*U*v1 - N_x*N_y
      E_L(3,4) =           -TWO*U*b1
      E_L(3,5) =           -TWO*U*b3

      E_L(4,1) =  UdotN*N_y-TWO*V*b4
      E_L(4,2) =            TWO*V*u1 - N_y*N_x
      E_L(4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(4,4) =           -TWO*V*b1
      E_L(4,5) =           -TWO*V*b3

      E_L(5,1) = -f1
      E_L(5,2) =  ZERO
      E_L(5,3) =  ZERO
      E_L(5,4) =  ZERO
      E_L(5,5) =  ONE

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
        Lambda(5) = UdotN
      endif

  END Subroutine Eigensystem_2D_2S
!
!-----------------------------------------------------------------------
!
  Subroutine Eigensystem_3D_2S (NV, N_x, N_y, N_z,                   &
                                r, U, V, W, h, f1, Gamma, D_Gamma,   &
                                E_R, E_L, Lambda)

  integer  :: NV
  REALTYPE :: N_x, N_y, N_z, r, U, V, W, h, f1, Gamma, D_Gamma
  REALTYPE :: gm1, q, c, e, X, Y, UdotN, M_x, M_y, M_z
  REALTYPE :: b1, b2, b3, b4, c1, u1, v1, w1
 
  REALTYPE, dimension(NV,NV)           :: E_R, E_L 
  REALTYPE, dimension(NV)   , OPTIONAL :: Lambda

#if defined (DEBUG)
  if (N_x**2 + N_y**2 + N_z**2 /= ONE) then
    write (6 ,*) 'STOP!  The normal was not normalized to 1.'
    write (99,*) 'STOP!  The normal was not normalized to 1.'
    STOP 
  endif
#endif

      gm1 = Gamma-ONE

      q = HALF*(U*U + V*V + W*W)
      c = SQRT(ABS(gm1*(h - q)))
      e = c**2/(Gamma*gm1)

      X  = D_Gamma*e ; Y = X/gm1

      UdotN = U*N_x + V*N_y + W*N_z

!     if (UdotN == ZERO) UdotN = EPSILON(ONE)
      if (ABS(UdotN) <= EPSILON(ONE)) UdotN = EPSILON(ONE)
  
      M_x = N_x/UdotN ; M_y = N_y/UdotN ; M_z = N_z/UdotN

      b1 = HALF*gm1/c**2 ; b2 = q*b1 ; b3 = Y*b1 ; b4 = b2 - f1*b3
      c1 = HALF/c        ; u1 = U*b1 ; v1 = V*b1 ; w1 = W*b1

      E_R(1,1) = ONE
      E_R(1,2) = ONE
      E_R(1,3) = M_x
      E_R(1,4) = M_y
      E_R(1,5) = M_z
      E_R(1,6) = ZERO

      E_R(2,1) = U    + c*N_x
      E_R(2,2) = U    - c*N_x
      E_R(2,3) = ONE
      E_R(2,4) = ZERO
      E_R(2,5) = ZERO
      E_R(2,6) = ZERO

      E_R(3,1) = V    + c*N_y
      E_R(3,2) = V    - c*N_y
      E_R(3,3) = ZERO
      E_R(3,4) = ONE
      E_R(3,5) = ZERO
      E_R(3,6) = ZERO

      E_R(4,1) = W    + c*N_z
      E_R(4,2) = W    - c*N_z
      E_R(4,3) = ZERO
      E_R(4,4) = ZERO
      E_R(4,5) = ONE
      E_R(4,6) = ZERO

      E_R(5,1) =  h     + c*UdotN
      E_R(5,2) =  h     - c*UdotN
      E_R(5,3) = -q*M_x + U
      E_R(5,4) = -q*M_y + V
      E_R(5,5) = -q*M_z + W
      E_R(5,6) =        - Y

      E_R(6,1) = f1
      E_R(6,2) = f1
      E_R(6,3) = f1*M_x
      E_R(6,4) = f1*M_y
      E_R(6,5) = f1*M_z
      E_R(6,6) = ONE

      E_L(1,1) =     - c1*UdotN + b4
      E_L(1,2) =       c1*N_x   - u1
      E_L(1,3) =       c1*N_y   - v1
      E_L(1,4) =       c1*N_z   - w1
      E_L(1,5) =                  b1
      E_L(1,6) =                  b3

      E_L(2,1) =       c1*UdotN + b4
      E_L(2,2) =     - c1*N_x   - u1
      E_L(2,3) =     - c1*N_y   - v1
      E_L(2,4) =     - c1*N_z   - w1
      E_L(2,5) =                  b1
      E_L(2,6) =                  b3

      E_L(3,1) =  UdotN*N_x-TWO*U*b4
      E_L(3,2) =            TWO*U*u1 - N_x*N_x + ONE
      E_L(3,3) =            TWO*U*v1 - N_x*N_y
      E_L(3,4) =            TWO*U*w1 - N_x*N_z
      E_L(3,5) =           -TWO*U*b1
      E_L(3,6) =           -TWO*U*b3

      E_L(4,1) =  UdotN*N_y-TWO*V*b4
      E_L(4,2) =            TWO*V*u1 - N_y*N_x
      E_L(4,3) =            TWO*V*v1 - N_y*N_y + ONE
      E_L(4,4) =            TWO*V*w1 - N_y*N_z
      E_L(4,5) =           -TWO*V*b1
      E_L(4,6) =           -TWO*V*b3

      E_L(5,1) =  UdotN*N_z-TWO*W*b4
      E_L(5,2) =            TWO*W*u1 - N_z*N_x
      E_L(5,3) =            TWO*W*v1 - N_z*N_y 
      E_L(5,4) =            TWO*W*w1 - N_z*N_z + ONE
      E_L(5,5) =           -TWO*W*b1
      E_L(5,6) =           -TWO*W*b3

      E_L(6,1) = -f1
      E_L(6,2) =  ZERO
      E_L(6,3) =  ZERO
      E_L(6,4) =  ZERO
      E_L(6,5) =  ZERO
      E_L(6,6) =  ONE

      if (PRESENT(Lambda)) then
        Lambda(1) = UdotN+c
        Lambda(2) = UdotN-c
        Lambda(3) = UdotN
        Lambda(4) = UdotN
        Lambda(5) = UdotN
        Lambda(6) = UdotN
      endif

  END Subroutine Eigensystem_3D_2S 
#endif



#if 0
!
!-----------------------------------------------------------------------
!
  Subroutine PS_MPI_Euler_Eigenvalue (Index, Lambda_p, Lambda, MPI_Comm_Type)

  integer                             :: Index
  REALTYPE, dimension(:)              :: Lambda_p, Lambda

  integer ,                  OPTIONAL :: MPI_Comm_Type
#if defined (PARALLEL_MPI)
  integer                             :: MPI_Communicator
  integer                             :: MPI_Communicator_1D
  integer                             :: N_Dims, Status
  logical , dimension(:), ALLOCATABLE :: Sub_Grid
#endif

  Lambda = Lambda_p

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) RETURN

                                MPI_Communicator = MPI_Comm_World
    if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

    call MPI_TOPO_Test (MPI_Communicator, Status, MPI_Error_Status)

    if (Status == MPI_CART) then
      call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

      ALLOCATE (Sub_Grid(N_Dims))

      Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

      call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                     MPI_Error_Status)

      call MPI_ALLREDUCE (Lambda_p, Lambda, SIZE(Lambda),                 &
                          MPI_REALTYPE, MPI_MAX, MPI_Communicator_1D,     &
                          MPI_Error_Status)

      DEALLOCATE (Sub_Grid)

      call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)
    else
      call MPI_ALLREDUCE (Lambda_p, Lambda, SIZE(Lambda),                 &
                          MPI_REALTYPE, MPI_MAX, MPI_Communicator,        &
                          MPI_Error_Status)
    endif
#endif

  END Subroutine PS_MPI_Euler_Eigenvalue
#endif

END Module Euler_Eigensystem 
