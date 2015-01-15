#include "pseudopack.h"

! ----------------------------------------------------------------------
! FUNCTION NAME: Characeristic_Decomposition
! AUTHOR       : Wai-Sun Don
! DESCRIPTION  : Characteristic Boundary Condition based on the
!                Linearized Euler Equations.
!                  U_t + Div (F(U)) = 1/Re Div (Grad(F_v(U))
!
! Q      Assumed-shape Array of size (3)     for 1D problem, 
!                                    (:,4)   for 2D problem,
!                                    (:,:,5) for 3D problem,
!          For the 2D problem, the first dimension's size is the number
!            of grid points along one dimension boundary line.
!
!          For the 3D problem, the first and second dimension's sizes are 
!            the number of grid points on the two dimension boundary surface.
!
!          The last dimension's size is the number of conservative variables.
!           =3 for one   dimensional  (Rho, RhoU, E)
!           =4 for two   dimensional  (Rho, RhoU, RhoV, E)
!           =5 for three dimensional  (Rho, RhoU, RhoV, RhoW, E)
! 
! Q0     Assumed-shape Array of size (3)     for 1D problem, 
!                                    (:,4)   for 2D problem,
!                                    (:,:,5) for 3D problem,
!        The linearized quantity of Q
!    
! Normal Assumed-shape Array of size (1)     for 1D problem, 
!                                    (:,2)   for 2D problem,
!                                    (:,:,3) for 3D problem,
!        Unit outward normal vector
!        OPTIONAL
!
! gamma  Specific heat ratio
!        OPTIONAL
!
!*********************************************************************

MODULE Characteristic_Decomposition

implicit NONE

INTERFACE OPERATOR (.dot.)
  MODULE PROCEDURE   Dot_Product_1D
END INTERFACE

INTERFACE OPERATOR (.cross.)
  MODULE PROCEDURE Cross_Product
END INTERFACE

INTERFACE Characteristic_Boundary
  MODULE PROCEDURE Characteristic_BC_Point
  MODULE PROCEDURE Characteristic_BC_Line
  MODULE PROCEDURE Characteristic_BC_Surface
  MODULE PROCEDURE Characteristic_BC_Line_1
  MODULE PROCEDURE Characteristic_BC_Surface_1
END INTERFACE

INTERFACE Characteristic_Variable
  MODULE PROCEDURE Characteristic_Variable_1D
  MODULE PROCEDURE Characteristic_Variable_2D
  MODULE PROCEDURE Characteristic_Variable_3D
END INTERFACE

INTERFACE Similarity_Transform
  MODULE PROCEDURE Similarity_Transform_2D
END INTERFACE

INTERFACE Similarity_Transform_Inverse
  MODULE PROCEDURE Similarity_Transform_Inverse_2D
END INTERFACE

INTERFACE Euler_Jacobian
  MODULE PROCEDURE Euler_Jacobian_2D
END INTERFACE

INTERFACE Viscous_Jacobian
  MODULE PROCEDURE Viscous_Jacobian_2D
END INTERFACE

INTERFACE Viscous_Variable
  MODULE PROCEDURE Viscous_Variable_2D
END INTERFACE

integer                           :: Case_3D

REALTYPE, dimension(3), parameter :: i = (/ ONE , ZERO, ZERO /)
REALTYPE, dimension(3), parameter :: j = (/ ZERO, ONE , ZERO /)
REALTYPE, dimension(3), parameter :: k = (/ ZERO, ZERO, ONE  /)

PRIVATE
PUBLIC  :: Similarity_Transform
PUBLIC  :: Similarity_Transform_Inverse
PUBLIC  :: Eigenvalue
PUBLIC  :: Euler_Jacobian
PUBLIC  :: Characteristic
PUBLIC  :: Characteristic_Variable
PUBLIC  :: Viscous_Jacobian
PUBLIC  :: Viscous_Variable
PUBLIC  :: Characteristic_Boundary

CONTAINS

!******************************************************
  FUNCTION   Dot_Product_1D (A, B) RESULT (C)

  REALTYPE, dimension(:), INTENT(IN) :: A, B
  REALTYPE                           :: C

  C = DOT_PRODUCT(A, B)

  END FUNCTION Dot_Product_1D

!******************************************************
  FUNCTION Cross_Product (A, B) RESULT (C)

  REALTYPE, dimension(:), INTENT(IN) :: A, B
  REALTYPE, dimension(SIZE(A))       :: C

  SELECT CASE (SIZE(A))
    CASE DEFAULT
      C(1) =   A(1)*B(2) - A(2)*B(1)
      C(2) =   ZERO

    CASE (3)
      C(1) =  (A(2)*B(3) - A(3)*B(2))
      C(2) = -(A(1)*B(3) - A(3)*B(1))
      C(3) =  (A(1)*B(2) - A(2)*B(1))

  END SELECT

  END FUNCTION Cross_Product
!
!=========================================================================
!
  FUNCTION Sound_Speed (Q, gamma) RESULT (C)

  REALTYPE, dimension(:)           :: Q 
  REALTYPE              , OPTIONAL :: gamma

  integer                          :: V
  REALTYPE                         :: g, gm1
  REALTYPE                         :: D, E, P, C
  REALTYPE, dimension(SIZE(Q)-2)   :: M, U

  V = SIZE(Q)

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  gm1 = g - ONE

  D = Q(1) ; M = Q(2:V-1) ; E = Q(V)

  U = M/D
  P = gm1*(HALF*D*(U .dot. U) - (M .dot. U) + E)
  C = SQRT(g*P/D)

  END FUNCTION Sound_Speed
!
!=========================================================================
!
  FUNCTION Similarity_Transform_2D (Q, Normal, gamma) RESULT (S)

  REALTYPE, dimension(:)                   :: Q
  REALTYPE, dimension(SIZE(Q),SIZE(Q))     :: S
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  integer                          :: V
  REALTYPE                         :: g, gm1, a
  REALTYPE                         :: D, E, P, C, H
  REALTYPE, dimension(SIZE(Q)-2)   :: M, U, N, K

  V = SIZE(Q)

  if (PRESENT(Normal)) then ; N = Normal ; else ; N = ZERO ; N(1) = -ONE ; endif

  K(1) = N(2) ; K(2) =  -N(1)

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  gm1 = g - ONE

  D = Q(1) ; M = Q(2:V-1) ; E = Q(V)

  U = Q(2:V-1)/Q(1)
  C = Sound_Speed (Q, g)
  P = gm1*(HALF*D*(U .dot. U) - (M .dot. U) + E)
  H = (E+P)/D

  a = HALF/C

  S(1,1) = ONE 
  S(1,2) = ZERO 
  S(1,3) = a
  S(1,4) = a

  S(2,1) =    U(1) 
  S(2,2) =           K(1)
  S(2,3) = a*(U(1)-C*N(1)) 
  S(2,4) = a*(U(1)+C*N(1))

  S(3,1) =    U(2) 
  S(3,2) =           K(2)
  S(3,3) = a*(U(2)-C*N(2)) 
  S(3,4) = a*(U(2)+C*N(2))

  S(4,1) =     HALF*(U .dot. U) 
  S(4,2) =          (U .dot. K)
  S(4,3) = a*(H - C*(U .dot. N)) 
  S(4,4) = a*(H + C*(U .dot. N))

  END FUNCTION Similarity_Transform_2D 
!
!=========================================================================
!
  FUNCTION Similarity_Transform_Inverse_2D (Q, Normal, gamma) RESULT (S_Inverse)

  REALTYPE, dimension(:)                   :: Q
  REALTYPE, dimension(SIZE(Q),SIZE(Q))     :: S_Inverse
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  integer                          :: V
  REALTYPE                         :: C, g, gm1, U2, t, s
  REALTYPE, dimension(SIZE(Q)-2)   :: U, N, K

  V = SIZE(Q)

  if (PRESENT(Normal)) then ; N = Normal ; else ; N = ZERO ; N(1) = -ONE ; endif

  K(1) = N(2) ; K(2) =  -N(1)

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  gm1 = g - ONE

  U = Q(2:V-1)/Q(1)
  C = Sound_Speed (Q, g)

  U2 = HALF*(U .dot. U) ; t = gm1/C**2 ; s = gm1/C

  S_Inverse(1,1) =          ONE - t*U2
  S_Inverse(1,2) =                t*U(1)
  S_Inverse(1,3) =                t*U(2)
  S_Inverse(1,4) =              - t

  S_Inverse(2,1) = -(U .dot. K)
  S_Inverse(2,2) =   K(1)
  S_Inverse(2,3) =   K(2)
  S_Inverse(2,4) = ZERO

  S_Inverse(3,1) =  (U .dot. N) + s*U2 
  S_Inverse(3,2) = -(N(1)       + s*U(1))
  S_Inverse(3,3) = -(N(2)       + s*U(2))
  S_Inverse(3,4) =                s

  S_Inverse(4,1) = -(U .dot. N) + s*U2 
  S_Inverse(4,2) =  (N(1)       - s*U(1))
  S_Inverse(4,3) =  (N(2)       - s*U(2))
  S_Inverse(4,4) =                s

  END FUNCTION Similarity_Transform_Inverse_2D
!
!=========================================================================
!
  FUNCTION Euler_Jacobian_2D (Q, Normal, gamma) RESULT (J)

  REALTYPE, dimension(:)                   :: Q
  REALTYPE, dimension(SIZE(Q),SIZE(Q))     :: J
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  integer                          :: V
  REALTYPE                         :: C, g, gm1, gm2, U2, t
  REALTYPE, dimension(SIZE(Q)-2)   :: U, N, K

  V = SIZE(Q)

  if (PRESENT(Normal)) then ; N = Normal ; else ; N = ZERO ; N(1) = -ONE ; endif

  K(1) = N(2) ; K(2) =  -N(1)

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  gm1 = g - ONE ; gm2 = g - TWO

  U = Q(2:V-1)/Q(1)
  C = Sound_Speed (Q, g)

  U2 = HALF*(U .dot. U) ; t = gm1/C**2

  J(1,1) = ZERO
  J(1,2) = N(1)
  J(1,3) = N(2)
  J(1,4) = ZERO

  J(2,1) =  gm1*U2  *N(1) - U(1)*(U .dot. N)
  J(2,2) = -gm2*U(1)*N(1) +      (U .dot. N)
  J(2,3) = -gm2*U(2)*N(1) +      (U .dot. K)
  J(2,4) =  gm1*     N(1)

  J(3,1) =  gm1*U2  *N(2) - U(2)*(U .dot. N)
  J(3,2) = -gm2*U(1)*N(2) -      (U .dot. K)
  J(3,3) = -gm2*U(2)*N(2) +      (U .dot. N)
  J(3,4) =  gm1*     N(2)

  J(4,1) =                    (gm2*U2-ONE/t)*(U .dot. N)
  J(4,2) = (U2 + ONE/t)*N(1) -      gm1*U(1)*(U .dot. N) 
  J(4,3) = (U2 + ONE/t)*N(2) -      gm1*U(2)*(U .dot. N) 
  J(4,4) =                    (gm1+ONE)     *(U .dot. N)    

  END FUNCTION Euler_Jacobian_2D
!
!=========================================================================
!
  FUNCTION Viscous_Jacobian_2D (Q, lambda, mu, Kappa, Pr, Normal, gamma) &
                                RESULT (J)

  REALTYPE, dimension(:)                   :: Q
  REALTYPE, dimension(SIZE(Q),SIZE(Q),2)   :: J
  REALTYPE                                 :: lambda, mu, Kappa, Pr
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  integer                          :: V
  REALTYPE                         :: D, C, g, gm1, U2, t, c22, c32, gkPr
  REALTYPE, dimension(SIZE(Q)-2)   :: U, N, K

  REALTYPE, dimension(SIZE(Q),SIZE(Q)) :: J_xx, J_xy, J_yy

  V = SIZE(Q)

  if (PRESENT(Normal)) then ; N = Normal ; else ; N = ZERO ; N(1) = -ONE ; endif

  K(1) = N(2) ; K(2) =  -N(1)

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  gm1 = g - ONE 

  D = Q(1) 

  U = Q(2:V-1)/Q(1)
  C = Sound_Speed (Q, g)

  U2 = HALF*(U .dot. U) ; t = gm1/C**2

  c22  = HALF*(lambda +     mu)
  c32  = HALF*(lambda + TWO*mu)
  gkPr = g*Kappa/Pr

  J_xx = ZERO ; J_xy = ZERO ; J_yy = ZERO

  J_xx(2,1) = -TWO*c32*U(1)
  J_xx(2,2) =  TWO*c32

  J_xx(3,1) = -mu*U(2)
  J_xx(3,3) =  mu

  J_xx(4,1) = -TWO*c22*U(1)**2 - TWO*mu*U2 + (U2-ONE/(t*g))*gkPr
  J_xx(4,2) = (TWO*c32-gkPr)*U(1)
  J_xx(4,3) = (     mu-gkPr)*U(2)
  J_xx(4,4) = gkPr

  J_xx = J_xx/D

  J_xy(2,1) = -U(2)
  J_xy(2,3) = ONE

  J_xy(3,1) = -U(1)
  J_xy(3,2) = ONE

  J_xy(4,1) = -TWO*U(1)*U(2)
  J_xy(4,2) =  U(2)
  J_xy(4,3) =  U(1)

  J_xy = J_xy/D*c22

  J_yy(2,1) = -mu*U(1)
  J_yy(2,2) =  mu

  J_yy(3,1) = -TWO*c32*U(2)
  J_yy(3,3) =  TWO*c32

  J_yy(4,1) = -TWO*c22*U(2)**2 - TWO*mu*U2 + (U2-ONE/(t*g))*gkPr
  J_yy(4,2) = (     mu-gkPr)*U(1)
  J_yy(4,2) = (TWO*c32-gkPr)*U(2)
  J_yy(4,4) = gkPr

  J_yy = J_yy/D

  J(:,:,1) = J_xx*N(1) + J_xy*N(2)
  J(:,:,2) = J_xy*N(1) + J_yy*N(2)

  END FUNCTION Viscous_Jacobian_2D

#if 1
!
!=========================================================================
!
  FUNCTION Viscous_Variable_2D (Q, R_x, R_y,                          &
                                lambda, mu, Kappa, Pr, Normal, gamma) &
                                RESULT (R_v)

  REALTYPE, dimension(:)                   :: Q, R_x, R_y
  REALTYPE, dimension(SIZE(Q))             :: R_v
  REALTYPE                                 :: lambda, mu, Kappa, Pr
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  REALTYPE, dimension(SIZE(Q),SIZE(Q),2) :: J
  REALTYPE, dimension(SIZE(Q),SIZE(Q))   :: S, S_Inverse, B_x, B_y

  J         = Viscous_Jacobian (Q, lambda, mu, Kappa, Pr, Normal, gamma)

  S         = Similarity_Transform         (Q, Normal, gamma)
  S_Inverse = Similarity_Transform_Inverse (Q, Normal, gamma)

  B_x = MATMUL(J(:,:,1),S)
  B_x = MATMUL(S_Inverse, B_x)

  B_y = MATMUL(J(:,:,2),S)
  B_y = MATMUL(S_Inverse, B_y)

  R_v = MATMUL(B_x, R_x) + MATMUL(B_y, R_y)

  END FUNCTION Viscous_Variable_2D
#if 0
!
!=========================================================================
!
  FUNCTION Viscous_Variable_3D (Q, R_x, R_y, R_z,                     &
                                lambda, mu, Kappa, Pr, Normal, gamma) &
                                RESULT (R_v)

  REALTYPE, dimension(:)                   :: Q, R_x, R_y, R_z
  REALTYPE, dimension(SIZE(Q))             :: R_v
  REALTYPE                                 :: lambda, mu, Kappa, Pr
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  REALTYPE, dimension(SIZE(Q),SIZE(Q),SIZE(Q),3) :: J
  REALTYPE, dimension(SIZE(Q),SIZE(Q),SIZE(Q))   :: S, S_Inverse, B_x, B_y, B_z

  J         = Viscous_Jacobian (Q, lambda, mu, Kappa, Pr, Normal, gamma)

  S         = Similarity_Transform         (Q, Normal, gamma)
  S_Inverse = Similarity_Transform_Inverse (Q, Normal, gamma)

  B_x = MATMUL(J(:,:,:,1),S)
  B_x = MATMUL(S_Inverse, B_x)

  B_y = MATMUL(J(:,:,:,2),S)
  B_y = MATMUL(S_Inverse, B_y)

  B_z = MATMUL(J(:,:,:,3),S)
  B_z = MATMUL(S_Inverse, B_z)

  R_v = MATMUL(B_x, R_x) + MATMUL(B_y, R_y) + MATMUL(B_z, R_z)

  END FUNCTION Viscous_Variable_3D
#endif
!
!=========================================================================
!
  FUNCTION Eigenvalue (Q, Normal, gamma) RESULT (Lambda)

  REALTYPE, dimension(:)                   :: Q
  REALTYPE, dimension(SIZE(Q))             :: Lambda
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  integer                          :: V
  REALTYPE                         :: C, g
  REALTYPE, dimension(SIZE(Q)-2)   :: U, N

  V = SIZE(Q)

  if (PRESENT(Normal)) then ; N = Normal ; else ; N = ZERO ; N(1) = -ONE ; endif

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  U = Q(2:V-1)/Q(1)
  C = Sound_Speed (Q, g)

  Lambda(1:V-2) =  U .dot. N
  Lambda(  V-1) = (U .dot. N) - C
  Lambda(  V  ) = (U .dot. N) + C

  END FUNCTION Eigenvalue
!
!=========================================================================
!
  FUNCTION Characteristic (Q, Q_0, Normal, gamma) RESULT (R)

  REALTYPE, dimension(:)                   :: Q, Q_0
  REALTYPE, dimension(SIZE(Q))             :: R
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  REALTYPE, dimension(SIZE(Q),SIZE(Q)) :: S_Inverse

  S_Inverse = Similarity_Transform_Inverse (     Q_0, Normal, gamma)

  R = MATMUL(S_Inverse, Q)

  END FUNCTION Characteristic 
!
!=========================================================================
!
  Subroutine Characteristic_BC_Point (Q, Q_0, Normal, gamma) 

  REALTYPE, dimension(:)                   :: Q, Q_0
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  REALTYPE, dimension(SIZE(Q))         :: Lambda, R_Inflow, R_Outflow, R
  REALTYPE, dimension(SIZE(Q),SIZE(Q)) :: S, S_Inverse

  S         = Similarity_Transform         (Q_0, Normal, gamma)
  S_Inverse = Similarity_Transform_Inverse (Q_0, Normal, gamma)
  Lambda    = Eigenvalue                   (Q_0, Normal, gamma)

  R_Inflow  = MATMUL(S_Inverse, Q_0)
  R_Outflow = MATMUL(S_Inverse, Q  )

  WHERE (Lambda <  ZERO) R =       R_Inflow
  WHERE (Lambda == ZERO) R = HALF*(R_Inflow + R_Outflow)
  WHERE (Lambda >  ZERO) R =                  R_Outflow

  Q = MATMUL(S, R)

  END Subroutine Characteristic_BC_Point 
!
!=========================================================================
!
  Subroutine Characteristic_BC_Line (Q, Q_0, Normal, gamma)

  REALTYPE, dimension(:,:)          , TARGET :: Q, Q_0
  REALTYPE, dimension(:,:), OPTIONAL, TARGET :: Normal
  REALTYPE, dimension(:)  , OPTIONAL         :: gamma

  integer                         :: i
  integer                         :: N1, NV, NN
  REALTYPE, dimension(:), POINTER :: p_Q, p_Q_0, p_N

  N1 = SIZE(Q, DIM=1) ; NV = SIZE(Q, DIM=2) 
  NN = SIZE(Normal,DIM=2)

#if defined (DEBUG)
  if (PRESENT(Normal) .AND. (NN /= NV-2)) &
    call Dimension_Check ('Characteristic_BC_Line', NN, NV-2)
#endif

  Nullify (p_Q, p_Q_0, p_N)

  do i = 1, N1

    p_Q   => Q  (i,:)
    p_Q_0 => Q_0(i,:)

    if (PRESENT(Normal)) then
      p_N     => Normal(i,:)

      if (PRESENT(gamma)) then
        call Characteristic_BC_Point (p_Q, p_Q_0, p_N, gamma=gamma(i))
      else
        call Characteristic_BC_Point (p_Q, p_Q_0, p_N                )
      endif
    else
      if (PRESENT(gamma)) then
        call Characteristic_BC_Point (p_Q, p_Q_0,      gamma=gamma(i))
      else
        call Characteristic_BC_Point (p_Q, p_Q_0                     )
      endif
    endif

  enddo

  Nullify (p_Q, p_Q_0, p_N)

  END Subroutine Characteristic_BC_Line
!
!=========================================================================
!
  Subroutine Characteristic_BC_Surface (Q, Q_0, Normal, gamma)

  REALTYPE, dimension(:,:,:)          , TARGET :: Q, Q_0
  REALTYPE, dimension(:,:,:), OPTIONAL, TARGET :: Normal
  REALTYPE, dimension(:,:)  , OPTIONAL         :: gamma

  integer                           :: i
  integer                           :: N1, NV, NN
  REALTYPE, dimension(:,:), POINTER :: p_Q, p_Q_0, p_N

  N1 = SIZE(Q, DIM=1) ; NV = SIZE(Q, DIM=3) 
  NN = SIZE(Normal,DIM=3)

#if defined (DEBUG)
  if (PRESENT(Normal) .AND. (NN /= NV-2)) &
    call Dimension_Check ('Characteristic_BC_Surface', NN, NV-2)
#endif

  Nullify (p_Q, p_Q_0, p_N)

  do i = 1, N1

    p_Q   => Q  (i,:,:)
    p_Q_0 => Q_0(i,:,:)

    if (PRESENT(Normal)) then
      p_N     => Normal(i,:,:)
      
      if (PRESENT(gamma)) then
        call Characteristic_BC_Line  (p_Q, p_Q_0, p_N, gamma=gamma(i,:))
      else
        call Characteristic_BC_Line  (p_Q, p_Q_0, p_N                  )
      endif
    else
      if (PRESENT(gamma)) then
        call Characteristic_BC_Line  (p_Q, p_Q_0,      gamma=gamma(i,:))
      else
        call Characteristic_BC_Line  (p_Q, p_Q_0                       )
      endif
    endif

  enddo

  Nullify (p_Q, p_Q_0, p_N)

  END Subroutine Characteristic_BC_Surface
!
!=========================================================================
!
  Subroutine Characteristic_BC_Line_1 (Q, Q_0, gamma, Normal)

  REALTYPE, dimension(:,:)          , TARGET :: Q, Q_0
  REALTYPE                                   :: gamma
  REALTYPE, dimension(:,:), OPTIONAL, TARGET :: Normal

  integer                         :: i
  integer                         :: N1, NV, NN
  REALTYPE, dimension(:), POINTER :: p_Q, p_Q_0, p_N

  N1 = SIZE(Q, DIM=1) ; NV = SIZE(Q, DIM=2) 
  NN = SIZE(Normal,DIM=2)

#if defined (DEBUG)
  if (PRESENT(Normal) .AND. (NN /= NV-2)) &
    call Dimension_Check ('Characteristic_BC_Line_1', NN, NV-2)
#endif

  Nullify (p_Q, p_Q_0, p_N)

  do i = 1, N1

    p_Q   => Q  (i,:)
    p_Q_0 => Q_0(i,:)

    if (PRESENT(Normal)) then
      p_N     => Normal(i,:)
      
      call Characteristic_BC_Point (p_Q, p_Q_0, p_N, gamma)
    else
      call Characteristic_BC_Point (p_Q, p_Q_0,      gamma=gamma)
    endif

  enddo

  Nullify (p_Q, p_Q_0, p_N)

  END Subroutine Characteristic_BC_Line_1
!
!=========================================================================
!
  Subroutine Characteristic_BC_Surface_1 (Q, Q_0, gamma, Normal)

  REALTYPE, dimension(:,:,:)          , TARGET :: Q, Q_0
  REALTYPE                                     :: gamma
  REALTYPE, dimension(:,:,:), OPTIONAL, TARGET :: Normal

  integer                           :: i
  integer                           :: N1, NV, NN
  REALTYPE, dimension(:,:), POINTER :: p_Q, p_Q_0, p_N

  N1 = SIZE(Q, DIM=1) ; NV = SIZE(Q, DIM=3) 
  NN = SIZE(Normal,DIM=3)

#if defined (DEBUG)
  if (PRESENT(Normal) .AND. (NN /= NV-2)) &
    call Dimension_Check ('Characteristic_BC_Surface_1', NN, NV-2)
#endif

  Nullify (p_Q, p_Q_0, p_N)

  do i = 1, N1

    p_Q   => Q  (i,:,:)
    p_Q_0 => Q_0(i,:,:)

    if (PRESENT(Normal)) then
      p_N     => Normal(i,:,:)
      
      call Characteristic_BC_Line_1 (p_Q, p_Q_0, gamma, p_N)
    else
      call Characteristic_BC_Line_1 (p_Q, p_Q_0, gamma     )
    endif

  enddo

  Nullify (p_Q, p_Q_0, p_N)

  END Subroutine Characteristic_BC_Surface_1
!
!=========================================================================
!
  Subroutine Characteristic_Variable_0D (NN, NV, Q, Q_0, R, Normal, gamma)

  integer                                      :: NN, NV
  REALTYPE, dimension(NN,NV)          , TARGET :: Q, Q_0, R
  REALTYPE, dimension(NN,NV), OPTIONAL, TARGET :: Normal
  REALTYPE, dimension(NN)   , OPTIONAL         :: gamma

  integer                         :: i
  REALTYPE, dimension(:), POINTER :: p_Q, p_Q_0, p_N

  Nullify (p_Q, p_Q_0, p_N)

  do i = 1, NN

    p_Q   => Q  (i,:)
    p_Q_0 => Q_0(i,:)

    if (PRESENT(Normal)) then
      p_N     => Normal(i,:)
      
      if (PRESENT(gamma)) then
        R(i,:) = Characteristic (p_Q, p_Q_0, p_N, gamma=gamma(i))
      else
        R(i,:) = Characteristic (p_Q, p_Q_0, p_N                )
      endif
    else
      if (PRESENT(gamma)) then
        R(i,:) = Characteristic (p_Q, p_Q_0,      gamma=gamma(i))
      else
        R(i,:) = Characteristic (p_Q, p_Q_0                     )
      endif
    endif

  enddo

  Nullify (p_Q, p_Q_0, p_N)

  END Subroutine Characteristic_Variable_0D
!
!=========================================================================
!
  Subroutine Characteristic_Variable_1D (Q, Q_0, R, Normal, gamma)

  REALTYPE, dimension(:,:)           :: Q, Q_0, R
  REALTYPE, dimension(:,:), OPTIONAL :: Normal
  REALTYPE, dimension(:)  , OPTIONAL :: gamma

  integer  :: N1, NV, NN

  N1 = SIZE(Q, DIM=1) ; NV = SIZE(Q, DIM=2)
  NN = SIZE(Normal,DIM=2)

#if defined (DEBUG)
  if (PRESENT(Normal) .AND. (NN /= NV-2)) &
    call Dimension_Check ('Characteristic_Variable_1D', NN, NV-2)
#endif

  call Characteristic_Variable_0D (N1, NV, Q, Q_0, R, Normal, gamma)

  END Subroutine Characteristic_Variable_1D
!
!=========================================================================
!
  Subroutine Characteristic_Variable_2D (Q, Q_0, R, Normal, gamma)

  REALTYPE, dimension(:,:,:)          , TARGET :: Q, Q_0, R
  REALTYPE, dimension(:,:,:), OPTIONAL, TARGET :: Normal
  REALTYPE, dimension(:,:)  , OPTIONAL, TARGET :: gamma

  integer  :: N1, N2, NV, NN

  N1 = SIZE(Q, DIM=1) ; N2 = SIZE(Q, DIM=2) ; NV = SIZE(Q, DIM=3)
  NN = SIZE(Normal,DIM=3)

#if defined (DEBUG)
  if (PRESENT(Normal) .AND. (NN /= NV-2)) &
    call Dimension_Check ('Characteristic_Variable_2D', NN, NV-2)
#endif

  call Characteristic_Variable_0D (N1*N2, NV, Q, Q_0, Normal, gamma)

  END Subroutine Characteristic_Variable_2D
!
!=========================================================================
!
  Subroutine Characteristic_Variable_3D (Q, Q_0, R, Normal, gamma)

  REALTYPE, dimension(:,:,:,:)          , TARGET :: Q, Q_0, R
  REALTYPE, dimension(:,:,:,:), OPTIONAL, TARGET :: Normal
  REALTYPE, dimension(:,:,:)  , OPTIONAL, TARGET :: gamma

  integer  :: N1, N2, N3, NV, NN

  N1 = SIZE(Q, DIM=1) ; N2 = SIZE(Q, DIM=2) ; N3 = SIZE(Q, DIM=3)
  NV = SIZE(Q, DIM=4) 
  NN = SIZE(Normal,DIM=4)

#if defined (DEBUG)
  if (PRESENT(Normal) .AND. (NN /= NV-2)) &
    call Dimension_Check ('Characteristic_Variable_3D', NN, NV-2)
#endif

  call Characteristic_Variable_0D (N1*N2*N3, NV, Q, Q_0, Normal, gamma)

  END Subroutine Characteristic_Variable_3D
#endif
!
!=========================================================================
!
  Subroutine Dimension_Check (Subroutine_Name, N1, N2)

  integer          :: N1, N2
  Character(LEN=*) :: Subroutine_Name

  write (6 ,100) TRIM(Subroutine_Name), N1, N2
  write (99,100) TRIM(Subroutine_Name), N1, N2
  STOP

 100 format (1x,72('*')/  &
             1x,'FATAL : ',(a) /&
             1x,'        Incompatible dimensions of Q and the Normal vector'/&
             1x,'        Normal dimension     = ',i1,' should be equal to' /&
             1x,'        Q      dimension - 2 = ',i1 / &
             1x,72('*')/)

  END Subroutine Dimension_Check


!***************************************************************************!
!                    OBSOLETE AND OLD CODE                                  !
!***************************************************************************!

#if 0
!
!=========================================================================
!
  FUNCTION Viscous_Variable_2D (Q, R_x, R_y,                      &
                                lambda, mu, Kappa, Pr, Normal, gamma) &
                                RESULT (R_v)

  REALTYPE, dimension(:)                   :: Q, R_x, R_y
  REALTYPE, dimension(SIZE(Q))             :: R_v
  REALTYPE                                 :: lambda, mu, Kappa, Pr
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  integer                          :: V
  REALTYPE                         :: D, C, g, gm1
  REALTYPE                         :: gkPr, c00, c22, c32, c44, c55
  REALTYPE, dimension(SIZE(Q)-2)   :: N, K

  REALTYPE, dimension(SIZE(Q),SIZE(Q),2) :: J
  REALTYPE, dimension(SIZE(Q),SIZE(Q))   :: S, S_Inverse, B_x, B_y

  V = SIZE(Q)

  if (PRESENT(Normal)) then ; N = Normal ; else ; N = ZERO ; N(1) = -ONE ; endif

  K(1) = N(2) ; K(2) =  -N(1)

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  gm1 = g - ONE 

  D = Q(1) 
  C = Sound_Speed (Q, g)

  gkPr =  g*Kappa/Pr

  c00  =  gm1/g
  c22  =  HALF*(lambda +     mu)
  c32  =  HALF*(lambda + TWO*mu)
  c44  =  c32 + HALF*gkPr*c00
  c55  = -c32 + HALF*gkPr*c00

  B_x(1,1) =  gkPr/g          *N(1)/D
  B_x(1,1) =  ZERO            *N(1)/D
  B_x(1,3) = -HALF*gkPr*c00/C *N(1)/D
  B_x(1,4) = -HALF*gkPr*c00/C *N(1)/D

  B_x(2,1) =  ZERO            
  B_x(2,2) =  mu              *N(1)/D
  B_x(2,3) = -c22             *N(2)/D
  B_x(2,4) =  c22             *N(2)/D

  B_x(3,1) = -C*gkPr/g        *N(1)/D
  B_x(3,2) = -TWO*c22         *N(2)/D
  B_x(3,3) =  c44             *N(1)/D
  B_x(3,4) =  c55             *N(1)/D

  B_x(4,1) = -C*gkPr/g        *N(1)/D
  B_x(4,2) =  TWO*c22         *N(2)/D
  B_x(4,3) =  c55             *N(1)/D
  B_x(4,4) =  c44             *N(1)/D

  B_y(1,1) =  gkPr/g          *N(2)/D
  B_y(1,1) =  ZERO            *N(2)/D
  B_y(1,3) = -HALF*gkPr*c00/C *N(2)/D
  B_y(1,4) = -HALF*gkPr*c00/C *N(2)/D

  B_y(2,1) =  ZERO            
  B_y(2,2) =  mu              *N(2)/D
  B_y(2,3) =  c22             *N(1)/D
  B_y(2,4) = -c22             *N(1)/D

  B_y(3,1) = -C*gkPr/g        *N(2)/D
  B_y(3,2) =  TWO*c22         *N(1)/D
  B_y(3,3) =  c44             *N(2)/D
  B_y(3,4) =  c55             *N(2)/D

  B_y(4,1) = -C*gkPr/g        *N(2)/D
  B_y(4,2) = -TWO*c22         *N(1)/D
  B_y(4,3) =  c55             *N(2)/D
  B_y(4,4) =  c44             *N(2)/D

  R_v(1) = (gkPr*c00)/(TWO*C)*(                                &
             (TWO*C/(g*c00)*R_x(1)-R_x(3)+R_x(4))*N(1) +       &
             (TWO*C/(g*c00)*R_y(1)-R_y(3)+R_y(4))*N(2) ) /D

  R_v(2) = ((mu*R_x(2)*N(1)-c22*(R_x(3)-R_x(4))*K(1)) +        &
            (mu*R_y(2)*N(2)-c22*(R_y(3)-R_y(4))*K(2)) ) /D

  R_v(3) = ( (( (-C*gkPr/g)*R_x(1)+c44*R_x(3)+c55*R_x(4))*N(1) -  &
                TWO*c22*R_x(2)*K(1)) +                            &
             (( (-C*gkPr/g)*R_y(1)+c44*R_y(3)+c55*R_y(4))*N(2) -  &
                TWO*c22*R_y(2)*K(2)) ) /D

  R_v(4) = ( (( (-C*gkPr/g)*R_x(1)+c44*R_x(3)+c55*R_x(4))*N(1) +  &
                TWO*c22*R_x(2)*K(1)) +                            &
             (( (-C*gkPr/g)*R_y(1)+c44*R_y(3)+c55*R_y(4))*N(2) +  &
                TWO*c22*R_y(2)*K(2)) ) /D

  END FUNCTION Viscous_Variable_2D
!
!=========================================================================
!
  FUNCTION Eigenvalue (Q, Normal, gamma) RESULT (Lambda)

  REALTYPE, dimension(:)                   :: Q
  REALTYPE, dimension(SIZE(Q))             :: Lambda
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  integer                          :: V
  REALTYPE                         :: C, g
  REALTYPE, dimension(SIZE(Q)-2)   :: U, N

  V = SIZE(Q)

  if (PRESENT(Normal)) then ; N = Normal ; else ; N = ZERO ; N(1) = -ONE ; endif

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  U = Q(2:V-1)/Q(1)
  C = Sound_Speed (Q, g)

  Lambda(1:V-2) =  U .dot. N
  Lambda(  V-1) = (U .dot. N) - C
  Lambda(  V  ) = (U .dot. N) + C

  END FUNCTION Eigenvalue
!
!=========================================================================
!
  FUNCTION Characteristic (Q, Q_0, Normal, gamma) RESULT (R)

  REALTYPE, dimension(:)                   :: Q, Q_0
  REALTYPE, dimension(SIZE(Q))             :: R
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  integer                          :: V
  integer                          :: N_x_i, N_x_j, N_x_k
  REALTYPE                         :: g, gm1
  REALTYPE                         :: D, E, P, C
  REALTYPE, dimension(SIZE(Q)-2)   :: U_0, M, N, Cross

  V = SIZE(Q)

  if (PRESENT(Normal)) then ; N = Normal ; else ; N = ZERO ; N(1) = -ONE ; endif

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  gm1 = g - ONE 

  U_0 = Q_0(2:V-1)/Q_0(1)

  D = Q(1) ; M = Q(2:V-1) ; E = Q(V)

  P = gm1*(HALF*D*(U_0 .dot. U_0) - (M .dot. U_0) + E)
  C = Sound_Speed (Q_0, g)

  R(1)   =        D*C*C           - P
  R(V-1) = -((M - D*U_0) .dot. N) + P/C
  R(V)   =  ((M - D*U_0) .dot. N) + P/C

  Cross  =   (M - D*U_0) .cross. N

  SELECT CASE (V)
    CASE (4)
      R(2) =  Cross(1)

    CASE (5)
      N_x_i = 1 ; if (SUM(N .cross. i) == ZERO) N_x_i = 0
      N_x_j = 1 ; if (SUM(N .cross. j) == ZERO) N_x_j = 0
      N_x_k = 1 ; if (SUM(N .cross. k) == ZERO) N_x_k = 0

      SELECT CASE (N_x_i)
        CASE DEFAULT
          R(2) = Cross(1) 

          SELECT CASE (N_x_j)
            CASE DEFAULT
              R(3) = Cross(2) ; Case_3D = 12 ; RETURN

            CASE (0)
              if (N_x_k == 1) then
                R(3) = Cross(3) ; Case_3D = 13 ; RETURN
              endif

          END SELECT

        CASE (0)
          SELECT CASE (N_x_j)
            CASE DEFAULT
              R(2) = Cross(2) 
  
            CASE (0)
              if (N_x_k == 1) then
                R(3) = Cross(3) ; Case_3D = 23 ; RETURN
              endif

          END SELECT

      END SELECT

      write (6 ,100) N_x_i, N_x_j, N_x_k
      write (99,100) N_x_i, N_x_j, N_x_k

  END SELECT

 100 format (1x,72('*')/  &
             1x,'FATAL : Characteristic_Variable_3D'/ &
             1x,'        Unable to set the characteristic variable R(2), R(3)'/&
             1x,'        The normal vector N is PARALLEL to two axis'/&
             1x,'        N cross i = ',G7.2 / &
             1x,'        N cross j = ',G7.2 / &
             1x,'        N cross k = ',G7.2 / &
             1x,72('*')/)

  END FUNCTION Characteristic
!
!=========================================================================
!
  Subroutine Characteristic_BC_Point (Q, Q_0, Normal, gamma) 

  REALTYPE, dimension(:)                   :: Q, Q_0
  REALTYPE, dimension(SIZE(Q)-2), OPTIONAL :: Normal
  REALTYPE                      , OPTIONAL :: gamma

  integer                          :: V
  REALTYPE                         :: g, gm1, C
  REALTYPE, dimension(SIZE(Q)-2)   :: U, N
  REALTYPE, dimension(SIZE(Q))     :: Lambda, R_Inflow, R_Outflow, F, Z

  V = SIZE(Q)

  if (PRESENT(Normal)) then ; N = Normal ; else ; N = ZERO ; N(1) = -ONE ; endif

  if (PRESENT(gamma)) then ; g = gamma ; else ; g = 1.40d0 ; endif 

  gm1 = g - ONE 

  U         = Q_0(2:V-1)/Q_0(1)
  C         = Sound_Speed    (Q_0,              g)
  Lambda    = Eigenvalue     (Q_0,      Normal, g)
  R_Inflow  = Characteristic (Q_0, Q_0, Normal, g)
  R_Outflow = Characteristic (Q  , Q_0, Normal, g)

  WHERE (Lambda <  ZERO) F =       R_Inflow
  WHERE (Lambda == ZERO) F = HALF*(R_Inflow + R_Outflow)
  WHERE (Lambda >  ZERO) F =                  R_Outflow

  SELECT CASE (V)
    CASE DEFAULT
      Z(1) = C*(F(3)+F(2))/TWO
      Z(2) =   (F(3)-F(2))/TWO
      Z(3) =    F(1)
  
      Q(1) = (Z(1) + Z(3))/(C*C)
      Q(2) =  Q(1)*U(1) + Z(2)*N(1)
      Q(3) =  Z(1)/gm1 + (Q(2:2) - HALF*Q(1)*U .dot. U)

    CASE (4)
      Z(1) = C*(F(4)+F(3))/TWO
      Z(2) =   (F(4)-F(3))/TWO
      Z(3) =    F(2)
      Z(4) =    F(1)
  
      Q(1) = (Z(1) + Z(4))/(C*C)
      Q(2) =  Q(1)*U(1) + (Z(2)*N(1) + Z(3)*N(2)) 
      Q(3) =  Q(1)*U(2) + (Z(2)*N(2) - Z(3)*N(1)) 
      Q(4) =  Z(1)/gm1 + (Q(2:3) - HALF*Q(1)*U .dot. U)

    CASE (5)
      Z(1) = C*(F(5)+F(4))/TWO
      Z(2) =   (F(5)-F(4))/TWO     
      Z(3) =    F(2)               
      Z(4) =    F(3)
      Z(5) =    F(1)

      Q(1) = (Z(1) + Z(5))/(C*C)

      SELECT CASE (Case_3D)
        CASE (12)
          Q(2) =  Q(1)*U(1) + (Z(2)* N(1)                  -  &
                               Z(3)* N(1)*N(2)       /N(3) -  &
                               Z(4)*(N(2)**2+N(3)**2)/N(3) )
          Q(3) =  Q(1)*U(2) + (Z(2)* N(2)                  +  &
                               Z(3)*(N(1)**2+N(3)**2)/N(3) +  &
                               Z(4)* N(1)*N(2)       /N(3) )
          Q(4) =  Q(1)*U(3) + (Z(2)* N(3)                  -  &
                               Z(3)* N(2)                  +  &
                               Z(4)* N(1)                  )

        CASE (13)
          Q(2) =  Q(1)*U(1) + (Z(2)* N(1)                  +  &
                               Z(3)* N(1)*N(3)       /N(2) +  &
                               Z(4)*(N(2)**2+N(3)**2)/N(2) )
          Q(3) =  Q(1)*U(2) + (Z(2)* N(2)                  +  &
                               Z(3)* N(3)                  -  &
                               Z(4)* N(1)                  )
          Q(4) =  Q(1)*U(3) + (Z(2)* N(3)                  -  &
                               Z(3)*(N(1)**2+N(2)**2)/N(2) -  &
                               Z(4)* N(1)*N(3)       /N(2) )

        CASE (23)
          Q(2) =  Q(1)*U(1) + (Z(2)* N(1)                  -  &
                               Z(3)* N(3)                  +  &
                               Z(4)* N(2)                  )
          Q(3) =  Q(1)*U(2) + (Z(2)* N(2)                  -  &
                               Z(3)* N(2)*N(3)       /N(1) -  &
                               Z(4)*(N(1)**2+N(3)**2)/N(1) )
          Q(4) =  Q(1)*U(3) + (Z(2)* N(3)                  +  &
                               Z(3)*(N(1)**2+N(2)**2)/N(1) +  &
                               Z(4)* N(2)*N(3)       /N(1) )
          
      END SELECT

      Q(5) =  Z(1)/gm1 + (Q(2:4) - HALF*Q(1)*U .dot. U)

  END SELECT

  END Subroutine Characteristic_BC_Point 
!
!=========================================================================
!
  Subroutine Characteristic_BC_Line (Q, Q_0, Normal, gamma) 

  REALTYPE, dimension(:,:)          , TARGET :: Q, Q_0
  REALTYPE, dimension(:,:), OPTIONAL, TARGET :: Normal
  REALTYPE                , OPTIONAL         :: gamma

  integer                         :: i
  REALTYPE, dimension(:), POINTER :: p_Q, p_Q_0, p_N

  if (PRESENT(Normal)) then
    if (SIZE(Normal, DIM=2) /= SIZE(Q, DIM=2)-2) then
      write (6 ,100) SIZE(Normal, DIM=2), SIZE(Q, DIM=2)-2
      write (99,100) SIZE(Normal, DIM=2), SIZE(Q, DIM=2)-2
      STOP
    endif
  endif

  Nullify (p_Q, p_Q_0, p_N)

  do i = 1, SIZE(Q, DIM=1)

    p_Q   => Q  (i,:)
    p_Q_0 => Q_0(i,:)

    if (PRESENT(Normal)) then
      p_N => Normal(i,:)

      call Characteristic_BC_Point (p_Q, p_Q_0, p_N, gamma) 
    else
      call Characteristic_BC_Point (p_Q, p_Q_0,      gamma=gamma)
    endif

  enddo

  Nullify (p_Q, p_Q_0, p_N)

 100 format (1x,72('*')/  &
             1x,'FATAL : Characteristic_BC_Line    '/&
             1x,'        Incompatible dimensions of Q and the normal vector'/&
             1x,'        Normal dimension     = ',i1,' should be equal to' /&
             1x,'        Q      dimension - 2 = ',i1 / &
             1x,72('*')/)

  END Subroutine Characteristic_BC_Line 
!
!=========================================================================
!
  Subroutine Characteristic_BC_Surface (Q, Q_0, Normal, gamma)

  REALTYPE, dimension(:,:,:)          , TARGET :: Q, Q_0
  REALTYPE, dimension(:,:,:), OPTIONAL, TARGET :: Normal
  REALTYPE                  , OPTIONAL         :: gamma

  integer                           :: i
  REALTYPE, dimension(:,:), POINTER :: p_Q, p_Q_0, p_N

  if (PRESENT(Normal)) then
    if (SIZE(Normal, DIM=3) /= SIZE(Q, DIM=3)-2) then
      write (6 ,100) SIZE(Normal, DIM=3), SIZE(Q, DIM=3)-2
      write (99,100) SIZE(Normal, DIM=3), SIZE(Q, DIM=3)-2
      STOP
    endif
  endif

  Nullify (p_Q, p_Q_0, p_N)

  do i = 1, SIZE(Q, DIM=1)

    p_Q   => Q  (i,:,:)
    p_Q_0 => Q_0(i,:,:)

    if (PRESENT(Normal)) then
      p_N => Normal(i,:,:)

      call Characteristic_BC_Line (p_Q, p_Q_0, p_N, gamma)
    else
      call Characteristic_BC_Line (p_Q, p_Q_0,      gamma=gamma)
    endif

  enddo

  Nullify (p_Q, p_Q_0, p_N)

 100 format (1x,72('*')/  &
             1x,'FATAL : Characteristic_BC_Surface '/&
             1x,'        Incompatible dimensions of Q and the normal vector'/&
             1x,'        Normal dimension     = ',i1,' should be equal to' /&
             1x,'        Q      dimension - 2 = ',i1 / &
             1x,72('*')/)

  END Subroutine Characteristic_BC_Surface
!
!=========================================================================
!
  Subroutine Characteristic_Variable_1D (Q, Q_0, R, Normal, gamma)

  REALTYPE, dimension(:,:)          , TARGET :: Q, Q_0, R
  REALTYPE, dimension(:,:), OPTIONAL, TARGET :: Normal
  REALTYPE, dimension(:)  , OPTIONAL, TARGET :: gamma

  integer                         :: i
  REALTYPE,               POINTER :: p_gamma
  REALTYPE, dimension(:), POINTER :: p_Q, p_Q_0, p_N

  if (PRESENT(Normal)) then
    if (SIZE(Normal, DIM=2) /= SIZE(Q, DIM=2)-2) then
      write (6 ,100) SIZE(Normal, DIM=2), SIZE(Q, DIM=2)-2
      write (99,100) SIZE(Normal, DIM=2), SIZE(Q, DIM=2)-2
      STOP
    endif
  endif

  Nullify (p_Q, p_Q_0, p_N, p_gamma)

  do i = 1, SIZE(Q, DIM=1)

    p_Q   => Q  (i,:)
    p_Q_0 => Q_0(i,:)

    if (PRESENT(Normal)) then
      p_N     => Normal(i,:)
      
      if (PRESENT(gamma)) then
        R(i,:) = Characteristic (p_Q, p_Q_0, p_N, gamma(i))
      else
        p_gamma => gamma(i)

        R(i,:) = Characteristic (p_Q, p_Q_0, p_N, p_gamma)
      endif

    else
      if (PRESENT(gamma)) then
        R(i,:) = Characteristic (p_Q, p_Q_0,      gamma=gamma(i))
      else
        p_gamma => gamma(i)

        R(i,:) = Characteristic (p_Q, p_Q_0,      gamma=p_gamma)
      endif
    endif

  enddo

  Nullify (p_Q, p_Q_0, p_N, p_gamma)

 100 format (1x,72('*')/  &
             1x,'FATAL : Characteristic_Variable_1D    '/&
             1x,'        Incompatible dimensions of Q and the normal vector'/&
             1x,'        Normal dimension     = ',i1,' should be equal to' /&
             1x,'        Q      dimension - 2 = ',i1 / &
             1x,72('*')/)

  END Subroutine Characteristic_Variable_1D
!
!=========================================================================
!
  Subroutine Characteristic_Variable_2D (Q, Q_0, R, Normal, gamma)

  REALTYPE, dimension(:,:,:)          , TARGET :: Q, Q_0, R
  REALTYPE, dimension(:,:,:), OPTIONAL, TARGET :: Normal
  REALTYPE, dimension(:,:)  , OPTIONAL, TARGET :: gamma

  integer                         :: i, j
  REALTYPE,               POINTER :: p_gamma
  REALTYPE, dimension(:), POINTER :: p_Q, p_Q_0, p_N

  if (PRESENT(Normal)) then
    if (SIZE(Normal, DIM=3) /= SIZE(Q, DIM=3)-2) then
      write (6 ,100) SIZE(Normal, DIM=3), SIZE(Q, DIM=3)-2
      write (99,100) SIZE(Normal, DIM=3), SIZE(Q, DIM=3)-2
      STOP
    endif
  endif

  Nullify (p_Q, p_Q_0, p_N, p_gamma)

  do i = 1, SIZE(Q, DIM=1)
    do j = 1, SIZE(Q, DIM=2)

      p_Q   => Q  (i,j,:)
      p_Q_0 => Q_0(i,j,:)

      if (PRESENT(Normal)) then
        p_N     => Normal(i,j,:)
      
        if (PRESENT(gamma)) then
          R(i,j,:) = Characteristic (p_Q, p_Q_0, p_N, gamma(i,j))
        else
          p_gamma => gamma(i,j)

          R(i,j,:) = Characteristic (p_Q, p_Q_0, p_N, p_gamma)
        endif

      else
        if (PRESENT(gamma)) then
          R(i,j,:) = Characteristic (p_Q, p_Q_0,      gamma=gamma(i,j)) 
        else
          p_gamma => gamma(i,j)

          R(i,j,:) = Characteristic (p_Q, p_Q_0,      gamma=p_gamma)
        endif
      endif

    enddo
  enddo

  Nullify (p_Q, p_Q_0, p_N, p_gamma)

 100 format (1x,72('*')/  &
             1x,'FATAL : Characteristic_Variable_2D    '/&
             1x,'        Incompatible dimensions of Q and the normal vector'/&
             1x,'        Normal dimension     = ',i1,' should be equal to' /&
             1x,'        Q      dimension - 2 = ',i1 / &
             1x,72('*')/)

  END Subroutine Characteristic_Variable_2D
!
!=========================================================================
!
  Subroutine Characteristic_Variable_3D (Q, Q_0, R, Normal, gamma)

  REALTYPE, dimension(:,:,:,:)          , TARGET :: Q, Q_0, R
  REALTYPE, dimension(:,:,:,:), OPTIONAL, TARGET :: Normal
  REALTYPE, dimension(:,:,:)  , OPTIONAL, TARGET :: gamma

  integer                         :: i, j, k
  REALTYPE,               POINTER :: p_gamma
  REALTYPE, dimension(:), POINTER :: p_Q, p_Q_0, p_N

  if (PRESENT(Normal)) then
    if (SIZE(Normal, DIM=4) /= SIZE(Q, DIM=4)-2) then
      write (6 ,100) SIZE(Normal, DIM=4), SIZE(Q, DIM=4)-2
      write (99,100) SIZE(Normal, DIM=4), SIZE(Q, DIM=4)-2
      STOP
    endif
  endif

  Nullify (p_Q, p_Q_0, p_N, p_gamma)

  do k = 1, SIZE(Q, DIM=3)
    do j = 1, SIZE(Q, DIM=2)
      do i = 1, SIZE(Q, DIM=1)

        p_Q   => Q  (i,j,k,:)
        p_Q_0 => Q_0(i,j,k,:)

        if (PRESENT(Normal)) then
          p_N     => Normal(i,j,k,:)
      
          if (PRESENT(gamma)) then
            R(i,j,k,:) = Characteristic (p_Q, p_Q_0, p_N, gamma(i,j,k))
          else
            p_gamma => gamma(i,j,k)

            R(i,j,k,:) = Characteristic (p_Q, p_Q_0, p_N, p_gamma)
          endif

        else
          if (PRESENT(gamma)) then
            R(i,j,k,:) = Characteristic (p_Q, p_Q_0,      gamma=gamma(i,j,k)) 
          else
            p_gamma => gamma(i,j,k)

            R(i,j,k,:) = Characteristic (p_Q, p_Q_0,      gamma=p_gamma)
          endif
        endif

      enddo
    enddo
  enddo

  Nullify (p_Q, p_Q_0, p_N, p_gamma)

 100 format (1x,72('*')/  &
             1x,'FATAL : Characteristic_Variable_3D    '/&
             1x,'        Incompatible dimensions of Q and the normal vector'/&
             1x,'        Normal dimension     = ',i1,' should be equal to' /&
             1x,'        Q      dimension - 2 = ',i1 / &
             1x,72('*')/)

  END Subroutine Characteristic_Variable_3D
#endif
!***************************************************************************!
!                                                                           !
!***************************************************************************!

END MODULE Characteristic_Decomposition
