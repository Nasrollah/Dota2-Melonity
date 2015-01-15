#include "pseudopack.h"

!***************************************************************************
! MODULE     : Interpolation
!
! SUBROUTINE : PS_Get_Operator_I_N_Size (N_Pnts, M, I_N_Size)
!              PS_Setup_Operator_I_N    (I_N, y, M, Grid_Index, &
!                                             Property, Domain, Mapping, Error)
!              PS_Setup_Operator_I_N    (I_N, y, M, N_Pnts,     &
!                                             Property, Domain, Mapping, Error)
!              PS_Interpolate           (I_N, f, g, Even_Or_Odd, Error)
!              PS_Interpolation_Weight_Output (I_N, lid)
!
! Function   : Perform Interpolation Setup and Operations.
!
!              Given the collocation grid and its functional value pair
!                    (x_i, f_i), i = 1,...,N_Pnts
!              and location y_j in [x(1),x(N_Pnts)], j=1,..,M, 
!              the subroutine compute the value g(y_j), j=1,..M. 
!           
!              That is,  g(y_j) = SUM (i=1,N_Pnts) I_N(i,j)*f(i)
!
!              
! DEBUG             Test the integrity of the input data 
!
! SINGLE_PRECISION  Single precision code instead of double precision by
!                   default
!
! N_Pnts Integer                                            [Input]
!          Total number of Collocation points x(i), i=0,..,N_Pnts-1
!            used for the Interpolation Weight Setup
!            before taking symmetry into account.
!
! M      Integer                                            [Input]
!          Extend of the Interpolated stencils y(j), j=1,..,M
!
! I_N_Size Integer                                          [Output]
!          Size of the Interpolation Operator I_N
!
! y      REALTYPE                                           [Input]
!          Location where the function is interpolated at.
!
! I_N    Integer of Size I_N_Size                           [Input]
!          Interpolation Operator
!
! Even_Or_Odd Integer
!          Even or Odd property of the interpolation function
!
! Property  TYPE(PS_Property)                               [Input, OPTIONAL]
!             Property of the Interpolation scheme
!               [Method, Point_Type, Algorithm, Symmetry]
!
! Domain    TYPE(PS_Domain)                                 [Input, OPTIONAL]
!             Domain of the Interpolation scheme  
!               [x0, x1, Map_F] for the Polynomial basis
!               [        Map_F] for the Trigonometric basis  [DEFAULT]
!
! Mapping   TYPE(PS_Mapping)                                [Input, OPTIONAL]
!             Mapping of the Interpolation scheme
!               [Map, Manual, alpha, beta]
!
! Grid_Index TYPE(PS_Grid_Index)                             [Input]
!             Grid Index of the Interpolation scheme
!               [N=N_Pnts, LDY]
!
! Error  Logical                                             [Output, OPTIONAL]
!             Error Flag
!
! lid    Integer                                             [Input, OPTIONAL]
!          Error Output Unit
!
! f      REALTYPE of Size N1                                 [Input]
!          Interpolation Function
!          Use subroutine
!            PS_Get_Dimension (N1, Property, Grid_Index)
!          to compute N1
!          N1 = N_Pnts for non-symmetric function
!          N1 = N_Pnts/2+1 for symmetric function
! 
! g      REALTYPE of Size M                                   [Output]
!          Interpolated Function
! 
! lid    Integer
!
! AUTHOR      :: Wai Sun Don, PhD
! DATE        :: 3-14-1999
! VERSION     :: 2.0
!***************************************************************************

MODULE Interpolation
  USE PseudoPack_Type

  USE PseudoPack_Grid   , ONLY: PS_Collocation_Grid
  USE PseudoPack_Utility, ONLY: PS_Get_Dimension

  USE  Chebyshev_Mapping
  USE    Fourier_Mapping

  USE Lagrangian_Weights

implicit NONE

INTERFACE PS_Setup_Operator_I_N
  MODULE PROCEDURE PS_Setup_Operator_I_N_0
  MODULE PROCEDURE PS_Setup_Operator_I_N_1
  MODULE PROCEDURE PS_Setup_Operator_I_N_2
  MODULE PROCEDURE PS_Setup_Operator_I_N_3
END INTERFACE

PRIVATE

PUBLIC  :: PS_Property
PUBLIC  :: PS_Grid_Index
PUBLIC  :: PS_Domain
PUBLIC  :: PS_Mapping

PUBLIC  :: PS_Get_Operator_I_N_Size
PUBLIC  :: PS_Setup_Operator_I_N
PUBLIC  :: PS_Interpolation_Weight_Output

PUBLIC  :: PS_Lagrangian_Weight
PUBLIC  :: PS_Lagrangian_Weight_Output

PUBLIC  :: PS_Interpolate

CONTAINS

  Subroutine PS_Get_Operator_I_N_Size (N_Pnts, M, I_N_Size)

  integer  :: N_Pnts, M, I_N_Size

  I_N_Size = HEADER_SIZE + N_Pnts + M + M + N_Pnts*M

  END Subroutine PS_Get_Operator_I_N_Size
!
! ---------------------------------------------------------------------------
!
  Subroutine PS_Interpolate (I_N, f, g, Even_Or_Odd, Error)
                               
  REALTYPE, dimension(*) :: I_N
  REALTYPE, dimension(*) :: f, g

  integer ,     OPTIONAL :: Even_Or_Odd
  logical ,     OPTIONAL :: Error

  integer  :: N_Pnts, M, N, i_n_s, i, j, N1, IEO
  REALTYPE :: g_j

  if (PRESENT(Error)) Error = .FALSE.

  N_Pnts = I_N(_I_N_PNTS)
  M      = I_N(_I_N_VECTORS)

  N1 = N_Pnts-1

  call PS_Get_Dimension (I_N, N)

  IEO = 0 ; if (PRESENT(Even_Or_Odd)) IEO = Even_Or_Odd

  if (N == N_Pnts) IEO = 0

#if defined (DEBUG)
  if (PRESENT(Even_Or_Odd)) then
    if (ABS(IEO) >  1) then 
      write (6,100) IEO ; write (99,100) IEO 

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif

 100 format (1x,60('*')                                             / &
             1x,'FATAL ERROR : PS_Interpolate '                     / &
             15x,'Incorrect Specification of Even_Or_Odd Argument'  / &
             15x,'              Even_Or_Odd = ',i5                  / &
             1x,60('*'))
#endif

  i_n_s = HEADER_SIZE + N_Pnts + M + M + 1

  call Interpolate_f_to_g (f, g, I_N(i_n_s))

  CONTAINS

    Subroutine Interpolate_f_To_g (f, g, W)

    REALTYPE, dimension(0:N-1)  :: f
    REALTYPE, dimension(M)      :: g
    REALTYPE, dimension(0:N1,M) :: W

!$OMP PARALLEL DO
    do j = 1,M

      g_j = SUM(W(0:N-1,j)*f(0:N-1))

      if (IEO /= 0) then
        do i = N,N1
          g_j = g_j + IEO*W(i,j)*f(N1-i)
        enddo
      endif

      g(j) = g_j
    enddo
!$OMP END PARALLEL DO

    END Subroutine Interpolate_f_To_g 

  END Subroutine PS_Interpolate
!
! ---------------------------------------------------------------------------
!
  Subroutine PS_Setup_Operator_I_N_0 (I_N, y, M, Grid_Index, &
                                      Property, Domain, Mapping, Error)

  integer                :: M
  REALTYPE, dimension(*) :: I_N  
  REALTYPE, dimension(M) :: y

  TYPE (PS_Grid_Index)           :: Grid_Index
  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping

  logical             , OPTIONAL :: Error

  integer  :: N_Pnts
  integer  :: Method, Point_Type, Symmetry, Map, Manual, Map_F
  REALTYPE :: x0, x1, alpha, beta

  call Set_Default_Data (Property, Method, Point_Type, Symmetry,   &
                         Domain  , x0, x1, Map_F,                  &
                         Mapping , Map, Manual, alpha, beta)

  N_Pnts = Grid_Index%N

  call Set_I_N_Parameters (I_N, Method, Point_Type, Symmetry, N_Pnts, M, &
                                x0, x1, Map_F, Map, Manual, alpha, beta)

  SELECT CASE (Method)
    CASE (0)
      call    Fourier_Weight (I_N, N_Pnts, y, M, Error)

    CASE (1,2)
      call Polynomial_Weight (I_N, N_Pnts, y, M, Error)

  END SELECT

  END Subroutine PS_Setup_Operator_I_N_0
!
! ---------------------------------------------------------------------------
!
  Subroutine PS_Setup_Operator_I_N_1 (I_N, y, M, N_Pnts, &
                                       Property, Domain, Mapping, Error)

  integer                :: M
  REALTYPE, dimension(*) :: I_N
  REALTYPE, dimension(M) :: y

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping

  logical             , OPTIONAL :: Error

  integer  :: N_Pnts
  integer  :: Method, Point_Type, Symmetry, Map_F, Map, Manual
  REALTYPE :: x0, x1, alpha, beta

  call Set_Default_Data (Property, Method, Point_Type, Symmetry,   &
                         Domain  , x0, x1, Map_F ,                 &
                         Mapping , Map, Manual, alpha, beta)

  call Set_I_N_Parameters (I_N, Method, Point_Type, Symmetry, N_Pnts, M, &
                                x0, x1, Map_F, Map, Manual, alpha, beta)

  SELECT CASE (Method)
    CASE (0)
      call    Fourier_Weight (I_N, N_Pnts, y, M, Error)

    CASE (1,2)
      call Polynomial_Weight (I_N, N_Pnts, y, M, Error)

  END SELECT

  END Subroutine PS_Setup_Operator_I_N_1
!
! ---------------------------------------------------------------------------
!
  Subroutine PS_Setup_Operator_I_N_2 (I_N, y, z, M, Grid_Index, &
                                       Property, Domain, Mapping, Error)

  integer                :: M
  REALTYPE, dimension(*) :: I_N
  REALTYPE, dimension(M) :: y, z

  TYPE (PS_Grid_Index)           :: Grid_Index
  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping

  logical             , OPTIONAL :: Error

  integer  :: N_Pnts
  integer  :: Method, Point_Type, Symmetry, Map_F, Map, Manual
  REALTYPE :: x0, x1, alpha, beta
  integer  :: j

  call Set_Default_Data (Property, Method, Point_Type, Symmetry,   &
                         Domain  , x0, x1, Map_F ,                 &
                         Mapping , Map, Manual, alpha, beta)

  N_Pnts = Grid_Index%N

  call Set_I_N_Parameters (I_N, Method, Point_Type, Symmetry, N_Pnts, M, &
                                x0, x1, Map_F, Map, Manual, alpha, beta)

  SELECT CASE (Method)
    CASE (0)
!$OMP PARALLEL DO
      do j = 1,M
        z(j) = F_Inverse_Mapping (I_N, y(j))
      enddo
!$OMP END PARALLEL DO

    CASE (1,2)
!$OMP PARALLEL DO
      do j = 1,M
        z(j) = C_Inverse_Mapping (I_N, y(j))
      enddo
!$OMP END PARALLEL DO

  END SELECT

  END Subroutine PS_Setup_Operator_I_N_2
!
! ---------------------------------------------------------------------------
!
  Subroutine PS_Setup_Operator_I_N_3 (I_N, y, z, M, N_Pnts, &
                                       Property, Domain, Mapping, Error)

  integer                :: M
  REALTYPE, dimension(*) :: I_N
  REALTYPE, dimension(M) :: y, z

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping

  logical             , OPTIONAL :: Error

  integer  :: N_Pnts
  integer  :: Method, Point_Type, Symmetry, Map_F, Map, Manual
  REALTYPE :: x0, x1, alpha, beta
  integer  :: j

  call Set_Default_Data (Property, Method, Point_Type, Symmetry,   &
                         Domain  , x0, x1, Map_F ,                 &
                         Mapping , Map, Manual, alpha, beta)

  call Set_I_N_Parameters (I_N, Method, Point_Type, Symmetry, N_Pnts, M, &
                                x0, x1, Map_F, Map, Manual, alpha, beta)

  SELECT CASE (Method)
    CASE (0)
!$OMP PARALLEL DO
      do j = 1,M
        z(j) = F_Inverse_Mapping (I_N, y(j))
      enddo
!$OMP END PARALLEL DO

    CASE (1,2)
!$OMP PARALLEL DO
      do j = 1,M
        z(j) = C_Inverse_Mapping (I_N, y(j))
      enddo
!$OMP END PARALLEL DO

  END SELECT

  END Subroutine PS_Setup_Operator_I_N_3
!
! -----------------------------------------------------------------------------
!
  Subroutine Set_Default_Data (Property, Method, Point_Type, Symmetry,   &
                               Domain  , x0, x1, Map_F,                  &
                               Mapping , Map, Manual, alpha, beta)

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping

  integer  :: Method, Point_Type, Symmetry, Map_F, Map, Manual
  REALTYPE :: x0, x1, alpha, beta

  Logical  :: L_Property, L_Domain, L_Mapping

  L_Property   = PRESENT(Property  )
  L_Domain     = PRESENT(Domain    )
  L_Mapping    = PRESENT(Mapping   )

  if (L_Property) then
    Method         = Property%Method
    Point_Type     = Property%Point_Type
    Symmetry       = Property%Symmetry
  else
#if defined (DEFAULT_FOURIER)
    Method = 0 ; Point_Type = 1 ; Symmetry = 0
#else
    Method = 1 ; Point_Type = 1 ; Symmetry = 0
#endif
  endif

  if (L_Domain) then
    x0             = Domain%x0
    x1             = Domain%x1

    Map_F          = Domain%Map_F
  else
    SELECT CASE (Method)
      CASE (0)
        x0 = ZERO ; x1 = TWO*PI ; Map_F = 0

      CASE (1,2)
        x0 = -ONE ; x1 = ONE    ; Map_F = 0
    END SELECT
  endif
  
  if (L_Mapping) then
    Map            = Mapping%Map
    Manual         = Mapping%Manual
    alpha          = Mapping%alpha
    beta           = Mapping%beta
  else
    Map = 0 ; Manual = 0 ; alpha = ZERO ; beta = ZERO 
  endif
    
  if (Method == 0) Manual = Map_F

  END Subroutine Set_Default_Data 
!
! -----------------------------------------------------------------------------
!
  Subroutine Set_I_N_Parameters (I_N, Method, Point_Type, Symmetry, N_Pnts, M, &
                                      x0, x1, Map_F, Map, Manual, alpha, beta)

  REALTYPE, dimension(*) :: I_N

  integer  :: N_Pnts, M
  integer  :: Method, Point_Type, Symmetry, Map_F, Map, Manual
  REALTYPE :: x0, x1, alpha, beta

  integer  :: i, N_Threads = 1
#if defined (PARALLEL_OPENMP)
  integer  :: OMP_GET_NUM_THREADS

!$OMP PARALLEL DO
  do i = 1,2
    N_Threads = OMP_GET_NUM_THREADS()
  enddo
!$OMP END PARALLEL DO
#endif

  I_N(_I_CASE)           = 2
  I_N(_I_METHOD)         = Method
  I_N(_I_POINT_TYPE)     = Point_Type
  I_N(_I_SYMMETRY)       = Symmetry
  I_N(_I_N_PNTS)         = N_Pnts
  I_N(_I_X_LEFT)         = x0
  I_N(_I_X_RIGHT)        = x1
  I_N(_I_MAP_F)          = Map_F
  I_N(_I_MAP)            = Map
  I_N(_I_MANUAL)         = Manual
  I_N(_I_ALPHA)          = alpha
  I_N(_I_BETA)           = beta
  I_N(_I_N_VECTORS)      = M
  I_N(_I_NUMBER_THREADS) = N_Threads

  END Subroutine Set_I_N_Parameters
!
! -----------------------------------------------------------------------------
!
  Subroutine    Fourier_Weight (I_N, N_Pnts, y, M, Error)

  integer                         :: N_Pnts, M
  REALTYPE, dimension(*)          :: I_N
  REALTYPE, dimension(0:N_Pnts-1) :: x
  REALTYPE, dimension(M)          :: y, z

  logical ,              OPTIONAL :: Error

  integer  :: N, j, i_n_s
  integer  :: Method, Point_Type, Map_F
  REALTYPE :: x0, x1

#if defined (DEBUG)
  x0         = OPERATOR_X_LEFT(I_N)
  x1         = OPERATOR_X_RIGHT(I_N)

  call Data_Integrity (x0, x1, M, y, Error)
#endif

  Method     = OPERATOR_METHOD(I_N)
  Point_Type = OPERATOR_POINT_TYPE(I_N)

  x = ZERO ; x0 = ZERO ; x1 = TWOPI ; Map_F = 0
  call PS_Collocation_Grid (x, N_Pnts, x0, x1, Map_F,  &
                               Method, Point_Type, Full_x=.TRUE.)

!$OMP PARALLEL DO
  do j = 1,M
    z(j) = F_Inverse_Mapping (I_N, y(j))
  enddo
!$OMP END PARALLEL DO

  N = N_Pnts - 1

  i_n_s = HEADER_SIZE + 1 ; I_N(i_n_s:i_n_s+N  ) = x
  i_n_s = i_n_s + N_Pnts  ; I_N(i_n_s:i_n_s+M-1) = y
  i_n_s = i_n_s + M       ; I_N(i_n_s:i_n_s+M-1) = z
  i_n_s = i_n_s + M       ; call Get_Fourier_Weight (I_N(i_n_s))

  CONTAINS

    Subroutine Get_Fourier_Weight (W)

    REALTYPE, dimension(0:N,M) :: W

    SELECT CASE (Point_Type)
      CASE (1,3)
!$OMP PARALLEL DO
        do j = 1,M
          call Compute_Fourier_Lobatto_Weight (z(j), N+1, x, W(0,j))
        enddo
!$OMP END PARALLEL DO

      CASE (2)
!$OMP PARALLEL DO
        do j = 1,M
          call Compute_Fourier_Radau_Weight   (z(j), N+1, x, W(0,j))
        enddo
!$OMP END PARALLEL DO

    END SELECT

    END Subroutine Get_Fourier_Weight

  END Subroutine    Fourier_Weight
!
! -------------------------------------------------------------------
!
  Subroutine Compute_Fourier_Lobatto_Weight (z, N, x, W)

  integer  :: N, i
  REALTYPE :: z, Angle

  REALTYPE, dimension(0:N-1) :: x, W
  
!$OMP PARALLEL DO
  do i = 0,N-1
    Angle = HALF*(z-x(i))

    if (z == x(i)) then
      W(i) = ONE
    else
      W(i) = SIN(N*Angle)/TAN(Angle)/N
    endif
  enddo
!$OMP END PARALLEL DO

  END Subroutine Compute_Fourier_Lobatto_Weight 
!
! -------------------------------------------------------------------
!
  Subroutine Compute_Fourier_Radau_Weight (z, N, x, W)

  integer  :: N, i
  REALTYPE :: z, Angle

  REALTYPE, dimension(0:N-1) :: x, W
  
!$OMP PARALLEL DO
  do i = 0,N-1
    Angle = HALF*(z-x(i))

    if (z == x(i)) then
      W(i) = ONE
    else
      W(i) = SIN(N*Angle)/SIN(Angle)/N
    endif
  enddo
!$OMP END PARALLEL DO

  END Subroutine Compute_Fourier_Radau_Weight 
!
! -----------------------------------------------------------------------------
!
  Subroutine Polynomial_Weight (I_N, N_Pnts, y, M, Error)

  integer                         :: N_Pnts, M
  REALTYPE, dimension(*)          :: I_N
  REALTYPE, dimension(0:N_Pnts-1) :: x
  REALTYPE, dimension(M)          :: y, z

  logical ,              OPTIONAL :: Error

  integer  :: N, j, i_n_s, Order
  integer  :: Method, Point_Type, Map_F
  REALTYPE :: x0, x1

#if defined (DEBUG)
  x0         = OPERATOR_X_LEFT(I_N)
  x1         = OPERATOR_X_RIGHT(I_N)

  call Data_Integrity (x0, x1, M, y, Error)
#endif

  Method     = OPERATOR_METHOD(I_N)
  Point_Type = OPERATOR_POINT_TYPE(I_N)

  x = ZERO ; x0 = -ONE ; x1 = ONE ; Map_F = 0 ; Order = 0
  call PS_Collocation_Grid (x, N_Pnts, x0, x1, Map_F, &
                               Method, Point_Type, Full_x=.TRUE.)

!$OMP PARALLEL DO
  do j = 1,M
    z(j) = C_Inverse_Mapping (I_N, y(j))
  enddo
!$OMP END PARALLEL DO

  N = N_Pnts - 1

  i_n_s = HEADER_SIZE + 1 ; I_N(i_n_s:i_n_s+N  ) = x
  i_n_s = i_n_s + N_Pnts  ; I_N(i_n_s:i_n_s+M-1) = y
  i_n_s = i_n_s + M       ; I_N(i_n_s:i_n_s+M-1) = z
  i_n_s = i_n_s + M       ; call Get_Lagrangian_Weight (I_N(i_n_s))

  CONTAINS

    Subroutine Get_Lagrangian_Weight (W)

    integer  :: l
    REALTYPE, dimension(0:N,M)   :: W

!$OMP PARALLEL DO
    do j = 1,M
      l = j
      call PS_Lagrangian_Weight (z(l), Order, N, x, N+1, W(:,l:l), Error=Error)
    enddo
!$OMP END PARALLEL DO

    END Subroutine Get_Lagrangian_Weight

  END Subroutine Polynomial_Weight
!
! -----------------------------------------------------------------------------
!
  Subroutine PS_Interpolation_Weight_Output (I_N, lid)

  REALTYPE, dimension(*)   :: I_N

  integer  :: N_Pnts, N, M, j, i_n_s, Order, ldy, lid

  N_Pnts = OPERATOR_N_PNTS(I_N)
  M      = OPERATOR_N_VECTORS(I_N)

  N = N_Pnts-1 

  call Get_Lagrangian_Weight

  CONTAINS

    Subroutine Get_Lagrangian_Weight 

    REALTYPE, dimension(0:N) :: x
    REALTYPE, dimension(M)   :: y

    i_n_s = HEADER_SIZE + 1 ; x = I_N(i_n_s:i_n_s+N  )
    i_n_s = i_n_s + N_Pnts  ; y = I_N(i_n_s:i_n_s+M-1) 
    i_n_s = i_n_s + M       
    i_n_s = i_n_s + M       ; call Output_Lagrangian_Weight (x, y, I_N(i_n_s))

    END Subroutine Get_Lagrangian_Weight 

    Subroutine Output_Lagrangian_Weight (x, y, W)

    REALTYPE, dimension(0:N)     :: x
    REALTYPE, dimension(M)       :: y
    REALTYPE, dimension(0:N,M)   :: W

    Order = 0 ; ldy = SIZE(W, DIM=1)

!$OMP PARALLEL DO
    do j = 1,M
      call PS_Lagrangian_Weight_Output (y(j), Order, N, x, ldy, W(:,j:j), lid)
    enddo
!$OMP END PARALLEL DO

    END Subroutine Output_Lagrangian_Weight 

  END Subroutine PS_Interpolation_Weight_Output 

#if defined (DEBUG)
!
! -----------------------------------------------------------------------------
!
  Subroutine Data_Integrity (x0, x1, M, y, Error)

  integer  :: M, j
  REALTYPE :: x0, x1, z0, z1

  REALTYPE, dimension(M) :: y

  logical ,              OPTIONAL :: Error

  if (PRESENT(Error)) Error = .FALSE.

  z0 = x0 ; z1 = x1
  if (x0 >  x1) then ; z0 = x1 ; z1 = x0 ; endif

  do j = 1,M
    if (((z0-y(j)) >  EPSILON(ONE)) .OR. ((y(j)-z1) >  EPSILON(ONE))) then
      write (6,100) z0, z1, j, y(j) ; write (99,100) z0, z1, j, y(j)

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  enddo

 100 format (1x,60('*')                                             / &
             1x,'FATAL ERROR : Polynomial_Weight '                  / &
             15x,'Domain Bound Error : y is outside the bound of x' / &
             15x,'                     x = (',f10.5,',',f10.5,')'   / &
             15x,'                  j, y = (',i10  ,',',f10.5,')'   / &
             1x,60('*'))

  END Subroutine Data_Integrity
#endif

END MODULE Interpolation
