#include "pseudopack.h"


MODULE Grid_Generation
  USE Grid_Auxillary_Subroutine

  USE   Fourier_Grid_Point
  USE Chebyshev_Grid_Point
  USE  Legendre_Grid_Point

  USE   Fourier_Mapping
  USE Chebyshev_Mapping

implicit NONE

PRIVATE

PUBLIC  ::   Grid_And_Mapping
PUBLIC  :: F_Grid_And_Mapping
PUBLIC  :: C_Grid_And_Mapping
PUBLIC  :: L_Grid_And_Mapping

CONTAINS

! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! FUNCTION NAME: Grid_And_Mapping
! DESCRIPTION  : For All Methods (Fourier, Chebyshev and Legendre)
!                Computes various things associated with Mapping including
!                the mapped grid points x and its Transformation Metrics.
!
!  if the first elements of x is coded with -99, then regardless of the
!    SYMMETRY setting, the collocations points x are returned in
!    full domain rather than half domain even if SYMMETRY=1.
!
!  N+1     = Total  number of points regardless Symmetry (Chebyshev/Legendre)
!  N       = Total  number of points regardless Symmetry (Fourier)
!  N0+1    = Actual number of points taking into account of symmetry
! 
!  Point_Type = 1 Lobatto Point
!             = 2 Radau   Point
!             = 3 Gauss   Point
!
!  x_Left  = Left/Bottom  Physical Domain Points
!  x_Right = Right/Top    Physical Domain Points
!  x       = Computational Grid points (Input)
!            Mapped (if Mapping used) Phsyical Grid points (Output)
!  A       = coefficients of the Transformation Metrics
!            that make up the Mapped k-th Derivative
!
!  if Map  < 0, use user defined Mapping functions.
!  if Map  = 0, No mapping is used
!  if Map  > 0, use built in Mapping functions.
!  if  1 < ABS(Map) <= 10, it is a mapping for FINITE        domain.
!  if 10 < ABS(Map) <= 20, it is a mapping for SEMI-INFINITE domain.
!  if 20 < ABS(Map) <= 30, it is a mapping for INFINITE      domain.
!  if Linear Mapping is applied, then Map is set to be 99.
!
!  For Chebyshev and Legendre Method (Method=1,2) and Tal-Ezer Map (Map=1) ONLY
!    Manual = 0, the value of alpha is computed from epsilon
!    Manual = 1, the value of alpha is input
!
!  For (Fourier Method)  Manual = Map_F from the calling statments.
!    Manual  = 0, [x_Left, x_Right] are defaulted to [0, 2 Pi]
!            =-1, Map [0,2 Pi] to [x_Left, x_Right] in Compuational Space.
!            = 1, Map [0,2 Pi] to [x_Left, x_Right] in Physical Space.
!
!    In Computatonal Space, (x_Left-x_Right = 360 degree = 2*PI).
!      In the other words, it is a simple translation.
!      As  Input, x_Left and x_Right are in degrees.
!      As Output, x_Left and x_Right are in Radians.
!
!    In Physical     Space, [0, 2 Pi] <-> [x_Left, x_Right].
!      x_Left and x_Right can be any real number.
!
!    For Symmetry case, only the first half of x are returned.
!       x(0) = x_Left, x(N/2)=(x_Left+x_Right)/2
!
!    if x(0) = -99, then full range of x are returned regardless of
!      the setting of the parameter Symmetry.
!
! ----------------------------------------------------------------------
  Subroutine Grid_And_Mapping (Method, Point_Type,              &
                               N, N0, Map, Manual, alpha, beta, &
                               x_Left, x_Right, x, A)

  integer  :: Method, Point_Type, N, N0, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right
  integer  :: N2
  logical  :: Full_x = .TRUE.

  REALTYPE, dimension(:)           :: x
  REALTYPE, dimension(*), OPTIONAL :: A

#ifdef DEBUG
  ASSERT(N>0, "Grid_And_Mapping:N<=0")

  if (Method == 0) then
    if ((Point_Type == 1) .AND. (MOD(N,2) == 1)) then
      write (6 ,*) 'Grid_And_Mapping : Fourier Lobatto => N must be Even' 
      write (99,*) 'Grid_And_Mapping : Fourier Lobatto => N must be Even' 
      STOP
    endif
    if ((Point_Type == 2) .AND. (MOD(N,2) == 0)) then
      write (6 ,*) 'Grid_And_Mapping : Fourier Radau => N must be Odd' 
      write (99,*) 'Grid_And_Mapping : Fourier Radau => N must be Odd' 
      STOP
    endif
  endif
#endif

#if 0
  write (6,*) 'Full_x', Method, x(1)
  Full_x = .FALSE. 
  write (6,*) 'Full_x', Full_x
  if (x(1) == -99) Full_x = .TRUE. 
  write (6,*) 'Full_x', Full_x
! Full_x = (x(1) == -99)
  write (6,*) 'Method', Method
#else
  Full_x = (x(1) == -99)
#endif

  SELECT CASE (Method)
    CASE (0)
      SELECT CASE (Point_Type)
        CASE (1)
          if (Full_x .OR. (N-1 == N0)) then
            call F_GL_Point_Full (N, x) ; N2 = N-1
          else
            call F_GL_Point_Half (N, x) ; N2 = N0
          endif

        CASE (2)
            call F_GR_Point_Full (N, x) ; N2 = N-1

        CASE (3)
          if (Full_x .OR. (N-1 == N0)) then
            call F_GL_Point_Full (N, x) ; N2 = N-1
          else
            call F_GL_Point_Half (N, x) ; N2 = N0
          endif

      END SELECT

      call F_Grid_And_Mapping (N,N0,N2, Map, Manual, alpha, beta, &
                                        x_Left,x_Right, x, A)

   CASE (1)
      SELECT CASE (Point_Type)
        CASE (1)
          if (Full_x .OR. (N == N0)) then
            call C_GL_Point_Full (N, x) ; N2 = N
          else
            call C_GL_Point_Half (N, x) ; N2 = N0
          endif

        CASE (2)
            call C_GR_Point_Full (N, x) ; N2 = N

        CASE (3)
          if (Full_x .OR. (N == N0)) then
            call C_Ga_Point_Full (N, x) ; N2 = N
          else
            call C_Ga_Point_Half (N, x) ; N2 = N0
          endif

      END SELECT

      call C_Grid_And_Mapping (N,N0,N2, Map, Manual, alpha, beta, &
                                        x_Left,x_Right, x, A)

   CASE (2)
      SELECT CASE (Point_Type)
        CASE (1)
          if (Full_x .OR. (N == N0)) then
            call L_GL_Point_Full (N, x) ; N2 = N
          else
            call L_GL_Point_Half (N, x) ; N2 = N0
          endif

        CASE (2)
            call L_GR_Point_Full (N, x) ; N2 = N

        CASE (3)
          if (Full_x .OR. (N == N0)) then
            call L_Ga_Point_Full (N, x) ; N2 = N
          else
            call L_Ga_Point_Half (N, x) ; N2 = N0
          endif

      END SELECT

      call L_Grid_And_Mapping (N,N0,N2, Map, Manual, alpha, beta, &
                                        x_Left,x_Right, x, A)

    END SELECT

  END Subroutine Grid_And_Mapping 

! ----------------------------------------------------------------------
! FUNCTION NAME: F_Grid_And_Mapping
! DESCRIPTION  : For Fourier Methods.
!                Computes various things associated with Mapping including
!                the mapped grid points x and its Transformation Metrics.
!
!  N0+1    = Acutual number of points taking in account of symmetry
!  x       = Phsyical Grid points (Mapped if possible)
!  A       = coefficients of the Different unMapped Derivatives in the
!            linear combinations that make up the Mapped k-th Derivative
!
!  Map     = 0 No Mapping
!          = 1 Do the mapping with parameters alpha and beta.
!
!  Map_F   =-1, Map [0,2 Pi] to [x_Left, x_Right] in Compuational Space.
!          = 0, Defaulted to [0,2 Pi], x_Left and x_Right are not used.
!          = 1, Map [0,2 Pi] to [x_Left, x_Right] in Physical Space.
!
!    In Computatonal Space, (x_Left-x_Right = 360 degree = 2*PI).
!      In the other words, it is a simple translation.
!      As  Input, x_Left and x_Right are in degrees.
!      As Output, x_Left and x_Right are in Radians.
!
!    In Physical     Space, [0, 2 Pi] <-> [x_Left, x_Right].
!      x_Left and x_Right can be any real number.
!
! ----------------------------------------------------------------------
  Subroutine F_Grid_And_Mapping (N,N0,N2, Map, Map_F,  alpha, beta, &
                                          x_Left,x_Right, x, A)

  integer  :: N, N0, N2, Map, Map_F
  REALTYPE :: alpha, beta, x_Left, x_Right
  REALTYPE :: Slope, x_Temp

  REALTYPE, dimension(0:N2)             :: x
  REALTYPE, dimension(0:N0,*), OPTIONAL :: A

  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2, &
                               D_Map_D_x_3, D_Map_D_x_4

#ifdef DEBUG
  logical  :: Bug

  ASSERT(N>0, "F_Grid_And_Map:N<=0")
  Bug =  (Map_F ==  0) .OR.  (Map_F == 1) .OR.  &
        ((Map_F == -1) .AND. (ABS(int(x_Right-x_Left)) == 360))
  ASSERT(Bug, "F_Grid_And_Map:x_Right-x_Left<>360")
#endif

  SELECT CASE (Map_F)
    CASE (-1)
        beta  = beta   *D_2_R 
      x_Left  = x_Left *D_2_R 
      x_Right = x_Right*D_2_R

    CASE (0)   
        beta  = beta   *D_2_R 
      x_Left  = ZERO
      x_Right = TWOPI
 
  END SELECT

  call F_Mapping (N, N0, N2, Map, Map_F,  alpha, beta, x_Left, x_Right, &
                  x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)
 
  if (PRESENT(A)) &
    call Setup_Metric_A (N0, MAX_DERIVATIVE, D_Map_D_x_1, D_Map_D_x_2,  &
                                             D_Map_D_x_3, D_Map_D_x_4, A)

  Slope = TWOPI/(x_Right-x_Left)

  if (ABS(Slope) /= ONE) then
    if (PRESENT(A)) call Apply_Linear_Map (N0, MAX_DERIVATIVE, Slope, A)
   
    if (Map == 0) Map = 99
  else
    if (Slope == -ONE) then
      x_Temp  = x_Right ; x_Right = x_Left ; x_Left  = x_Temp
    endif
  endif

  call Map_Computational_2_Physical (N2+1,   ZERO,  TWOPI, x, &
                                           x_Left,x_Right, x)

  END Subroutine F_Grid_And_Mapping 

! ----------------------------------------------------------------------
! FUNCTION NAME: C_Grid_And_Mapping
! DESCRIPTION  : For Chebyshev AND Legendre Method.
!                Computes various things associated with Mapping including
!                the mapped grid points x AND its Transformation Metrics.
!
! See Grid_And_Mapping for description of parameters.
!
! ----------------------------------------------------------------------
  Subroutine C_Grid_And_Mapping (N,N0,N2, Map, Manual, alpha, beta, &
                                          x_Left,x_Right, x, A)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right
  REALTYPE :: Slope

  REALTYPE, dimension(0:N2)             :: x
  REALTYPE, dimension(0:N0,*), OPTIONAL :: A

  REALTYPE, dimension(0:N0) :: D_Map_D_x_1, D_Map_D_x_2, &
                               D_Map_D_x_3, D_Map_D_x_4

#ifdef DEBUG
  ASSERT(N>0, "Grid_And_Map:N<=0")
#endif
 
  call C_Mapping (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
                  x, D_Map_D_x_1, D_Map_D_x_2, D_Map_D_x_3, D_Map_D_x_4)
 
  if (PRESENT(A)) &
    call Setup_Metric_A (N0, MAX_DERIVATIVE, D_Map_D_x_1, D_Map_D_x_2,  &
                                             D_Map_D_x_3, D_Map_D_x_4, A)


  SELECT CASE (ABS(Map))
    CASE DEFAULT
      Slope = TWO/(x_Right-x_Left)

      if (PRESENT(A)) call Apply_Linear_Map (N0, MAX_DERIVATIVE, Slope, A)

      if ((Slope /= ONE) .AND. (Map == 0)) Map = 99

      call Map_Computational_2_Physical (N2+1,   -ONE,    ONE, x, &
                                               x_Left,x_Right, x)

    CASE (11:20)
      x_Left  =  x(N) ; x_Right = x(0)

    CASE (21:30)
      x_Left  = -x(0) ; x_Right = x(0)

  END SELECT

#if 0
  if ((ABS(Map) >= 11) .AND. (ABS(Map) <= 20)) then
    x_Left  =  x(N) ; x_Right = x(0) ; RETURN
  endif

  if ((ABS(Map) >= 21) .AND. (ABS(Map) <= 30)) then
    x_Left  = -x(0) ; x_Right = x(0) ; RETURN
  endif

  Slope = TWO/(x_Right-x_Left)

  if (PRESENT(A)) call Apply_Linear_Map (N0, MAX_DERIVATIVE, Slope, A)

  if ((Slope /= ONE) .AND. (Map == 0)) Map = 99

#if defined (TESTING)
  if (ABS(Slope) /= ONE) then
    call Apply_Linear_Map (N0, MAX_DERIVATIVE, Slope, A)

    if (Map == 0) Map = 99
  else
    if (Slope == -ONE) then
      x_Temp  = x_Right ; x_Right = x_Left ; x_Left  = x_Temp
    endif

    if (x_Right == ONE) RETURN
  endif
#endif

  call Map_Computational_2_Physical (N2+1,   -ONE,    ONE, x, &
                                           x_Left,x_Right, x)
#endif

  END Subroutine C_Grid_And_Mapping 
!
! ----------------------------------------------------------------------
!
  Subroutine L_Grid_And_Mapping (N,N0,N2, Map, Manual, alpha, beta, &
                                          x_Left,x_Right, x, A)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2)             :: x
  REALTYPE, dimension(0:N0,*), OPTIONAL :: A

  if (PRESENT(A)) then
    call C_Grid_And_Mapping (N,N0,N2, Map, Manual, alpha, beta, &
                                      x_Left,x_Right, x, A)
  else
    call C_Grid_And_Mapping (N,N0,N2, Map, Manual, alpha, beta, &
                                      x_Left,x_Right, x)
  endif

  END Subroutine L_Grid_And_Mapping 

END MODULE Grid_Generation
