#include "pseudopack.h"


MODULE EOD_Auxillary_Subroutine
  USE Gauss_Quadrature

implicit NONE

PRIVATE

PUBLIC  :: Diagonal_Element_By_Sum_EOD
PUBLIC  :: Diagonal_Element_By_MXM_EOD
PUBLIC  :: Adjust_EOD_Matrix
PUBLIC  :: Slope_Check
PUBLIC  :: Apply_Mapping_EOD
PUBLIC  :: Smoothing_Matrix

CONTAINS

  Subroutine Diagonal_Element_By_Sum_EOD (N, M, x, D_Even, D_Odd)

  integer  :: N, M

  REALTYPE, dimension(0:N/2,0:N/2,M) :: D_Even, D_Odd
  REALTYPE, dimension(0:N/2)         :: x

  integer  :: k, j, N0, N2
  REALTYPE :: Sum

  N2  = N/2

  do k = 0,N2

    Sum = ZERO
    do j = 0,N2
      if (j /= k) Sum = Sum + D_Even(k,j,M)
    enddo

    D_Even(k,k,M) = -Sum
  enddo

                     N0 = N2
  if (MOD(N,2) == 0) N0 = N2-1

#if 1
  do k = 0,N0

    Sum = ZERO
    do j = 0,N2
      if (j /= k) Sum = Sum + D_Odd(k,j,M)*x(j)
    enddo

    if (M == 1) Sum = Sum-ONE

    D_Odd (k,k,M) = -Sum/x(k)
  enddo

  if (MOD(N,2) == 0) D_Odd (N2,N2,M) = ZERO
#else
  do k = 0,N0

    Sum = ZERO
    do j = 0,N0
      if (j /= k) Sum = Sum + D_Odd(k,j,M)*x(j)
    enddo

    D_Odd (k,k,M) = -Sum
  enddo

  if (M == 1) then
    do k = 0,N0
      D_Odd (k,k,M) = (ONE+D_Odd (k,k,M))/x(k)
    enddo
  else
    do k = 0,N0
      D_Odd (k,k,M) =  D_Odd (k,k,M) /x(k)
    enddo
  endif

  if (MOD(N,2) == 0) then
    k = N2
      D_Odd (k,k,M) = ZERO
  endif
#endif

  END Subroutine Diagonal_Element_By_Sum_EOD
!
!--------------------------------------------------------------------------
!
  Subroutine Diagonal_Element_By_MXM_EOD (N, M,    D_Even, D_Odd)

  integer  :: N, M

  REALTYPE, dimension(0:N/2,0:N/2,M) :: D_Even, D_Odd

  integer  :: j, k, N2
  REALTYPE :: D_OE, D_EO, D_kk, D_kl

  N2 = N/2

  do k = 0,N2

    D_kk = ZERO
    D_kl = ZERO
    do j = 0,N2
      D_OE = D_Odd (k,j,M-1)*D_Even(j,k,1)
      D_EO = D_Even(k,j,M-1)*D_Odd (j,k,1)

      D_kk = D_kk + (D_OE + D_EO)
      D_kl = D_kl + (D_OE - D_EO)
    enddo
    D_kk = HALF*D_kk
    D_kl = HALF*D_kl

    D_Even(k,k,M) = D_kk + D_kl
    D_Odd (k,k,M) = D_kk - D_kl
  enddo

  END Subroutine Diagonal_Element_By_MXM_EOD 
!
!--------------------------------------------------------------------------
!
  Subroutine Adjust_EOD_Matrix (N2, D_Even, D_Odd)

  integer  :: N2

  REALTYPE, dimension(0:N2,0:N2) :: D_Even, D_Odd

  D_Even = HALF*D_Even
  D_Odd  = HALF*D_Odd 

  END Subroutine Adjust_EOD_Matrix 
! ----------------------------------------------------------------------
! FUNCTION NAME: Slope_Check
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa  (Rio)
! DESCRIPTION  : Check to make sure that the transformation metric
!                aka slope is an even function.
!
!  N        = number of grid points.
!  Slope    = Mapping Derivative
!  EOD_Code = 0 (OK, do it here)
!           = 1 (Apply the mapping slope at the end of EO_Diff_)
!
! ----------------------------------------------------------------------
  Subroutine Slope_Check (N, Slope, EOD_Code)

  integer  ::  N, EOD_Code
  REALTYPE :: Odd_Sum, Odd_Sum_Tolerance

  REALTYPE, dimension(0:*) :: Slope

  integer  :: i, N2

  Odd_Sum_Tolerance = HUNDRED*EPSILON(ONE)

  N2 = N/2

  Odd_Sum = ZERO
  do i = 1,N2
    Odd_Sum = Odd_Sum + ABS(Slope(i)-Slope(N-i))
  enddo

                                    EOD_Code = 0
  if (Odd_Sum >  Odd_Sum_Tolerance) EOD_Code = 1

#if defined(DEBUG)
  if (EOD_Code == 1) then
    write (6 ,*) 'Caution : Non-Even Mapping Metric, Odd_Sum = ',Odd_Sum
    write (99,*) 'Caution : Non-Even Mapping Metric, Odd_Sum = ',Odd_Sum
  endif
#endif

  END Subroutine Slope_Check 

! ----------------------------------------------------------------------
! FUNCTION NAME: Apply_Mapping_EOD
! DESCRIPTION  : Apply transformation metric of mapping A to
!                Differentiation Matrix D up to Number of Derivative.
!
! ----------------------------------------------------------------------
  Subroutine Apply_Mapping_EOD (Map, N0, N2, Number_Of_Derivatives, &
                                A, D_Even, D_odd)

  integer  :: N0, N2, Number_Of_Derivatives, Map

  REALTYPE, dimension(0:N0,*)                          :: A
  REALTYPE, dimension(0:N2,0:N2,Number_Of_Derivatives) :: D_Even, D_Odd

  integer  :: j

  if (Map == 0) RETURN

  if (Number_Of_Derivatives == 4) then
    do j = 0,N2
      D_Even(:,j,4) = A(0:N2,10)*D_Even(:,j,4) + A(0:N2,9)*D_Even(:,j,3) + &
                      A(0:N2,8 )*D_Even(:,j,2) + A(0:N2,7)*D_Even(:,j,1)
      D_Odd (:,j,4) = A(0:N2,10)*D_Odd (:,j,4) + A(0:N2,9)*D_Odd (:,j,3) + &
                      A(0:N2,8 )*D_Odd (:,j,2) + A(0:N2,7)*D_Odd (:,j,1)
    enddo
  endif

  if (Number_Of_Derivatives >= 3) then
    do j = 0,N2
      D_Even(:,j,3) = A(0:N2,6)*D_Even(:,j,3) + A(0:N2,5)*D_Even(:,j,2) + &
                      A(0:N2,4)*D_Even(:,j,1)
      D_Odd (:,j,3) = A(0:N2,6)*D_Odd (:,j,3) + A(0:N2,5)*D_Odd (:,j,2) + &
                      A(0:N2,4)*D_Odd (:,j,1)
    enddo
  endif

  if (Number_Of_Derivatives >= 2) then
    do j = 0,N2
      D_Even(:,j,2) = A(0:N2,3)*D_Even(:,j,2) + A(0:N2,2)*D_Even(:,j,1) 
      D_Odd (:,j,2) = A(0:N2,3)*D_Odd (:,j,2) + A(0:N2,2)*D_Odd (:,j,1)
    enddo
  endif

  if (Number_Of_Derivatives >= 1) then
    do j = 0,N2
      D_Even(:,j,1) = A(0:N2,1)*D_Even(:,j,1)
      D_Odd (:,j,1) = A(0:N2,1)*D_Odd (:,j,1)
    enddo
  endif

  END Subroutine Apply_Mapping_EOD

! ----------------------------------------------------------------------
! FUNCTION NAME: Smoothing_Matrix
! DESCRIPTION  : Construct the Legendre Filter Matrix corresponding
!                to the array of Filter coefficients Sigma for all points type.
!
!                    Method = 1, Chebyshev
!                    Method = 2, Legendre
!
!                Point_Type = 1, Lobatto
!                Point_Type = 2, Radua
!                Point_Type = 3, Gauss
!
!                if S_Type = 0, make S_Even
!                   S_Type = 1, make S_Odd
!                   S_Type = 2, make regular S.
!
!                S is of size (0:N0,0:N1)
!                N+1 is total number of grid points.
!                Sigma is the filter of size (0:N)
!
! Kind = 1 ! Legendre, w(x) = 1, (-1,1)
! Kind = 2 ! Chebyshev of the first kind, w(x) = 1/sqrt(1-x*x), (-1,1)
! Kind = 3 ! Chebyshev of the second Kind, w(x) = sqrt(1-x*x),  (-1,1)
! Kind = 4 ! Hermite, w(x) = exp(-x*x), (-inf,inf)
! Kind = 5 ! Jacobi, w(x) = (1-x)**alpha*(1+x)**beta, (-1,1), alpha, beta > -1
! Kind = 6 ! Generalized Laguerre, w(x) = exp(-x)*x**alpha, (0, inf), alpha>-1
!
! ----------------------------------------------------------------------

  Subroutine Smoothing_Matrix (Method, Point_Type, S_Type, N, N0, N1, Sigma, S)

  integer  :: Method, Point_Type, N, N0, N1, S_Type
  integer  :: Kind, Nfix, IERR
  logical  :: Symmetry, Even, Odd
  REALTYPE :: gamma, Ratio_1, Ratio_2

  REALTYPE, dimension(0:N0,0:N1) :: S
  REALTYPE, dimension(0:N)       :: Sigma

  REALTYPE, dimension(2)   :: qxfix
  REALTYPE, dimension(0:N) :: x, w, Work, Scale
  REALTYPE, dimension(0:N) :: Prev, Now, Next

  integer  :: i, k, j, N2

  Even     =       (S_Type == 0)
  Odd      =       (S_Type == 1)
  Symmetry = .NOT. (S_Type == 2)

  SELECT CASE (Method)
    CASE (1)
      Kind = 2

    CASE (2)
      Kind = 1

  END SELECT

  SELECT CASE (Point_Type)
    CASE (1)
      Nfix =  2 ; qxfix(1) = -ONE ; qxfix(2) =  ONE

    CASE (2)
      Nfix =  1 ; qxfix(1) = -ONE

      Symmetry = .FALSE.

    CASE (3)
      Nfix =  0

  END SELECT

  call GaussQ (Kind, N+1, ZERO, ZERO, Nfix, qxfix, work, x, w, IERR)

  if (IERR /= 0) then
    write (6 ,*) "Smoothing_Matrix: Gauss Quadrature failed. Quitting."
    write (99,*) "Smoothing_Matrix: Gauss Quadrature failed. Quitting."
    STOP
  endif

  if (Point_Type /= 2) then
    N2 = N/2

    if (MOD(N,2) == 0) x(N2) = ZERO

    do k = 0,N2
      x(N-k) = -x(k)
      w(N-k) =  w(k)
    enddo
  endif

  x = -x

  do k = 0,N
    gamma = ONE/(k+HALF) ; Scale(k) = Sigma(k)/gamma
  enddo

  if (Point_Type == 1) then
    gamma = TWO/N        ; Scale(N) = Sigma(N)/gamma
  endif

  Prev = ONE ; Now = x ; S = ZERO

  do i = 0,N
    if ((.NOT. Symmetry) .OR. (MOD(i,2) == S_Type)) then
      do j = 0,N1
        S(:,j) = S(:,j) + Scale(i)*Prev(0:N0)*Prev(j)
      enddo
    endif

    Ratio_1 = (2*i+3)/(i+TWO)
    Ratio_2 = (  i+1)/(i+TWO)

    Next = Ratio_1*x*Now - Ratio_2*Prev
    Prev = Now
    Now  = Next
  enddo

  do j = 0,N1
    S(:,j) = w(j)*S(:,j)
  enddo

  if (Symmetry) then
    S(:,0:N0) = TWO*S(:,0:N0)

    if (MOD(N,2) == 0) then
      if (Even) then
        S(:,N0) = HALF*S(:,N0)
      else
        S(:,N0) = ZERO
     endif
   endif

  endif

  END Subroutine Smoothing_Matrix

END MODULE EOD_Auxillary_Subroutine
