#include "pseudopack.h"


MODULE MXM_Auxillary_Subroutine
  USE Gauss_Quadrature

implicit NONE

PRIVATE 

PUBLIC  :: Diagonal_Element_By_Sum_MXM
PUBLIC  :: Diagonal_Element_By_MXM_MXM
PUBLIC  :: Centro_Symmetry_Property
PUBLIC  :: Flip_Up
PUBLIC  :: Flip_Down
PUBLIC  :: Apply_Mapping
PUBLIC  :: Smoothing_Matrix

CONTAINS

! ------------------------------------------------------------------------
  Subroutine Diagonal_Element_By_Sum_MXM (Symmetry, N, O_Of_D, D)

  integer  :: N, O_Of_D
  logical  :: Symmetry

  REALTYPE, dimension(0:N,0:N,O_Of_D) :: D

  integer  :: j, k, M, N2
  REALTYPE :: Sum

  M = O_Of_D

  if (Symmetry) then ; N2 = N/2 ; else ; N2 = N ; endif

  do k = 0,N2

    Sum = ZERO
    do j = 0,k-1
      Sum = Sum - D(k,j,M)
    enddo
    do j = k+1,N
      Sum = Sum - D(k,j,M)
    enddo

    D(k,k,M) = Sum
  enddo

  END Subroutine Diagonal_Element_By_Sum_MXM

! ------------------------------------------------------------------------
  Subroutine Diagonal_Element_By_MXM_MXM (Symmetry, N, O_Of_D, D)

  integer  :: N, O_Of_D
  logical  :: Symmetry

  REALTYPE, dimension(0:N,0:N,O_Of_D) :: D

  integer  :: j, k, M, N2
  REALTYPE :: Sum, Factor

  M = O_Of_D

  if (Symmetry) then

        N2 = N/2
    Factor = (-1)**(M-1)

    if (MOD(N,2) .ne. 0) then

      do k = 0,N2

        Sum = ZERO
        do j = 0,N2
          Sum = Sum + D(k,j,1)*D(j,k,M-1) + Factor*D(k,N-j,1)*D(j,N-k,M-1)
        enddo

        D(k,k,M) = Sum
      enddo

    else

      do k = 0,N2

        Sum = ZERO
        do j = 0,N2-1
          Sum = Sum + D(k,j,1)*D(j,k,M-1) + Factor*D(k,N-j,1)*D(j,N-k,M-1)
        enddo

        D(k,k,M) = Sum + Factor*D(k,N2,1)*D(N2,N-k,M-1)
      enddo

    endif

  else

    do k = 0,N

      Sum = ZERO
      do j = 0,N
        Sum = Sum + D(k,j,1)*D(j,k,M-1)
      enddo

      D(k,k,M) = Sum
    enddo

  endif

  END Subroutine Diagonal_Element_By_MXM_MXM

! ------------------------------------------------------------------------

  Subroutine Centro_Symmetry_Property (Symmetry, Direction, N, O_Of_D, D)

  integer  :: N, O_Of_D, Direction
  logical  :: Symmetry

  REALTYPE, dimension(0:N,0:N,O_Of_D) :: D

  integer  :: M

  if (.NOT. Symmetry) RETURN

  do M = 1, O_Of_D
    if (Direction ==  1) call Flip_Up   (M, N, D(0,0,M))
    if (Direction == -1) call Flip_Down (M, N, D(0,0,M))
  enddo

  END Subroutine Centro_Symmetry_Property

! ------------------------------------------------------------------------
  Subroutine Flip_Up   (O_Of_D, N, D)

  integer  :: O_Of_D, N

  REALTYPE, dimension(0:N,0:N) :: D

  integer  :: j, k, N2
  REALTYPE :: Factor

  N2 = N/2

  Factor = (-1)**O_Of_D

  if (MOD(N,2) == 0) then
    do k = 0,N2-1
      do j = 0,N
        D(k,j) = Factor*D(N-k,N-j)
      enddo
    enddo

      do j = 0,N2-1
        D(N2,j) = Factor*D(N2,N-j)
      enddo
  else
    do k = 0,N2
      do j = 0,N
        D(k,j) = Factor*D(N-k,N-j)
      enddo
    enddo
  endif

  END Subroutine Flip_Up

! ------------------------------------------------------------------------
  Subroutine Flip_Down (O_Of_D, N, D)

  integer  :: O_Of_D, N

  REALTYPE, dimension(0:N,0:N) :: D

  integer  :: j, k, N2
  REALTYPE :: Factor

  N2 = N/2

  Factor = (-1)**O_Of_D

  if (MOD(N,2) == 0) then
    do k = 0,N2-1
      do j = 0,N
        D(N-k,N-j) = Factor*D(k,j)
      enddo
    enddo

      do j = 0,N2-1
        D(N2,N-j) = Factor*D(N2,j)
      enddo
  else
    do k = 0,N2
      do j = 0,N
        D(N-k,N-j) = Factor*D(k,j)
      enddo
    enddo
  endif

  END Subroutine Flip_Down

! ----------------------------------------------------------------------
! FUNCTION NAME: Apply_Mapping
! DESCRIPTION  : Apply transformation metric of mapping A to
!                Differentiation Matrix D up to Number of Derivative.
!
! ----------------------------------------------------------------------
  Subroutine Apply_Mapping (Map, N, Number_Of_Derivatives, A, D)
                                    

  integer  :: N, Number_Of_Derivatives, Map

  REALTYPE, dimension(0:N,*)                         :: A
  REALTYPE, dimension(0:N,0:N,Number_Of_Derivatives) :: D

  integer  :: j

  if (Map == 0) RETURN

  if (Number_Of_Derivatives == 4) then
    do j = 0,N
      D(:,j,4) = A(:,10)*D(:,j,4) + A(:,9)*D(:,j,3) + &
                 A(:,8 )*D(:,j,2) + A(:,7)*D(:,j,1)
    enddo
  endif

  if (Number_Of_Derivatives >= 3) then
    do j = 0,N
      D(:,j,3) = A(:,6)*D(:,j,3) + A(:,5)*D(:,j,2) + &
                 A(:,4)*D(:,j,1)
    enddo
  endif

  if (Number_Of_Derivatives >= 2) then
    do j = 0,N
      D(:,j,2) = A(:,3)*D(:,j,2) + A(:,2)*D(:,j,1) 
    enddo
  endif

  if (Number_Of_Derivatives >= 1) then
    do j = 0,N
      D(:,j,1) = A(:,1)*D(:,j,1)
    enddo
  endif

  END Subroutine Apply_Mapping 

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

  integer  :: i, j, k, N2

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

  call gaussq (Kind, N+1, ZERO, ZERO, Nfix, qxfix, work, x, w, IERR)

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

END MODULE MXM_Auxillary_Subroutine
