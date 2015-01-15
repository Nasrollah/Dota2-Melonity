#include "pseudopack.h"


!***************************************************************************
! MODULE     : Lagrangian_Weights 
!
! SUBROUTINE : Lagrangian_Weight         
!
! Function   : Compute the Lagrangian weights w(0:N,k) at a given point z
!              from a given set of N+1 data points x(j),j=0,..,N for the 
!              first k=0,..,Order derivatives.
!
! Name       : Lagrangian_Weight_Output  (Subroutine)
!              Formatted print out the Lagrangian weight w(0:N,0:Order)
!              and its sum for Order <= 4. 
!  
! C-Preprocessor :
! 
! DEBUG             Test the integrity of the input data 
! SINGLE_PRECISION  Single precision code instead of double precision by
!                     default
!
! z      REALTYPE                                           [Input]
!          Grid point at which the Interpolation and/or derivative 
!          to be approximated.  
!
! Order  the order of derivatives desired                   [Input]
!
! N+1    Integer                                            [Input]
!          Extend of the stencils of x(j), j=0,...,N
!
! x      REALTYPE Array of size (N+1)                       [Input]
!          Distinct grid points.  
!          Strictly increasing or decreasing
!
! ldy    Integer                                            [Input]
!          Leading dimension of w(ldy,0:Order)
!
! w      REALTYPE Array of size (ldy,Order+1)               [Output]
!          Lagrangian weight of the first k=0,..,Order derivatives.
!
! Method Integer (OPTIONAL)                                 [Input/Output]
!          1   Use Formula according to Case 0 and Case 1 (see below)
!          2   Use Formula according to Case 0 and Case 2 (see below)
!          3   Use Formula according to Case 0 and Case 3 (see below)
!          4   Use Formula Case 2                         (see below)
!          5   Use Formula Case 3                         (see below)
!          Default = 2
!
! To approximate the m_th derivative of f(z) based on the functional value
!  f(x_0),...,f(x_N),
!     f^m(z) = sum (j=0,N) w(j,m) f(x_j)  , m = 0,..,Order
!
!
! Case 0:
!   If z IS NOT one of the grid point x(j),
!     Uses the formula
!        w(j,m) = [1/(z-x_j)][P^m(x)/c_j - m w(j,m-1)]
!     and recursion formula
!        P_l^m(z) = m P_(l-1)^(m-1)(z) + (z-x_l) P_(l-1)^m(z)
!
! Case 1:
!   If z IS one of the grid point x(j), then recursion relationship
!       w(j  ,1) = (z-x(j))*a(j_z)/a(j)                   if j <> j_z
!       w(j_z,1) = -Sum (j=0,..,m, j <> j_z) w(j,1)
!
!     where a(j) = Sum (l=0,N, l<>j) (x(j)-x(l)).
!
!     for m = 2,min(N,Order)
!       w(j  ,m) = m*[w(j_z,m-1)*w(j,1)-w(j,m-1)/(z-x(j)] if j <> j_z
!       w(j_z,m) = -Sum (j=0,..,m, j <> j_z) w(j,m)
!
!     where j_z is the index such that z = x(j_z)
!
! Case 2:
!   If z IS or IS NOT one of the grid point x(j),
!     See Bengt Fornberg,
!       A practical guide to Pseudospectral methods.
!       Generation of finite difference formulae on
!         arbitrary spaced grids, Math. Comp. 51 (1988), pp. 699-706
!       w(i,j,k) =  [(z-x(j)  )*w(i,j-1  ,k)-k*w(i  ,j-1,k-1)]/(x(i)-x(j)),i<>j
!       w(j,j,k) = -[(z-x(j-1))*w(j-1,j-1,k)-k*w(j-1,j-1,k-1)]*a(j-1)/a(j)
!
!     where a(j) = -Sum (l=0,..,j-1) (x(j)-x(l))
!
!       w(i,k) = w(i,N,k)
!
!   It is rewritten and extended the range of applicable N>700 than the
!     the original formula by computing a(j-1)/a(j) 
!     with EXP and LOG rather than quotient.
!
! Case 3:
!   If z IS or IS NOT one of the grid point x(j),
!     using the original code from Bengt Fornberg's.
!
! AUTHOR      :: Wai Sun Don, PhD
! DATE        :: 3-14-1999
! VERSION     :: 2.0
!***************************************************************************


MODULE Lagrangian_Weights

implicit NONE 

#if defined (SINGLE_PRECISION)
integer , parameter :: N_Limit = 100
REALTYPE, parameter :: Factor  = 1.0e10
#else
integer , parameter :: N_Limit = 700
REALTYPE, parameter :: Factor  = 1.0d7
#endif

PRIVATE 

PUBLIC   :: PS_Lagrangian_Weight
PUBLIC   :: PS_Lagrangian_Weight_Output

CONTAINS

  Subroutine PS_Lagrangian_Weight   (z, Order, N, x, ldy, w, Method, Error)

  integer , INTENT(IN) :: Order, N, ldy
  REALTYPE, INTENT(IN) :: z

  integer , OPTIONAL   :: Method
  logical , OPTIONAL   :: Error

  REALTYPE, dimension(0:N)            , INTENT(IN)  :: x
  REALTYPE, dimension(0:ldy-1,0:Order), INTENT(OUT) :: w

  integer  :: j_z
  logical  :: Grid_Point

#if defined (DEBUG)
  call Check_Input 
#endif

#if 1
  call Is_z_At_Grid_Point 
#endif

  w(0:N,0:Order) = ZERO

  if (PRESENT(Method)) then
    if ((Method == 1) .AND. (N > N_Limit)) Method = 4

    SELECT CASE (Method)
      CASE (1)
        if (Grid_Point) then
          call Lagrangian_Weight_0  (z, Order, N, x, ldy, w, Error)
        else
          call Lagrangian_Weight_1  (z, Order, N, x, ldy, w, Error)
        endif

      CASE (2)
        if (Grid_Point) then
          call Lagrangian_Weight_0  (z, Order, N, x, ldy, w, Error)
        else
          call Lagrangian_Weight_2  (z, Order, N, x, ldy, w, Error)
        endif

      CASE (3)
        if (Grid_Point) then
          call Lagrangian_Weight_0  (z, Order, N, x, ldy, w, Error)
        else
          call Lagrangian_Weight_3  (z, Order, N, x, ldy, w, Error)
        endif

      CASE (4)
          call Lagrangian_Weight_2  (z, Order, N, x, ldy, w, Error)

      CASE (5)
          call Lagrangian_Weight_3  (z, Order, N, x, ldy, w, Error)

    END SELECT 
  else
          call Lagrangian_Weight_2  (z, Order, N, x, ldy, w, Error)
  endif

  CONTAINS

#if 1
    Subroutine Is_z_At_Grid_Point 

    integer              :: j

    Grid_Point = .FALSE.

    do j = 0,N
      if (z <= x(j)) then
        j_z = j 
        if (z == x(j)) then ; Grid_Point = .TRUE. ; EXIT ; endif
      endif
    enddo
 
    END Subroutine Is_z_At_Grid_Point 
#endif

#if defined (DEBUG)
!---------------------------------------------------------------------
    Subroutine Check_Input      

    if (PRESENT(Error)) Error = .FALSE.

    if (PRESENT(Method)) then
      if ((Method <  0) .OR. (Method >  5)) then
        write (6,140) Method ; write (99,140) Method
 
        if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
      endif
    endif

    if (N <  1) then 
      write (6,100) N        ; write (99,100) N

      if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
    endif

    if (N <  Order) then 
      write (6,110) N, Order ; write (99,110) N, Order

      if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
    endif

    if (N >  ldy-1) then
      write (6,120) N, ldy-1 ; write (99,120) N, ldy-1 

      if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
    endif

    if (Order <  0) then 
      write (6,130) Order    ; write (99,130) Order

      if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
    endif

 100 format (1x,'Lagrangian_Weight : N<',i8,'< 1')
 110 format (1x,'Lagrangian_Weight : N=',i8,'< Order=',i2)
 120 format (1x,'Lagrangian_Weight : N=',i8,'> ldy-1=',i5)
 130 format (1x,'Lagrangian_Weight : Order=',i2,'< 0')
 140 format (1x,'Lagrangian_Weight : Method=',i2,'< 0 or > 4')

    END Subroutine Check_Input 
#endif

  END Subroutine PS_Lagrangian_Weight 

!---------------------------------------------------------------------

  Subroutine Lagrangian_Weight_0 (z, Order, N, x, ldy, w, Error)

  integer , INTENT(IN) :: Order, N, ldy
  REALTYPE, INTENT(IN) :: z

  REALTYPE, dimension(0:N)            , INTENT(IN)  :: x
  REALTYPE, dimension(0:ldy-1,0:Order), INTENT(OUT) :: w

  logical , OPTIONAL   :: Error

  integer                  :: j, k, j_z
  logical                  :: Grid_Point
  REALTYPE, dimension(0:N) :: a, b

  call Get_Ratio_a_0 (N, x, a, Error)

#if 1
  call Is_z_At_Grid_Point 
#endif

  w(j_z,0) = ONE

  if (MIN(N,Order) <  1) RETURN

  do j = 0,N
    if (j /= j_z) b(j) = ONE/(z-x(j))
  enddo
  b(j_z) = ZERO

  w(0:N,1) = b*a/a(j_z)

  w(j_z,1) = ZERO
  w(j_z,1) = -SUM(w(0:N,1))

  do k = 2,MIN(N,Order)
    w(0:N,k) = k*(w(j_z,k-1)*w(0:N,1)-b*w(0:N,k-1))

    w(j_z,k) = ZERO
    w(j_z,k) = -SUM(w(0:N,k))
  enddo

#if 1
  CONTAINS

    Subroutine Is_z_At_Grid_Point 

    integer              :: j

    Grid_Point = .FALSE.

    do j = 0,N
      if (z <= x(j)) then
        j_z = j 
        if (z == x(j)) then ; Grid_Point = .TRUE. ; EXIT ; endif
      endif
    enddo
 
    END Subroutine Is_z_At_Grid_Point 
#endif

  END Subroutine Lagrangian_Weight_0 

!---------------------------------------------------------------------
  Subroutine Lagrangian_Weight_1 (z, Order, N, x, ldy, w, Error)

  integer , INTENT(IN) :: Order, N, ldy
  REALTYPE, INTENT(IN) :: z

  REALTYPE, dimension(0:N)            , INTENT(IN)  :: x
  REALTYPE, dimension(0:ldy-1,0:Order), INTENT(OUT) :: w

  logical , OPTIONAL   :: Error

  integer                  :: l, k
  REALTYPE, dimension(0:N) :: a, b

  call Get_Ratio_a_0 (N, x, a, Error)

  b = ONE/(z-x)

    w(0,0) = (z-x(0))
  do l = 1,N
    w(l,0) = (z-x(l))*w(l-1,0)
  enddo

  if (MIN(N,Order) == 0) then
    w(0:N,0) = w(N,0)*a*b

    RETURN
  endif

  do k = 1,MIN(N,Order)
    if (k == 1) w(0,1) = ONE

    do l = 1,N
      w(l,k) = k*w(l-1,k-1) + (z-x(l))*w(l-1,k)
    enddo
  enddo

    w(0:N,0) = w(N,0)*a*b

  do k = 1,MIN(N,Order)
    w(0:N,k) = (w(N,k)*a-k*w(0:N,k-1))*b
  enddo

  END Subroutine Lagrangian_Weight_1 

!---------------------------------------------------------------------
  Subroutine Lagrangian_Weight_2 (z, Order, N, x, ldy, w, Error)

  integer , INTENT(IN) :: Order, N, ldy
  REALTYPE, INTENT(IN) :: z

  REALTYPE, dimension(0:N)            , INTENT(IN)  :: x
  REALTYPE, dimension(0:ldy-1,0:Order), INTENT(OUT) :: w

  logical , OPTIONAL   :: Error

  integer                  :: i, j, k
  REALTYPE                 :: x_ij
  REALTYPE, dimension(0:N) :: a, b

  call Get_Ratio_a_2 (N, x, a, Error)

  b = z-x

  w(0,0) = ONE

  do j = 1,N
    do i = 0,j-1

      x_ij = ONE/(x(i)-x(j))

      if (i == j-1) then
        do k = MIN(j,Order),1,-1
          w(j,k) = a(j)*(b(j-1)*w(j-1,k) + k*w(j-1,k-1))
        enddo

          w(j,0) = a(j)* b(j-1)*w(j-1,0)
      endif

        do k = MIN(j,Order),1,-1
          w(i,k) = x_ij*(b(j)*w(i,k) + k*w(i,k-1))
        enddo

          w(i,0) = x_ij* b(j)*w(i,0)

    enddo
  enddo
       
  END Subroutine Lagrangian_Weight_2 

!---------------------------------------------------------------------
  Subroutine Lagrangian_Weight_3 (z, Order, N, x, ldy, w, Error)

  integer , INTENT(IN) :: Order, N, ldy
  REALTYPE, INTENT(IN) :: z

  REALTYPE, dimension(0:N)            , INTENT(IN)  :: x
  REALTYPE, dimension(0:ldy-1,0:Order), INTENT(OUT) :: w

  logical , OPTIONAL   :: Error

  integer                  :: i, j, k
  REALTYPE                 :: c1,c2,c3,c4,c5
  REALTYPE, dimension(0:N) :: a

  call Get_Ratio_a_3 (N, x, a, Error)

  w(0,0) = ONE

  c1 = ONE ; c4 = x(0)-z

  do i = 1,N

    c2 = ONE ; c5 = c4 ; c4 = x(i)-z

    do j = 0,i-1
      c3 = x(i)-x(j) ; c2 = c2*c3

      if (j == i-1) then
        do k = MIN(i,Order),1,-1
          w(i,k) = -c1*(c5*w(i-1,k)-k*w(i-1,k-1))/c2
        enddo
          w(i,0) = -c1* c5*w(i-1,0)              /c2
      endif

        do k = MIN(i,Order),1,-1
          w(j,k) = (c4*w(j,k)-k*w(j,k-1))/c3
        enddo
          w(j,0) =  c4*w(j,0)            /c3
      enddo

    c1 = c2
  enddo

  END Subroutine Lagrangian_Weight_3 

!---------------------------------------------------------------------
  Subroutine Get_Ratio_a_0 (N, x, a, Error)

  integer , INTENT(IN) :: N

  REALTYPE, dimension(0:N), INTENT(IN) :: x
  REALTYPE, dimension(0:N)             :: a

  logical , OPTIONAL   :: Error

  integer              :: j, l
  REALTYPE             :: Product

  if (PRESENT(Error)) Error = .FALSE.

  if (N <= N_Limit) then
    do j = 0,N

      Product = ONE
      do l = 0,N
        if (l /= j) Product = Product * (x(j)-x(l))
      enddo

      if (ABS(Product) <  Factor*TINY(ONE)) then 
        write (6,100) Factor*TINY(ONE) ; write (99,100) Factor*TINY(ONE) 

        if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
      endif

      a(j) = ONE/Product
    enddo
  else
    do j = 0,N

      Product = ZERO
      do l = 0,N
        if (l /= j) Product = Product + LOG(ABS(x(j)-x(l)))
      enddo

      if (EXP(Product) <  Factor*TINY(ONE)) then 
        write (6,100) Factor*TINY(ONE) ; write (99,100) Factor*TINY(ONE) 

        if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
      endif

      a(j) = (-1)**(N-j)/EXP(Product)
    enddo
  endif

 100 format (1x,'FATAL ERROR : Get_Ratio_a_0 '/                            &
             1x,'              Product of Grid Spacing < ',g14.7)

  END Subroutine Get_Ratio_a_0

!---------------------------------------------------------------------
  Subroutine Get_Ratio_a_2 (N, x, a, Error)

  integer , INTENT(IN) :: N
#if 0
  REALTYPE             :: a_sign
#endif

  REALTYPE, dimension(0:N), INTENT(IN) :: x
  REALTYPE, dimension(0:N)             :: a

  logical , OPTIONAL   :: Error

  integer              :: j, l
  REALTYPE             :: Product

  if (PRESENT(Error)) Error = .FALSE.

#if 1
  a(0) = ONE

  do j = 1,N

    Product = ONE
    do l = 0,j-1
      Product = Product * (x(j)-x(l))
    enddo

    if (ABS(Product) < Factor*TINY(ONE)) then 
      write (6,100) Factor*TINY(ONE) ; write (99,100) Factor*TINY(ONE) 

      if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
    endif

    a(j) = Product
  enddo

  a(N:1:-1) = a(N-1:0:-1)/a(N:1:-1)
#else
  a(0) = ZERO

  do j = 1,N

    Product = ZERO
    do l = 0,j-1
      Product = Product + LOG(ABS(x(j)-x(l)))
    enddo

    if (EXP(Product) < Factor*TINY(ONE)) then 
      write (6,100) Factor*TINY(ONE) ; write (99,100) Factor*TINY(ONE) 

      if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
    endif

    a(j) = Product
  enddo

     a_sign = SIGN(ONE,x(1)-x(0))

  a(N:1:-1) = a_sign*EXP(a(N-1:0:-1)-a(N:1:-1))

#endif

 100 format (1x,'FATAL ERROR : Get_Ratio_a_2'/                         &
             1x,'              Product of Grid Spacing < ',g14.7)

  END Subroutine Get_Ratio_a_2

!---------------------------------------------------------------------
  Subroutine Get_Ratio_a_3 (N, x, a, Error)

  integer , INTENT(IN) :: N

  REALTYPE, dimension(0:N), INTENT(IN) :: x
  REALTYPE, dimension(0:N)             :: a

  logical , OPTIONAL   :: Error
#if 0
  integer              :: j, l
  REALTYPE             :: Product
#endif

  if (PRESENT(Error)) Error = .FALSE.

#if 0
  do j = 0,N

    Product = ONE
    do l = 0,N
      if (l /= j) Product = Product * (x(j)-x(l))
    enddo

    if (ABS(Product) < Factor*TINY(ONE)) then
      write (6,100) Factor*TINY(ONE) ; write (99,100) Factor*TINY(ONE) 

      if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
    endif

  enddo
#else
  if (N >  N_Limit) then
    write (6,110) N ; write (99,110) N 

    if (PRESENT(Error)) then ; Error = .TRUE. ; else ; STOP ; endif
  endif
#endif

 100 format (1x,'FATAL ERROR : Get_Ratio_a_3'/                         &
             1x,'              Product of Grid Spacing < ',g14.7)
 110 format (1x,'FATAL ERROR : Get_Ratio_a_3'/                         &
             1x,'              N = ',i8,' > N_Limit= ',i8)

  END Subroutine Get_Ratio_a_3

!---------------------------------------------------------------------
  Subroutine PS_Lagrangian_Weight_Output (z, Order, N, x, ldy, w, lid, Method)

  integer , INTENT(IN) :: Order, N, ldy, lid
  REALTYPE, INTENT(IN) :: z

  integer , INTENT(IN), OPTIONAL :: Method

  REALTYPE, dimension(0:N)            , INTENT(IN)  :: x
  REALTYPE, dimension(0:ldy-1,0:Order), INTENT(IN)  :: w

  integer                      :: i, k, i_Max
  REALTYPE, dimension(0:Order) :: Total_Sum

  if (PRESENT(Method)) then
    write (lid,99) Method, z, Order
  else
    write (lid,99) 2     , z, Order
  endif

  if (Order == 0) write (lid,200)
  if (Order == 1) write (lid,201)
  if (Order == 2) write (lid,202)
  if (Order == 3) write (lid,203)
  if (Order == 4) write (lid,204)

  do i = 0,N
    write (lid,101) i, x(i), (w(i,k), k=0,Order)
  enddo
  write (lid,*)

  Total_Sum(0) = SUM(w(0:N,0))
  do k = 1,Order

    i_Max = 0
    do i = 0,N
      if (ABS(w(i,k)) >  ABS(w(i_Max,k))) i_Max = i
    enddo

    Total_Sum(k) = SUM(w(0:N,k))/(ABS(w(i_Max,k)+EPSILON(ONE)))
  enddo

  write (lid,102) (Total_Sum(k), k=0,Order)
  write (lid,103)
        
  99  format (1x/1x,'Method = ',i1,', z = ',g20.10,' Order = ',i4/)
 101 format (i4    ,     6(e12.4,:))
 102 format (' Sum',12x, 5(e12.4,:))
 103 format (1x,78('-'))

 200  format (2x,'i',5x,'x(i)',7x,'w(i,0)')
 201  format (2x,'i',5x,'x(i)',7x,'w(i,0)',6x,'w(i,1)')
 202  format (2x,'i',5x,'x(i)',7x,'w(i,0)',6x,'w(i,1)',6x,'w(i,2)')
 203  format (2x,'i',5x,'x(i)',7x,'w(i,0)',6x,'w(i,1)',6x,'w(i,2)',    &
                               6x,'w(i,3)')
 204  format (2x,'i',5x,'x(i)',7x,'w(i,0)',6x,'w(i,1)',6x,'w(i,2)',    &
                               6x,'w(i,3)',6x,'w(i,4)')

  END Subroutine PS_Lagrangian_Weight_Output 

END MODULE Lagrangian_Weights
