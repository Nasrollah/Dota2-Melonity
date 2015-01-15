#if defined (LOCAL_GAUSS_QUADRATURE)

#include "pseudopack.h"


! ---------------------------------------------------------------------
! This set of routines computes the nodes t(j) and weights
!   w(j) for gaussian-type quadrature rules with pre-assigned nodes.
! These are used when one wishes to approximate
!   integral (from a to b)  f(x) w(x) dx
!
!      n
!     by     sum w  f(t )
!     j=1         j    j
!
! (Note w(x) and w(j) have no connection with each other.)
! Here w(x) is one of six possible non-negative weight
! functions (listed below), and f(x) is the function to be integrated.  
!
! Gaussian quadrature is particularly useful on infinite intervals 
!  (with appropriate weight Functions), since then other techniques often fail.
!
! Associated with each weight Function w(x) is a set of orthogonal polynomials.
! The nodes t(j) are just the zeroes of the proper n-th degree polynomial.
!
! Input parameters (all real numbers are in REALTYPE)
!
! Kind     an integer between 1 and 6 giving the type of quadrature rule:
!
! Kind = 1:  Legendre quadrature, w(x) = 1 on (-1, 1)
! Kind = 2:  Chebyshev quadrature of the first Kind
!              w(x) = 1/sqrt(1 - x*x) on (-1, +1)
! Kind = 3:  Chebyshev quadrature of the second Kind
!              w(x) = sqrt(1 - x*x) on (-1, 1)
! Kind = 4:  Hermite quadrature, w(x) = exp(-x*x) on (-infinity, +infinity)
! Kind = 5:  Jacobi quadrature, w(x) = (1-x)**alpha * (1+x)** beta on (-1, 1),
!              alpha, beta .gt. -1.
!            Note: Kind=2 and 3 are a special case of this.
! Kind = 6:  Generalized Laguerre quadrature, w(x) = exp(-x)* x**alpha 
!              on (0, +infinity), alpha .gt. -1
!
! n        the number of points used for the quadrature rule
! alpha    Real parameter used only for gauss-Jacobi and gauss-
!            Laguerre quadrature (otherwise use 0.d0).
! beta     real parameter used only for gauss-Jacobi quadrature--
!            (otherwise use 0.d0)
! kpts     (integer) normally 0, unless the left or right end-point
!             (or both) of the interval is required to be a
!             node (this is called gauss-radau or gauss-lobatto
!             quadrature).  then kpts is the number of fixed
!             endpoints (1 or 2).
! End_Points  real array of length 2.  contains the values of
!               any fixed endpoints, if kpts = 1 or 2.
! b       real scratch array of length n
!
! Output parameters (both REALTYPE arrays of length n)
!
! t will contain the desired nodes.
! w will contain the desired weights w(j).
!
! Underflow may sometimes occur, but is harmless.
!
!     References
! 1.  Golub, G. H., and Welsch, J. H., "Calculation of Gaussian
!     Quadrature Rules," Mathematics of Computation 23 (April,
!     1969), pp. 221-230.
! 2.  Golub, G. H., "Some Modified Matrix Eigenvalue Problems,"
!     Siam Review 15 (April, 1973), pp. 318-334 (section 7).
! 3.  Stroud and Secrest, Gaussian quadrature formulas, prentice-hall
!     Englewood Cliffs, N.J., 1966.
!
! original version 20 jan 1975 from stanford
! modified 21 dec 1983 by eric grosse
!   imtql2 => gausq2
!   hex constant => d1mach (from core library)
!   compute pi using datan
!   removed accuracy claims, description of method
!   added single precision version
! ---------------------------------------------------------------------

MODULE Gauss_Quadrature
#if defined (LOCAL_GAMMA)
  USE Gamma_Function
#endif

implicit NONE

PRIVATE

PUBLIC  :: GaussQ

CONTAINS

#ifndef Gauss_F77

  Subroutine GaussQ (Kind, n, alpha, beta, kpts, End_Points, b, t, w, IERR)

  integer  :: Kind, n, kpts, IERR
  REALTYPE :: Mu_Zero, t1, alpha, beta, kappa

  REALTYPE, dimension(2) :: End_Points
  REALTYPE, dimension(n) :: b, t, w
 
  call Class (Kind, n, alpha, beta, b, t, Mu_Zero)
 
! The matrix of coefficients is assumed to be symmetric.
! The array t contains the     diagonal elements
! The array b contains the off-diagonal elements.
! Make appropriate changes in the lower right 2 by 2 submatrix.
 
  SELECT CASE (kpts)
    CASE (1)   ! only t(n) must be changed 
      t(n) = Solve(End_Points(1),n,t,b)*b(n-1)**2+End_Points(1)
 
    CASE (2)   ! t(n) and b(n-1) must be recomputed
      kappa =  Solve(End_Points(1),n,t,b)
         t1 = ((End_Points(1)-End_Points(2))/(Solve(End_Points(2),n,t,b)-kappa))

     b(n-1) = SQRT(t1)
     t(n)   = End_Points(1) + kappa*t1

  END SELECT
 
! Note that the indices of the elements of b run from 1 to n-1
!   and thus the value of b(n) is arbitrary.
! Now compute the eigenvalues of the symmetric tridiagonal matrix,
!   which has been modified as necessary.
! The method used is a QL-type method with origin shifting
 
  w(1) = ONE ; w(2:n) = ZERO

  call GaussQ_2 (n, t, b, w, IERR)

  w = Mu_Zero*w*w
 
  END Subroutine GaussQ 

! ---------------------------------------------------------------------
! This procedure performs elimination to Solve for the
!   n-th component of the solution delta to the equation
!
!     (jn - shift*identity)*delta  = en,
!
! where en is the vector of all zeroes except for 1 in the n-th position.
!
! The matrix jn is symmetric tridiagonal, with diagonal elements a(i), 
!   off-diagonal elements b(i).  
! 
! This equation must be Solved to obtain the appropriate changes in the lower
! 2 by 2 submatrix of coefficients for orthogonal polynomials.
! ---------------------------------------------------------------------
 
  REALTYPE FUNCTION Solve (shift, n, a, b)

  integer  :: n, i
  REALTYPE :: shift, alpha

  REALTYPE, dimension(n) :: a, b
 
  alpha = a(1) - shift

  do i = 2, n-1
    alpha = a(i) - shift - b(i-1)**2/alpha
  enddo

  Solve = ONE/alpha

  END FUNCTION Solve 
! ---------------------------------------------------------------------
!
! This procedure supplies the coefficients a(j), b(j) of the recurrence relation
!
!     b p (x) = (x - a ) p   (x) - b   p   (x)
!      j j    j   j-1   j-1 j-2
!
!  for the various Classical (normalized) orthogonal polynomials,
!  and the zero-th moment
!
!     Mu_Zero = integral w(x) dx
!
!  of the given polynomial's weight FUNCTION w(x).  
!  Since the polynomials are orthonormalized, the tridiagonal matrix is
!  guaranteed to be symmetric.
!
!  The input parameter alpha is used only for Laguerre and
!  Jacobi polynomials, and the parameter beta is used only for
!  Jacobi polynomials.  
!  The Laguerre and Jacobi polynomials require the Gamma Function.
! ---------------------------------------------------------------------

  Subroutine Class (Kind, n, alpha, beta, b, a, Mu_Zero)

  integer  :: Kind, n, i
  REALTYPE :: Mu_Zero, alpha, beta
  REALTYPE :: abi, b2a2, ab
  REALTYPE :: a1_d, b1_d, ab_d

  REALTYPE, dimension(n) :: a, b 

  SELECT CASE (Kind)
! Legendre polynomials P(x) on (-1, +1), w(x) = 1.
    CASE (1)   
      Mu_Zero = TWO

      a = ZERO
      do i = 1, n
        b(i) = i/SQRT(4*i*i - ONE)
      enddo
 
! Chebyshev polynomials of the first Kind T(x) on (-1,1), 
!   w(x) = 1/sqrt(1 - x*x)
    CASE (2)   
      Mu_Zero = PI

      a = ZERO
      b(1) = SQRT(HALF) ; b(2:n) = ZERO

! Chebyshev polynomials of the second Kind u(x) on (-1,1),
!   w(x) = sqrt(1 - x*x)
 
    CASE (3)
      Mu_Zero = PI/TWO

      a = ZERO
      b = HALF
 
! Hermite polynomials h(x) on (-infinity,+infinity), w(x) = exp(-x**2)
    CASE (4) 
      Mu_Zero = SQRT(PI)

      a = ZERO
     do i = 1, n
       b(i) = SQRT(HALF*i)
     enddo
 
! Jacobi polynomials p(alpha, beta)(x) on (-1,1), 
!   w(x)=(1-x)**alpha*(1+x)**beta, alpha>-1, beta>-1
    CASE (5) 
       ab  = alpha + beta
       abi = TWO + ab

      a1_d = alpha+ONE
      b1_d =  beta+ONE
      ab_d = a1_d+b1_d

      Mu_Zero = TWO**(ab+ONE)*GAMMA(a1_d)*GAMMA(b1_d)/GAMMA(ab_d)

      b2a2 = beta*beta - alpha*alpha

      a(1) = (beta - alpha)/abi
      b(1) = (TWO/abi)*SQRT((alpha+ONE)*(beta+ONE)/(abi+ONE))
      do i = 2, n
         abi = 2*i + ab

        a(i) = b2a2/((abi-TWO)*abi)
        b(i) = (TWO/abi)*SQRT(i*(alpha+i)*(beta+i)*(ab+i)/(abi*abi-ONE))
      enddo

! Laguerre polynomials l(alpha)(x) on (0, infinity), 
!   w(x)=exp(-x)*x**alpha, alpha > -1
    CASE (6) 
        a1_d = alpha+ONE

      Mu_Zero = GAMMA(a1_d)

      do i = 1, n
        a(i) = 2*i-1 + alpha
        b(i) = SQRT(i*(i+alpha))
      enddo

  END SELECT

  END Subroutine Class 

! ---------------------------------------------------------------------
! This Subroutine is a translation of an algol procedure,
! num. math. 12, 377-383(1968) by martin and wilkinson,
! as modified in num. math. 15, 450(1970) by dubrulle.
! handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
! This is a modified version of the 'eispack' routine imtql2.
!
! This Subroutine finds the eigenvalues and first components of the
! eigenvectors of a symmetric tridiagonal matrix by the implicit QL method.
!
! On input:
!
!   n is the order of the matrix;
!
!   d contains the diagonal elements of the input matrix;
!
!   e contains the subdiagonal elements of the input matrix
!     in its first n-1 positions.  e(n) is arbitrary;
!
!   z contains the first row of the identity matrix.
!
! On output:
!
!   d contains the eigenvalues in ascending order.  if an
!     error exit is made, the eigenvalues are correct but
!     unordered for indices 1, 2, ..., IERR-1;
!
!   e has been destroyed;
!
!   z contains the first components of the orthonormal eigenvectors
!     of the symmetric tridiagonal matrix.  
!
!     If an error exit is made,
!     z contains the eigenvectors associated with the stored eigenvalues.
!
!  IERR is set to
!   0  Normal RETURN,
!   j  If the j-th eigenvalue has not been determined after 100 iterations.
!
! ---------------------------------------------------------------------
 
  Subroutine GaussQ_2 (n, d, e, z, IERR)

  integer, parameter :: Max_Iteration = 100

  integer  :: i, j, k, l, m, n, Iteration, IERR
  REALTYPE :: b, c, f, g, p, r, s

  REALTYPE, dimension(n) :: d, e, z

  IERR = 0
  if (n == 1) RETURN
 
  e(n) = ZERO
  do l = 1, n

    Iteration = 0
    do while (Iteration <= Max_Iteration)

! Look for small sub-diagonal element 
      do m = l, n-1
        if (ABS(e(m)) <= EPSILON(ONE)*(ABS(d(m))+ABS(d(m+1)))) EXIT
      enddo

      p = d(l)
      if (m == l) EXIT

      Iteration = Iteration + 1

! Form shift 
      g = (d(l+1) - p)/(TWO*e(l))
      r = SQRT(g*g+ONE)
      g = d(m) - p + e(l)/(g + SIGN(r, g))

      s = ONE ; c = ONE

      p = ZERO
 
! For i=m-1 step -1 until l 
      do i = m-1,l,-1
        f = s*e(i) ; b = c*e(i)

        if (ABS(f) >= ABS(g)) then
          c = g/f ; r = SQRT(c*c+ONE)

          e(i+1) = f*r ; s = ONE/r ; c = c*s
        else
          s = f/g ; r = SQRT(s*s+ONE)

          e(i+1) = g*r ; c = ONE/r ; s = s*c
        endif

        g = d(i+1) - p
        r = (d(i) - g)*s + TWO*c*b
        p = s*r

        d(i+1) = g + p
        g = c*r - b

! Form first component of vector 
        f = z(i+1) ; z(i+1) = s*z(i) + c*f ; z(i) = c*z(i) - s*f
      enddo
 
      d(l) = d(l) - p
      e(l) = g
      e(m) = ZERO
    enddo

    if (Iteration >  Max_Iteration) then ; IERR = 1 ; RETURN ; endif
  enddo
 
! Order eigenvalues and eigenvectors ::::::::::
  do i = 1, n-1

! Find the minimum of d at k
    k = i ; p = d(i)
    do j = i+1, n
      if (d(j) <  p) then ;  k = j ; p = d(j) ; endif
    enddo
 
    if (k == i) CYCLE

! Swap d and z at index i with index k
    d(k) = d(i) ; d(i) = p

    p = z(i) ; z(i) = z(k) ; z(k) = p
  enddo

  END Subroutine GaussQ_2 

#else

      Subroutine GaussQ (Kind, n, alpha, beta, kpts, End_Points,  &
                         b, t, w, IERR)

      integer    Kind, n, kpts, IERR
      REALTYPE   b(n), t(n), w(n), End_Points(2), Mu_Zero, t1, &
                 Solve, alpha, beta, kappa
 
      call class (Kind, n, alpha, beta, b, t, Mu_Zero)
 
      if (kpts .eq. 1) then
        t(n) = Solve(End_Points(1),n,t,b)*b(n-1)**2+End_Points(1)
      endif
 
      if (kpts .eq. 2) then
         kappa =  Solve(End_Points(1),n,t,b)
            t1 = ((End_Points(1)-End_Points(2))/ &
                 (Solve(End_Points(2),n,t,b)-kappa))

        b(n-1) = SQRT(t1)
        t(n)   = End_Points(1) + kappa*t1
      endif
 
      w(1) = ONE
      do i = 2, n
        w(i) = ZERO
      enddo
 
      call GaussQ_2 (n, t, b, w, IERR)

      do i = 1, n
        w(i) = Mu_Zero * w(i) * w(i)
      enddo
 
      return
      end
!
! ======================================================================
!
      REALTYPE Function Solve(shift, n, a, b)

      integer  n, i
      REALTYPE shift, a(n), b(n), alpha
 
      alpha = a(1) - shift

      do i = 2, n-1
         alpha = a(i) - shift - b(i-1)**2/alpha
      enddo

      Solve = ONE/alpha

      return
      end
!
! ======================================================================
!
      Subroutine class (Kind, n, alpha, beta, b, a, Mu_Zero)

      integer   Kind, n
      REALTYPE  a(n), b(n), Mu_Zero, alpha, beta
      REALTYPE  abi, b2a2, ab

      REALTYPE GAMMA, a1_d, b1_d, ab_d
      external GAMMA

! Kind = 1:  Legendre polynomials p(x)
!            on (-1, +1), w(x) = 1.
 
      if (Kind .eq. 1) then
        Mu_Zero = TWO

        do i = 1, n
          a(i) = ZERO
          b(i) = i/SQRT(4*i*i - ONE)
        enddo
      endif
 
! Kind = 2:  Chebyshev polynomials of the first Kind t(x) on
!            (-1, +1), w(x) = 1 / sqrt(1 - x*x)
 
      if (Kind .eq. 2) then
        Mu_Zero = PI

        do i = 1, n
          a(i) = ZERO
          b(i) = HALF
        enddo
          b(1) = SQRT(HALF)
      endif
 
! Kind = 3:  Chebyshev polynomials of the second Kind u(x) on
!            (-1, +1), w(x) = sqrt(1 - x*x)
 
      if (Kind .eq. 3) then
        Mu_Zero = PI/TWO

        do i = 1, n
          a(i) = ZERO
          b(i) = HALF
        enddo
      endif
 
! Kind = 4:  Hermite polynomials h(x) on 
!            (-infinity,+infinity), w(x) = exp(-x**2)
 
      if (Kind .eq. 4) then
        Mu_Zero = SQRT(PI)

        do i = 1, n
          a(i) = ZERO
          b(i) = SQRT(HALF*i)
        enddo
      endif
 
! Kind = 5:  Jacobi polynomials p(alpha, beta)(x) on
!            (-1, +1), w(x)=(1-x)**alpha*(1+x)**beta, alpha>-1, beta>-1
 
      if (Kind .eq. 5) then
             ab = alpha + beta
            abi = TWO + ab

             a1_d = alpha+ONE
             b1_d =  beta+ONE
             ab_d = a1_d+b1_d

        Mu_Zero = TWO**(ab+ONE)*GAMMA(a1_d)*GAMMA(b1_d)/GAMMA(ab_d)

        b2a2 = beta*beta - alpha*alpha

        a(1) = (beta - alpha)/abi
        b(1) = (TWO/abi)*SQRT((alpha+ONE)*(beta+ONE)/(abi+ONE))
        do i = 2, n
           abi = 2*i + ab

          a(i) = b2a2/((abi-TWO)*abi)
          b(i) = (TWO/abi)*SQRT(i*(alpha+i)*(beta+i)*(ab+i)/ &
                                (abi*abi-ONE))
        enddo
      endif
 
! Kind = 6:  Laguerre polynomials l(alpha)(x) on
!            (0, +infinity), w(x)=exp(-x)*x**alpha, alpha > -1
 
      if (Kind .eq. 6) then
           a1_d = alpha+ONE

        Mu_Zero = GAMMA(a1_d)

        do i = 1, n
          a(i) = 2*i-1 + alpha
          b(i) = SQRT(i*(i+alpha))
        enddo
      endif

      return
      end
!
! ======================================================================
!
      Subroutine GaussQ_2 (n, d, e, z, IERR)

#if 0
      integer   i, j, k, l, m, n, ii, IERR
      REALTYPE  d(n), e(n), z(n), b, c, f, g, p, r, s, machep
      REALTYPE  d1mach
 
      machep=d1mach(4)
 
      IERR = 0
      if (n .eq. 1) return
 
      e(n) = ZERO
      do l = 1, n
         j = 0

!     :::::::::: look for small sub-diagonal element ::::::::::
  105    do m = l, n
           if (m .eq. n) go to 120
           if (ABS(e(m)) .le. machep*(ABS(d(m))+ABS(d(m+1)))) goto 120
         enddo
  120    continue

         p = d(l)

         if (m .eq. l) then
           IERR = 0
         else
           if (j .eq. 100) then
             IERR = 1
             return
           endif

           j = j + 1

!     :::::::::: form shift ::::::::::
           g = HALF*(d(l+1) - p)/e(l)
           r = SQRT(g*g + ONE)

           g = d(m) - p + e(l)/(g+SIGN(r,g))
           s = ONE
           c = ONE
           p = ZERO
 
!     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
           do ii = 1, m - l
             i = m - ii
             f = s*e(i)
             b = c*e(i)

             if (ABS(f) .lt. ABS(g)) then
               s = f/g
               r = SQRT(s*s + ONE)

               e(i+1) = g*r
                    c = ONE/r
                    s = s*c
             else
               c = g/f
               r = SQRT(c*c + ONE)
 
               e(i+1) = f*r
                    s = ONE/r
                    c = c*s
             endif

             g =  d(i+1)    - p
             r = (d(i)-g)*s + TWO*c*b
             p = s*r

             d(i+1) = g   + p
             g      = c*r - b

!     :::::::::: form first component of vector ::::::::::
                  f = z(i+1)
             z(i+1) = s*z(i) + c*f
             z(i  ) = c*z(i) - s*f
           enddo
 
           d(l) = d(l) - p
           e(l) = g
           e(m) = ZERO

           go to 105
         endif
      enddo
 
!     :::::::::: order eigenvalues and eigenvectors ::::::::::
      do i = 1, n-1
        k = i
        p = d(i)
 
! Find the minimum d at k
        do j = i+1, n
          if (d(j) .lt. p) then
            k = j
            p = d(j)
          endif
        enddo

! Swap d and z at index i with index k
        if (k .ne. i) then
          d(k) = d(i)
          d(i) = p

          p    = z(i)
          z(i) = z(k)
          z(k) = p
        endif
      enddo

#else

      integer i, j, k, l, m, n, ii, mml, IERR
      REALTYPE d(n), e(n), z(n), b, c, f, g, p, r, s, machep

      machep=EPSILON(ONE)
 
      IERR = 0
      if (n .eq. 1) return
 
      e(n) = ZERO
      do 240 l = 1, n
         j = 0
!     :::::::::: look for small sub-diagonal element ::::::::::
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if (ABS(e(m)) .le. machep * (ABS(d(m)) + ABS(d(m+1)))) go to 120 
  110    continue
!
  120    p = d(l)
         if (m .eq. l) go to 240

         if (j .eq. 100) then
           IERR = 1
           return
         endif

         j = j + 1
!     :::::::::: form shift ::::::::::
         g = (d(l+1) - p) / (TWO * e(l))
         r = SQRT(g*g+ONE)
         g = d(m) - p + e(l) / (g + SIGN(r, g))
         s = ONE
         c = ONE
         p = ZERO
         mml = m - l
 
!     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
         do ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (ABS(f) .lt. ABS(g)) go to 150
            c = g / f
            r = SQRT(c*c+ONE)
            e(i+1) = f * r
            s = ONE / r
            c = c * s
            go to 160
  150       s = f / g
            r = SQRT(s*s+ONE)
            e(i+1) = g * r
            c = ONE / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + TWO * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
!     :::::::::: form first component of vector ::::::::::
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
            z(i) = c * z(i) - s * f
         enddo
 
         d(l) = d(l) - p
         e(l) = g
         e(m) = ZERO
         go to 105
  240 continue
 
!     :::::::::: order eigenvalues and eigenvectors ::::::::::
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
 
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
 
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         p = z(i)
         z(i) = z(k)
         z(k) = p
  300 continue
 
#endif

      return
      end
#endif

END MODULE Gauss_Quadrature

#else

MODULE Gauss_Quadrature

implicit NONE

END MODULE Gauss_Quadrature

#endif
