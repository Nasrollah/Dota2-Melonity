#include "pseudopack.h"


MODULE Euler_Projection

implicit NONE

INTERFACE PS_Euler_Right_Projection
  MODULE PROCEDURE Euler_Right_Projection
END INTERFACE

INTERFACE PS_Euler_Left_Projection
  MODULE PROCEDURE Euler_Left_Projection
END INTERFACE

PRIVATE
PUBLIC  :: PS_Euler_Right_Projection
PUBLIC  :: PS_Euler_Left_Projection

CONTAINS

  Subroutine Euler_Left_Projection (Order, N0,N5, i0,i1, NV, m, E, F, G)

  integer  :: Order, N0,N5, i0,i1, NV, m

  REALTYPE, dimension(N0:N5,NV,NV)                       :: E
  REALTYPE, dimension(N0:N5,NV,2)                        :: F
  REALTYPE, dimension(N0:N5,-(Order-1)/2:(Order-1)/2,2)  :: G

  integer  :: i, k, M_s

  M_s = (Order-1)/2

  SELECT CASE (NV)
    CASE (1)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2)
        enddo
      enddo

    CASE (2)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1) + E(i,m, 2)*F(i+k  , 2,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2) + E(i,m, 2)*F(i-k+1, 2,2) 
        enddo
      enddo

    CASE (3)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1) + E(i,m, 2)*F(i+k  , 2,1) +  &
                     E(i,m, 3)*F(i+k  , 3,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2) + E(i,m, 2)*F(i-k+1, 2,2) +  &
                     E(i,m, 3)*F(i-k+1, 3,2)
        enddo
      enddo

    CASE (4)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1) + E(i,m, 2)*F(i+k  , 2,1) +  &
                     E(i,m, 3)*F(i+k  , 3,1) + E(i,m, 4)*F(i+k  , 4,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2) + E(i,m, 2)*F(i-k+1, 2,2) +  &
                     E(i,m, 3)*F(i-k+1, 3,2) + E(i,m, 4)*F(i-k+1, 4,2) 
        enddo
      enddo

    CASE (5)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1) + E(i,m, 2)*F(i+k  , 2,1) +  &
                     E(i,m, 3)*F(i+k  , 3,1) + E(i,m, 4)*F(i+k  , 4,1) +  &
                     E(i,m, 5)*F(i+k  , 5,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2) + E(i,m, 2)*F(i-k+1, 2,2) +  &
                     E(i,m, 3)*F(i-k+1, 3,2) + E(i,m, 4)*F(i-k+1, 4,2) +  &
                     E(i,m, 5)*F(i-k+1, 5,2) 
        enddo
      enddo

    CASE (6)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1) + E(i,m, 2)*F(i+k  , 2,1) +  &
                     E(i,m, 3)*F(i+k  , 3,1) + E(i,m, 4)*F(i+k  , 4,1) +  &
                     E(i,m, 5)*F(i+k  , 5,1) + E(i,m, 6)*F(i+k  , 6,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2) + E(i,m, 2)*F(i-k+1, 2,2) +  &
                     E(i,m, 3)*F(i-k+1, 3,2) + E(i,m, 4)*F(i-k+1, 4,2) +  &
                     E(i,m, 5)*F(i-k+1, 5,2) + E(i,m, 6)*F(i-k+1, 6,2) 
        enddo
      enddo

    CASE (7)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1) + E(i,m, 2)*F(i+k  , 2,1) +  &
                     E(i,m, 3)*F(i+k  , 3,1) + E(i,m, 4)*F(i+k  , 4,1) +  &
                     E(i,m, 5)*F(i+k  , 5,1) + E(i,m, 6)*F(i+k  , 6,1) +  &
                     E(i,m, 7)*F(i+k  , 7,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2) + E(i,m, 2)*F(i-k+1, 2,2) +  &
                     E(i,m, 3)*F(i-k+1, 3,2) + E(i,m, 4)*F(i-k+1, 4,2) +  &
                     E(i,m, 5)*F(i-k+1, 5,2) + E(i,m, 6)*F(i-k+1, 6,2) +  &
                     E(i,m, 7)*F(i-k+1, 7,2) 
        enddo
      enddo

    CASE (8)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1) + E(i,m, 2)*F(i+k  , 2,1) +  &
                     E(i,m, 3)*F(i+k  , 3,1) + E(i,m, 4)*F(i+k  , 4,1) +  &
                     E(i,m, 5)*F(i+k  , 5,1) + E(i,m, 6)*F(i+k  , 6,1) +  &
                     E(i,m, 7)*F(i+k  , 7,1) + E(i,m, 8)*F(i+k  , 8,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2) + E(i,m, 2)*F(i-k+1, 2,2) +  &
                     E(i,m, 3)*F(i-k+1, 3,2) + E(i,m, 4)*F(i-k+1, 4,2) +  &
                     E(i,m, 5)*F(i-k+1, 5,2) + E(i,m, 6)*F(i-k+1, 6,2) +  &
                     E(i,m, 7)*F(i-k+1, 7,2) + E(i,m, 8)*F(i-k+1, 8,2) 
        enddo
      enddo

    CASE (9)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1) + E(i,m, 2)*F(i+k  , 2,1) +  &
                     E(i,m, 3)*F(i+k  , 3,1) + E(i,m, 4)*F(i+k  , 4,1) +  &
                     E(i,m, 5)*F(i+k  , 5,1) + E(i,m, 6)*F(i+k  , 6,1) +  &
                     E(i,m, 7)*F(i+k  , 7,1) + E(i,m, 8)*F(i+k  , 8,1) +  &
                     E(i,m, 9)*F(i+k  , 9,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2) + E(i,m, 2)*F(i-k+1, 2,2) +  &
                     E(i,m, 3)*F(i-k+1, 3,2) + E(i,m, 4)*F(i-k+1, 4,2) +  &
                     E(i,m, 5)*F(i-k+1, 5,2) + E(i,m, 6)*F(i-k+1, 6,2) +  &
                     E(i,m, 7)*F(i-k+1, 7,2) + E(i,m, 8)*F(i-k+1, 8,2) +  &
                     E(i,m, 9)*F(i-k+1, 9,2) 
        enddo
      enddo

    CASE (10)
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = E(i,m, 1)*F(i+k  , 1,1) + E(i,m, 2)*F(i+k  , 2,1) +  &
                     E(i,m, 3)*F(i+k  , 3,1) + E(i,m, 4)*F(i+k  , 4,1) +  &
                     E(i,m, 5)*F(i+k  , 5,1) + E(i,m, 6)*F(i+k  , 6,1) +  &
                     E(i,m, 7)*F(i+k  , 7,1) + E(i,m, 8)*F(i+k  , 8,1) +  &
                     E(i,m, 9)*F(i+k  , 9,1) + E(i,m,10)*F(i+k  ,10,1)
          G(i,k,2) = E(i,m, 1)*F(i-k+1, 1,2) + E(i,m, 2)*F(i-k+1, 2,2) +  &
                     E(i,m, 3)*F(i-k+1, 3,2) + E(i,m, 4)*F(i-k+1, 4,2) +  &
                     E(i,m, 5)*F(i-k+1, 5,2) + E(i,m, 6)*F(i-k+1, 6,2) +  &
                     E(i,m, 7)*F(i-k+1, 7,2) + E(i,m, 8)*F(i-k+1, 8,2) +  &
                     E(i,m, 9)*F(i-k+1, 9,2) + E(i,m,10)*F(i-k+1,10,2) 
        enddo
      enddo

  END SELECT

  if (NV >  10) then
      do k = -M_s,M_s
        do i = i0,i1
          G(i,k,1) = G(i,k,1) + SUM(E(i,m,11:NV)*F(i+k  ,11:NV,1))
          G(i,k,2) = G(i,k,2) + SUM(E(i,m,11:NV)*F(i-k+1,11:NV,2))
        enddo
      enddo
  endif

  G(N0:i0-1,:,1) = ZERO ; G(i1+1:N5,:,1) = ZERO
  G(N0:i0-1,:,2) = ZERO ; G(i1+1:N5,:,2) = ZERO

! WHERE (ABS(G) <  1.0d-12) G = ZERO

  END Subroutine Euler_Left_Projection
!
! ======================================================================
!
  Subroutine Euler_Right_Projection (N0,N5, i0,i1, NV, E, F, G)

  integer  :: N0,N5, i0,i1, NV

  REALTYPE, dimension(N0:N5,NV,NV) :: E
  REALTYPE, dimension(N0:N5,NV)    :: F, G

  integer  :: i, m

  SELECT CASE (NV)
    CASE (1)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1)
        enddo
      enddo

    CASE (2)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1) + E(i,m, 2)*F(i, 2) 
        enddo
      enddo

    CASE (3)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1) + E(i,m, 2)*F(i, 2) +  &
                   E(i,m, 3)*F(i, 3)
        enddo
      enddo

    CASE (4)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1) + E(i,m, 2)*F(i, 2) +  &
                   E(i,m, 3)*F(i, 3) + E(i,m, 4)*F(i, 4) 
        enddo
      enddo

    CASE (5)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1) + E(i,m, 2)*F(i, 2) +  &
                   E(i,m, 3)*F(i, 3) + E(i,m, 4)*F(i, 4) +  &
                   E(i,m, 5)*F(i, 5) 
        enddo
      enddo

    CASE (6)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1) + E(i,m, 2)*F(i, 2) +  &
                   E(i,m, 3)*F(i, 3) + E(i,m, 4)*F(i, 4) +  &
                   E(i,m, 5)*F(i, 5) + E(i,m, 6)*F(i, 6) 
        enddo
      enddo

    CASE (7)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1) + E(i,m, 2)*F(i, 2) +  &
                   E(i,m, 3)*F(i, 3) + E(i,m, 4)*F(i, 4) +  &
                   E(i,m, 5)*F(i, 5) + E(i,m, 6)*F(i, 6) +  &
                   E(i,m, 7)*F(i, 7) 
        enddo
      enddo

    CASE (8)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1) + E(i,m, 2)*F(i, 2) +  &
                   E(i,m, 3)*F(i, 3) + E(i,m, 4)*F(i, 4) +  &
                   E(i,m, 5)*F(i, 5) + E(i,m, 6)*F(i, 6) +  &
                   E(i,m, 7)*F(i, 7) + E(i,m, 8)*F(i, 8) 
        enddo
      enddo

    CASE (9)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1) + E(i,m, 2)*F(i, 2) +  &
                   E(i,m, 3)*F(i, 3) + E(i,m, 4)*F(i, 4) +  &
                   E(i,m, 5)*F(i, 5) + E(i,m, 6)*F(i, 6) +  &
                   E(i,m, 7)*F(i, 7) + E(i,m, 8)*F(i, 8) +  &
                   E(i,m, 9)*F(i, 9) 
        enddo
      enddo

    CASE (10)
      do m = 1,NV
        do i = i0,i1
          G(i,m) = E(i,m, 1)*F(i, 1) + E(i,m, 2)*F(i, 2) +  &
                   E(i,m, 3)*F(i, 3) + E(i,m, 4)*F(i, 4) +  &
                   E(i,m, 5)*F(i, 5) + E(i,m, 6)*F(i, 6) +  &
                   E(i,m, 7)*F(i, 7) + E(i,m, 8)*F(i, 8) +  &
                   E(i,m, 9)*F(i, 9) + E(i,m,10)*F(i,10) 
        enddo
      enddo

  END SELECT

  if (NV >  10) then
      do m = 11,NV
        do i = i0,i1
          G(i,m) = G(i,m) + SUM(E(i,m,11:NV)*F(i,11:NV))
        enddo
      enddo
  endif

  G(N0:i0-1,:) = ZERO ; G(i1+1:N5,:) = ZERO

!  WHERE (ABS(G) <  1.0d-12) G = ZERO

  END Subroutine Euler_Right_Projection

END MODULE Euler_Projection
