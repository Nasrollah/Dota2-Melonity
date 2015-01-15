#include "pseudopack.h"


MODULE FFT_2D_Module

implicit NONE

INTERFACE PS_RCFT_2D
  MODULE PROCEDURE  RCFT_2D_2
  MODULE PROCEDURE  RCFT_2D_3
  MODULE PROCEDURE  RCFT_2D_4
END INTERFACE

INTERFACE PS_CRFT_2D
  MODULE PROCEDURE  CRFT_2D_2
  MODULE PROCEDURE  CRFT_2D_3
  MODULE PROCEDURE  CRFT_2D_4
END INTERFACE

PRIVATE

PUBLIC  :: PS_RCFT_2D
PUBLIC  :: PS_CRFT_2D

CONTAINS
!
! -------------------------------------------------------------------
!
! x    real array of size (0:N1-1, 0:N2-1)    (0:31,0:1)  (N1 = 32 N2 = 2)
! y complex array of size (0:N1/2, 0:N2-1)    (0:16,0:1)

  Subroutine RCFT_2D_2 (x, y, ISign, Scale)

  REALTYPE, dimension(0:,0:) :: x
  CMPLTYPE, dimension(0:,0:) :: y

  integer            :: ISign
  REALTYPE, OPTIONAL :: Scale

  CMPLTYPE :: Sum, Arg_1, Arg_2, Arg_i, Arg_j, I
  integer  :: N1, N2, M1, M2
  integer  :: j1, j2
  integer  :: k1, k2

  N1 = SIZE(x,DIM=1) ; N2 = SIZE(x,DIM=2)
  M1 = N1/2+1        ; M2 = N2

  I = CMPLX(ZERO,ONE)

  Arg_1 = EXP(ISign*TWO*PI*I/N1)
  Arg_2 = EXP(ISign*TWO*PI*I/N2)

  do k2 = 0,M2-1
    do k1 = 0,M1-1
      
      Sum = CMPLX(ZERO, ZERO)
      do j2 = 0,N2-1
        Arg_j = Arg_2**(j2*k2)
        do j1 = 0,N1-1
          Arg_i = Arg_1**(j1*k1)

          Sum = Sum + x(j1,j2)*Arg_i*Arg_j
        enddo
      enddo

      y(k1,k2) = Sum
    enddo
  enddo

  if (PRESENT(Scale)) y = y*Scale
  
  END Subroutine RCFT_2D_2
!
! -------------------------------------------------------------------
!
  Subroutine RCFT_2D_3 (x, y, Index, ISign, Scale)

  REALTYPE, dimension(0:,0:,0:), TARGET  :: x
  CMPLTYPE, dimension(0:,0:,0:), TARGET  :: y

  REALTYPE, dimension(:,:)     , POINTER :: p_x
  CMPLTYPE, dimension(:,:)     , POINTER :: p_y

  integer            :: Index, ISign
  REALTYPE, OPTIONAL :: Scale

  integer            :: k

    SELECT CASE (Index)
      CASE (1)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call RCFT_2D_2 (x(k,:,:), y(k,:,:), ISign, Scale)
        enddo

      CASE (2)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call RCFT_2D_2 (x(:,k,:), y(:,k,:), ISign, Scale)
        enddo
        
      CASE (3)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call RCFT_2D_2 (x(:,:,k), y(:,:,k), ISign, Scale)
        enddo
        
    END SELECT

  END Subroutine RCFT_2D_3
!
! -------------------------------------------------------------------
!
  Subroutine RCFT_2D_4 (x, y, Index, ISign, Scale)

  REALTYPE, dimension(0:,0:,0:,:), TARGET  :: x
  CMPLTYPE, dimension(0:,0:,0:,:), TARGET  :: y

  REALTYPE, dimension(:,:)       , POINTER :: p_x
  CMPLTYPE, dimension(:,:)       , POINTER :: p_y

  integer            :: Index, ISign
  REALTYPE, OPTIONAL :: Scale

  integer            :: k, n

    SELECT CASE (Index)
      CASE (1)
        do n = LBOUND(x,DIM=4),UBOUND(x,DIM=4)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call RCFT_2D_2 (x(k,:,:,n), y(k,:,:,n), ISign, Scale)
        enddo
        enddo

      CASE (2)
        do n = LBOUND(x,DIM=4),UBOUND(x,DIM=4)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call RCFT_2D_2 (x(:,k,:,n), y(:,k,:,n), ISign, Scale)
        enddo
        enddo
        
      CASE (3)
        do n = LBOUND(x,DIM=4),UBOUND(x,DIM=4)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call RCFT_2D_2 (x(:,:,k,n), y(:,:,k,n), ISign, Scale)
        enddo
        enddo
        
    END SELECT

  END Subroutine RCFT_2D_4
!
! -------------------------------------------------------------------
!
! x  complex array of size (0:N1/2,0:N2-1)   (0:16,0:1)  N1 = 33, N2 = 2
! y     real array of size (0:N1-1,0:N2-1)   (0:31,0:1)

  Subroutine CRFT_2D_2 (x, y, ISign, Scale)

  CMPLTYPE, dimension(0:,0:) :: x
  REALTYPE, dimension(0:,0:) :: y

  integer            :: ISign
  REALTYPE, OPTIONAL :: Scale

  CMPLTYPE :: Sum, Arg_1, Arg_2, Arg_i, Arg_j, I, Sum_k
  integer  :: N1, N2, M1, M2
  integer  :: j1, j2
  integer  :: k1, k2

  N1 = 2*(SIZE(x,DIM=1)-1) ; N2 = SIZE(x,DIM=2)
  M1 = N1/2+1              ; M2 = N2

  I = CMPLX(ZERO,ONE)

  Arg_1 = EXP(-ISign*TWO*PI*I/N1)
  Arg_2 = EXP(-ISign*TWO*PI*I/N2)

  do j2 = 0,N2-1
    do j1 = 0,N1-1
      
      Sum = CMPLX(ZERO,ZERO)
      do k2 = 0,M2-1
        Arg_j = Arg_2**(j2*k2)

        Sum_k = (x(0,k2) + x(M1-1,k2)*(-1)**(j1))/TWO
        do k1 = 1,M1-2
          Arg_i = Arg_1**(j1*k1)

          Sum_k = Sum_k + x(k1,k2)*Arg_i 
        enddo
        Sum = Sum + Sum_k*Arg_j
      enddo

      y(j1,j2) = TWO*Sum

    enddo
  enddo

  if (PRESENT(Scale)) y = y*Scale
  
  END Subroutine CRFT_2D_2
!
! -------------------------------------------------------------------
!
  Subroutine CRFT_2D_3 (x, y, Index, ISign, Scale)

  CMPLTYPE, dimension(0:,0:,0:), TARGET  :: x
  REALTYPE, dimension(0:,0:,0:), TARGET  :: y

  CMPLTYPE, dimension(:,:)     , POINTER :: p_x
  REALTYPE, dimension(:,:)     , POINTER :: p_y

  integer            :: Index, ISign
  REALTYPE, OPTIONAL :: Scale

  integer            :: k

    SELECT CASE (Index)
      CASE (1)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call CRFT_2D_2 (x(k,:,:), y(k,:,:), ISign, Scale)
        enddo

      CASE (2)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call CRFT_2D_2 (x(:,k,:), y(:,k,:), ISign, Scale)
        enddo
        
      CASE (3)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call CRFT_2D_2 (x(:,:,k), y(:,:,k), ISign, Scale)
        enddo
        
    END SELECT

  END Subroutine CRFT_2D_3
!
! -------------------------------------------------------------------
!
  Subroutine CRFT_2D_4 (x, y, Index, ISign, Scale)

  CMPLTYPE, dimension(0:,0:,0:,:), TARGET  :: x
  REALTYPE, dimension(0:,0:,0:,:), TARGET  :: y

  CMPLTYPE, dimension(:,:)       , POINTER :: p_x
  REALTYPE, dimension(:,:)       , POINTER :: p_y

  integer            :: Index, ISign
  REALTYPE, OPTIONAL :: Scale

  integer            :: k, n

    SELECT CASE (Index)
      CASE (1)
        do n = LBOUND(x,DIM=4),UBOUND(x,DIM=4)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call CRFT_2D_2 (x(k,:,:,n), y(k,:,:,n), ISign, Scale)
        enddo
        enddo

      CASE (2)
        do n = LBOUND(x,DIM=4),UBOUND(x,DIM=4)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call CRFT_2D_2 (x(:,k,:,n), y(:,k,:,n), ISign, Scale)
        enddo
        enddo
        
      CASE (3)
        do n = LBOUND(x,DIM=4),UBOUND(x,DIM=4)
        do k = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
          call CRFT_2D_2 (x(:,:,k,n), y(:,:,k,n), ISign, Scale)
        enddo
        enddo
        
    END SELECT

  END Subroutine CRFT_2D_4

END MODULE FFT_2D_Module
