!----------
#if 1
!----------

#include "pseudopack.h"


MODULE Specific_Heat_Gamma

IMPLICIT NONE

INTERFACE PS_Specific_Heat_Gamma
  MODULE PROCEDURE Specific_Heat_Gamma_0D 
  MODULE PROCEDURE Specific_Heat_Gamma_1D 
  MODULE PROCEDURE Specific_Heat_Gamma_2D
  MODULE PROCEDURE Specific_Heat_Gamma_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Specific_Heat_Gamma
PUBLIC :: PS_Gamma_Prime

CONTAINS
!
!==========================================================================
!
  Subroutine Specific_Heat_Gamma_0D (Formulation, Q, Gamma_Gas, Gamma, NF, NS)

  integer  :: Formulation, NF, NS

  REALTYPE, dimension(NF+NS)       :: Q
  REALTYPE, dimension(:,:)         :: Gamma_Gas
  REALTYPE                         :: Gamma

  integer  :: n

  REALTYPE, dimension(NS) :: f
  REALTYPE                :: f_Cp, f_Cv

  SELECT CASE (Formulation)
    CASE DEFAULT
      do n = 1, NS-1
       f(n) = Q(n+NF)/Q(1)
     enddo

     f(1:NS-1) = MIN(MAX(ZERO,f(1:NS-1)),ONE)
     f(NS    ) = ONE -    SUM(f(1:NS-1))

        f_Cp = SUM(f(1:NS)*Gamma_Gas(1:NS,2))
        f_Cv = SUM(f(1:NS)*Gamma_Gas(1:NS,3))

      Gamma = f_Cp/f_Cv

    CASE (2)
      Gamma = Gamma_Gas(1,1)

  END SELECT

  END Subroutine Specific_Heat_Gamma_0D
!
!==========================================================================
!
  Subroutine Specific_Heat_Gamma_1D (Formulation, Q, Gamma_Gas, Gamma)

  integer  :: Formulation

  REALTYPE, dimension(:,:)           :: Q
  REALTYPE, dimension(:,:)           :: Gamma_Gas
  REALTYPE, dimension(:)             :: Gamma

  integer  :: NF, NV, NS, n, i

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2)-2) :: f
  REALTYPE, dimension(SIZE(Q,DIM=1)) :: f_Cp, f_Cv

  SELECT CASE (Formulation)
    CASE DEFAULT
      NF = 3 ; NV = SIZE(Q,DIM=2) ; NS = NV-NF+1

      do n = 1, NS-1
        f(:,n) = Q(:,n+NF)/Q(:,1)
      enddo

      f(:,1:NS-1) = MIN(MAX(ZERO,f(:,1:NS-1)),ONE)
      f(:,NS    ) = ONE -    SUM(f(:,1:NS-1), DIM=2)

      do i = 1,SIZE(Q,DIM=1)
        f_Cp(i) = SUM(f(i,1:NS)*Gamma_Gas(1:NS,2))
        f_Cv(i) = SUM(f(i,1:NS)*Gamma_Gas(1:NS,3))
      enddo

      Gamma = f_Cp/f_Cv

    CASE (2)
      Gamma = Gamma_Gas(1,1)

  END SELECT

  END Subroutine Specific_Heat_Gamma_1D 
!
!==========================================================================
!
  Subroutine Specific_Heat_Gamma_2D (Formulation, Q, Gamma_Gas, Gamma)

  integer  :: Formulation

  REALTYPE, dimension(:,:,:) :: Q
  REALTYPE, dimension(:,:)   :: Gamma_Gas
  REALTYPE, dimension(:,:)   :: Gamma

  integer  :: NF, NV, NS, n, i, j

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3)-3) :: f
  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2))                 :: f_Cp, f_Cv

  SELECT CASE (Formulation)
    CASE DEFAULT
      NF = 4 ; NV = SIZE(Q,DIM=3) ; NS = NV-NF+1

      do n = 1, NS-1
        f(:,:,n) = Q(:,:,n+NF)/Q(:,:,1)
      enddo

      f(:,:,1:NS-1) = MIN(MAX(ZERO,f(:,:,1:NS-1)),ONE)
      f(:,:,NS    ) = ONE -    SUM(f(:,:,1:NS-1), DIM=3)

      do j = 1,SIZE(Q,DIM=2)
        do i = 1,SIZE(Q,DIM=1)
          f_Cp(i,j) = SUM(f(i,j,1:NS)*Gamma_Gas(1:NS,2))
          f_Cv(i,j) = SUM(f(i,j,1:NS)*Gamma_Gas(1:NS,3))
        enddo
      enddo

      Gamma = f_Cp/f_Cv

    CASE (2)
      Gamma = Gamma_Gas(1,1)

  END SELECT

  END Subroutine Specific_Heat_Gamma_2D 
!
!==========================================================================
!
  Subroutine Specific_Heat_Gamma_3D (Formulation, Q, Gamma_Gas, Gamma)

  integer  :: Formulation

  REALTYPE, dimension(:,:,:,:) :: Q
  REALTYPE, dimension(:,:)     :: Gamma_Gas
  REALTYPE, dimension(:,:,:)   :: Gamma

  integer  :: NF, NV, NS, n, i, j, k

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3),          &
                                                  SIZE(Q,DIM=4)-4) :: f
  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3))   :: f_Cp, f_Cv

  SELECT CASE (Formulation)
    CASE DEFAULT
      NF = 5 ; NV = SIZE(Q,DIM=4) ; NS = NV-NF+1

      do n = 1, NS-1
        f(:,:,:,n) = Q(:,:,:,n+NF)/Q(:,:,:,1)
      enddo

      f(:,:,:,1:NS-1) = MIN(MAX(ZERO,f(:,:,:,1:NS-1)),ONE)
      f(:,:,:,NS    ) = ONE -    SUM(f(:,:,:,1:NS-1), DIM=4)

      do k = 1,SIZE(Q,DIM=3)
        do j = 1,SIZE(Q,DIM=2)
          do i = 1,SIZE(Q,DIM=1)
            f_Cp(i,j,k) = SUM(f(i,j,k,1:NS)*Gamma_Gas(1:NS,2))
            f_Cv(i,j,k) = SUM(f(i,j,k,1:NS)*Gamma_Gas(1:NS,3))
          enddo
        enddo
      enddo

      Gamma = f_Cp/f_Cv

    CASE (2)
      Gamma = Gamma_Gas(1,1)

  END SELECT

  END Subroutine Specific_Heat_Gamma_3D 
!
!==========================================================================
!
  REALTYPE FUNCTION PS_Gamma_Prime (Formulation, n, f, Gamma_Gas) 

  integer                  :: Formulation, n
  REALTYPE, dimension(:)   :: f
  REALTYPE, dimension(:,:) :: Gamma_Gas
  REALTYPE                 :: Gamma_Prime

  REALTYPE :: f_Cp, F_Cv

  SELECT CASE (Formulation)
    CASE DEFAULT
      f_Cp = SUM(f*Gamma_Gas(:,2))
      f_Cv = SUM(f*Gamma_Gas(:,3))

      Gamma_Prime = (Gamma_Gas(n,2)*f_Cv - Gamma_Gas(n,3)*f_Cp)/f_Cv**2

    CASE (2)
      Gamma_Prime = ZERO

  END SELECT 

  PS_Gamma_Prime = Gamma_Prime

  END FUNCTION PS_Gamma_Prime

END MODULE Specific_Heat_Gamma

!----------
#else
!----------

#include "pseudopack.h"


MODULE Specific_Heat_Gamma

IMPLICIT NONE

INTERFACE PS_Specific_Heat_Gamma
  MODULE PROCEDURE Specific_Heat_Gamma_0D 
  MODULE PROCEDURE Specific_Heat_Gamma_1D 
  MODULE PROCEDURE Specific_Heat_Gamma_2D
  MODULE PROCEDURE Specific_Heat_Gamma_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Specific_Heat_Gamma
PUBLIC :: PS_Gamma_Prime

CONTAINS
!
!==========================================================================
!
  Subroutine Specific_Heat_Gamma_0D (Formulation, Q, Gamma_Gas, Gamma, NF, NS)

  integer  :: Formulation, NF, NS

  REALTYPE, dimension(NF+NS)       :: Q
  REALTYPE, dimension(:,:)         :: Gamma_Gas
  REALTYPE                         :: Gamma

  integer  :: n

  REALTYPE, dimension(NS) :: f
  REALTYPE                :: f_Cp, f_Cv

  do n = 1, NS-1
    f(n) = Q(n+NF)/Q(1)
  enddo

  f(1:NS-1) = MIN(MAX(ZERO,f(1:NS-1)),ONE)
  f(NS    ) = ONE -    SUM(f(1:NS-1))

  SELECT CASE (Formulation)
    CASE (0)
      Gamma = Gamma

    CASE (1)
        f_Cp = SUM(f(1:NS)*Gamma_Gas(1:NS,2))
        f_Cv = SUM(f(1:NS)*Gamma_Gas(1:NS,3))

      Gamma = f_Cp/f_Cv

  END SELECT

  END Subroutine Specific_Heat_Gamma_0D
!
!==========================================================================
!
  Subroutine Specific_Heat_Gamma_1D (Formulation, Q, Gamma_Gas, Gamma)

  integer  :: Formulation

  REALTYPE, dimension(:,:)           :: Q
  REALTYPE, dimension(:,:)           :: Gamma_Gas
  REALTYPE, dimension(:)             :: Gamma

  integer  :: NF, NV, NS, n, i

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2)-2) :: f
  REALTYPE, dimension(SIZE(Q,DIM=1)) :: f_Cp, f_Cv

  NF = 3 ; NV = SIZE(Q,DIM=2) ; NS = NV-NF+1

  do n = 1, NS-1
    f(:,n) = Q(:,n+NF)/Q(:,1)
  enddo

  f(:,1:NS-1) = MIN(MAX(ZERO,f(:,1:NS-1)),ONE)
  f(:,NS    ) = ONE -    SUM(f(:,1:NS-1), DIM=2)

  SELECT CASE (Formulation)
    CASE DEFAULT
      do i = 1,SIZE(Q,DIM=1)
        f_Cp(i) = SUM(f(i,1:NS)*Gamma_Gas(1:NS,2))
        f_Cv(i) = SUM(f(i,1:NS)*Gamma_Gas(1:NS,3))
      enddo

      Gamma = f_Cp/f_Cv

    CASE (2)
      Gamma = Gamma

  END SELECT

  END Subroutine Specific_Heat_Gamma_1D 
!
!==========================================================================
!
  Subroutine Specific_Heat_Gamma_2D (Formulation, Q, Gamma_Gas, Gamma)

  integer  :: Formulation

  REALTYPE, dimension(:,:,:) :: Q
  REALTYPE, dimension(:,:)   :: Gamma_Gas
  REALTYPE, dimension(:,:)   :: Gamma

  integer  :: NF, NV, NS, n, i, j

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3)-3) :: f
  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2))                 :: f_Cp, f_Cv

  NF = 4 ; NV = SIZE(Q,DIM=3) ; NS = NV-NF+1

  do n = 1, NS-1
    f(:,:,n) = Q(:,:,n+NF)/Q(:,:,1)
  enddo

  f(:,:,1:NS-1) = MIN(MAX(ZERO,f(:,:,1:NS-1)),ONE)
  f(:,:,NS    ) = ONE -    SUM(f(:,:,1:NS-1), DIM=3)

  SELECT CASE (Formulation)
    CASE DEFAULT
      do j = 1,SIZE(Q,DIM=2)
        do i = 1,SIZE(Q,DIM=1)
          f_Cp(i,j) = SUM(f(i,j,1:NS)*Gamma_Gas(1:NS,2))
          f_Cv(i,j) = SUM(f(i,j,1:NS)*Gamma_Gas(1:NS,3))
        enddo
      enddo

      Gamma = f_Cp/f_Cv

    CASE (2)
      Gamma = Gamma

  END SELECT

  END Subroutine Specific_Heat_Gamma_2D 
!
!==========================================================================
!
  Subroutine Specific_Heat_Gamma_3D (Formulation, Q, Gamma_Gas, Gamma)

  integer  :: Formulation

  REALTYPE, dimension(:,:,:,:) :: Q
  REALTYPE, dimension(:,:)     :: Gamma_Gas
  REALTYPE, dimension(:,:,:)   :: Gamma

  integer  :: NF, NV, NS, n, i, j, k

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3),          &
                                                  SIZE(Q,DIM=4)-4) :: f
  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3))   :: f_Cp, f_Cv

  NF = 5 ; NV = SIZE(Q,DIM=4) ; NS = NV-NF+1

  do n = 1, NS-1
    f(:,:,:,n) = Q(:,:,:,n+NF)/Q(:,:,:,1)
  enddo

  f(:,:,:,1:NS-1) = MIN(MAX(ZERO,f(:,:,:,1:NS-1)),ONE)
  f(:,:,:,NS    ) = ONE -    SUM(f(:,:,:,1:NS-1), DIM=4)

  SELECT CASE (Formulation)
    CASE DEFAULT
      do k = 1,SIZE(Q,DIM=3)
        do j = 1,SIZE(Q,DIM=2)
          do i = 1,SIZE(Q,DIM=1)
            f_Cp(i,j,k) = SUM(f(i,j,k,1:NS)*Gamma_Gas(1:NS,2))
            f_Cv(i,j,k) = SUM(f(i,j,k,1:NS)*Gamma_Gas(1:NS,3))
          enddo
        enddo
      enddo

      Gamma = f_Cp/f_Cv

    CASE (2)
      Gamma = Gamma

  END SELECT

  END Subroutine Specific_Heat_Gamma_3D 
!
!==========================================================================
!
  REALTYPE FUNCTION PS_Gamma_Prime (Formulation, n, f, Gamma_Gas) 

  integer                  :: Formulation, n
  REALTYPE, dimension(:)   :: f
  REALTYPE, dimension(:,:) :: Gamma_Gas
  REALTYPE                 :: Gamma_Prime

  REALTYPE :: f_Cp, F_Cv

  SELECT CASE (Formulation)
    CASE DEFAULT
      f_Cp = SUM(f*Gamma_Gas(:,2))
      f_Cv = SUM(f*Gamma_Gas(:,3))

      Gamma_Prime = (Gamma_Gas(n,2)*f_Cv - Gamma_Gas(n,3)*f_Cp)/f_Cv**2

    CASE (2)
      Gamma_Prime = ZERO

  END SELECT 

  PS_Gamma_Prime = Gamma_Prime

  END FUNCTION PS_Gamma_Prime

END MODULE Specific_Heat_Gamma

!----------
#endif
!----------
