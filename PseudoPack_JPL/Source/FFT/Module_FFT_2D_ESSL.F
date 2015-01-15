#include "pseudopack.h"


MODULE Two_Dimensional_FFT_ESSL

implicit NONE

#if defined (ESSL)
INTERFACE PS_Real_2_Complex_Transform
  MODULE PROCEDURE Real_2_Complex_Transform_2
  MODULE PROCEDURE Real_2_Complex_Transform_3
  MODULE PROCEDURE Real_2_Complex_Transform_4
END INTERFACE

INTERFACE PS_Complex_2_Real_Transform
  MODULE PROCEDURE Complex_2_Real_Transform_2
  MODULE PROCEDURE Complex_2_Real_Transform_3
  MODULE PROCEDURE Complex_2_Real_Transform_4
END INTERFACE

PRIVATE

PUBLIC :: PS_Real_2_Complex_Transform 
PUBLIC :: PS_Complex_2_Real_Transform 

CONTAINS
!
! ----------------------------------------------------------------------
!
  Subroutine Real_2_Complex_Transform_2 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  REALTYPE, dimension(:,:), TARGET  :: Q
  CMPLTYPE, dimension(:,:), TARGET  :: a

  REALTYPE, dimension(:,:), POINTER :: x
  CMPLTYPE, dimension(:,:), POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling

  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  x => Q(:,:)

  Nx = SIZE(x, DIM=1) ; Ny = SIZE(x, DIM=2)

  call Get_RCFT2_Storage_Size (Nx, Ny, N_WorkSpace, N_Aux2)

  Allocate (WorkSpace(N_WorkSpace), Aux2(N_Aux2), Aux3(1))

  x => Q(:,:) ; y => a(:,:)

  call Real_2_Complex_Transform (INIT, x, y, Nx, Ny, ISign, Scale, &
                                        WorkSpace, Aux2, Aux3)

  DEALLOCATE (WorkSpace, Aux2, Aux3)

  END Subroutine Real_2_Complex_Transform_2
!
! ----------------------------------------------------------------------
!
  Subroutine Real_2_Complex_Transform_3 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  REALTYPE, dimension(:,:,:), TARGET  :: Q
  CMPLTYPE, dimension(:,:,:), TARGET  :: a

  REALTYPE, dimension(:,:)  , POINTER :: x
  CMPLTYPE, dimension(:,:)  , POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling

  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  x => Q(:,:,1)

  Nx = SIZE(x, DIM=1) ; Ny = SIZE(x, DIM=2)

  call Get_RCFT2_Storage_Size (Nx, Ny, N_WorkSpace, N_Aux2)

  Allocate (WorkSpace(N_WorkSpace), Aux2(N_Aux2), Aux3(1))

    do n = 1,NV
      x => Q(:,:,n) ; y => a(:,:,n)

      call Real_2_Complex_Transform (INIT, x, y, Nx, Ny, ISign, Scale, &
                                            WorkSpace, Aux2, Aux3)
    enddo

  DEALLOCATE (WorkSpace, Aux2, Aux3)

  END Subroutine Real_2_Complex_Transform_3
!
! ----------------------------------------------------------------------
!
  Subroutine Real_2_Complex_Transform_4 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  REALTYPE, dimension(:,:,:,:), TARGET  :: Q
  CMPLTYPE, dimension(:,:,:,:), TARGET  :: a

  REALTYPE, dimension(:,:)    , POINTER :: x
  CMPLTYPE, dimension(:,:)    , POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling
  
  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  SELECT CASE (Index)
    CASE (1) 
      x => Q(1,:,:,1) 

    CASE (2) 
      x => Q(:,1,:,1) 

    CASE (3) 
      x => Q(:,:,1,1) 

  END SELECT

  Nx = SIZE(x, DIM=1) ; Ny = SIZE(x, DIM=2) 

  call Get_RCFT2_Storage_Size (Nx, Ny, N_WorkSpace, N_Aux2)

  Allocate (WorkSpace(N_WorkSpace), Aux2(N_Aux2), Aux3(1))
  
  do n = 1,NV
    do l = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
      SELECT CASE (Index)
        CASE (1) 
          x => Q(l,:,:,n) ; y => a(l,:,:,n)

        CASE (2) 
          x => Q(:,l,:,n) ; y => a(:,l,:,n)

        CASE (3) 
          x => Q(:,:,l,n) ; y => a(:,:,l,n)

      END SELECT

      call Real_2_Complex_Transform (INIT, x, y, Nx, Ny, ISign, Scale, &
                                            WorkSpace, Aux2, Aux3)
    enddo
  enddo

  DEALLOCATE (WorkSpace, Aux2, Aux3)

  END Subroutine Real_2_Complex_Transform_4
!
! ----------------------------------------------------------------------
!
  Subroutine Get_RCFT2_Storage_Size (Nx, Ny, N_WorkSpace, N_Aux2)

  implicit NONE

  integer  :: Nx, Ny, N_Workspace, N_Aux2

#if defined (SINGLE_PRECISION)
                            N_Workspace = 45000
  if (MAX(Nx/2,Ny) >  8192) N_WorkSpace = 40000+0.82*Nx+1.14*Ny

                                      N_Aux2 = 20000
  if ((Nx > 16384) .AND. (Ny <  252)) N_Aux2 = 20000+0.57*Nx
  if (Ny >= 252) N_Aux2 = N_Aux2 + (Ny+256)*(1.14+MIN(64,1+Nx/2))

#else
                            N_Workspace = 42000
  if (MAX(Nx/2,Ny) >  2048) N_WorkSpace = 40000+1.64*Nx+2.28*Ny

                                      N_Aux2 = 20000
  if ((Nx >  4096) .AND. (Ny <  252)) N_Aux2 = 20000+1.14*Nx
  if (Ny >= 252) N_Aux2 = N_Aux2 + (2*Ny+256)*(2.28+MIN(64,1+Nx/2))
#endif

  END Subroutine Get_RCFT2_Storage_Size 
!
! ----------------------------------------------------------------------
!
  Subroutine Real_2_Complex_Transform (INIT, x, y, Nx, Ny, ISign, Scale,  &
                                              WorkSpace, Aux2, Aux3)

  implicit NONE

  integer  :: INIT, Nx, Ny, ISign
  REALTYPE :: Scale

  REALTYPE, dimension(:,:) :: x
  CMPLTYPE, dimension(:,:) :: y
  REALTYPE, dimension(:)   :: WorkSpace, Aux2, Aux3

  if (INIT == 1) then
    call RCFT2 (INIT, x, SIZE(x,DIM=1),  y, SIZE(y,DIM=1),  Nx, Ny, &
                      ISign, Scale, WorkSpace, SIZE(WorkSpace),     &
                      Aux2, SIZE(Aux2), Aux3, SIZE(Aux3))

    INIT = 0
  endif

    call RCFT2 (INIT, x, SIZE(x,DIM=1),  y, SIZE(y,DIM=1),  Nx, Ny, &
                      ISign, Scale, WorkSpace, SIZE(WorkSpace),     &
                      Aux2, SIZE(Aux2), Aux3, SIZE(Aux3))

! WHERE (ABS(y) <= THOUSAND*EPSILON(ONE)) y = ZERO

 END Subroutine Real_2_Complex_Transform 
!
! ----------------------------------------------------------------------
!
  Subroutine Complex_2_Real_Transform_2 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  CMPLTYPE, dimension(:,:)    , TARGET  :: Q
  REALTYPE, dimension(:,:)    , TARGET  :: a

  CMPLTYPE, dimension(:,:)    , POINTER :: x
  REALTYPE, dimension(:,:)    , POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling

  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  x => Q(:,:)

  Nx = SIZE(x, DIM=1) ; Ny = SIZE(x, DIM=2)

  call Get_CRFT2_Storage_Size (Nx, Ny, N_WorkSpace, N_Aux2)

  Allocate (WorkSpace(N_WorkSpace), Aux2(N_Aux2), Aux3(1))

  x => Q(:,:) ; y => a(:,:)

  call Complex_2_Real_Transform (INIT, x, y, Nx, Ny, ISign, Scale, &
                                        WorkSpace, Aux2, Aux3)

  DEALLOCATE (WorkSpace, Aux2, Aux3)

  END Subroutine Complex_2_Real_Transform_2
!
! ----------------------------------------------------------------------
!
  Subroutine Complex_2_Real_Transform_3 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  CMPLTYPE, dimension(:,:,:)  , TARGET  :: Q
  REALTYPE, dimension(:,:,:)  , TARGET  :: a

  CMPLTYPE, dimension(:,:)    , POINTER :: x
  REALTYPE, dimension(:,:)    , POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling

  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  x => Q(:,:,1)

  Nx = SIZE(x, DIM=1) ; Ny = SIZE(x, DIM=2)

  call Get_CRFT2_Storage_Size (Nx, Ny, N_WorkSpace, N_Aux2)

  Allocate (WorkSpace(N_WorkSpace), Aux2(N_Aux2), Aux3(1))

    do n = 1,NV
      x => Q(:,:,n) ; y => a(:,:,n)

      call Complex_2_Real_Transform (INIT, x, y, Nx, Ny, ISign, Scale, &
                                            WorkSpace, Aux2, Aux3)
    enddo

  DEALLOCATE (WorkSpace, Aux2, Aux3)

  END Subroutine Complex_2_Real_Transform_3
!
! ----------------------------------------------------------------------
!
  Subroutine Complex_2_Real_Transform_4 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  CMPLTYPE, dimension(:,:,:,:), TARGET  :: Q
  REALTYPE, dimension(:,:,:,:), TARGET  :: a

  CMPLTYPE, dimension(:,:)    , POINTER :: x
  REALTYPE, dimension(:,:)    , POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling
  
  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  SELECT CASE (Index)
    CASE (1) 
      x => Q(1,:,:,1) 

    CASE (2) 
      x => Q(:,1,:,1) 

    CASE (3) 
      x => Q(:,:,1,1) 

  END SELECT

  Nx = SIZE(x, DIM=1) ; Ny = SIZE(x, DIM=2) 

  call Get_CRFT2_Storage_Size (Nx, Ny, N_WorkSpace, N_Aux2)

  Allocate (WorkSpace(N_WorkSpace), Aux2(N_Aux2), Aux3(1))
  
  do n = 1,NV
    do l = LBOUND(x,DIM=Index),UBOUND(x,DIM=Index)
      SELECT CASE (Index)
        CASE (1) 
          x => Q(l,:,:,n) ; y => a(l,:,:,n)

        CASE (2) 
          x => Q(:,l,:,n) ; y => a(:,l,:,n)

        CASE (3) 
          x => Q(:,:,l,n) ; y => a(:,:,l,n)

      END SELECT

      call Complex_2_Real_Transform (INIT, x, y, Nx, Ny, ISign, Scale, &
                                            WorkSpace, Aux2, Aux3)
    enddo
  enddo

  DEALLOCATE (WorkSpace, Aux2, Aux3)

  END Subroutine Complex_2_Real_Transform_4
!
! ----------------------------------------------------------------------
!
  Subroutine Get_CRFT2_Storage_Size (Nx, Ny, N_WorkSpace, N_Aux2)

  implicit NONE

  integer  :: Nx, Ny, N_Workspace, N_Aux2

#if defined (SINGLE_PRECISION)
                            N_Workspace = 45000
  if (MAX(Nx/2,Ny) >  8192) N_WorkSpace = 40000+0.82*Nx+1.14*Ny

                                      N_Aux2 = 20000
  if ((Nx > 16384) .AND. (Ny <  252)) N_Aux2 = 20000+0.57*Nx
  if (Ny >= 252) N_Aux2 = N_Aux2 + (Ny+256)*(1.14+MIN(64,1+Nx/2))

#else
                            N_Workspace = 42000
  if (MAX(Nx/2,Ny) >  2048) N_WorkSpace = 40000+1.64*Nx+2.28*Ny

                                      N_Aux2 = 20000
  if ((Nx >  4096) .AND. (Ny <  252)) N_Aux2 = 20000+1.14*Nx
  if (Ny >= 252) N_Aux2 = N_Aux2 + (2*Ny+256)*(2.28+MIN(64,1+Nx/2))
#endif

  END Subroutine Get_CRFT2_Storage_Size 
!
! ----------------------------------------------------------------------
!
  Subroutine Complex_2_Real_Transform (INIT, x, y, Nx, Ny, ISign, Scale,  &
                                              WorkSpace, Aux2, Aux3)

  implicit NONE

  integer  :: INIT, Nx, Ny, ISign
  REALTYPE :: Scale

  CMPLTYPE, dimension(:,:) :: x
  REALTYPE, dimension(:,:) :: y
  REALTYPE, dimension(:)   :: WorkSpace, Aux2, Aux3

  if (INIT == 1) then
    call CRFT2 (INIT, x, SIZE(x,DIM=1),  y, SIZE(y,DIM=1),  Nx, Ny, &
                      ISign, Scale, WorkSpace, SIZE(WorkSpace),     &
                      Aux2, SIZE(Aux2), Aux3, SIZE(Aux3))

    INIT = 0
  endif

    call CRFT2 (INIT, x, SIZE(x,DIM=1),  y, SIZE(y,DIM=1),  Nx, Ny, &
                      ISign, Scale, WorkSpace, SIZE(WorkSpace),     &
                      Aux2, SIZE(Aux2), Aux3, SIZE(Aux3))

! WHERE (ABS(y) <= THOUSAND*EPSILON(ONE)) y = ZERO

 END Subroutine Complex_2_Real_Transform 

#endif

END MODULE Two_Dimensional_FFT_ESSL

