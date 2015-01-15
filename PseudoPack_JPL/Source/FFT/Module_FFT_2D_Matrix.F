#include "pseudopack.h"

MODULE Two_Dimensional_FFT_MATRIX
 
  USE FFT_2D_Module

implicit NONE

INTERFACE PS_Real_2_Complex_Transform
  MODULE PROCEDURE Real_2_Complex_Transform_5
  MODULE PROCEDURE Real_2_Complex_Transform_6
  MODULE PROCEDURE Real_2_Complex_Transform_7
END INTERFACE

INTERFACE PS_Complex_2_Real_Transform
  MODULE PROCEDURE Complex_2_Real_Transform_5
  MODULE PROCEDURE Complex_2_Real_Transform_6
  MODULE PROCEDURE Complex_2_Real_Transform_7
END INTERFACE

PRIVATE

PUBLIC :: PS_Real_2_Complex_Transform 
PUBLIC :: PS_Complex_2_Real_Transform 

CONTAINS
!
! ----------------------------------------------------------------------
!
  Subroutine Real_2_Complex_Transform_5 (Index, NV, Q, a,                &
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

  x => Q(:,:) ; y => a(:,:)

  call PS_RCFT_2D (x, y, ISign, Scale)

  END Subroutine Real_2_Complex_Transform_5
!
! ----------------------------------------------------------------------
!
  Subroutine Real_2_Complex_Transform_6 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  REALTYPE, dimension(:,:,:), TARGET  :: Q
  CMPLTYPE, dimension(:,:,:), TARGET  :: a

  REALTYPE, dimension(:,:,:), POINTER :: x
  CMPLTYPE, dimension(:,:,:), POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling

  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  x => Q(:,:,:) ; y => a(:,:,:)

  call PS_RCFT_2D (x, y, Index, ISign, Scale)

  END Subroutine Real_2_Complex_Transform_6
!
! ----------------------------------------------------------------------
!
  Subroutine Real_2_Complex_Transform_7 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  REALTYPE, dimension(:,:,:,:), TARGET  :: Q
  CMPLTYPE, dimension(:,:,:,:), TARGET  :: a

  REALTYPE, dimension(:,:,:,:), POINTER :: x
  CMPLTYPE, dimension(:,:,:,:), POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling
  
  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  x => Q(:,:,:,:) ; y => a(:,:,:,:)

  call PS_RCFT_2D (x, y, Index, ISign, Scale)

  END Subroutine Real_2_Complex_Transform_7
!
! ----------------------------------------------------------------------
!
  Subroutine Complex_2_Real_Transform_5 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  CMPLTYPE, dimension(:,:), TARGET  :: Q
  REALTYPE, dimension(:,:), TARGET  :: a

  CMPLTYPE, dimension(:,:), POINTER :: x
  REALTYPE, dimension(:,:), POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling

  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  x => Q(:,:) ; y => a(:,:)

  call PS_CRFT_2D (x, y, ISign, Scale)

  END Subroutine Complex_2_Real_Transform_5
!
! ----------------------------------------------------------------------
!
  Subroutine Complex_2_Real_Transform_6 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  CMPLTYPE, dimension(:,:,:), TARGET  :: Q
  REALTYPE, dimension(:,:,:), TARGET  :: a

  CMPLTYPE, dimension(:,:,:), POINTER :: x
  REALTYPE, dimension(:,:,:), POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling

  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  x => Q(:,:,:) ; y => a(:,:,:)

  call PS_CRFT_2D (x, y, Index, ISign, Scale)

  END Subroutine Complex_2_Real_Transform_6
!
! ----------------------------------------------------------------------
!
  Subroutine Complex_2_Real_Transform_7 (Index, NV, Q, a,                &
                                                Initialize, Sign, Scaling)

  CMPLTYPE, dimension(:,:,:,:), TARGET  :: Q
  REALTYPE, dimension(:,:,:,:), TARGET  :: a

  CMPLTYPE, dimension(:,:,:,:), POINTER :: x
  REALTYPE, dimension(:,:,:,:), POINTER :: y

  integer , OPTIONAL :: Initialize, Sign
  REALTYPE, OPTIONAL :: Scaling
  
  integer  :: Index, NV
  integer  :: INIT, ISign, Nx, Ny, N_Workspace, N_Aux2, l, n
  REALTYPE :: Scale

  REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace, Aux2, Aux3

  INIT  =   1 ; if (PRESENT(Initialize)) INIT  = Initialize
  ISign =   1 ; if (PRESENT(Sign))       ISign = Sign
  Scale = ONE ; if (PRESENT(Scaling))    Scale = Scaling

  x => Q(:,:,:,:) ; y => a(:,:,:,:)

  call PS_CRFT_2D (x, y, Index, ISign, Scale)

  END Subroutine Complex_2_Real_Transform_7

END MODULE Two_Dimensional_FFT_MATRIX
