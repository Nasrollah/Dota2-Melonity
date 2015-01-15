#include "pseudopack.h"


#if 1
! ----------------------------------------------------------------------
! FUNCTION NAME: User_Defined_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_f_Map_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE User_Defined_Mapping
  USE Grid_Auxillary_Subroutine, ONLY: Map_Computational_2_Physical, &
                                       Map_Physical_2_Computational

  USE User_Defined_Map_Function

implicit NONE

PRIVATE

PUBLIC :: User_Defined_Map
PUBLIC ::         Map_User_Defined
PUBLIC :: Inverse_Map_User_Defined

CONTAINS

  Subroutine User_Defined_Map                                        &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_f_Map_x_1, D_f_Map_x_2, D_f_Map_x_3, D_f_Map_x_4)

  integer  :: i
  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_f_Map_x_1, D_f_Map_x_2,   &
                               D_f_Map_x_3, D_f_Map_x_4

  do i = 0,N0
    call D_f_Map_x (x(i), x_Left, x_Right, alpha, beta,  &
                    D_f_Map_x_1(i), D_f_Map_x_2(i),      &
                    D_f_Map_x_3(i), D_f_Map_x_4(i), Map)
  enddo

  if (Map == 0) RETURN

  do i = 0,N2
    x(i) = f_Map_x (x(i), x_Left, x_Right, alpha, beta)
  enddo

  END Subroutine User_Defined_Map

! ----------------------------------------------------------------------
! FUNCTION NAME: Map_User_Defined
! DESCRIPTION  : Map the Chebyshev Point to the 
!                  User Defined Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_User_Defined (D, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: D

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right
  REALTYPE :: xx

  Map     = OPERATOR_MAP(D)
  alpha   = OPERATOR_ALPHA(D)
  beta    = OPERATOR_BETA(D)
  x_Left  = OPERATOR_X_LEFT(D)
  x_Right = OPERATOR_X_RIGHT(D)

  if (Map == 0) then
    Map_User_Defined = x ; RETURN
  endif

  xx = f_Map_x (x, x_Left, x_Right, alpha, beta)

  call Map_Computational_2_Physical (1,   -ONE,    ONE, xx, &
                                        x_Left,x_Right, xx)

  Map_User_Defined = xx

  END FUNCTION Map_User_Defined 

! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_User_Defined
! DESCRIPTION  : Map the User Defined Mapped Physical Point back to 
!                 Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_User_Defined  (D, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: D

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right
  REALTYPE :: xx

  Map     = OPERATOR_MAP(D)
  alpha   = OPERATOR_ALPHA(D)
  beta    = OPERATOR_BETA(D)
  x_Left  = OPERATOR_X_LEFT(D)
  x_Right = OPERATOR_X_RIGHT(D)

#if defined (DEBUG)
  if (((x_Left <  x_Right) .AND. ((x <  x_Left) .OR. (x >  x_Right))) .OR. &
      ((x_Left >  x_Right) .AND. ((x >  x_Left) .OR. (x <  x_Right)))) then
    write (6 ,100) x, x_Left, x_Right ; write (99,100) x, x_Left, x_Right ; STOP
  endif

 100 format (1x,'Inverse_Map_User_Defined : '           / &
             1x,'   x is out of range [x_Left,x_Right]' / &
             1x,'   x, x_Left, x_Right = ',3(1x,f14.7))
#endif

  if (Map == 0) then
    Inverse_Map_User_Defined = x ; RETURN
  endif

  call Map_Physical_2_Computational (1, x_Left, x_Right,  x, &
                                          -ONE,     ONE, xx)

  xx = Inverse_f_Map_x (x, x_Left, x_Right, alpha, beta)

  Inverse_Map_User_Defined = xx

  END FUNCTION Inverse_Map_User_Defined 

END MODULE User_Defined_Mapping

#else

! ----------------------------------------------------------------------
! FUNCTION NAME: User_Defined_Map
! DESCRIPTION  : Compute the Mapped grid Points and its four transformation
!                metrics D_f_Map_x_k, k=1,2,3,4.
! ----------------------------------------------------------------------

MODULE User_Defined_Mapping
  USE Grid_Auxillary_Subroutine, ONLY: Map_Computational_2_Physical, &
                                       Map_Physical_2_Computational


implicit NONE

PRIVATE

PUBLIC :: User_Defined_Map
PUBLIC ::         Map_User_Defined
PUBLIC :: Inverse_Map_User_Defined

CONTAINS

  Subroutine User_Defined_Map                                        &
             (N, N0, N2, Map, Manual, alpha, beta, x_Left, x_Right,  &
              x, D_f_Map_x_1, D_f_Map_x_2, D_f_Map_x_3, D_f_Map_x_4)

  integer  :: N, N0, N2, Map, Manual
  REALTYPE :: alpha, beta, x_Left, x_Right

  REALTYPE, dimension(0:N2) :: x
  REALTYPE, dimension(0:N0) :: D_f_Map_x_1, D_f_Map_x_2,   &
                               D_f_Map_x_3, D_f_Map_x_4

  integer  :: i

  if (Map == 0) RETURN
  
  do i = 0,N0
    D_f_Map_x_1(i) = ONE
    D_f_Map_x_2(i) = ZERO
    D_f_Map_x_3(i) = ZERO
    D_f_Map_x_4(i) = ZERO
  enddo

  x = x

  END Subroutine User_Defined_Map

! ----------------------------------------------------------------------
! FUNCTION NAME: Map_User_Defined
! DESCRIPTION  : Map the Chebyshev Point to the 
!                  User Defined Mapped Physical Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Map_User_Defined (D, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: D

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right
  REALTYPE :: xx

  Map     = OPERATOR_MAP(D)
  alpha   = OPERATOR_ALPHA(D)
  beta    = OPERATOR_BETA(D)
  x_Left  = OPERATOR_X_LEFT(D)
  x_Right = OPERATOR_X_RIGHT(D)

  if (Map == 0) then
    Map_User_Defined = x ; RETURN
  endif

  xx = x

  call Map_Computational_2_Physical (1,   -ONE,    ONE, xx, &
                                        x_Left,x_Right, xx)

  Map_User_Defined = xx

  END FUNCTION Map_User_Defined 

! ----------------------------------------------------------------------
! FUNCTION NAME: Inverse_Map_User_Defined
! DESCRIPTION  : Map the User Defined Mapped Physical Point back to 
!                 Chebyshev Point
! ----------------------------------------------------------------------
  REALTYPE FUNCTION Inverse_Map_User_Defined  (D, x)

  REALTYPE               :: x
  REALTYPE, dimension(*) :: D

  integer  :: Map
  REALTYPE :: alpha, beta, x_Left, x_Right
  REALTYPE :: xx

  Map     = OPERATOR_MAP(D)
  alpha   = OPERATOR_ALPHA(D)
  beta    = OPERATOR_BETA(D)
  x_Left  = OPERATOR_X_LEFT(D)
  x_Right = OPERATOR_X_RIGHT(D)

#if defined (DEBUG)
  if (((x_Left <  x_Right) .AND. ((x <  x_Left) .OR. (x >  x_Right))) .OR. &
      ((x_Left >  x_Right) .AND. ((x >  x_Left) .OR. (x <  x_Right)))) then
    write (6 ,100) x, x_Left, x_Right ; write (99,100) x, x_Left, x_Right ; STOP
  endif

 100 format (1x,'Inverse_Map_User_Defined : x is out of range [x_Left,x_Right]'/ &
             1x,'   x, x_Left, x_Right = ',3(1x,f14.7))
#endif

  if (Map == 0) then
    Inverse_Map_User_Defined = x ; RETURN
  endif

  call Map_Physical_2_Computational (1, x_Left, x_Right,  x, &
                                          -ONE,     ONE, xx)

  xx = x

  Inverse_Map_User_Defined = xx

  END FUNCTION Inverse_Map_User_Defined 

END MODULE User_Defined_Mapping
#endif
