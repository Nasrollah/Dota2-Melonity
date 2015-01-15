#include "pseudopack.h"


! --------------------------------------------------------
! Date    : 2-28-05
! Version : 1.0
! 
! Note    : Should be correct for both even and odd MR and even and odd 
!             number of grid points.
!           Used the same coefficients beta and renamed the beta to a_x .
!           Two dimensional array is added.
!
! Date    : 3-5-05
! Version : 1.1
!
! Note    : Fix up the wrong order of coefficients for the Odd order MR.
!         : Enforce sysmetry for the Even order MR.
!         : Use index shift to simplify the Even and Odd order MR coding.
!         : Replace the Check_Parameter_Even and Check_Parameter_Odd with
!             Check_Parameter that will deal with both Even and Odd MR.
!         : Same fixes applied to 2D.
!
! Date    : 3-7-05
! Version : 1.2
!           Improve performance by removing the SUM with an explicit expression
!           without loops.  
!           Take advantage of the symmetry of the Even order MR.
!           Take advantage of the special form for the computing 
!           SUM(Beta*F) of the Odd order MR at ii+1
!           Similar fixes applied to 2D
!
! Date    : 3-8-05
! Version : 1.3
!           Added 3D MR
!
! Date    : 3-15-05
! Version : 1.3
!           Moved the argument Scale in the argument list to the end
!           Corrected MR_Scale from integer to REALTYPE
!
! Date    : 5-3-05
! Version : 1.3a
!           Complete overhaul of the subroutine arguments list to avoid
!           compiler errors on the SUN workstation.   List of arguments
!           are modified to reflect a more efficient usage of the routines. 
!
!           Q is of size (N0:N5,M0:M5,K0:K5,NV) with 3 Direction (1,2,3) and
!             NV number of Index (1,..,NV).  
! --------------------------------------------------------

Module Multi_Resolution

IMPLICIT NONE   

Interface PS_MR_Analysis
  Module Procedure Multi_Resolution_0

  Module Procedure Multi_Resolution_1a
  Module Procedure Multi_Resolution_2a
  Module Procedure Multi_Resolution_3a

  Module Procedure Multi_Resolution_1D
  Module Procedure Multi_Resolution_2D
  Module Procedure Multi_Resolution_3D
END Interface

integer                :: MR_Activation = 0
integer , dimension(3) :: MR_Activate   = 0

integer  :: MR_Default     = 0

integer  :: MR_Order       = 2
integer  :: MR_Index       = 1
integer  :: MR_Level       = 1
REALTYPE :: MR_Tolerance   = 1.0d-3
REALTYPE :: MR_Jump_Factor = ONE
integer  :: MR_i0          = 0
integer  :: MR_i1          = 0

integer  :: MR_Initialization = 1
integer  :: MR_Output_Init    = 1

character(LEN=256) :: MR_Flag_Filename
character(LEN=256) :: MR_Input_Filename

integer , dimension(0:3) :: MR_WENO
REALTYPE, dimension(0:3) :: MR_Ratio

PRIVATE

PUBLIC  :: PS_WENO_MR_Index_Update
PUBLIC  :: PS_MR_Analysis
PUBLIC  :: PS_MR_Adjust_Ghost_Index 

PUBLIC  :: MR_Activation, MR_Activate, MR_Default
PUBLIC  :: MR_Order, MR_Index, MR_Level, MR_Jump_Factor, MR_Tolerance
PUBLIC  :: MR_i0, MR_i1
PUBLIC  :: MR_Initialization
PUBLIC  :: MR_Input_Filename, MR_Output_Init, MR_Flag_Filename
PUBLIC  :: MR_WENO, MR_Ratio

CONTAINS
!
! ======================================================================
!
  Subroutine Multi_Resolution_0  (Order_x, N0, N5, i0, i1,                    &
                                                        Q, Flag,              &
                                  Tolerance, Jump_Factor, Scale, Default_Huge )

  integer            :: Order_x, N0, N5, i0, i1
  REALTYPE, OPTIONAL :: Tolerance, Jump_Factor
  integer , OPTIONAL :: Default_Huge

  REALTYPE, dimension(N0:N5)           :: Q
  integer , dimension(N0:N5)           :: Flag
  REALTYPE, dimension(N0:N5), OPTIONAL :: Scale

  call Multi_Resolution_1D (Order_x, N0, N5, i0, i1,                          &
                                    1, Q             , Flag,                  &
                            Tolerance, Jump_Factor, Scale, Default_Huge       )

  END Subroutine Multi_Resolution_0
!
! ======================================================================
!
  Subroutine Multi_Resolution_1a (Order_x, N0, N5, i0, i1,                    &
                                  NV, Direction, Index, Q, Flag,              &
                                  Tolerance, Jump_Factor, Scale, Default_Huge )

  integer            :: Order_x, N0, N5, i0, i1
  integer            :: NV, Direction
  integer            :: Index
  REALTYPE, OPTIONAL :: Tolerance, Jump_Factor
  integer , OPTIONAL :: Default_Huge

  REALTYPE, dimension(N0:N5,NV)        :: Q
  integer , dimension(N0:N5)           :: Flag
  REALTYPE, dimension(N0:N5), OPTIONAL :: Scale

  call Multi_Resolution_1D (Order_x, N0, N5, i0, i1,                          &
                                    1, Q(:,Index)    , Flag,                  &
                            Tolerance, Jump_Factor, Scale, Default_Huge       )

  END Subroutine Multi_Resolution_1a
!
! ======================================================================
!
  Subroutine Multi_Resolution_2a (Order_x, N0, N5, i0, i1,                    &
                                  Order_y, M0, M5, j0, j1,                    &
                                  NV, Direction, Index, Q, Flag,              &
                                  Tolerance, Jump_Factor, Scale, Default_Huge )

  integer            :: Order_x, N0, N5, i0, i1
  integer            :: Order_y, M0, M5, j0, j1
  integer            :: NV, Direction
  integer            :: Index
  REALTYPE, OPTIONAL :: Tolerance, Jump_Factor
  integer , OPTIONAL :: Default_Huge

  REALTYPE, dimension(N0:N5,M0:M5,NV)        :: Q
  integer , dimension(N0:N5,M0:M5)           :: Flag
  REALTYPE, dimension(N0:N5,M0:M5), OPTIONAL :: Scale

  call Multi_Resolution_2D (Order_x, N0, N5, i0, i1,                          &
                            Order_y, M0, M5, j0, j1,                          &
                            Direction, Q(:,:,Index)  , Flag,                  &
                            Tolerance, Jump_Factor, Scale, Default_Huge       )

  END Subroutine Multi_Resolution_2a
!
! ======================================================================
!
  Subroutine Multi_Resolution_3a (Order_x, N0, N5, i0, i1,                    &
                                  Order_y, M0, M5, j0, j1,                    &
                                  Order_z, K0, K5, l0, l1,                    &
                                  NV, Direction, Index, Q, Flag,              &
                                  Tolerance, Jump_Factor, Scale, Default_Huge )

  integer            :: Order_x, N0, N5, i0, i1
  integer            :: Order_y, M0, M5, j0, j1
  integer            :: Order_z, K0, K5, l0, l1
  integer            :: NV, Direction
  integer            :: Index
  REALTYPE, OPTIONAL :: Tolerance, Jump_Factor
  integer , OPTIONAL :: Default_Huge

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)        :: Q
  integer , dimension(N0:N5,M0:M5,K0:K5)           :: Flag
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5), OPTIONAL :: Scale

  call Multi_Resolution_3D (Order_x, N0, N5, i0, i1,                          &
                            Order_y, M0, M5, j0, j1,                          &
                            Order_z, K0, K5, l0, l1,                          &
                            Direction, Q(:,:,:,Index), Flag,                  &
                            Tolerance, Jump_Factor, Scale, Default_Huge       )

  END Subroutine Multi_Resolution_3a
!
! ======================================================================
!
  Subroutine Multi_Resolution_1D (Order_x, N0, N5, i0, i1,                &
                                  Direction, Q, Flag,                     &
                                  Tolerance, Jump_Factor,                 &
                                  Scale, Default_Huge                     )

  integer            :: Order_x, N0, N5, i0, i1
  integer            :: Direction
  REALTYPE, OPTIONAL :: Tolerance, Jump_Factor
  integer , OPTIONAL :: Default_Huge

  REALTYPE, dimension(N0:N5)           :: Q
  integer , dimension(N0:N5)           :: Flag
  REALTYPE, dimension(N0:N5), OPTIONAL :: Scale

  REALTYPE, dimension(N0:N5)                        :: dd
  REALTYPE, dimension(-(Order_x+1)/2:(Order_x+1)/2) :: a_x
  REALTYPE, dimension(:), ALLOCATABLE               :: Q_1, Q_2, Q_3

  integer  :: P, L, R
  integer  :: i, N2, NN, ii, i2
  REALTYPE :: Tolerance_Level, Factor

                            Tolerance_Level = 1.0d-3
  if (PRESENT(Tolerance  )) Tolerance_Level = Tolerance 

                                     Factor = ONE
  if (PRESENT(Jump_Factor))          Factor = Jump_Factor

  Tolerance_Level = Factor*Tolerance_Level

  dd = ZERO ; if (PRESENT(Default_Huge)) dd = HUGE(ONE)

  call MR_Coefficients (Order_x, a_x)

  NN = i1-i0+1 ; N2 = NN/2 ; if (MOD(NN,2) == 1) N2 = N2+1

  P = Order_x+1 ; L = P/2 ; R = Order_x-L

  i2 = 0 ; if (MOD(Order_x,2) == 1) i2 = 1

  ALLOCATE (Q_1(-L:N2+R), Q_2(0:N2-1), Q_3(0:N2-1))

#if defined (DEBUG)
  call Check_Parameter (Order_x, N0, N5, i0, i1, N2, R)
#endif

  do i = -L, N2-1+R+i2
    ii = i0 + 2*i

    Q_1(i) = HALF*(Q(ii) + Q(ii+1))
  enddo

#if 1
  call Direct_Sum_1D_x (N2, Order_x, L, R, a_x, Q_1, Q_2, Q_3)

  do i = 0, N2-1
    ii = i0 + 2*i

    dd(ii  ) = Q(ii  ) - Q_2(i)
    dd(ii+1) = Q(ii+1) - Q_3(i)
  enddo
#else
  do i = 0, N2-1
    ii = i0 + 2*i

    dd(ii  ) = Q(ii  ) - SUM(a_x(-L: R: 1)*Q_1(i-L   :i+R   ))
    dd(ii+1) = Q(ii+1) - SUM(a_x( R:-L:-1)*Q_1(i-L+i2:i+R+i2))
  enddo
#endif

  Flag = 1 ; WHERE (ABS(dd) <= Tolerance_Level) Flag = 0

  Flag(N0  :i0-1) = 0 
  Flag(i1+1:N5  ) = 0

  if (PRESENT(Scale)) Scale = dd

  DEALLOCATE (Q_1, Q_2, Q_3)

  END Subroutine Multi_Resolution_1D
!
! ======================================================================
!
  Subroutine Multi_Resolution_2D (Order_x, N0, N5, i0, i1,                &
                                  Order_y, M0, M5, j0, j1,                &
                                  Direction, Q, Flag,                     &
                                  Tolerance, Jump_Factor,                 &
                                  Scale, Default_Huge                     )

  integer            :: Order_x, N0, N5, i0, i1
  integer            :: Order_y, M0, M5, j0, j1
  integer            :: Direction
  REALTYPE, OPTIONAL :: Tolerance, Jump_Factor
  integer , OPTIONAL :: Default_Huge

  REALTYPE, dimension(N0:N5,M0:M5)           :: Q
  integer , dimension(N0:N5,M0:M5)           :: Flag
  REALTYPE, dimension(N0:N5,M0:M5), OPTIONAL :: Scale

  REALTYPE, dimension(N0:N5,M0:M5)                  :: dd
  REALTYPE, dimension(-(Order_x+1)/2:(Order_x+1)/2) :: a_x
  REALTYPE, dimension(-(Order_y+1)/2:(Order_y+1)/2) :: a_y
  REALTYPE, dimension(:,:), ALLOCATABLE             :: Q_1, Q_2, Q_3

  integer  :: P, L, R
  integer  :: i, N2, NN, ii, i2
  integer  :: j, M2, MM, jj, j2
  REALTYPE :: Tolerance_Level, Factor

                            Tolerance_Level = 1.0d-3
  if (PRESENT(Tolerance  )) Tolerance_Level = Tolerance 

                                     Factor = ONE
  if (PRESENT(Jump_Factor))          Factor = Jump_Factor

  Tolerance_Level = Factor*Tolerance_Level

  dd = ZERO ; if (PRESENT(Default_Huge)) dd = HUGE(ONE)

  SELECT CASE (Direction)
    CASE (1)
      call MR_Coefficients (Order_x, a_x)
#if 0
      NN = i1-i0+1 ; N2 = NN/2 ; if (MOD(NN,2) == 1) N2 = N2+1

      P = Order_x+1 ; L = P/2 ; R = Order_x-L

      i2 = 0 ; if (MOD(Order_x,2) == 1) i2 = 1

!     ALLOCATE (Q_1(-L:N2+R+1,j0:j1), Q_2(0:N2-1,j0:j1), Q_3(0:N2-1,j0:j1))
      ALLOCATE (Q_1(-L:N2+R,j0:j1), Q_2(0:N2-1,j0:j1), Q_3(0:N2-1,j0:j1))

!     do i = -L, N2  +R+i2
      do i = -L, N2-1+R+i2
        ii = i0 + 2*i

        Q_1(i,:) = HALF*(Q(ii,j0:j1) + Q(ii+1,j0:j1))
write (6,*) i, ii, Q(ii,j0), Q(ii+1,j0)
      enddo

!     do i = 0, N2
      do i = 0, N2-1
        ii = i0 + 2*i

        do j = j0, j1
          dd(ii  ,j) = Q(ii  ,j) - SUM(a_x(-L: R: 1)*Q_1(i-L   :i+R   ,j))
          dd(ii+1,j) = Q(ii+1,j) - SUM(a_x( R:-L:-1)*Q_1(i-L+i2:i+R+i2,j))
        enddo
      enddo

write (6,*) 'N0, N5, i0, i1, Order : ', N0, N5, i0, i1, Order_x 
write (6,*) 'P, L, R               : ', P, L, R
write (6,*) 'NN, N2, i2, ii        : ', NN, N2, i2, ii

write (8,*) 'N0, N5, i0, i1, Order : ', N0, N5, i0, i1, Order_x 
write (8,*) 'P, L, R               : ', P, L, R
write (8,*) 'NN, N2, i2, ii        : ', NN, N2, i2, ii

#else

      NN = i1-i0+1 ; N2 = NN/2 ; if (MOD(NN,2) == 1) N2 = N2+1

      P = Order_x+1 ; L = P/2 ; R = Order_x-L

      i2 = 0 ; if (MOD(Order_x,2) == 1) i2 = 1

      ALLOCATE (Q_1(-L:N2+R,j0:j1), Q_2(0:N2-1,j0:j1), Q_3(0:N2-1,j0:j1))

#if defined (DEBUG)
      call Check_Parameter (Order_x, N0, N5, i0, i1, N2, R)
#endif

      do i = -L, N2-1+R+i2
        ii = i0 + 2*i

        Q_1(i,:) = HALF*(Q(ii,j0:j1) + Q(ii+1,j0:j1))
      enddo

#if 1
      call Direct_Sum_2D_x (N2, Order_x, L, R, a_x, Q_1, Q_2, Q_3)

      do i = 0, N2-1
        ii = i0 + 2*i

          dd(ii  ,j0:j1) = Q(ii  ,j0:j1) - Q_2(i,:)
          dd(ii+1,j0:j1) = Q(ii+1,j0:j1) - Q_3(i,:)
      enddo
#else
      do i = 0, N2-1
        ii = i0 + 2*i

        do j = j0, j1
          dd(ii  ,j) = Q(ii  ,j) - SUM(a_x(-L: R: 1)*Q_1(i-L   :i+R   ,j))
          dd(ii+1,j) = Q(ii+1,j) - SUM(a_x( R:-L:-1)*Q_1(i-L+i2:i+R+i2,j))
        enddo
      enddo
#endif

#endif

      Flag = 1 ; WHERE (ABS(dd) <= Tolerance_Level) Flag = 0

#if 0
     write (6,100) 'Q    ',Q(:,j0)
     write (8,100) 'Q    ',Q(:,j0)
     write (6,100) 'Q_1  ',Q_1(:,j0)
     write (8,100) 'Q_1  ',Q_1(:,j0)
     write (6,100) 'dd   ',ABS(dd(:,j0))
     write (8,100) 'dd   ',ABS(dd(:,j0))
     write (6,100) 'Flag ',Flag(:,j0)
     write (8,100) 'Flag ',Flag(:,j0)
     STOP
100 format (1x,a8,100(1x,g10.4,:))
#endif

      Flag(N0  :i0-1,j0:j1) = 0 
      Flag(i1+1:N5  ,j0:j1) = 0

    CASE (2)
      call MR_Coefficients (Order_y, a_y)

      MM = j1-j0+1 ; M2 = MM/2 ; if (MOD(MM,2) == 1) M2 = M2+1

      P = Order_y+1 ; L = P/2 ; R = Order_y-L

      j2 = 0 ; if (MOD(Order_y,2) == 1) j2 = 1
       
      ALLOCATE (Q_1(i0:i1,-L:M2+R), Q_2(i0:i1,0:M2-1), Q_3(i0:i1,0:M2-1))

#if defined (DEBUG)
      call Check_Parameter (Order_y, M0, M5, j0, j1, M2, R)
#endif

      do j = -L, M2-1+R+j2
        jj = j0 + 2*j

        Q_1(:,j) = HALF*(Q(i0:i1,jj) + Q(i0:i1,jj+1))
      enddo

#if 1
      call Direct_Sum_2D_y (M2, Order_y, L, R, a_y, Q_1, Q_2, Q_3)

      do j = 0, M2-1
        jj = j0 + 2*j

          dd(i0:i1,jj  ) = Q(i0:i1,jj  ) - Q_2(:,j)
          dd(i0:i1,jj+1) = Q(i0:i1,jj+1) - Q_3(:,j)
      enddo
#else
      do j = 0, M2-1
        jj = j0 + 2*j

        do i = i0, i1
          dd(i,jj  ) = Q(i,jj  ) - SUM(a_y(-L: R: 1)*Q_1(i,j-L   :j+R   ))
          dd(i,jj+1) = Q(i,jj+1) - SUM(a_y( R:-L:-1)*Q_1(i,j-L+j2:j+R+j2))
        enddo
      enddo
#endif

      Flag = 1 ; WHERE (ABS(dd) <= Tolerance_Level) Flag = 0

      Flag(i0:i1,M0  :j0-1) = 0 
      Flag(i0:i1,j1+1:M5  ) = 0

  END SELECT

  if (PRESENT(Scale)) Scale = dd

  DEALLOCATE (Q_1, Q_2, Q_3)

  END Subroutine Multi_Resolution_2D
!
! ======================================================================
!
  Subroutine Multi_Resolution_3D (Order_x, N0, N5, i0, i1,                &
                                  Order_y, M0, M5, j0, j1,                &
                                  Order_z, K0, K5, l0, l1,                &
                                  Direction, Q, Flag,                     &
                                  Tolerance, Jump_Factor,                 &
                                  Scale, Default_Huge                     )

  integer            :: Order_x, N0, N5, i0, i1
  integer            :: Order_y, M0, M5, j0, j1
  integer            :: Order_z, K0, K5, l0, l1
  integer            :: Direction
  REALTYPE, OPTIONAL :: Tolerance, Jump_Factor
  integer , OPTIONAL :: Default_Huge

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)           :: Q
  integer , dimension(N0:N5,M0:M5,K0:K5)           :: Flag
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5), OPTIONAL :: Scale

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)            :: dd
  REALTYPE, dimension(-(Order_x+1)/2:(Order_x+1)/2) :: a_x
  REALTYPE, dimension(-(Order_y+1)/2:(Order_y+1)/2) :: a_y
  REALTYPE, dimension(-(Order_z+1)/2:(Order_z+1)/2) :: a_z
  REALTYPE, dimension(:,:,:), ALLOCATABLE           :: Q_1, Q_2, Q_3

  integer  :: P, L, R
  integer  :: i, N2, NN, ii, i2
  integer  :: j, M2, MM, jj, j2
  integer  :: k, K2, KK, ll, l2
  REALTYPE :: Tolerance_Level, Factor

                            Tolerance_Level = 1.0d-3
  if (PRESENT(Tolerance  )) Tolerance_Level = Tolerance 

                                     Factor = ONE
  if (PRESENT(Jump_Factor))          Factor = Jump_Factor

  Tolerance_Level = Factor*Tolerance_Level

  dd = ZERO ; if (PRESENT(Default_Huge)) dd = HUGE(ONE)

  SELECT CASE (Direction)
    CASE (1)
      call MR_Coefficients (Order_x, a_x)

      NN = i1-i0+1 ; N2 = NN/2 ; if (MOD(NN,2) == 1) N2 = N2+1

      P = Order_x+1 ; L = P/2 ; R = Order_x-L

      i2 = 0 ; if (MOD(Order_x,2) == 1) i2 = 1

      ALLOCATE (Q_1(-L:N2+R,j0:j1,l0:l1), Q_2(0:N2-1,j0:j1,l0:l1),  &
                                          Q_3(0:N2-1,j0:j1,l0:l1))

#if defined (DEBUG)
      call Check_Parameter (Order_x, N0, N5, i0, i1, N2, R)
#endif

      do i = -L, N2-1+R+i2
        ii = i0 + 2*i

        Q_1(i,:,:) = HALF*(Q(ii,j0:j1,l0:l1) + Q(ii+1,j0:j1,l0:l1))
      enddo

#if 1
      call Direct_Sum_3D_x (N2, Order_x, L, R, a_x, Q_1, Q_2, Q_3)

      do i = 0, N2-1
        ii = i0 + 2*i

          dd(ii  ,j0:j1,l0:l1) = Q(ii  ,j0:j1,l0:l1) - Q_2(i,:,:)
          dd(ii+1,j0:j1,l0:l1) = Q(ii+1,j0:j1,l0:l1) - Q_3(i,:,:)
      enddo
#else
      do i = 0, N2-1
        ii = i0 + 2*i

        do k = l0, l1
          do j = j0, j1
            dd(ii  ,j,k) = Q(ii  ,j,k) -  &
                             SUM(a_x(-L: R: 1)*Q_1(i-L   :i+R   ,j,k))
            dd(ii+1,j,k) = Q(ii+1,j,k) -  &
                             SUM(a_x( R:-L:-1)*Q_1(i-L+i2:i+R+i2,j,k))
          enddo
        enddo
      enddo
#endif

      Flag = 1 ; WHERE (ABS(dd) <= Tolerance_Level) Flag = 0

      Flag(N0  :i0-1,j0:j1,l0:l1) = 0 
      Flag(i1+1:N5  ,j0:j1,l0:l1) = 0

    CASE (2)
      call MR_Coefficients (Order_y, a_y)

      MM = j1-j0+1 ; M2 = MM/2 ; if (MOD(MM,2) == 1) M2 = M2+1

      P = Order_y+1 ; L = P/2 ; R = Order_y-L

      j2 = 0 ; if (MOD(Order_y,2) == 1) j2 = 1
       
      ALLOCATE (Q_1(i0:i1,-L:M2+R,l0:l1), Q_2(i0:i1,0:M2-1,l0:l1), &
                                          Q_3(i0:i1,0:M2-1,l0:l1))

#if defined (DEBUG)
      call Check_Parameter (Order_y, M0, M5, j0, j1, M2, R)
#endif

      do j = -L, M2-1+R+j2
        jj = j0 + 2*j

        Q_1(:,j,:) = HALF*(Q(i0:i1,jj,l0:l1) + Q(i0:i1,jj+1,l0:l1))
      enddo

#if 1
      call Direct_Sum_3D_y (M2, Order_y, L, R, a_y, Q_1, Q_2, Q_3)

      do j = 0, M2-1
        jj = j0 + 2*j

          dd(i0:i1,jj  ,l0:l1) = Q(i0:i1,jj  ,l0:l1) - Q_2(:,j,:)
          dd(i0:i1,jj+1,l0:l1) = Q(i0:i1,jj+1,l0:l1) - Q_3(:,j,:)
      enddo
#else
      do j = 0, M2-1
        jj = j0 + 2*j

        do k = l0, l1
          do i = i0, i1
            dd(i,jj  ,k) = Q(i,jj  ,k) - &
                             SUM(a_y(-L: R: 1)*Q_1(i,j-L   :j+R   ,k))
            dd(i,jj+1,k) = Q(i,jj+1,k) - &
                             SUM(a_y( R:-L:-1)*Q_1(i,j-L+j2:j+R+j2,k))
          enddo
        enddo
      enddo
#endif

      Flag = 1 ; WHERE (ABS(dd) <= Tolerance_Level) Flag = 0

      Flag(i0:i1,M0  :j0-1,l0:l1) = 0 
      Flag(i0:i1,j1+1:M5  ,l0:l1) = 0

    CASE (3)
      call MR_Coefficients (Order_z, a_z)

      KK = l1-l0+1 ; K2 = KK/2 ; if (MOD(KK,2) == 1) K2 = K2+1

      P = Order_z+1 ; L = P/2 ; R = Order_z-L

      l2 = 0 ; if (MOD(Order_z,2) == 1) l2 = 1

      ALLOCATE (Q_1(i0:i1,j0:j1,-L:K2+R), Q_2(i0:i1,j0:j1,0:K2-1),  &
                                          Q_3(i0:i1,j0:j1,0:K2-1))

#if defined (DEBUG)
      call Check_Parameter (Order_z, K0, K5, l0, l1, K2, R)
#endif

      do k = -L, K2-1+R+l2
        ll = l0 + 2*k

        Q_1(:,:,k) = HALF*(Q(i0:i1,j0:j1,ll) + Q(i0:i1,j0:j1,ll+1))
      enddo

#if 1
      call Direct_Sum_3D_z (K2, Order_z, L, R, a_z, Q_1, Q_2, Q_3)

      do k = 0, K2-1
        ll = l0 + 2*k

          dd(i0:i1,j0:j1,ll  ) = Q(i0:i1,j0:j1,ll  ) - Q_2(:,:,k)
          dd(i0:i1,j0:j1,ll+1) = Q(i0:i1,j0:j1,ll+1) - Q_3(:,:,k)
      enddo
#else
      do k = 0, K2-1
        ll = l0 + 2*k

        do j = j0, j1
          do i = i0, i1
            dd(i,j,ll  ) = Q(i,j,ll  ) -  &
                             SUM(a_z(-L: R: 1)*Q_1(i,j,k-L   :k+R   ))
            dd(i,j,ll+1) = Q(i,j,ll+1) -  &
                             SUM(a_z( R:-L:-1)*Q_1(i,j,k-L+l2:k+R+l2))
          enddo
        enddo
      enddo
#endif

      Flag = 1 ; WHERE (ABS(dd) <= Tolerance_Level) Flag = 0

      Flag(i0:i1,j0:j1,K0  :l0-1) = 0 
      Flag(i0:i1,j0:j1,l1+1:K5  ) = 0

  END SELECT

  if (PRESENT(Scale)) Scale = dd

  DEALLOCATE (Q_1, Q_2, Q_3)

  END Subroutine Multi_Resolution_3D
!
! ======================================================================
!
  Subroutine MR_Coefficients (Order, Beta)
  
  integer                                       :: Order
  REALTYPE, dimension(-(Order+1)/2:(Order+1)/2) :: Beta

  REALTYPE                                      :: Factor

  Beta = ZERO

  SELECT CASE (Order)
    CASE (1)
      Factor   =  4.0d0

      Beta(-1) =  1.0d0
      Beta( 0) =  3.0d0

    CASE DEFAULT
      Factor   =  8.0d0

      Beta(-1) =  1.0d0
      Beta( 0) =  8.0d0

    CASE (3)
      Factor   =  64.0d0

      Beta(-2) = - 3.0d0
      Beta(-1) =  17.0d0
      Beta( 0) =  55.0d0
      Beta( 1) = - 5.0d0

    CASE (4)
      Factor   = 128.0d0

      Beta(-2) = - 3.0d0
      Beta(-1) =  22.0d0
      Beta( 0) = 128.0d0

    CASE (5)
      Factor   =  512.0d0

      Beta(-3) =    5.0d0
      Beta(-2) = - 37.0d0
      Beta(-1) =  138.0d0
      Beta( 0) =  462.0d0
      Beta( 1) = - 63.0d0
      Beta( 2) =    7.0d0

    CASE (6)
      Factor   = 1024.0d0

      Beta(-3) =    5.0d0
      Beta(-2) = - 44.0d0
      Beta(-1) =  201.0d0
      Beta( 0) = 1024.0d0

    CASE (7)
      Factor   =  16384.0d0

      Beta(-4) = -   35.0d0
      Beta(-3) =    325.0d0
      Beta(-2) = - 1439.0d0
      Beta(-1) =   4441.0d0
      Beta( 0) =  15159.0d0
      Beta( 1) = - 2481.0d0
      Beta( 2) =    459.0d0
      Beta( 3) = -   45.0d0

    CASE (8)
      Factor   =  32768.0d0

      Beta(-4) = -   35.0d0
      Beta(-3) =    370.0d0
      Beta(-2) = - 1898.0d0
      Beta(-1) =   6922.0d0
      Beta( 0) =  32768.0d0

    CASE (9)
      Factor   =  131072.0d0

      Beta(-5) =      63.0d0
      Beta(-4) = -   707.0d0
      Beta(-3) =    3748.0d0
      Beta(-2) = - 12884.0d0
      Beta(-1) =   35626.0d0
      Beta( 0) =  123134.0d0
      Beta( 1) = - 22396.0d0
      Beta( 2) =    5324.0d0
      Beta( 3) = -   913.0d0
      Beta( 4) =      77.0d0

    CASE (10)
      Factor   =  262144.0d0

      Beta(-5) =      63.0d0
      Beta(-4) = -   784.0d0
      Beta(-3) =    4661.0d0
      Beta(-2) = - 18208.0d0
      Beta(-1) =   58022.0d0
      Beta( 0) =  262144.0d0

    CASE (11)
      Factor   =  2097152.0d0

      Beta(-6) = -    231.0d0
      Beta(-5) =     3045.0d0
      Beta(-4) = -  18977.0d0
      Beta(-3) =    75403.0d0
      Beta(-2) = - 221894.0d0
      Beta(-1) =   570898.0d0
      Beta( 0) =  1990430.0d0
      Beta( 1) = - 387946.0d0
      Beta( 2) =   107549.0d0
      Beta( 3) = -  24583.0d0
      Beta( 4) =     3731.0d0
      Beta( 5) = -    273.0d0  

  END SELECT

  if (MOD(Order,2) == 0) Beta(1:Order/2) = -Beta(-1:-Order/2:-1)

  Beta = Beta/Factor

  END Subroutine MR_Coefficients 
!
! ======================================================================
!
  Subroutine Direct_Sum_1D_x (N2, Order, L, R, A, B, C, D)  

  integer  :: N2, Order, L, R

  REALTYPE, dimension(-(Order+1)/2:(Order+1)/2) :: A 
  REALTYPE, dimension(-L:N2+R)                  :: B
  REALTYPE, dimension(:)                        :: C, D

  integer :: i, N

  N = N2-1

  SELECT CASE (Order)
    CASE (1)
      C = A(0)*B(0:N) + A(-1)*B(-1:N-1) 

      D = C + A(-1)      *(B(1:N+1)-B(-1:N-1)) 

    CASE (2)
      C = A(0)*B(0:N) + A(-1)*(B(-1:N-1)-B(1:N+1))

    CASE (3)
      C = A(0)*B(0:N) + A(-1)*B(-1:N-1) + A(-2)*B(-2:N-2) + &
                        A( 1)*B( 1:N+1)
      
      D = C + A(-2)      *(B(2:N+2)-B(-2:N-2))            + &
             (A(-1)-A(1))*(B(1:N+1)-B(-1:N-1)) 

    CASE (4)
      C = A(0)*B(0:N) + A(-1)*(B(-1:N-1)-B(1:N+1))        + &
                        A(-2)*(B(-2:N-2)-B(2:N+2)) 

    CASE (5)
      C = A(0)*B(0:N) + A(-1)*B(-1:N-1) + A(-2)*B(-2:N-2) + &
                        A(-3)*B(-3:N-3)                   + &
                        A( 1)*B( 1:N+1) + A( 2)*B( 2:N+2) 
      
      D = C + A(-3)      *(B(3:N+3)-B(-3:N-3))            + &
             (A(-2)-A(2))*(B(2:N+2)-B(-2:N-2))            + &
             (A(-1)-A(1))*(B(1:N+1)-B(-1:N-1)) 

    CASE (6)
      C = A(0)*B(0:N) + A(-1)*(B(-1:N-1)-B(1:N+1))        + &
                        A(-2)*(B(-2:N-2)-B(2:N+2))        + & 
                        A(-3)*(B(-3:N-3)-B(3:N+3)) 

    CASE (7)
      C = A(0)*B(0:N) + A(-1)*B(-1:N-1) + A(-2)*B(-2:N-2) + &
                        A(-3)*B(-3:N-3) + A(-4)*B(-4:N-4) + &
                        A( 1)*B( 1:N+1) + A( 2)*B( 2:N+2) + &
                        A( 3)*B( 3:N+3) 

      D = C + A(-4)      *(B(4:N+4)-B(-4:N-4)) + &
             (A(-3)-A(3))*(B(3:N+3)-B(-3:N-3)) + & 
             (A(-2)-A(2))*(B(2:N+2)-B(-2:N-2)) + &
             (A(-1)-A(1))*(B(1:N+1)-B(-1:N-1)) 

    CASE (8)
      C = A(0)*B(0:N) + A(-1)*(B(-1:N-1)-B(1:N+1)) + &
                        A(-2)*(B(-2:N-2)-B(2:N+2)) + & 
                        A(-3)*(B(-3:N-3)-B(3:N+3)) + &
                        A(-4)*(B(-4:N-4)-B(4:N+4)) 

    CASE (9)
      C = A(0)*B(0:N) + A(-1)*B(-1:N-1) + A(-2)*B(-2:N-2) + &
                        A(-3)*B(-3:N-3) + A(-4)*B(-4:N-4) + &
                        A(-5)*B(-5:N-5)                   + &
                        A( 1)*B( 1:N+1) + A( 2)*B( 2:N+2) + &
                        A( 3)*B( 3:N+3) + A( 4)*B( 4:N+4) 

      D = C + A(-5)      *(B(5:N+5)-B(-5:N-5)) + &
             (A(-4)-A(4))*(B(4:N+4)-B(-4:N-4)) + & 
             (A(-3)-A(3))*(B(3:N+3)-B(-3:N-3)) + & 
             (A(-2)-A(2))*(B(2:N+2)-B(-2:N-2)) + &
             (A(-1)-A(1))*(B(1:N+1)-B(-1:N-1)) 

    CASE (10)
      C = A(0)*B(0:N) + A(-1)*(B(-1:N-1)-B(1:N+1)) + &
                        A(-2)*(B(-2:N-2)-B(2:N+2)) + & 
                        A(-3)*(B(-3:N-3)-B(3:N+3)) + &
                        A(-4)*(B(-4:N-4)-B(4:N+4)) + &
                        A(-5)*(B(-5:N-5)-B(5:N+5)) 

    CASE (11)
      C = A(0)*B(0:N) + A(-1)*B(-1:N-1) + A(-2)*B(-2:N-2) + &
                        A(-3)*B(-3:N-3) + A(-4)*B(-4:N-4) + &
                        A(-5)*B(-5:N-5) + A(-6)*B(-6:N-6) + &
                        A( 1)*B( 1:N+1) + A( 2)*B( 2:N+2) + &
                        A( 3)*B( 3:N+3) + A( 4)*B( 4:N+4) + &
                        A( 5)*B( 5:N+5) 

      D = C + A(-6)      *(B(6:N+6)-B(-6:N-6)) + &
             (A(-5)-A(5))*(B(5:N+5)-B(-5:N-5)) + & 
             (A(-4)-A(4))*(B(4:N+4)-B(-4:N-4)) + & 
             (A(-3)-A(3))*(B(3:N+3)-B(-3:N-3)) + & 
             (A(-2)-A(2))*(B(2:N+2)-B(-2:N-2)) + &
             (A(-1)-A(1))*(B(1:N+1)-B(-1:N-1)) 

  END SELECT

  if (MOD(Order,2) /= 0) RETURN

  D = -C + TWO*B(0:N)

  END Subroutine Direct_Sum_1D_x
!
! ======================================================================
!
  Subroutine Direct_Sum_2D_x (N2, Order, L, R, A, B, C, D)  

  integer  :: N2, Order, L, R

  REALTYPE, dimension(-(Order+1)/2:(Order+1)/2) :: A 
  REALTYPE, dimension(-L:,:)                    :: B
  REALTYPE, dimension(:,:)                      :: C, D

  integer :: i, N

  N = N2-1

  SELECT CASE (Order)
    CASE (1)
      C = A(0)*B(0:N,:) + A(-1)*B(-1:N-1,:) 

      D = C + A(-1)      *(B(1:N+1,:)-B(-1:N-1,:)) 

    CASE (2)
      C = A(0)*B(0:N,:) + A(-1)*(B(-1:N-1,:)-B(1:N+1,:))

    CASE (3)
      C = A(0)*B(0:N,:) + A(-1)*B(-1:N-1,:) + A(-2)*B(-2:N-2,:) + &
                          A( 1)*B( 1:N+1,:)
      
      D = C + A(-2)      *(B(2:N+2,:)-B(-2:N-2,:))              + &
             (A(-1)-A(1))*(B(1:N+1,:)-B(-1:N-1,:)) 

    CASE (4)
      C = A(0)*B(0:N,:) + A(-1)*(B(-1:N-1,:)-B(1:N+1,:))        + &
                          A(-2)*(B(-2:N-2,:)-B(2:N+2,:)) 

    CASE (5)
      C = A(0)*B(0:N,:) + A(-1)*B(-1:N-1,:) + A(-2)*B(-2:N-2,:) + &
                          A(-3)*B(-3:N-3,:)                     + &
                          A( 1)*B( 1:N+1,:) + A( 2)*B( 2:N+2,:) 
      
      D = C + A(-3)      *(B(3:N+3,:)-B(-3:N-3,:))              + &
             (A(-2)-A(2))*(B(2:N+2,:)-B(-2:N-2,:))              + &
             (A(-1)-A(1))*(B(1:N+1,:)-B(-1:N-1,:)) 

    CASE (6)
      C = A(0)*B(0:N,:) + A(-1)*(B(-1:N-1,:)-B(1:N+1,:))        + &
                          A(-2)*(B(-2:N-2,:)-B(2:N+2,:))        + & 
                          A(-3)*(B(-3:N-3,:)-B(3:N+3,:)) 

    CASE (7)
      C = A(0)*B(0:N,:) + A(-1)*B(-1:N-1,:) + A(-2)*B(-2:N-2,:) + &
                          A(-3)*B(-3:N-3,:) + A(-4)*B(-4:N-4,:) + &
                          A( 1)*B( 1:N+1,:) + A( 2)*B( 2:N+2,:) + &
                          A( 3)*B( 3:N+3,:) 

      D = C + A(-4)      *(B(4:N+4,:)-B(-4:N-4,:)) + &
             (A(-3)-A(3))*(B(3:N+3,:)-B(-3:N-3,:)) + & 
             (A(-2)-A(2))*(B(2:N+2,:)-B(-2:N-2,:)) + &
             (A(-1)-A(1))*(B(1:N+1,:)-B(-1:N-1,:)) 

    CASE (8)
      C = A(0)*B(0:N,:) + A(-1)*(B(-1:N-1,:)-B(1:N+1,:)) + &
                          A(-2)*(B(-2:N-2,:)-B(2:N+2,:)) + & 
                          A(-3)*(B(-3:N-3,:)-B(3:N+3,:)) + &
                          A(-4)*(B(-4:N-4,:)-B(4:N+4,:)) 

    CASE (9)
      C = A(0)*B(0:N,:) + A(-1)*B(-1:N-1,:) + A(-2)*B(-2:N-2,:) + &
                          A(-3)*B(-3:N-3,:) + A(-4)*B(-4:N-4,:) + &
                          A(-5)*B(-5:N-5,:)                     + &
                          A( 1)*B( 1:N+1,:) + A( 2)*B( 2:N+2,:) + &
                          A( 3)*B( 3:N+3,:) + A( 4)*B( 4:N+4,:) 

      D = C + A(-5)      *(B(5:N+5,:)-B(-5:N-5,:)) + &
             (A(-4)-A(4))*(B(4:N+4,:)-B(-4:N-4,:)) + & 
             (A(-3)-A(3))*(B(3:N+3,:)-B(-3:N-3,:)) + & 
             (A(-2)-A(2))*(B(2:N+2,:)-B(-2:N-2,:)) + &
             (A(-1)-A(1))*(B(1:N+1,:)-B(-1:N-1,:)) 

    CASE (10)
      C = A(0)*B(0:N,:) + A(-1)*(B(-1:N-1,:)-B(1:N+1,:)) + &
                          A(-2)*(B(-2:N-2,:)-B(2:N+2,:)) + & 
                          A(-3)*(B(-3:N-3,:)-B(3:N+3,:)) + &
                          A(-4)*(B(-4:N-4,:)-B(4:N+4,:)) + &
                          A(-5)*(B(-5:N-5,:)-B(5:N+5,:)) 

    CASE (11)
      C = A(0)*B(0:N,:) + A(-1)*B(-1:N-1,:) + A(-2)*B(-2:N-2,:) + &
                          A(-3)*B(-3:N-3,:) + A(-4)*B(-4:N-4,:) + &
                          A(-5)*B(-5:N-5,:) + A(-6)*B(-6:N-6,:) + &
                          A( 1)*B( 1:N+1,:) + A( 2)*B( 2:N+2,:) + &
                          A( 3)*B( 3:N+3,:) + A( 4)*B( 4:N+4,:) + &
                          A( 5)*B( 5:N+5,:) 

      D = C + A(-6)      *(B(6:N+6,:)-B(-6:N-6,:)) + &
             (A(-5)-A(5))*(B(5:N+5,:)-B(-5:N-5,:)) + & 
             (A(-4)-A(4))*(B(4:N+4,:)-B(-4:N-4,:)) + & 
             (A(-3)-A(3))*(B(3:N+3,:)-B(-3:N-3,:)) + & 
             (A(-2)-A(2))*(B(2:N+2,:)-B(-2:N-2,:)) + &
             (A(-1)-A(1))*(B(1:N+1,:)-B(-1:N-1,:)) 

  END SELECT

  if (MOD(Order,2) /= 0) RETURN

  D = -C + TWO*B(0:N,:)

  END Subroutine Direct_Sum_2D_x
!
! ======================================================================
!
  Subroutine Direct_Sum_2D_y (N2, Order, L, R, A, B, C, D)  

  integer  :: N2, Order, L, R

  REALTYPE, dimension(-(Order+1)/2:(Order+1)/2) :: A 
  REALTYPE, dimension(:,-L:)                    :: B
  REALTYPE, dimension(:,:)                      :: C, D

  integer :: N

  N = N2-1

  SELECT CASE (Order)
    CASE (1)
      C = A(0)*B(:,0:N) + A(-1)*B(:,-1:N-1) 

      D = C + A(-1)      *(B(:,1:N+1)-B(:,-1:N-1)) 

    CASE (2)
      C = A(0)*B(:,0:N) + A(-1)*(B(:,-1:N-1)-B(:,1:N+1))

    CASE (3)
      C = A(0)*B(:,0:N) + A(-1)*B(:,-1:N-1) + A(-2)*B(:,-2:N-2) + &
                          A( 1)*B(:, 1:N+1)
      
      D = C + A(-2)      *(B(:,2:N+2)-B(:,-2:N-2))              + &
             (A(-1)-A(1))*(B(:,1:N+1)-B(:,-1:N-1)) 

    CASE (4)
      C = A(0)*B(:,0:N) + A(-1)*(B(:,-1:N-1)-B(:,1:N+1))        + &
                          A(-2)*(B(:,-2:N-2)-B(:,2:N+2)) 

    CASE (5)
      C = A(0)*B(:,0:N) + A(-1)*B(:,-1:N-1) + A(-2)*B(:,-2:N-2) + &
                          A(-3)*B(:,-3:N-3)                     + &
                          A( 1)*B(:, 1:N+1) + A( 2)*B(:, 2:N+2) 
      
      D = C + A(-3)      *(B(:,3:N+3)-B(:,-3:N-3))              + &
             (A(-2)-A(2))*(B(:,2:N+2)-B(:,-2:N-2))              + &
             (A(-1)-A(1))*(B(:,1:N+1)-B(:,-1:N-1)) 

    CASE (6)
      C = A(0)*B(:,0:N) + A(-1)*(B(:,-1:N-1)-B(:,1:N+1))        + &
                          A(-2)*(B(:,-2:N-2)-B(:,2:N+2))        + & 
                          A(-3)*(B(:,-3:N-3)-B(:,3:N+3)) 

    CASE (7)
      C = A(0)*B(:,0:N) + A(-1)*B(:,-1:N-1) + A(-2)*B(:,-2:N-2) + &
                          A(-3)*B(:,-3:N-3) + A(-4)*B(:,-4:N-4) + &
                          A( 1)*B(:, 1:N+1) + A( 2)*B(:, 2:N+2) + &
                          A( 3)*B(:, 3:N+3) 

      D = C + A(-4)      *(B(:,4:N+4)-B(:,-4:N-4)) + &
             (A(-3)-A(3))*(B(:,3:N+3)-B(:,-3:N-3)) + & 
             (A(-2)-A(2))*(B(:,2:N+2)-B(:,-2:N-2)) + &
             (A(-1)-A(1))*(B(:,1:N+1)-B(:,-1:N-1)) 

    CASE (8)
      C = A(0)*B(:,0:N) + A(-1)*(B(:,-1:N-1)-B(:,1:N+1)) + &
                          A(-2)*(B(:,-2:N-2)-B(:,2:N+2)) + & 
                          A(-3)*(B(:,-3:N-3)-B(:,3:N+3)) + &
                          A(-4)*(B(:,-4:N-4)-B(:,4:N+4)) 

    CASE (9)
      C = A(0)*B(:,0:N) + A(-1)*B(:,-1:N-1) + A(-2)*B(:,-2:N-2) + &
                          A(-3)*B(:,-3:N-3) + A(-4)*B(:,-4:N-4) + &
                          A(-5)*B(:,-5:N-5)                     + &
                          A( 1)*B(:, 1:N+1) + A( 2)*B(:, 2:N+2) + &
                          A( 3)*B(:, 3:N+3) + A( 4)*B(:, 4:N+4) 

      D = C + A(-5)      *(B(:,5:N+5)-B(:,-5:N-5)) + &
             (A(-4)-A(4))*(B(:,4:N+4)-B(:,-4:N-4)) + & 
             (A(-3)-A(3))*(B(:,3:N+3)-B(:,-3:N-3)) + & 
             (A(-2)-A(2))*(B(:,2:N+2)-B(:,-2:N-2)) + &
             (A(-1)-A(1))*(B(:,1:N+1)-B(:,-1:N-1)) 

    CASE (10)
      C = A(0)*B(:,0:N) + A(-1)*(B(:,-1:N-1)-B(:,1:N+1)) + &
                          A(-2)*(B(:,-2:N-2)-B(:,2:N+2)) + & 
                          A(-3)*(B(:,-3:N-3)-B(:,3:N+3)) + &
                          A(-4)*(B(:,-4:N-4)-B(:,4:N+4)) + &
                          A(-5)*(B(:,-5:N-5)-B(:,5:N+5)) 

    CASE (11)
      C = A(0)*B(:,0:N) + A(-1)*B(:,-1:N-1) + A(-2)*B(:,-2:N-2) + &
                          A(-3)*B(:,-3:N-3) + A(-4)*B(:,-4:N-4) + &
                          A(-5)*B(:,-5:N-5) + A(-6)*B(:,-6:N-6) + &
                          A( 1)*B(:, 1:N+1) + A( 2)*B(:, 2:N+2) + &
                          A( 3)*B(:, 3:N+3) + A( 4)*B(:, 4:N+4) + &
                          A( 5)*B(:, 5:N+5) 

      D = C + A(-6)      *(B(:,6:N+6)-B(:,-6:N-6)) + &
             (A(-5)-A(5))*(B(:,5:N+5)-B(:,-5:N-5)) + & 
             (A(-4)-A(4))*(B(:,4:N+4)-B(:,-4:N-4)) + & 
             (A(-3)-A(3))*(B(:,3:N+3)-B(:,-3:N-3)) + & 
             (A(-2)-A(2))*(B(:,2:N+2)-B(:,-2:N-2)) + &
             (A(-1)-A(1))*(B(:,1:N+1)-B(:,-1:N-1)) 

  END SELECT

  if (MOD(Order,2) /= 0) RETURN

  D = -C + TWO*B(:,0:N)

  END Subroutine Direct_Sum_2D_y
!
! ======================================================================
!
  Subroutine Direct_Sum_3D_x (N2, Order, L, R, A, B, C, D)  

  integer  :: N2, Order, L, R

  REALTYPE, dimension(-(Order+1)/2:(Order+1)/2) :: A 
  REALTYPE, dimension(-L:,:,:)                  :: B
  REALTYPE, dimension(  :,:,:)                  :: C, D

  integer :: i, N

  N = N2-1

  SELECT CASE (Order)
    CASE (1)
      C = A(0)*B(0:N,:,:) + A(-1)*B(-1:N-1,:,:) 

      D = C + A(-1)      *(B(1:N+1,:,:)-B(-1:N-1,:,:)) 

    CASE (2)
      C = A(0)*B(0:N,:,:) + A(-1)*(B(-1:N-1,:,:)-B(1:N+1,:,:))

    CASE (3)
      C = A(0)*B(0:N,:,:) + A(-1)*B(-1:N-1,:,:) + A(-2)*B(-2:N-2,:,:) + &
                            A( 1)*B( 1:N+1,:,:)
      
      D = C + A(-2)      *(B(2:N+2,:,:)-B(-2:N-2,:,:))                + &
             (A(-1)-A(1))*(B(1:N+1,:,:)-B(-1:N-1,:,:)) 

    CASE (4)
      C = A(0)*B(0:N,:,:) + A(-1)*(B(-1:N-1,:,:)-B(1:N+1,:,:))        + &
                            A(-2)*(B(-2:N-2,:,:)-B(2:N+2,:,:)) 

    CASE (5)
      C = A(0)*B(0:N,:,:) + A(-1)*B(-1:N-1,:,:) + A(-2)*B(-2:N-2,:,:) + &
                            A(-3)*B(-3:N-3,:,:)                       + &
                            A( 1)*B( 1:N+1,:,:) + A( 2)*B( 2:N+2,:,:) 
      
      D = C + A(-3)      *(B(3:N+3,:,:)-B(-3:N-3,:,:))                + &
             (A(-2)-A(2))*(B(2:N+2,:,:)-B(-2:N-2,:,:))                + &
             (A(-1)-A(1))*(B(1:N+1,:,:)-B(-1:N-1,:,:)) 

    CASE (6)
      C = A(0)*B(0:N,:,:) + A(-1)*(B(-1:N-1,:,:)-B(1:N+1,:,:))        + &
                            A(-2)*(B(-2:N-2,:,:)-B(2:N+2,:,:))        + & 
                            A(-3)*(B(-3:N-3,:,:)-B(3:N+3,:,:)) 

    CASE (7)
      C = A(0)*B(0:N,:,:) + A(-1)*B(-1:N-1,:,:) + A(-2)*B(-2:N-2,:,:) + &
                            A(-3)*B(-3:N-3,:,:) + A(-4)*B(-4:N-4,:,:) + &
                            A( 1)*B( 1:N+1,:,:) + A( 2)*B( 2:N+2,:,:) + &
                            A( 3)*B( 3:N+3,:,:) 

      D = C + A(-4)      *(B(4:N+4,:,:)-B(-4:N-4,:,:)) + &
             (A(-3)-A(3))*(B(3:N+3,:,:)-B(-3:N-3,:,:)) + & 
             (A(-2)-A(2))*(B(2:N+2,:,:)-B(-2:N-2,:,:)) + &
             (A(-1)-A(1))*(B(1:N+1,:,:)-B(-1:N-1,:,:)) 

    CASE (8)
      C = A(0)*B(0:N,:,:) + A(-1)*(B(-1:N-1,:,:)-B(1:N+1,:,:)) + &
                            A(-2)*(B(-2:N-2,:,:)-B(2:N+2,:,:)) + & 
                            A(-3)*(B(-3:N-3,:,:)-B(3:N+3,:,:)) + &
                            A(-4)*(B(-4:N-4,:,:)-B(4:N+4,:,:)) 

    CASE (9)
      C = A(0)*B(0:N,:,:) + A(-1)*B(-1:N-1,:,:) + A(-2)*B(-2:N-2,:,:) + &
                            A(-3)*B(-3:N-3,:,:) + A(-4)*B(-4:N-4,:,:) + &
                            A(-5)*B(-5:N-5,:,:)                       + &
                            A( 1)*B( 1:N+1,:,:) + A( 2)*B( 2:N+2,:,:) + &
                            A( 3)*B( 3:N+3,:,:) + A( 4)*B( 4:N+4,:,:) 

      D = C + A(-5)      *(B(5:N+5,:,:)-B(-5:N-5,:,:)) + &
             (A(-4)-A(4))*(B(4:N+4,:,:)-B(-4:N-4,:,:)) + & 
             (A(-3)-A(3))*(B(3:N+3,:,:)-B(-3:N-3,:,:)) + & 
             (A(-2)-A(2))*(B(2:N+2,:,:)-B(-2:N-2,:,:)) + &
             (A(-1)-A(1))*(B(1:N+1,:,:)-B(-1:N-1,:,:)) 

    CASE (10)
      C = A(0)*B(0:N,:,:) + A(-1)*(B(-1:N-1,:,:)-B(1:N+1,:,:)) + &
                            A(-2)*(B(-2:N-2,:,:)-B(2:N+2,:,:)) + & 
                            A(-3)*(B(-3:N-3,:,:)-B(3:N+3,:,:)) + &
                            A(-4)*(B(-4:N-4,:,:)-B(4:N+4,:,:)) + &
                            A(-5)*(B(-5:N-5,:,:)-B(5:N+5,:,:)) 

    CASE (11)
      C = A(0)*B(0:N,:,:) + A(-1)*B(-1:N-1,:,:) + A(-2)*B(-2:N-2,:,:) + &
                            A(-3)*B(-3:N-3,:,:) + A(-4)*B(-4:N-4,:,:) + &
                            A(-5)*B(-5:N-5,:,:) + A(-6)*B(-6:N-6,:,:) + &
                            A( 1)*B( 1:N+1,:,:) + A( 2)*B( 2:N+2,:,:) + &
                            A( 3)*B( 3:N+3,:,:) + A( 4)*B( 4:N+4,:,:) + &
                            A( 5)*B( 5:N+5,:,:) 

      D = C + A(-6)      *(B(6:N+6,:,:)-B(-6:N-6,:,:)) + &
             (A(-5)-A(5))*(B(5:N+5,:,:)-B(-5:N-5,:,:)) + & 
             (A(-4)-A(4))*(B(4:N+4,:,:)-B(-4:N-4,:,:)) + & 
             (A(-3)-A(3))*(B(3:N+3,:,:)-B(-3:N-3,:,:)) + & 
             (A(-2)-A(2))*(B(2:N+2,:,:)-B(-2:N-2,:,:)) + &
             (A(-1)-A(1))*(B(1:N+1,:,:)-B(-1:N-1,:,:)) 

  END SELECT

  if (MOD(Order,2) /= 0) RETURN

  D = -C + TWO*B(0:N,:,:)

  END Subroutine Direct_Sum_3D_x
!
! ======================================================================
!
  Subroutine Direct_Sum_3D_y (N2, Order, L, R, A, B, C, D)  

  integer  :: N2, Order, L, R

  REALTYPE, dimension(-(Order+1)/2:(Order+1)/2) :: A 
  REALTYPE, dimension(:,-L:,:)                  :: B
  REALTYPE, dimension(:,  :,:)                  :: C, D

  integer :: i, N

  N = N2-1

  SELECT CASE (Order)
    CASE (1)
      C = A(0)*B(:,0:N,:) + A(-1)*B(:,-1:N-1,:) 

      D = C + A(-1)      *(B(:,1:N+1,:)-B(:,-1:N-1,:)) 

    CASE (2)
      C = A(0)*B(:,0:N,:) + A(-1)*(B(:,-1:N-1,:)-B(:,1:N+1,:))

    CASE (3)
      C = A(0)*B(:,0:N,:) + A(-1)*B(:,-1:N-1,:) + A(-2)*B(:,-2:N-2,:) + &
                            A( 1)*B(:, 1:N+1,:)
      
      D = C + A(-2)      *(B(:,2:N+2,:)-B(:,-2:N-2,:))                + &
             (A(-1)-A(1))*(B(:,1:N+1,:)-B(:,-1:N-1,:)) 

    CASE (4)
      C = A(0)*B(:,0:N,:) + A(-1)*(B(:,-1:N-1,:)-B(:,1:N+1,:))        + &
                            A(-2)*(B(:,-2:N-2,:)-B(:,2:N+2,:)) 

    CASE (5)
      C = A(0)*B(:,0:N,:) + A(-1)*B(:,-1:N-1,:) + A(-2)*B(:,-2:N-2,:) + &
                            A(-3)*B(:,-3:N-3,:)                       + &
                            A( 1)*B(:, 1:N+1,:) + A( 2)*B(:, 2:N+2,:) 
      
      D = C + A(-3)      *(B(:,3:N+3,:)-B(:,-3:N-3,:))                + &
             (A(-2)-A(2))*(B(:,2:N+2,:)-B(:,-2:N-2,:))                + &
             (A(-1)-A(1))*(B(:,1:N+1,:)-B(:,-1:N-1,:)) 

    CASE (6)
      C = A(0)*B(:,0:N,:) + A(-1)*(B(:,-1:N-1,:)-B(:,1:N+1,:))        + &
                            A(-2)*(B(:,-2:N-2,:)-B(:,2:N+2,:))        + & 
                            A(-3)*(B(:,-3:N-3,:)-B(:,3:N+3,:)) 

    CASE (7)
      C = A(0)*B(:,0:N,:) + A(-1)*B(:,-1:N-1,:) + A(-2)*B(:,-2:N-2,:) + &
                            A(-3)*B(:,-3:N-3,:) + A(-4)*B(:,-4:N-4,:) + &
                            A( 1)*B(:, 1:N+1,:) + A( 2)*B(:, 2:N+2,:) + &
                            A( 3)*B(:, 3:N+3,:) 

      D = C + A(-4)      *(B(:,4:N+4,:)-B(:,-4:N-4,:)) + &
             (A(-3)-A(3))*(B(:,3:N+3,:)-B(:,-3:N-3,:)) + & 
             (A(-2)-A(2))*(B(:,2:N+2,:)-B(:,-2:N-2,:)) + &
             (A(-1)-A(1))*(B(:,1:N+1,:)-B(:,-1:N-1,:)) 

    CASE (8)
      C = A(0)*B(:,0:N,:) + A(-1)*(B(:,-1:N-1,:)-B(:,1:N+1,:)) + &
                            A(-2)*(B(:,-2:N-2,:)-B(:,2:N+2,:)) + & 
                            A(-3)*(B(:,-3:N-3,:)-B(:,3:N+3,:)) + &
                            A(-4)*(B(:,-4:N-4,:)-B(:,4:N+4,:)) 

    CASE (9)
      C = A(0)*B(:,0:N,:) + A(-1)*B(:,-1:N-1,:) + A(-2)*B(:,-2:N-2,:) + &
                            A(-3)*B(:,-3:N-3,:) + A(-4)*B(:,-4:N-4,:) + &
                            A(-5)*B(:,-5:N-5,:)                       + &
                            A( 1)*B(:, 1:N+1,:) + A( 2)*B(:, 2:N+2,:) + &
                            A( 3)*B(:, 3:N+3,:) + A( 4)*B(:, 4:N+4,:) 

      D = C + A(-5)      *(B(:,5:N+5,:)-B(:,-5:N-5,:)) + &
             (A(-4)-A(4))*(B(:,4:N+4,:)-B(:,-4:N-4,:)) + & 
             (A(-3)-A(3))*(B(:,3:N+3,:)-B(:,-3:N-3,:)) + & 
             (A(-2)-A(2))*(B(:,2:N+2,:)-B(:,-2:N-2,:)) + &
             (A(-1)-A(1))*(B(:,1:N+1,:)-B(:,-1:N-1,:)) 

    CASE (10)
      C = A(0)*B(:,0:N,:) + A(-1)*(B(:,-1:N-1,:)-B(:,1:N+1,:)) + &
                            A(-2)*(B(:,-2:N-2,:)-B(:,2:N+2,:)) + & 
                            A(-3)*(B(:,-3:N-3,:)-B(:,3:N+3,:)) + &
                            A(-4)*(B(:,-4:N-4,:)-B(:,4:N+4,:)) + &
                            A(-5)*(B(:,-5:N-5,:)-B(:,5:N+5,:)) 

    CASE (11)
      C = A(0)*B(:,0:N,:) + A(-1)*B(:,-1:N-1,:) + A(-2)*B(:,-2:N-2,:) + &
                            A(-3)*B(:,-3:N-3,:) + A(-4)*B(:,-4:N-4,:) + &
                            A(-5)*B(:,-5:N-5,:) + A(-6)*B(:,-6:N-6,:) + &
                            A( 1)*B(:, 1:N+1,:) + A( 2)*B(:, 2:N+2,:) + &
                            A( 3)*B(:, 3:N+3,:) + A( 4)*B(:, 4:N+4,:) + &
                            A( 5)*B(:, 5:N+5,:) 

      D = C + A(-6)      *(B(:,6:N+6,:)-B(:,-6:N-6,:)) + &
             (A(-5)-A(5))*(B(:,5:N+5,:)-B(:,-5:N-5,:)) + & 
             (A(-4)-A(4))*(B(:,4:N+4,:)-B(:,-4:N-4,:)) + & 
             (A(-3)-A(3))*(B(:,3:N+3,:)-B(:,-3:N-3,:)) + & 
             (A(-2)-A(2))*(B(:,2:N+2,:)-B(:,-2:N-2,:)) + &
             (A(-1)-A(1))*(B(:,1:N+1,:)-B(:,-1:N-1,:)) 

  END SELECT

  if (MOD(Order,2) /= 0) RETURN

  D = -C + TWO*B(:,0:N,:)

  END Subroutine Direct_Sum_3D_y
!
! ======================================================================
!
  Subroutine Direct_Sum_3D_z (N2, Order, L, R, A, B, C, D)  

  integer  :: N2, Order, L, R

  REALTYPE, dimension(-(Order+1)/2:(Order+1)/2) :: A 
  REALTYPE, dimension(:,:,-L:)                  :: B
  REALTYPE, dimension(:,:,  :)                  :: C, D

  integer :: i, N

  N = N2-1

  SELECT CASE (Order)
    CASE (1)
      C = A(0)*B(:,:,0:N) + A(-1)*B(:,:,-1:N-1) 

      D = C + A(-1)      *(B(:,:,1:N+1)-B(:,:,-1:N-1)) 

    CASE (2)
      C = A(0)*B(:,:,0:N) + A(-1)*(B(:,:,-1:N-1)-B(:,:,1:N+1))

    CASE (3)
      C = A(0)*B(:,:,0:N) + A(-1)*B(:,:,-1:N-1) + A(-2)*B(:,:,-2:N-2) + &
                            A( 1)*B(:,:, 1:N+1)
      
      D = C + A(-2)      *(B(:,:,2:N+2)-B(:,:,-2:N-2))                + &
             (A(-1)-A(1))*(B(:,:,1:N+1)-B(:,:,-1:N-1)) 

    CASE (4)
      C = A(0)*B(:,:,0:N) + A(-1)*(B(:,:,-1:N-1)-B(:,:,1:N+1))        + &
                            A(-2)*(B(:,:,-2:N-2)-B(:,:,2:N+2)) 

    CASE (5)
      C = A(0)*B(:,:,0:N) + A(-1)*B(:,:,-1:N-1) + A(-2)*B(:,:,-2:N-2) + &
                            A(-3)*B(:,:,-3:N-3)                       + &
                            A( 1)*B(:,:, 1:N+1) + A( 2)*B(:,:, 2:N+2) 
      
      D = C + A(-3)      *(B(:,:,3:N+3)-B(:,:,-3:N-3))                + &
             (A(-2)-A(2))*(B(:,:,2:N+2)-B(:,:,-2:N-2))                + &
             (A(-1)-A(1))*(B(:,:,1:N+1)-B(:,:,-1:N-1)) 

    CASE (6)
      C = A(0)*B(:,:,0:N) + A(-1)*(B(:,:,-1:N-1)-B(:,:,1:N+1))        + &
                            A(-2)*(B(:,:,-2:N-2)-B(:,:,2:N+2))        + & 
                            A(-3)*(B(:,:,-3:N-3)-B(:,:,3:N+3)) 

    CASE (7)
      C = A(0)*B(:,:,0:N) + A(-1)*B(:,:,-1:N-1) + A(-2)*B(:,:,-2:N-2) + &
                            A(-3)*B(:,:,-3:N-3) + A(-4)*B(:,:,-4:N-4) + &
                            A( 1)*B(:,:, 1:N+1) + A( 2)*B(:,:, 2:N+2) + &
                            A( 3)*B(:,:, 3:N+3) 

      D = C + A(-4)      *(B(:,:,4:N+4)-B(:,:,-4:N-4)) + &
             (A(-3)-A(3))*(B(:,:,3:N+3)-B(:,:,-3:N-3)) + & 
             (A(-2)-A(2))*(B(:,:,2:N+2)-B(:,:,-2:N-2)) + &
             (A(-1)-A(1))*(B(:,:,1:N+1)-B(:,:,-1:N-1)) 

    CASE (8)
      C = A(0)*B(:,:,0:N) + A(-1)*(B(:,:,-1:N-1)-B(:,:,1:N+1)) + &
                            A(-2)*(B(:,:,-2:N-2)-B(:,:,2:N+2)) + & 
                            A(-3)*(B(:,:,-3:N-3)-B(:,:,3:N+3)) + &
                            A(-4)*(B(:,:,-4:N-4)-B(:,:,4:N+4)) 

    CASE (9)
      C = A(0)*B(:,:,0:N) + A(-1)*B(:,:,-1:N-1) + A(-2)*B(:,:,-2:N-2) + &
                            A(-3)*B(:,:,-3:N-3) + A(-4)*B(:,:,-4:N-4) + &
                            A(-5)*B(:,:,-5:N-5)                       + &
                            A( 1)*B(:,:, 1:N+1) + A( 2)*B(:,:, 2:N+2) + &
                            A( 3)*B(:,:, 3:N+3) + A( 4)*B(:,:, 4:N+4) 

      D = C + A(-5)      *(B(:,:,5:N+5)-B(:,:,-5:N-5)) + &
             (A(-4)-A(4))*(B(:,:,4:N+4)-B(:,:,-4:N-4)) + & 
             (A(-3)-A(3))*(B(:,:,3:N+3)-B(:,:,-3:N-3)) + & 
             (A(-2)-A(2))*(B(:,:,2:N+2)-B(:,:,-2:N-2)) + &
             (A(-1)-A(1))*(B(:,:,1:N+1)-B(:,:,-1:N-1)) 

    CASE (10)
      C = A(0)*B(:,:,0:N) + A(-1)*(B(:,:,-1:N-1)-B(:,:,1:N+1)) + &
                            A(-2)*(B(:,:,-2:N-2)-B(:,:,2:N+2)) + & 
                            A(-3)*(B(:,:,-3:N-3)-B(:,:,3:N+3)) + &
                            A(-4)*(B(:,:,-4:N-4)-B(:,:,4:N+4)) + &
                            A(-5)*(B(:,:,-5:N-5)-B(:,:,5:N+5)) 

    CASE (11)
      C = A(0)*B(:,:,0:N) + A(-1)*B(:,:,-1:N-1) + A(-2)*B(:,:,-2:N-2) + &
                            A(-3)*B(:,:,-3:N-3) + A(-4)*B(:,:,-4:N-4) + &
                            A(-5)*B(:,:,-5:N-5) + A(-6)*B(:,:,-6:N-6) + &
                            A( 1)*B(:,:, 1:N+1) + A( 2)*B(:,:, 2:N+2) + &
                            A( 3)*B(:,:, 3:N+3) + A( 4)*B(:,:, 4:N+4) + &
                            A( 5)*B(:,:, 5:N+5) 

      D = C + A(-6)      *(B(:,:,6:N+6)-B(:,:,-6:N-6)) + &
             (A(-5)-A(5))*(B(:,:,5:N+5)-B(:,:,-5:N-5)) + & 
             (A(-4)-A(4))*(B(:,:,4:N+4)-B(:,:,-4:N-4)) + & 
             (A(-3)-A(3))*(B(:,:,3:N+3)-B(:,:,-3:N-3)) + & 
             (A(-2)-A(2))*(B(:,:,2:N+2)-B(:,:,-2:N-2)) + &
             (A(-1)-A(1))*(B(:,:,1:N+1)-B(:,:,-1:N-1)) 

  END SELECT

  if (MOD(Order,2) /= 0) RETURN

  D = -C + TWO*B(:,:,0:N)

  END Subroutine Direct_Sum_3D_z
#if defined (DEBUG)
!
! ======================================================================
!
  Subroutine Check_Parameter       (Order, N0, N5, i0, i1, N2, R)

  integer  :: Order, N0, N5, i0, i1, N2, R
  logical  :: Error

  integer  :: ii, L1

  ii = 0 ; if (MOD(Order,2) == 1) ii = 1

  Error = .FALSE.

  L1 = i0-Order
  if (N0 >  L1) then
    write (6,*) "N0 > i0-Order           ", N0, i0, Order
    Error = .TRUE.
  endif

  L1 = i0+2*(N2-1+R+ii)+1
  if (N5 <  L1) then
    write (6,*) "N5 < i0+2*(N2-1+R+ii)+1 ", N5, i0, N2, R, ii, L1
    Error = .TRUE.
  endif

  L1 = i0+2*(N2-1)+1
  if (N5 <  L1) then
    write (6,*) "N5 < i0+2*(N2-1)+1      ", N5, i0, N2, L1
    Error = .TRUE.
  endif

  if (Error) STOP

  END Subroutine Check_Parameter
#endif
!
! ==================================================================
!
  Subroutine PS_WENO_MR_Index_Update (WENO_Order, MR_Order, N0, N2, N3, N5, &
                                      Symmetry, N_Ghost)

  integer                            :: WENO_Order, MR_Order, N0, N2, N3, N5
  integer ,                 OPTIONAL :: Symmetry
  integer , dimension(0:1), OPTIONAL :: N_Ghost
  integer                            :: N_Ghost_Point

  if (MOD(MR_Order,2) == 0) then
    N_Ghost_Point = MAX(MR_Order  , (WENO_Order+1)/2)
  else
    N_Ghost_Point = MAX(MR_Order+1, (WENO_Order+1)/2)
  endif

  N0 = N2 - N_Ghost_Point
  N5 = N3 + N_Ghost_Point 

  if (MOD(N3-N2+1,2) == 1) N0 = N0-1
  if (MOD(N3-N2+1,2) == 1) N5 = N5+1

  if (PRESENT(Symmetry)) then
    if ((Symmetry /= 1) .AND. (MOD(N3-N2+1,2) == 1)) N0 = N0+1
  endif

  if (PRESENT(N_Ghost)) then
    N_Ghost(0) = N2-N0 ; N_Ghost(1) = N5-N3
  endif

  END Subroutine PS_WENO_MR_Index_Update
!
! ==================================================================
!
  Subroutine PS_MR_Adjust_Ghost_Index (MR_Order, N2, N3, N0, N5)

  integer  :: MR_Order, N2, N3, N0, N5
  integer  :: MR_Ghost

  MR_Ghost = MR_Order ; if (MOD(MR_ORDER,2) == 1) MR_Ghost = MR_Order+1

  N0 = MIN(N2 - MR_Ghost, N0) - 1
  N5 = MAX(N3 + MR_Ghost, N5) + 1 

  END Subroutine PS_MR_Adjust_Ghost_Index 
  
END Module Multi_Resolution
