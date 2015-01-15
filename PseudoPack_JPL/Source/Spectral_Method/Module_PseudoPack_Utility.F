#include "pseudopack.h"


MODULE PseudoPack_Utility
  USE PseudoPack_Type
  USE PS_IO_Unit

#if defined (LOCAL_GEMM)
  USE General_Matrix_Multiply
#endif

implicit NONE

INTERFACE PS_Get_Dimension
  MODULE PROCEDURE PS_Get_Dimension_0
  MODULE PROCEDURE PS_Get_Dimension_1
  MODULE PROCEDURE PS_Get_Dimension_2
  MODULE PROCEDURE PS_Get_Dimension_3
  MODULE PROCEDURE PS_Get_Dimension_4
END INTERFACE

PRIVATE

PUBLIC  :: PS_Property
PUBLIC  :: PS_Grid_Index

PUBLIC  :: PS_Get_Grid_Size

PUBLIC  :: PS_Get_Dimension

PUBLIC  :: PS_Get_Matrix_Size

PUBLIC  :: PS_Get_Matrix
PUBLIC  :: PS_Put_Matrix

PUBLIC  :: PS_Change_Filter
PUBLIC  :: PS_Filter_Switch

PUBLIC  :: PS_Chebyshev_Lobatto_a_k
PUBLIC  :: PS_Chebyshev_Lobatto_T_k
PUBLIC  :: PS_Interpolation_Sum_a_T

PUBLIC  :: PS_Index_Search

CONTAINS

  Subroutine PS_Get_Grid_Size (N_Pnts, Method, Point_Type, Symmetry, N)

  integer  :: N_Pnts, Method, Point_Type, Symmetry, N, NN

  NN = N_Pnts

  if (Method == 0) then
    if ((Point_Type /= 2) .AND. (MOD(NN,2) /= 0)) NN = NN-1
    if ((Point_Type == 2) .AND. (MOD(NN,2) /= 1)) NN = NN-1
  endif

  if (Symmetry == 1) then
    SELECT CASE (Method)
      CASE (0)                     ! Fourier Method
        N = _H(NN+1)

      CASE (1,2)                   ! Chebyshev/Legendre Method
        N = _H(NN)

    END SELECT
  else
    N = NN
  endif

  END Subroutine PS_Get_Grid_Size
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Dimension (D, N)

  integer                :: N
  REALTYPE, dimension(*) :: D

  integer  :: Method, Point_Type, Symmetry, N_Pnts

  Method     = OPERATOR_METHOD(D)
  Point_Type = OPERATOR_POINT_TYPE(D)
  Symmetry   = OPERATOR_SYMMETRY(D)
  N_Pnts     = OPERATOR_N_PNTS(D)

  N = N_Pnts

  if (Method == 0) then
    if ((Point_Type /= 2) .AND. (MOD(N,2) /= 0)) N = N-1
    if ((Point_Type == 2) .AND. (MOD(N,2) /= 1)) N = N-1
  endif

  if (Symmetry == 1) then
    SELECT CASE (Method)
      CASE (0)                     ! Fourier Method
        N = _H(N+1)

      CASE (1,2)                   ! Chebyshev/Legendre Method
        N = _H(N)

    END SELECT
  endif

  END Subroutine PS_Dimension
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Get_Dimension_0 (N, N_Pnts, Method, Point_Type, Symmetry)

  integer            :: N_Pnts, N
  integer            :: Method, Point_Type
  integer , OPTIONAL :: Symmetry
  logical            :: L_Symmetry

#if defined (DEBUG)
  integer            :: IO_Error

  IO_Error = 0

  call Data_Integrity (lid6) ; call Data_Integrity (lid99)

  if (IO_Error == 1) STOP
#endif

  N = N_Pnts

  if (Method == 0) then
    if ((Point_Type /= 2) .AND. (MOD(N,2) /= 0)) N = N-1
    if ((Point_Type == 2) .AND. (MOD(N,2) /= 1)) N = N-1
  endif

  if (Point_Type == 2) RETURN

  L_Symmetry = PRESENT(Symmetry)

  if (L_Symmetry) then
    if (Symmetry == 0) RETURN

    SELECT CASE (Method)
      CASE (0)                     ! Fourier Method
        N = _H(N+1)

      CASE (1,2)                   ! Chebyshev/Legendre Method
        N = _H(N)

    END SELECT
  endif

#if defined (DEBUG)
  CONTAINS

    Subroutine Data_Integrity (lid)

    integer   :: lid

    if (N_Pnts <= 0) then
      write (lid,*) 'Error : N_Pnts = ', N_Pnts,' <= 0'
      IO_Error = 1
    endif

    if ((Method <  0) .OR. (Method >  2)) then
      write (lid,*) 'Error : Method = ', Method
      write (lid,*) '        Method = 0, Fourier   Method'
      write (lid,*) '               = 1, Chebyshev Method'
      write (lid,*) '               = 2, Legendre  Method'
      write (lid,*) '----------------------------------------'
      IO_Error = 1
    endif

    if ((Point_Type <  1) .OR. (Point_Type >  3)) then
      write (lid,*) 'Error : Point_Type = ', Point_Type
      write (lid,*) '        Point_Type = 1, Lobatto Point'
      write (lid,*) '                   = 2, Radau Point'
      write (lid,*) '                   = 3, Gauss Point'
      write (lid,*) '----------------------------------------'
      IO_Error = 1
    endif

    if (Method == 0) then
      if ((Point_Type == 1) .AND. (MOD(N,2) /= 0)) then
        write (lid,*) 'Error : Method     = ', Method,' Fourier'
        write (lid,*) '        Point_Type = ', Point_Type,' Lobatto'
        write (lid,*) '                 N = ', N,' must be an EVEN number'
        write (lid,*) '----------------------------------------'
        IO_Error = 2
      endif

      if ((Point_Type == 2) .AND. (MOD(N,2) /= 1)) then
        write (lid,*) 'Error : Method     = ', Method,' Fourier'
        write (lid,*) '        Point_Type = ', Point_Type,' Radau'
        write (lid,*) '                 N = ', N,' must be an ODD  number'
        write (lid,*) '----------------------------------------'
        IO_Error = 2
      endif
    endif
   
    if (L_Symmetry) then
      if ((Symmetry <  0) .OR. (Symmetry >  1)) then
        write (lid,*) 'Error : Symmetry = ', Symmetry
        write (lid,*) '        Symmetry = 0, No Symmetry'
        write (lid,*) '                 = 1, Even/Odd Symmetric Data'
        write (lid,*) '----------------------------------------'
        IO_Error = 1
      endif
    endif

    if (IO_Error /= 0) then
      write (lid,*) 
      write (lid,*) 'Arguments Error in PS_Get_Dimension_0 Subroutine'
      write (lid,*) '==============================================='
    endif

    END Subroutine Data_Integrity
#endif

  END Subroutine PS_Get_Dimension_0
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Get_Dimension_1 (N, Property, Grid_Index)

  integer  :: N
  integer  :: Method, Point_Type, Symmetry

  TYPE (PS_Property)   :: Property
  TYPE (PS_Grid_Index) :: Grid_Index

  N = Grid_Index%N

  Method     = Property%Method
  Point_Type = Property%Point_Type
  Symmetry   = Property%Symmetry

  if (Method == 0) then
    if ((Point_Type /= 2) .AND. (MOD(N,2) /= 0)) N = N-1
    if ((Point_Type == 2) .AND. (MOD(N,2) /= 1)) N = N-1
  endif

  if (Point_Type == 2) RETURN
  if (Symmetry   == 0) RETURN

  SELECT CASE (Method)
    CASE (0)                     ! Fourier Method
      N = _H(N+1)

    CASE (1,2)                   ! Chebyshev/Legendre Method
      N = _H(N)

  END SELECT

  END Subroutine PS_Get_Dimension_1 
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Get_Dimension_2 (D_1, N_1)

  integer                :: N_1
  REALTYPE, dimension(*) :: D_1

  call PS_Dimension (D_1, N_1)

  END Subroutine PS_Get_Dimension_2
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Get_Dimension_3 (D_1, D_2, N_1, N_2, Ordered)

  integer                :: N_1, N_2
  REALTYPE, dimension(*) :: D_1, D_2
  integer , OPTIONAL     :: Ordered

  integer                :: Dimension_1, Dimension_2, N

  call PS_Dimension (D_1, N_1)
  call PS_Dimension (D_2, N_2)

  if (.NOT. PRESENT(Ordered)) RETURN
  if (Ordered == 0)           RETURN

  Dimension_1 = OPERATOR_RANK(D_1)
  Dimension_2 = OPERATOR_RANK(D_2)
 
#if defined (DEBUG)
  if (Dimension_1 == Dimension_2) then
    write (lid6,100) Dimension_1 ; write (lid99,100) Dimension_1 ; STOP
  endif
#endif

  if (Dimension_2 <  Dimension_1) then ; N = N_1 ; N_1 = N_2 ; N_2 = N ; endif

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_Dimension '/  &
             15x,'D_1 and D_2 are both operating on same dimension = ',i5/ &
             15x,'Cannot Reorder the extents for a 2D data array!'/ &
             1x,60('*'))
#endif

  END Subroutine PS_Get_Dimension_3
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Get_Dimension_4 (D_1, D_2, D_3, N_1, N_2, N_3, Ordered)

  integer                :: N_1, N_2, N_3
  REALTYPE, dimension(*) :: D_1, D_2, D_3
  integer , OPTIONAL     :: Ordered

  logical                :: Error
  integer                :: Dimension_1, Dimension_2, Dimension_3, N

  call PS_Dimension (D_1, N_1)
  call PS_Dimension (D_2, N_2)
  call PS_Dimension (D_3, N_3)

  if (.NOT. PRESENT(Ordered)) RETURN
  if (Ordered == 0)           RETURN

  Dimension_1 = OPERATOR_RANK(D_1)
  Dimension_2 = OPERATOR_RANK(D_2)
  Dimension_3 = OPERATOR_RANK(D_3)
 
#if defined (DEBUG)
  Error = .FALSE.

  if (Dimension_1 == Dimension_2) then
    write (lid6,100) Dimension_1 ; write (lid99,100) Dimension_1 
    Error = .TRUE.
  endif

  if (Dimension_1 == Dimension_3) then
    write (lid6,101) Dimension_1 ; write (lid99,101) Dimension_1 
    Error = .TRUE.
  endif

  if (Dimension_2 == Dimension_3) then
    write (lid6,102) Dimension_2 ; write (lid99,102) Dimension_2 
    Error = .TRUE.
  endif

  if (Error) STOP
#endif

  if (Dimension_2 <  Dimension_1) then ; N = N_2 ; N_2 = N_1 ; N_1 = N ; endif
  if (Dimension_3 <  Dimension_1) then ; N = N_3 ; N_3 = N_1 ; N_1 = N ; endif
  if (Dimension_2 <  Dimension_3) then ; N = N_3 ; N_3 = N_2 ; N_2 = N ; endif

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_Dimension '/  &
             15x,'D_1 and D_2 are both operating on same dimension = ',i5/ &
             15x,'Cannot find the extents of a three dimensional data array!'/ &
             1x,60('*'))
 101 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_Dimension '/  &
             15x,'D_1 and D_3 are both operating on same dimension = ',i5/ &
             15x,'Cannot find the extents of a three dimensional data array!'/ &
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_Dimension '/  &
             15x,'D_2 and D_3 are both operating on same dimension = ',i5/ &
             15x,'Cannot find the extents of a three dimensional data array!'/ &
             1x,60('*'))
#endif

  END Subroutine PS_Get_Dimension_4


!
! -----------------------------------------------------------------------
!
  Subroutine PS_Get_Matrix_Size (Operator, N, Error)

  integer            :: N
  logical , OPTIONAL :: Error

  REALTYPE, dimension(*) :: Operator

  integer  :: Method, Algorithm, N_Pnts, N1

  if (PRESENT(Error)) Error = .FALSE.

  Method     = INT(Operator(_I_METHOD))
  Algorithm  = INT(Operator(_I_ALGORITHM))
  N_Pnts     = INT(Operator(_I_N_PNTS))

  SELECT CASE (Algorithm)
    CASE (0)
      N = N_Pnts

    CASE (1)
      SELECT CASE (Method)
        CASE (0)
          N1 = N_Pnts+1 ; N = _H(N1)

        CASE (1,2)
          N1 = N_Pnts   ; N = _H(N1)
 
      END SELECT
         
    CASE DEFAULT
      write (lid6,100) ; write (lid99,100) 

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif

  END SELECT

 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_Matrix_Size '/  &
             15x,'Unable to Extract the Matrix Size from Operator.'/&
             15x,'ONLY Matrix-Matric Multiply (MXM) Algorithm = 0 or '/ &
             15x,'     Even-Odd Decomposition (EOD) Algorithm = 1 is valid.'/ &
             1x,60('*'))

  END Subroutine PS_Get_Matrix_Size
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Get_Matrix (Operator, D, D_Even, D_Odd, Order, N, Error)

  integer , OPTIONAL :: N
  integer , OPTIONAL :: Order
  logical , OPTIONAL :: Error

  REALTYPE, dimension(*)             :: Operator
  REALTYPE, dimension(:,:), OPTIONAL :: D, D_Even, D_Odd

  logical  :: L_D, L_D_Even, L_D_Odd

  integer  :: Algorithm, D_Or_S

  if (PRESENT(Error)) Error = .FALSE.

  L_D = PRESENT(D) ; L_D_Even = PRESENT(D_Even) ; L_D_Odd  = PRESENT(D_Odd)

#if 1
  if (.NOT. (L_D .OR. L_D_Even .OR. L_D_Odd)) then
#else
  if (.NOT. (L_D .OR. (L_D_Even .AND. L_D_Odd))) then
#endif
    write (lid6,102) ; write (lid99,102) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif

  D_Or_S    = INT(Operator(_I_CASE))
  Algorithm = INT(Operator(_I_ALGORITHM))

  SELECT CASE (D_Or_S)
    CASE (0)                  ! Smoothing
      SELECT CASE (Algorithm)
        CASE (0)
          call PS_Get_MXM_S (Operator, D, N, Error)

        CASE (1)
          call PS_Get_EOD_S (Operator, D_Even, D_Odd, N, Error)

        CASE DEFAULT
          write (lid6,100) ; write (lid99,100) 

          if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN  
                              else ; STOP ; endif

      END SELECT

    CASE (1)                  ! Differentiation
      SELECT CASE (Algorithm)
        CASE (0)
          call PS_Get_MXM_D (Operator, D, Order, N, Error)

        CASE (1)
          call PS_Get_EOD_D (Operator, D_Even, D_Odd, Order, N, Error)

        CASE DEFAULT
          write (lid6,101) ; write (lid99,101) 

          if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN  
                              else ; STOP ; endif

      END SELECT

  END SELECT

 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_Matrix '/  &
             15x,'Unable to Extract the Smoothing Matrix from Operator.'/&
             15x,'ONLY Matrix-Matric Multiply (MXM) Algorithm = 0 or '/ &
             15x,'     Even-Odd Decomposition (EOD) Algorithm = 1 is valid.'/ &
             1x,60('*'))
 101 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_Matrix '/  &
             15x,'Unable to Extract the Differentiation Matrix from Operator.'/&
             15x,'ONLY Matrix-Matric Multiply (MXM) Algorithm = 0 or '/ &
             15x,'     Even-Odd Decomposition (EOD) Algorithm = 1 is valid.'/ &
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_Matrix '/  &
             15x,'Missing Array Argument(s)!                                 '/&
             15x,'Must specify either the full matrix D (S) or               '/&
             15x,'Both Even or Odd Matrices D_Even (S_Even) and D_Odd (S_Odd)'/&
             15x,'Respectively.                                              '/&
             1x,60('*'))

  END Subroutine PS_Get_Matrix 
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Put_Matrix (Operator, D, D_Even, D_Odd, Order, Error)

  integer , OPTIONAL :: Order
  logical , OPTIONAL :: Error

  REALTYPE, dimension(*)             :: Operator
  REALTYPE, dimension(:,:), OPTIONAL :: D, D_Even, D_Odd

  logical  :: L_D, L_D_Even, L_D_Odd

  integer  :: Algorithm, D_Or_S

  if (PRESENT(Error)) Error = .FALSE.

  L_D = PRESENT(D) ; L_D_Even = PRESENT(D_Even) ; L_D_Odd  = PRESENT(D_Odd)

#if 1
  if (.NOT. (L_D .OR. L_D_Even .OR. L_D_Odd)) then
#else
  if (.NOT. (L_D .OR. (L_D_Even .AND. L_D_Odd))) then
#endif
    write (lid6,102) ; write (lid99,102) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif

  D_Or_S    = INT(Operator(_I_CASE))
  Algorithm = INT(Operator(_I_ALGORITHM))

  SELECT CASE (D_Or_S)
    CASE (0)                  ! Smoothing
      SELECT CASE (Algorithm)
        CASE (0)
          call PS_Put_MXM_S (Operator, D, Error)

        CASE (1)
          call PS_Put_EOD_S (Operator, D_Even, D_Odd, Error)

        CASE DEFAULT
          write (lid6,100) ; write (lid99,100) 

          if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN 
                              else ; STOP ; endif

      END SELECT

    CASE (1)                  ! Differentiation
      SELECT CASE (Algorithm)
        CASE (0)
          call PS_Put_MXM_D (Operator, D, Order, Error)

        CASE (1)
          call PS_Put_EOD_D (Operator, D_Even, D_Odd, Order, Error)

        CASE DEFAULT
          write (lid6,101) ; write (lid99,101) 

          if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN 
                              else ; STOP ; endif

      END SELECT

  END SELECT

 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_Matrix '/  &
             15x,'Unable to Restore the Smoothing Matrix from Operator.'/&
             15x,'ONLY Matrix-Matric Multiply (MXM) Algorithm = 0 or '/ &
             15x,'     Even-Odd Decomposition (EOD) Algorithm = 1 is valid.'/ &
             1x,60('*'))
 101 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_Matrix '/  &
             15x,'Unable to Restore the Differentiation Matrix from Operator.'/&
             15x,'ONLY Matrix-Matric Multiply (MXM) Algorithm = 0 or '/ &
             15x,'     Even-Odd Decomposition (EOD) Algorithm = 1 is valid.'/ &
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_Matrix '/  &
             15x,'Missing Array Argument(s)!                                 '/&
             15x,'Must specify either the full matrix D (S) or               '/&
             15x,'Both Even or Odd Matrices D_Even (S_Even) and D_Odd (S_Odd)'/&
             15x,'Respectively.                                              '/&
             1x,60('*'))


  END Subroutine PS_Put_Matrix 
!
! -------------------------------------------------------------------------
!
  Subroutine PS_Get_MXM_D (Operator, D, Order, N, Error)

  integer , OPTIONAL :: N
  integer , OPTIONAL :: Order
  logical            :: L_D, L_Order

  logical , OPTIONAL :: Error

  integer            :: Algorithm, Max_Order, N_Pnts, O_Of_D, N1
  integer            :: i_Offset, i_d_s, D_Size
  integer            :: M1, M2

  REALTYPE, dimension(*)             :: Operator
! REALTYPE, dimension(:,:), OPTIONAL :: D
  REALTYPE, dimension(:,:)           :: D

  if (PRESENT(Error)) Error = .FALSE.

#if defined (DEBUG)
  Algorithm  = INT(Operator(_I_ALGORITHM))

  if (Algorithm /= 0) then
    write (lid6,100) ; write (lid99,100)

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif

  L_Order = PRESENT(Order)
! L_D     = PRESENT(D) 
  L_D     = .TRUE.

  if (.NOT. (L_D .OR. PRESENT(N))) then
    write (lid6,200) ; write (lid99,200) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif

  N_Pnts     = INT(Operator(_I_N_PNTS))
  Max_Order  = INT(Operator(_I_MAX_ORDER_OF_D))

  O_Of_D = 1 ; if (L_Order) O_Of_D = Order

#if defined (DEBUG)
  if (L_Order) then
    if ((Order <  0) .OR. (Order >  Max_Order)) then
      write (lid6,101) Order, Max_Order ; write (lid99,101) Order, Max_Order 

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif
#endif

  N1       = N_Pnts
  D_Size   = MXM_FULL_SIZE(N1)

  i_Offset = (O_Of_D-1)*MXM_FULL_SIZE(N1)
  i_d_s    = MXM_MATRIX_START(N1) + i_Offset

#if defined (DEBUG)
  M1 = SIZE(D, DIM=1) ; M2 = SIZE(D, DIM=2)

  if ((M1 <  N1) .OR. (M2 <  N1)) then
    write (lid6,102) N1, N1, M1, M2 ; write (lid99,102) N1, N1, M1, M2 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif
  
  if (PRESENT(N)) N = N1

  if (L_D) call Get_B_From_A (N1, Operator(i_d_s), D)

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_MXM_D '/  &
             15x,'Unable to Extract the Differentiation Matrix from Operator.'/&
             15x,'ONLY Matrix-Matric Multiply (MXM) Algorithm = 0 is valid.'/ &
             1x,60('*'))
 101 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_MXM_D '/  &
             15x,'Incorrect Specification of Order of Differentiation'/ &
             15x,'1 <= Order = ',i1,' <= ',i1/&
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_MXM_D '/  &
             15x,'Insufficient storage was allocated to matrix D'/ &
             15x,'Required  Size for D = (',i4,',',i4,')'/ &
             15x,'Actual    Size for D = (',i4,',',i4,')'/ &
             1x,60('*'))
#endif
 200 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_MXM_D '/  &
             15x,'Missing Argument : D and N were missing.' / &
             1x,60('*'))

  END Subroutine PS_Get_MXM_D 
!
! -------------------------------------------------------------------------
!
  Subroutine PS_Put_MXM_D (Operator, D, Order, Error)

  integer , OPTIONAL :: Order
  logical            :: L_Order

  logical , OPTIONAL :: Error

  integer            :: Algorithm, Max_Order, N_Pnts, O_Of_D, N1
  integer            :: i_Offset, i_d_s, D_Size
  integer            :: M1, M2

  REALTYPE, dimension(*)             :: Operator
! REALTYPE, dimension(:,:), OPTIONAL :: D
  REALTYPE, dimension(:,:)           :: D

  if (PRESENT(Error)) Error = .FALSE.

#if defined (DEBUG)
  Algorithm  = INT(Operator(_I_ALGORITHM))

  if (Algorithm /= 0) then
    write (lid6,100) ; write (lid99,100)

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif

  N_Pnts     = INT(Operator(_I_N_PNTS))
  Max_Order  = INT(Operator(_I_MAX_ORDER_OF_D))

  L_Order = PRESENT(Order)

  O_Of_D = 1 ; if (L_Order) O_Of_D = Order

#if defined (DEBUG)
  if (L_Order) then
    if ((Order <  0) .OR. (Order >  Max_Order)) then
      write (lid6,101) Order, Max_Order ; write (lid99,101) Order, Max_Order 

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif
#endif

  N1       = N_Pnts
  D_Size   = MXM_FULL_SIZE(N1)

  i_Offset = (O_Of_D-1)*MXM_FULL_SIZE(N1)
  i_d_s    = MXM_MATRIX_START(N1) + i_Offset

#if defined (DEBUG)
  M1 = SIZE(D, DIM=1) ; M2 = SIZE(D, DIM=2)

  if ((M1 >  N1) .OR. (M2 >  N1)) then
    write (lid6,102) N1, N1, M1, M2 ; write (lid99,102) N1, N1, M1, M2 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif
  
  call Put_B_Into_A (N1, Operator(i_d_s), D)

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_MXM_D '/  &
             15x,'Unable to Restore the Differentiation Matrix from Operator.'/&
             15x,'ONLY Matrix-Matric Multiply (MXM) Algorithm = 0 is valid.'/ &
             1x,60('*'))
 101 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_MXM_D '/  &
             15x,'Incorrect Specification of Order of Differentiation'/ &
             15x,'1 <= Order = ',i1,' <= ',i1/&
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_MXM_D '/  &
             15x,'Inconsistent storage was allocated to matrix D'/ &
             15x,'Allocated Size for D = (',i4,',',i4,')'/ &
             15x,'Actual    Size for D = (',i4,',',i4,')'/ &
             1x,60('*'))
#endif

  END Subroutine PS_Put_MXM_D 
!
! -------------------------------------------------------------------------
!
  Subroutine PS_Get_MXM_S (Operator, S, N, Error)

  integer , OPTIONAL :: N
  logical , OPTIONAL :: Error
  logical            :: L_S 

  integer            :: Algorithm, N_Pnts, N1
  integer            :: i_d_s, S_Size
  integer            :: M1, M2

  REALTYPE, dimension(*)             :: Operator
! REALTYPE, dimension(:,:), OPTIONAL :: S
  REALTYPE, dimension(:,:)           :: S

  if (PRESENT(Error)) Error = .FALSE.

#if defined (DEBUG)
  Algorithm  = INT(Operator(_I_ALGORITHM))

  if (Algorithm /= 0) then
    write (lid6,100) ; write (lid99,100)

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif

! L_S = PRESENT(S)
  L_S = .TRUE.

  if (.NOT. (L_S .OR. PRESENT(N))) then
    write (lid6,200) ; write (lid99,200) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif

  N_Pnts     = INT(Operator(_I_N_PNTS))

  N1       = N_Pnts
  S_Size   = MXM_FULL_SIZE(N1)

  i_d_s    = MXM_MATRIX_START(N1)

#if defined (DEBUG)
  M1 = SIZE(S, DIM=1) ; M2 = SIZE(S, DIM=2)

  if ((M1 <  N1) .OR. (M2 <  N1)) then
    write (lid6,102) N1, N1, M1, M2 ; write (lid99,102) N1, N1, M1, M2 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif
  
  if (PRESENT(N)) N = N1

  if (L_S) call Get_B_From_A (N1, Operator(i_d_s), S)

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_MXM_S '/  &
             15x,'Unable to Extract the Smoothing Matrix from Operator.'/ &
             15x,'ONLY Matrix-Matric Multiply (MXM) Algorithm = 0 is valid.'/ &
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_MXM_S '/  &
             15x,'Insufficient Storage was allocated to matrix S'/ &
             15x,'Required  Size for S = (',i4,',',i4,')'/ &
             15x,'Actual    Size for S = (',i4,',',i4,')'/ &
             1x,60('*'))
#endif
 200 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_MXM_S '/  &
             15x,'Missing Argument : D and N were missing.' / &
             1x,60('*'))

  END Subroutine PS_Get_MXM_S 
!
! -------------------------------------------------------------------------
!
  Subroutine PS_Put_MXM_S (Operator, S, Error)

  logical , OPTIONAL :: Error

  integer            :: Algorithm, N_Pnts, N1
  integer            :: i_d_s, S_Size
  integer            :: M1, M2

  REALTYPE, dimension(*)             :: Operator
! REALTYPE, dimension(:,:), OPTIONAL :: S
  REALTYPE, dimension(:,:)           :: S

  if (PRESENT(Error)) Error = .FALSE.

#if defined (DEBUG)
  Algorithm  = INT(Operator(_I_ALGORITHM))

  if (Algorithm /= 0) then
    write (lid6,100) ; write (lid99,100)

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif

  N_Pnts     = INT(Operator(_I_N_PNTS))

  N1       = N_Pnts
  S_Size   = MXM_FULL_SIZE(N1)

  i_d_s    = MXM_MATRIX_START(N1)

#if defined (DEBUG)
  M1 = SIZE(S, DIM=1) ; M2 = SIZE(S, DIM=2)

  if ((M1 >  N1) .OR. (M2 >  N1)) then
    write (lid6,102) N1, N1, M1, M2 ; write (lid99,102) N1, N1, M1, M2 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif
  
  call Put_B_Into_A (N1, Operator(i_d_s), S)

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_MXM_S '/  &
             15x,'Unable to Extract the Smoothing Matrix from Operator.'/ &
             15x,'ONLY Matrix-Matric Multiply (MXM) Algorithm = 0 is valid.'/ &
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_MXM_S '/  &
             15x,'Insufficient Storage was allocated to matrix S'/ &
             15x,'Allocated Size for S = (',i4,',',i4,')'/ &
             15x,'Actual    Size for S = (',i4,',',i4,')'/ &
             1x,60('*'))
#endif

  END Subroutine PS_Put_MXM_S 
!
! --------------------------------------------------------------------
!
  Subroutine PS_Get_EOD_D (Operator, D_Even, D_Odd, Order, N, Error)

  integer , OPTIONAL :: N

  integer , OPTIONAL :: Order
  logical            :: L_Order

  logical , OPTIONAL :: Error

  integer            :: Method, Algorithm, Max_Order, N_Pnts, O_Of_D, N1, N2
  integer            :: i_Offset, i_e_s, i_o_s, D_Size
  integer            :: M1, M2

  logical            :: L_D_Even, L_D_Odd

  REALTYPE, dimension(*)             :: Operator
  REALTYPE, dimension(:,:), OPTIONAL :: D_Even, D_Odd

  if (PRESENT(Error)) Error = .FALSE.

#if defined (DEBUG)
  Algorithm  = INT(Operator(_I_ALGORITHM))

  if (Algorithm /= 1) then
    write (lid6,100) ; write (lid99,100)

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif

  L_Order  = PRESENT(Order)
  L_D_Even = PRESENT(D_Even)
  L_D_Odd  = PRESENT(D_Odd )

  if (.NOT. (L_D_Even .OR. L_D_Odd .OR. PRESENT(N))) then
    write (lid6,200) ; write (lid99,200) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif

  Method     = INT(Operator(_I_METHOD))   
  N_Pnts     = INT(Operator(_I_N_PNTS))
  Max_Order  = INT(Operator(_I_MAX_ORDER_OF_D))

  O_Of_D = 1 ; if (L_Order) O_Of_D = Order

#if defined (DEBUG)
  if (L_Order) then
    if ((Order <  0) .OR. (Order >  Max_Order)) then
      write (lid6,101) Order, Max_Order ; write (lid99,101) Order, Max_Order

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif
#endif

  SELECT CASE (Method)
    CASE (0)
      N1 = N_Pnts+1 ; N2 = _H(N1)

    CASE (1,2)
      N1 = N_Pnts   ; N2 = _H(N1)

  END SELECT

  D_Size   = EOD_HALF_SIZE(N1)

  i_Offset = (O_Of_D-1)*EOD_HALF_SIZE(N1)
  i_e_s    = EOD_MATRIX_START_EVEN(N1) + i_Offset
  i_o_s    = EOD_MATRIX_START_ODD(N1)  + i_Offset

#if defined (DEBUG)
  M1 = SIZE(D_Even, DIM=1) ; M2 = SIZE(D_Even, DIM=2)

  if (L_D_Even) then
    if ((M1 <  N2) .OR. (M2 <  N2)) then
      write (lid6,102) N2, N2, M1, M2 ; write (lid99,102) N2, N2, M1, M2

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif

  M1 = SIZE(D_Odd , DIM=1) ; M2 = SIZE(D_Odd , DIM=2)

  if (L_D_Odd ) then
    if ((M1 <  N2) .OR. (M2 <  N2)) then
      write (lid6,103) N2, N2, M1, M2 ; write (lid99,102) N2, N2, M1, M2

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif
#endif
  
  if (PRESENT(N)) N = N2

  if (L_D_Even) call Get_B_From_A (N2, Operator(i_e_s), D_Even)
  if (L_D_Odd ) call Get_B_From_A (N2, Operator(i_o_s), D_Odd )

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_EOD_D '/  &
             15x,'Unable to Extract the Even And/Or the Odd '/ &
             15x,'  Differentiation Matrix from Operator.'/&
             15x,'ONLY Even-Odd Decomposition (EOD) Algorithm = 1 is valid'/&
             1x,60('*'))
 101 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_EOD_D '/  &
             15x,'Incorrect Specification of Order of Differentiation'/ &
             15x,'1 <= Order = ',i1,' <= ',i1/&
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_EOD_D '/  &
             15x,'Insufficient storage was allocated to matrix D_Even'/ &
             15x,'Required  Size for D_Even = (',i4,',',i4,')'/ &
             15x,'Actual    Size for D_Even = (',i4,',',i4,')'/ &
             1x,60('*'))
 103 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_EOD_D '/  &
             15x,'Insufficient storage was allocated to matrix D_Odd '/ &
             15x,'Required  Size for D_Odd  = (',i4,',',i4,')'/ &
             15x,'Actual    Size for D_Odd  = (',i4,',',i4,')'/ &
             1x,60('*'))
#endif
 200 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_EOD_D '/  &
             15x,'Missing Argument : D_Even, D_Odd and N were missing!' /&
             1x,60('*'))

  END Subroutine PS_Get_EOD_D
!
! --------------------------------------------------------------------
!
  Subroutine PS_Put_EOD_D (Operator, D_Even, D_Odd, Order, Error)

  integer , OPTIONAL :: Order
  logical            :: L_Order

  logical , OPTIONAL :: Error

  integer            :: Method, Algorithm, Max_Order, N_Pnts, O_Of_D, N1, N2
  integer            :: i_Offset, i_e_s, i_o_s, D_Size
  integer            :: M1, M2

  logical            :: L_D_Even, L_D_Odd

  REALTYPE, dimension(*)             :: Operator
  REALTYPE, dimension(:,:), OPTIONAL :: D_Even, D_Odd

  if (PRESENT(Error)) Error = .FALSE.

#if defined (DEBUG)
  Algorithm  = INT(Operator(_I_ALGORITHM))

  if (Algorithm /= 1) then 
    write (lid6,100) ; write (lid99,100) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif

  L_Order  = PRESENT(Order)
  L_D_Even = PRESENT(D_Even)
  L_D_Odd  = PRESENT(D_Odd )

  if (.NOT. (L_D_Even .OR. L_D_Odd)) then
    write (lid6,200) ; write (lid99,200) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif

  Method     = INT(Operator(_I_METHOD))   
  N_Pnts     = INT(Operator(_I_N_PNTS))
  Max_Order  = INT(Operator(_I_MAX_ORDER_OF_D))

  O_Of_D = 1 ; if (L_Order) O_Of_D = Order

#if defined (DEBUG)
  if (L_Order) then
    if ((Order <  0) .OR. (Order >  Max_Order)) then
      write (lid6,101) Order, Max_Order ; write (lid99,101) Order, Max_Order 

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif
#endif

  SELECT CASE (Method)
    CASE (0)
      N1 = N_Pnts+1 ; N2 = _H(N1)

    CASE (1,2)
      N1 = N_Pnts   ; N2 = _H(N1)

  END SELECT

  D_Size   = EOD_HALF_SIZE(N1)

  i_Offset = (O_Of_D-1)*EOD_HALF_SIZE(N1)
  i_e_s    = EOD_MATRIX_START_EVEN(N1) + i_Offset
  i_o_s    = EOD_MATRIX_START_ODD(N1)  + i_Offset

#if defined (DEBUG)
  M1 = SIZE(D_Even, DIM=1) ; M2 = SIZE(D_Even, DIM=2)

  if (L_D_Even) then
    if ((M1 >  N2) .OR. (M2 >  N2)) then
      write (lid6,102) N2, N2, M1, M2 ; write (lid99,102) N2, N2, M1, M2

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif

  M1 = SIZE(D_Odd , DIM=1) ; M2 = SIZE(D_Odd , DIM=2)

  if (L_D_Odd ) then
    if ((M1 >  N2) .OR. (M2 >  N2)) then
      write (lid6,103) N2, N2, M1, M2 ; write (lid99,102) N2, N2, M1, M2

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif
#endif
  
  if (L_D_Even) call Put_B_Into_A (N2, Operator(i_e_s), D_Even)
  if (L_D_Odd ) call Put_B_Into_A (N2, Operator(i_o_s), D_Odd )

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_EOD_D '/  &
             15x,'Unable to Restore the Even And/Or the Odd '/ &
             15x,'  Differentiation Matrix from Operator.'/&
             15x,'ONLY Even-Odd Decomposition (EOD) Algorithm = 1 is valid'/&
             1x,60('*'))
 101 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_EOD_D '/  &
             15x,'Incorrect Specification of Order of Differentiation'/ &
             15x,'1 <= Order = ',i1,' <= ',i1/&
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_EOD_D '/  &
             15x,'Insufficient storage was allocated to matrix D_Even'/ &
             15x,'Allocated Size for D_Even = (',i4,',',i4,')'/ &
             15x,'Actual    Size for D_Even = (',i4,',',i4,')'/ &
             1x,60('*'))
 103 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_EOD_D '/  &
             15x,'Insufficient storage was allocated to matrix D_Odd '/ &
             15x,'Allocated Size for D_Odd  = (',i4,',',i4,')'/ &
             15x,'Actual    Size for D_Odd  = (',i4,',',i4,')'/ &
             1x,60('*'))
#endif
 200 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_EOD_D '/  &
             15x,'Missing Argument : Both D_Even and D_Odd were missing!' /&
             1x,60('*'))

  END Subroutine PS_Put_EOD_D 
!
! ----------------------------------------------------------------------
!
  Subroutine PS_Get_EOD_S (Operator, S_Even, S_Odd, N, Error)

  integer , OPTIONAL :: N

  logical , OPTIONAL :: Error

  integer            :: Method, Algorithm, N_Pnts, N1, N2
  integer            :: i_e_s, i_o_s, S_Size
  integer            :: M1, M2

  logical            :: L_S_Even, L_S_Odd

  REALTYPE, dimension(*)             :: Operator
  REALTYPE, dimension(:,:), OPTIONAL :: S_Even, S_Odd

  if (PRESENT(Error)) Error = .FALSE.

#if defined (DEBUG)
  Algorithm  = INT(Operator(_I_ALGORITHM))

  if (Algorithm /= 1) then
    write (lid6,100) ; write (lid99,100)

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif

  L_S_Even = PRESENT(S_Even)
  L_S_Odd  = PRESENT(S_Odd )

  if (.NOT. (L_S_Even .OR. L_S_Odd .OR. PRESENT(N))) then
    write (lid6,200) ; write (lid99,200) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif

  Method     = INT(Operator(_I_METHOD))   
  N_Pnts     = INT(Operator(_I_N_PNTS))

  SELECT CASE (Method)
    CASE (0)
      N1 = N_Pnts+1 ; N2 = _H(N1)

    CASE (1,2)
      N1 = N_Pnts   ; N2 = _H(N1)

  END SELECT

  S_Size   = EOD_HALF_SIZE(N1)

  i_e_s    = EOD_MATRIX_START_EVEN(N1)
  i_o_s    = EOD_MATRIX_START_EVEN(N1) + EOD_HALF_SIZE(N1)

#if defined (DEBUG)
  M1 = SIZE(S_Even, DIM=1) ; M2 = SIZE(S_Even, DIM=2)

  if (L_S_Even) then
    if ((M1 <  N2) .OR. (M2 <  N2)) then
      write (lid6,102) N2, N2, M1, M2 ; write (lid99,102) N2, N2, M1, M2

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif

  M1 = SIZE(S_Odd , DIM=1) ; M2 = SIZE(s_Odd , DIM=2)

  if (L_S_Odd ) then
    if ((M1 <  N2) .OR. (M2 <  N2)) then
      write (lid6,103) N2, N2, M1, M2 ; write (lid99,102) N2, N2, M1, M2

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif
#endif
  
  if (PRESENT(N)) N = N2

  if (L_S_Even) call Get_B_From_A (N2, Operator(i_e_s), S_Even)
  if (L_S_Odd ) call Get_B_From_A (N2, Operator(i_o_s), S_Odd )

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_EOD_S '/  &
             15x,'Unable to Extract the Even And/Or the Odd '/ &
             15x,'  Differentiation Matrix from Operator.'/&
             15x,'ONLY Even-Odd Decomposition (EOD) Algorithm = 1 is valid'/&
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_EOD_S '/  &
             15x,'Insufficient storage was allocated to matrix S_Even'/ &
             15x,'Required  Size for S_Even = (',i4,',',i4,')'/ &
             15x,'Actual    Size for S_Even = (',i4,',',i4,')'/ &
             1x,60('*'))
 103 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_EOD_S '/  &
             15x,'Insufficient storage was allocated to matrix S_Odd '/ &
             15x,'Required  Size for S_Odd  = (',i4,',',i4,')'/ &
             15x,'Actual    Size for S_Odd  = (',i4,',',i4,')'/ &
             1x,60('*'))
#endif
 200 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Get_EOD_S '/  &
             15x,'Missing Argument : D_Even, D_Odd and N were missing!' /&
             1x,60('*'))

  END Subroutine PS_Get_EOD_S 
!
! --------------------------------------------------------------------
!
  Subroutine PS_Put_EOD_S (Operator, S_Even, S_Odd, Error)

  logical , OPTIONAL :: Error

  integer            :: Method, Algorithm, N_Pnts, N1, N2
  integer            :: i_e_s, i_o_s, S_Size
  integer            :: M1, M2

  logical            :: L_S_Even, L_S_Odd

  REALTYPE, dimension(*)             :: Operator
  REALTYPE, dimension(:,:), OPTIONAL :: S_Even, S_Odd

  if (PRESENT(Error)) Error = .FALSE.

#if defined (DEBUG)
  Algorithm  = INT(Operator(_I_ALGORITHM))

  if (Algorithm /= 1) then 
    write (lid6,100) ; write (lid99,100) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif
#endif

  L_S_Even = PRESENT(S_Even)
  L_S_Odd  = PRESENT(S_Odd )

  if (.NOT. (L_S_Even .OR. L_S_Odd)) then
    write (lid6,200) ; write (lid99,200) 

    if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
  endif

  Method     = INT(Operator(_I_METHOD))   
  N_Pnts     = INT(Operator(_I_N_PNTS))

  SELECT CASE (Method)
    CASE (0)
      N1 = N_Pnts+1 ; N2 = _H(N1)

    CASE (1,2)
      N1 = N_Pnts   ; N2 = _H(N1)

  END SELECT

  S_Size   = EOD_HALF_SIZE(N1)

  i_e_s    = EOD_MATRIX_START_EVEN(N1) 
  i_o_s    = EOD_MATRIX_START_EVEN(N1) + EOD_HALF_SIZE(N1)

#if defined (DEBUG)
  M1 = SIZE(S_Even, DIM=1) ; M2 = SIZE(S_Even, DIM=2)

  if (L_S_Even) then
    if ((M1 >  N2) .OR. (M2 >  N2)) then
      write (lid6,102) N2, N2, M1, M2 ; write (lid99,102) N2, N2, M1, M2

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif

  M1 = SIZE(S_Odd , DIM=1) ; M2 = SIZE(s_Odd , DIM=2)

  if (L_S_Odd ) then
    if ((M1 >  N2) .OR. (M2 >  N2)) then
      write (lid6,103) N2, N2, M1, M2 ; write (lid99,102) N2, N2, M1, M2

      if (PRESENT(Error)) then ; Error = .TRUE. ; RETURN ; else ; STOP ; endif
    endif
  endif
#endif
  
  if (L_S_Even) call Put_B_Into_A (N2, Operator(i_e_s), S_Even)
  if (L_S_Odd ) call Put_B_Into_A (N2, Operator(i_o_s), S_Odd )

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_EOD_S '/  &
             15x,'Unable to Restore the Even And/Or the Odd '/ &
             15x,'  Differentiation Matrix from Operator.'/&
             15x,'ONLY Even-Odd Decomposition (EOD) Algorithm = 1 is valid'/&
             1x,60('*'))
 102 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_EOD_S '/  &
             15x,'Insufficient storage was allocated to matrix S_Even'/ &
             15x,'Allocated Size for S_Even = (',i4,',',i4,')'/ &
             15x,'Actual    Size for S_Even = (',i4,',',i4,')'/ &
             1x,60('*'))
 103 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_EOD_S '/  &
             15x,'Insufficient storage was allocated to matrix S_Odd '/ &
             15x,'Allocated Size for S_Odd  = (',i4,',',i4,')'/ &
             15x,'Actual    Size for S_Odd  = (',i4,',',i4,')'/ &
             1x,60('*'))
#endif
 200 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Put_EOD_S '/  &
             15x,'Missing Argument : Both D_Even and D_Odd were missing!' /&
             1x,60('*'))

  END Subroutine PS_Put_EOD_S 
!
! ------------------------------------------------------------
!
  Subroutine Get_B_From_A (N, A, B)

  integer                  :: N
  REALTYPE, dimension(N,N) :: A
  REALTYPE, dimension(:,:) :: B

  B(1:N,1:N) = A

  END Subroutine Get_B_From_A 
!
! ------------------------------------------------------------
!
  Subroutine Put_B_Into_A (N, A, B)

  integer                  :: N
  REALTYPE, dimension(N,N) :: A
  REALTYPE, dimension(:,:) :: B

  A(1:N,1:N) = B

  END Subroutine Put_B_Into_A


! ----------------------------------------------------------------------
! FUNCTION NAME: PS_Change_Filter 
! DESCRIPTION  : Replace the built-in filter inside Operator
!                with a user specified filter function Sigma.
!                Fourier and Chebyshev Transform Algorithm ONLY.
!
!  Smooth_1   = Turn on Filtering before differentiation for Chebyshev Method
!               Turn on Filtering for Fourier Differentiation
!               Turn on Filtering for Smoothing
!  Smooth_2   = Turn on Filtering aftter differentiation for Chebyshev Method
!               No used otherwise.
!  Operator   = Arrays holding filtering parameters 
!  Sigma      = 1D array containing new filter coefficients with length
!                 appropriated for the particular Method (Fourier or Chebyshev)
! ----------------------------------------------------------------------

  Subroutine PS_Change_Filter (Operator, Sigma)

  integer  :: N_Pnts, N, N1, Mode_N, i_s_s, i
  integer  :: Method, Algorithm, Symmetry

  REALTYPE, dimension(0:) :: Sigma
  REALTYPE, dimension(*)  :: Operator

  Algorithm  = INT(Operator(_I_ALGORITHM))
  Method     = INT(Operator(_I_METHOD))
  Symmetry   = INT(Operator(_I_SYMMETRY))
  N_Pnts     = INT(Operator(_I_N_PNTS))

#if defined (DEBUG)
  if (Algorithm /= 2) then ; write (lid6,100) ; write (lid99,100) ; STOP ; endif
#endif

  SELECT CASE (Method)
    CASE (0)
      N      = N_Pnts
      Mode_N = N/2

      if (Symmetry == 0) then
        N1 = N
        i_s_s = FFT_SIGMA_START_FULL(N1)
      else
        N1 = 2*(N/2)+1
        i_s_s = FFT_SIGMA_START_HALF(N1)
      endif

    CASE (1)
      N      = N_Pnts-1
      Mode_N = N

      if (Symmetry == 0) then
        N1 = N+1
        i_s_s = CFT_SIGMA_START_FULL(N1)
      else
        N1 = N+1
        i_s_s = CFT_SIGMA_START_HALF(N1)
      endif

    CASE DEFAULT
      write (lid6,101) ; write (lid99,101) ; STOP

  END SELECT

  do i = 0,Mode_N
    Operator(i_s_s+i) = Sigma(i)
  enddo

#if defined (DEBUG)
 100 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Change_Filter '/  &
             15x,'Unable to Modify the Filter in the Operator.'/ &
             15x,'  ONLY Transform Algorithm = 2 can be modified.'/ &
             1x,60('*'))
#endif
 101 format (1x,60('*')/  &
             1x,'FATAL ERROR : PS_Change_Filter '/  &
             15x,'Unable to Modify the Filter in the Operator.'/ &
             15x,'  ONLY Filter for the Fourier Method (0) or '/&
             15x,'  the Chebyshev Method (1) can be modified.'/ &
             1x,60('*'))

  END Subroutine PS_Change_Filter 

! ----------------------------------------------------------------------
! FUNCTION NAME: PS_Filter_Switch 
! DESCRIPTION  : Routines to turn on/off the filtering for Operator
! ----------------------------------------------------------------------
  Subroutine PS_Filter_Switch (Operator, Smooth_1, Smooth_2)

  integer , OPTIONAL :: Smooth_1, Smooth_2

  REALTYPE, dimension(*) :: Operator

  if (PRESENT(Smooth_1)) Operator(_I_SMOOTH_1) = Smooth_1
  if (PRESENT(Smooth_2)) Operator(_I_SMOOTH_2) = Smooth_2

  END Subroutine PS_Filter_Switch 
!
! ==================================================================
! Should pre-compute the cos(i*k*Pi/N) for multiple j
!
  Subroutine PS_Chebyshev_Lobatto_a_k (f, a)

  integer  :: N, i, k
  REALTYPE :: PiN

  REALTYPE, dimension(0:) :: f, a

  N = SIZE(f)-1

  PiN  =  PI/N

  do k = 0,N

    a(k) = HALF*(f(0)+f(N)*(-1)**k)
    do i = 1,N-1
      a(k) = a(k) + f(i)*COS(i*k*PiN)
    enddo
    a(k) = a(k)*TWO/N

  enddo

  a(0) = HALF*a(0)
  a(N) = HALF*a(N)

  END Subroutine PS_Chebyshev_Lobatto_a_k
!
! ==================================================================
!
  Subroutine PS_Chebyshev_Lobatto_T_k (x, T)

  integer  :: k, N

  REALTYPE, dimension(:)    :: x
  REALTYPE, dimension(:,0:) :: T

  N = SIZE(T,DIM=2)-1

  do k = 0,N
    T(:,k) = COS(k*ACOS(x))
  enddo

  END Subroutine PS_Chebyshev_Lobatto_T_k
!
! ==================================================================
!
  Subroutine PS_Interpolation_Sum_a_T (a, T, Q)

  REALTYPE, dimension(:,:) :: a, T, Q

#if 0
  Q = MATMUL(T,a)
#else
 call GEMM ('N', 'N', SIZE(T,DIM=1), SIZE(a,DIM=2), SIZE(T,DIM=2), &
            ONE, T, SIZE(T,DIM=1), a, SIZE(a,DIM=1), ZERO, Q, SIZE(Q,DIM=1))
#endif

  END Subroutine PS_Interpolation_Sum_a_T
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Index_Search (z, x, m0, m1, N, i0, i1, Wind_Velocity, Error)

  integer                    :: m0, m1, N, i0, i1
  REALTYPE                   :: z
  REALTYPE, dimension(m0:m1) :: x
  REALTYPE, dimension(m0:m1) :: y
  REALTYPE                   :: Tolerance = TEN*EPSILON(ONE)

  integer , parameter        :: lid6=6, lid99=99
  integer                    :: i, i2, i3, i4, i5, M, Grid_Direction
  logical                    :: z_is_y_i

  integer                , OPTIONAL :: Wind_Velocity
  logical                           :: L_Speed
  integer                           :: W_Speed, W_Direction

  logical                , OPTIONAL :: Error
  logical                           :: L_Error

  L_Speed = PRESENT(Wind_Velocity)
  L_Error = PRESENT(Error)

  M = ABS(m1-m0+1)

                      Grid_Direction =  1
  if (x(m0) >  x(m1)) Grid_Direction = -1

#if defined (DEBUG)
  call Data_Integrity
#endif

  if (L_Error .AND. Error) RETURN

                            y(m0:m1: 1) = x(m0:m1)
  if (Grid_Direction == -1) y(m1:m0:-1) = x(m0:m1)

  i2 = m0 ; i3 = m1 ; i4 =  1 

  do i = i2, i3, i4
    if (y(i) >= z) then ; i0 = i ; EXIT ; endif 
  enddo

                                 z_is_y_i = .FALSE.
  if (ABS(z-y(i0)) <= Tolerance) z_is_y_i = .TRUE.
 
! write (6,*) z, y(i0), z_is_y_i

    W_Speed = 0 ; W_Direction = 0
  if (L_Speed .AND. (Wind_Velocity /= 0)) then
    W_Direction = 1 ; if (Wind_Velocity <  ZERO) W_Direction = -1

    W_Speed = W_Direction*MIN(ABS(Wind_Velocity), N/2)
  endif

#if 0
  i2 = N/2 ; if (MOD(N,2) /= 0) i2 = i2+1

  i3 = i0-i2 ; i4 = i0+i2 ; i5 = i0
#else
  i2 = N/2 ; i3 = i0-i2 ; i4 = i0+i2 ; i5 = i0
#endif

  i0 = i3 + W_Speed
  i1 = i4 + W_Speed

  if (.NOT. z_is_y_i) then
    if (W_Direction ==  1) then   !  Upwinding 
      if (i0 >= i5) then ; i0 = i0-1 ; i1 = i1-1 ; endif
    endif

    if (W_Direction == -1) then   !  Downwinding
      if (i1 <= i5) then ; i0 = i0+1 ; i1 = i1+1 ; endif
    endif
  endif

  if (i0 < m0) then ; i0 = m0   ; i1 = m0+N ; endif
  if (i1 > m1) then ; i0 = m1-N ; i1 = m1   ; endif

  if (Grid_Direction == -1) then
    i2 = M-(N+1) ; i3 = i0 ; i0 = i2-i1 ; i1 = i2-i3 
  endif

#if defined (DEBUG)
  CONTAINS

    Subroutine Data_Integrity 

    logical :: Err

                 Err   = .FALSE.
    if (L_Error) Error = .FALSE.

    if (N >  M) then
      write (lid6,103) N, M ; write (lid99,103) N, M 
      if (L_Error) then ; Error = .TRUE. ; else ; Err = .TRUE. ; endif
    endif

    if (x(m0) == x(m1)) then
      write (lid6 ,102) x(m0) ; write (lid99,102) x(m0) 
      if (L_Error) then ; Error = .TRUE. ; else ; Err = .TRUE. ; endif
    endif

    SELECT CASE (Grid_Direction)
      CASE (1)
        if (z <  x(m0)) then
          write (lid6 ,100)  '<', z, x(m0) ; write (lid99,100)  '<', z, x(m0) 
          if (L_Error) then ; Error = .TRUE. ; else ; Err = .TRUE. ; endif
        endif
        if (z >  x(m1)) then
          write (lid6 ,101)  '>', z, x(m1) ; write (lid99,101)  '>', z, x(m1) 
          if (L_Error) then ; Error = .TRUE. ; else ; Err = .TRUE. ; endif
        endif

      CASE (-1)
        if (z >  x(m0)) then
          write (lid6 ,100)  '>', z, x(m0) ; write (lid99,100)  '>', z, x(m0) 
          if (L_Error) then ; Error = .TRUE. ; else ; Err = .TRUE. ; endif
        endif
        if (z <  x(m1)) then
          write (lid6 ,101)  '<', z, x(m1) ; write (lid99,101)  '<', z, x(m1) 
          if (L_Error) then ; Error = .TRUE. ; else ; Err = .TRUE. ; endif
        endif

    END SELECT

    if (Err) STOP

  100 format (1x,'FATAL ERROR : PS_Index_Search : z ',a,' x(m0)'/          &
              1x,'              z  = ',g14.7/                          &
              1x,'           x(m0) = ',g14.7)
  101 format (1x,'FATAL ERROR : PS_Index_Search : z ',a,' x(m1)'/          &
              1x,'              z  = ',g14.7/                          &
              1x,'           x(m1) = ',g14.7)
  102 format (1x,'FATAL ERROR : PS_Index_Search : x(m0) = x(m1)'/           &
              1x,'     x(m1)=x(m0) = ',g14.7)
  103 format (1x,'FATAL ERROR : PS_Index_Search : N > M'/                 &
              1x,'              N  = ',i5   /                          &
              1x,'              M  = ',i5   )

    END Subroutine Data_Integrity 
#endif

  END Subroutine PS_Index_Search 

END MODULE PseudoPack_Utility
