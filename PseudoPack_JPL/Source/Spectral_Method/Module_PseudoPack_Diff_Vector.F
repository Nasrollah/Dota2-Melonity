#include "pseudopack.h"
  
  
  
MODULE Differentiation_Vector
  
  USE PseudoPack_Setup     , ONLY: PS_Property, PS_Grid_Index, PS_Domain, &
                                   PS_Mapping
  USE PseudoPack_Allocation, ONLY: PS_Get_Operator_Size
  USE PseudoPack_Setup     , ONLY: PS_Setup_Operator
  USE PseudoPack_Utility   , ONLY: PS_Get_Matrix
  
implicit NONE

INTERFACE PS_BC_Differentiation_Vector
  MODULE PROCEDURE PS_BC_Differentiation_Vector_0
  MODULE PROCEDURE PS_BC_Differentiation_Vector_1
END INTERFACE
  
INTERFACE PS_Get_Differentiation_Vector
  MODULE PROCEDURE PS_Get_Differentiation_Vector_0
  MODULE PROCEDURE PS_Get_Differentiation_Vector_1
  MODULE PROCEDURE PS_Get_Differentiation_Vector_2
  MODULE PROCEDURE PS_Get_Differentiation_Vector_3
  MODULE PROCEDURE PS_Get_Differentiation_Vector_4
END INTERFACE
  
PRIVATE

PUBLIC  :: PS_BC_Differentiation_Vector 
PUBLIC  :: PS_Get_Differentiation_Vector 
  
CONTAINS
  
  Subroutine PS_BC_Differentiation_Vector_0 (N_Pnts, Order,                &
                                             Property, Grid_Index, Domain, &
                                             D_0, D_N)
  
  integer                               :: N_Pnts, Order, N, M_D
  REALTYPE, dimension(0:N_Pnts-1)       :: D_0, D_N
  REALTYPE, dimension(:)  , ALLOCATABLE :: D_Operator, x
  REALTYPE, dimension(:,:), ALLOCATABLE :: D_Matrix
  
  TYPE (PS_Property)   :: Property, D_Property
  TYPE (PS_Grid_Index) :: Grid_Index
  TYPE (PS_Domain)     :: Domain
  
  N = N_Pnts-1
  
  D_Property = Property 
  
  D_Property%Max_Order_Of_D = Order
  D_Property%Algorithm      = 0 
  D_Property%Symmetry       = 0
  
  call PS_Get_Operator_Size ('D', N_Pnts, M_D, D_Property)
  
    Allocate (D_Operator(M_D), x(0:N))
  
  call PS_Setup_Operator    (N_Pnts, D_Operator, x, &
                                     D_Property, Grid_Index, Domain)
  
! Allocate Storage for the Differentiation Matrix D_Matrix
  Allocate (D_Matrix(0:N,0:N))
  
! Get D_Matrix from the differentiation operator D
  call PS_Get_Matrix (D_Operator, D_Matrix, Order=Order)
  
! Get First and Last row of D_Matrix and store in D_0 and D_N
  D_0 = D_Matrix(0,:)
  D_N = D_Matrix(N,:)
  
! Deallocate storage of x, D and D_Matrix
  DeAllocate (x, D_Operator, D_Matrix)
  
  END Subroutine PS_BC_Differentiation_Vector_0
!
!------------------------------------------------------------------------------
!
  Subroutine PS_BC_Differentiation_Vector_1 (N_Pnts, Order,                &
                                             Property, Grid_Index, Domain, &
                                             Mapping,                      &
                                             D_0, D_N)
  
  integer                               :: N_Pnts, Order, N, M_D
  REALTYPE, dimension(0:N_Pnts-1)       :: D_0, D_N
  REALTYPE, dimension(:)  , ALLOCATABLE :: D_Operator, x
  REALTYPE, dimension(:,:), ALLOCATABLE :: D_Matrix
  
  TYPE (PS_Property)   :: Property  , D_Property
  TYPE (PS_Grid_Index) :: Grid_Index
  TYPE (PS_Domain)     :: Domain
  TYPE (PS_Mapping)    :: Mapping
  
  N = N_Pnts-1
  
  D_Property = Property
  
  D_Property%Max_Order_Of_D = Order
  D_Property%Algorithm      = 0
  D_Property%Symmetry       = 0
  
  call PS_Get_Operator_Size ('D', N_Pnts, M_D, D_Property)
  
    Allocate (D_Operator(M_D), x(0:N))
  
  call PS_Setup_Operator    (N_Pnts, D_Operator, x, &
                                     D_Property, Grid_Index, Domain, Mapping)
  
! Allocate Storage for the Differentiation Matrix D_Matrix
  Allocate (D_Matrix(0:N,0:N))
  
! Get D_Matrix from the differentiation operator D
  call PS_Get_Matrix (D_Operator, D_Matrix, Order=Order)
  
! Get First and Last row of D_Matrix and store in D_0 and D_N
  D_0 = D_Matrix(0,:)
  D_N = D_Matrix(N,:)
  
! Deallocate storage of x, D and D_Matrix
  DeAllocate (x, D_Operator, D_Matrix)
  
  END Subroutine PS_BC_Differentiation_Vector_1
!
!------------------------------------------------------------------------------
!
  Subroutine PS_Get_Differentiation_Vector_0 (N_Pnts, Order,                &
                                              Property, Grid_Index, Domain, &
                                              Vector_Index, D_Vector)
  
  integer                               :: N_Pnts, Order, Vector_Index, N, M_D
  REALTYPE, dimension(0:N_Pnts-1)       :: D_Vector
  REALTYPE, dimension(:)  , ALLOCATABLE :: D_Operator, x
  REALTYPE, dimension(:,:), ALLOCATABLE :: D_Matrix
  
  TYPE (PS_Property)   :: Property  , D_Property
  TYPE (PS_Grid_Index) :: Grid_Index
  TYPE (PS_Domain)     :: Domain
  
  N = N_Pnts-1
  
#if defined (DEBUG)
  call Data_Integrity_0 (N, Vector_Index)
#endif
  
  D_Property = Property 
  
  D_Property%Max_Order_Of_D = Order
  D_Property%Algorithm      = 0 
  D_Property%Symmetry       = 0
  
  call PS_Get_Operator_Size ('D', N_Pnts, M_D, D_Property)
  
    Allocate (D_Operator(M_D), x(0:N))
  
  call PS_Setup_Operator    (N_Pnts, D_Operator, x, &
                                     D_Property, Grid_Index, Domain)
  
! Allocate Storage for the Differentiation Matrix D_Matrix
  Allocate (D_Matrix(0:N,0:N))
  
! Get D_Matrix from the differentiation operator D
  call PS_Get_Matrix (D_Operator, D_Matrix, Order=Order)
  
! Get the row of vector at vector_Index of D_Matrix and store in D_vector
  D_Vector = D_Matrix(Vector_Index,:)
  
! Deallocate storage of x, D and D_Matrix
  DeAllocate (x, D_Operator, D_Matrix)
  
  END Subroutine PS_Get_Differentiation_Vector_0
!
!------------------------------------------------------------------------------
!
  Subroutine PS_Get_Differentiation_Vector_1 (N_Pnts, Order,                &
                                              Property, Grid_Index, Domain, &
                                              Mapping,                      &
                                              Vector_Index, D_Vector)
  
  integer                               :: N_Pnts, Order, Vector_Index, N, M_D
  REALTYPE, dimension(0:N_Pnts-1)       :: D_Vector
  REALTYPE, dimension(:)  , ALLOCATABLE :: D_Operator, x
  REALTYPE, dimension(:,:), ALLOCATABLE :: D_Matrix
  
  TYPE (PS_Property)   :: Property  , D_Property
  TYPE (PS_Grid_Index) :: Grid_Index
  TYPE (PS_Domain)     :: Domain
  TYPE (PS_Mapping)    :: Mapping   
  
  N = N_Pnts-1
  
#if defined (DEBUG)
  call Data_Integrity_0 (N, Vector_Index)
#endif
  
  D_Property = Property
  
  D_Property%Max_Order_Of_D = Order
  D_Property%Algorithm      = 0
  D_Property%Symmetry       = 0
  
  call PS_Get_Operator_Size ('D', N_Pnts, M_D, D_Property)
  
    Allocate (D_Operator(M_D), x(0:N))
  
  call PS_Setup_Operator    (N_Pnts, D_Operator, x, &
                                     D_Property, Grid_Index, Domain, Mapping)
  
! Allocate Storage for the Differentiation Matrix D_Matrix
  Allocate (D_Matrix(0:N,0:N))
  
! Get D_Matrix from the differentiation operator D
  call PS_Get_Matrix (D_Operator, D_Matrix, Order=Order)
  
! Get the row of vector at vector_Index of D_Matrix and store in D_vector
  D_Vector = D_Matrix(Vector_Index,:)
  
! Deallocate storage of x, D and D_Matrix
  DeAllocate (x, D_Operator, D_Matrix)
  
  END Subroutine PS_Get_Differentiation_Vector_1
!
!------------------------------------------------------------------------------
!
  Subroutine PS_Get_Differentiation_Vector_2 (N_Pnts, Order,                 &
                                              Property, Grid_Index, Domain,  &
                                              Mapping,                       &
                                              N_Index, Vector_Index,         &
                                              D_Vector)
  
  integer                                     :: N_Pnts, Order, N_Index, N, M_D, i
  integer , dimension(N_Index)                :: Vector_Index
  REALTYPE, dimension(0:N_Pnts-1,N_Index)     :: D_Vector
  REALTYPE, dimension(:)  , ALLOCATABLE       :: D_Operator, x
  REALTYPE, dimension(:,:), ALLOCATABLE       :: D_Matrix
  
  TYPE (PS_Property)   :: Property  , D_Property
  TYPE (PS_Grid_Index) :: Grid_Index, D_Grid_Index
  TYPE (PS_Domain)     :: Domain    , D_Domain
  TYPE (PS_Mapping)    :: Mapping   , D_Mapping
  
  N = N_Pnts-1
  
#if defined (DEBUG)
  call Data_Integrity_1 (N, N_Index, Vector_Index)
#endif
  
  D_Property   = Property 
  D_Grid_Index = Grid_Index
  D_Domain     = Domain
  D_Mapping    = Mapping
  
  D_Property%Max_Order_Of_D = Order
  D_Property%Algorithm      = 0 
  D_Property%Symmetry       = 0
  
  call PS_Get_Operator_Size ('D', N_Pnts, M_D, D_Property)
  
    Allocate (D_Operator(M_D), x(0:N))
  
  call PS_Setup_Operator    (N_Pnts, D_Operator, x,                      &
                                     D_Property, D_Grid_Index, D_Domain, &
                                     D_Mapping)
  
  Allocate (D_Matrix(0:N,0:N))
  
  call PS_Get_Matrix (D_Operator, D=D_Matrix, Order=Order)
  
  do i = 1,N_Index
    D_Vector(:,i) = D_Matrix(Vector_Index(i),:)
  enddo
  
  DeAllocate (x, D_Operator, D_Matrix)
  
  END Subroutine PS_Get_Differentiation_Vector_2
!
!------------------------------------------------------------------------------
!
  Subroutine PS_Get_Differentiation_Vector_3 (N_Pnts, Order,                 &
                                              Property, Grid_Index, Domain,  &
                                              Mapping,                       &
                                              N_Index, Vector_Index,         &
                                              D_Vector_Even, D_Vector_Odd)
  
  integer                                     :: N_Pnts, N, N2 
  integer                                     :: Order, N_Index, M_D, i
  integer , dimension(N_Index)                :: Vector_Index
  REALTYPE, dimension(0:(N_Pnts-1)/2,N_Index) :: D_Vector_Even, D_Vector_Odd
  REALTYPE, dimension(:)  , ALLOCATABLE       :: D_Operator, x
  REALTYPE, dimension(:,:), ALLOCATABLE       :: D_Matrix_Even, D_Matrix_Odd
  
  TYPE (PS_Property)   :: Property  , D_Property
  TYPE (PS_Grid_Index) :: Grid_Index, D_Grid_Index
  TYPE (PS_Domain)     :: Domain    , D_Domain
  TYPE (PS_Mapping)    :: Mapping   , D_Mapping
  
  N = N_Pnts-1 ; N2 = N/2 
  
#if defined (DEBUG)
  call Data_Integrity_1 (N, N_Index, Vector_Index)
#endif
  
  D_Property   = Property 
  D_Grid_Index = Grid_Index
  D_Domain     = Domain
  D_Mapping    = Mapping
  
  D_Property%Max_Order_Of_D = Order
  D_Property%Algorithm      = 1 
  D_Property%Symmetry       = 1
  
  call PS_Get_Operator_Size ('D', N_Pnts, M_D, D_Property)
  
    Allocate (D_Operator(M_D), x(0:N))
  
  call PS_Setup_Operator    (N_Pnts, D_Operator, x,                      &
                                     D_Property, D_Grid_Index, D_Domain, &
                                     D_Mapping)
  
  Allocate (D_Matrix_Even(0:N2,0:N2), D_Matrix_Odd(0:N2,0:N2))
  
  call PS_Get_Matrix (D_Operator, D_Even=D_Matrix_Even, D_Odd=D_Matrix_Odd, &
                                  Order=Order)
  
  do i = 1,N_Index
    D_Vector_Even(:,i) = D_Matrix_Even(Vector_Index(i),:)
    D_Vector_Odd (:,i) = D_Matrix_Odd (Vector_Index(i),:)
  enddo
  
  DeAllocate (x, D_Operator, D_Matrix_Even, D_Matrix_Odd)
  
  END Subroutine PS_Get_Differentiation_Vector_3
!
!------------------------------------------------------------------------------
!
  Subroutine PS_Get_Differentiation_Vector_4 (N_Pnts, Order,                 &
                                              Property, Grid_Index, Domain,  &
                                              Mapping,                       &
                                              N_Index, Vector_Index, Case,   &
                                              D_Vector, D_Vector_Even,       &
                                                        D_Vector_Odd)
  
  integer                                     :: N_Pnts, N, N2, Case 
  integer                                     :: Order, N_Index, M_D, i
  logical                                     :: L_D, L_D_Even, L_D_Odd
  
  integer , dimension(N_Index)          :: Vector_Index
  REALTYPE, dimension(:)  , ALLOCATABLE :: D_Operator, x
  REALTYPE, dimension(:,:), ALLOCATABLE :: D_Matrix
  
  REALTYPE, dimension(0: N_Pnts     ,N_Index), OPTIONAL :: D_Vector
  REALTYPE, dimension(0:(N_Pnts-1)/2,N_Index), OPTIONAL :: D_Vector_Even
  REALTYPE, dimension(0:(N_Pnts-1)/2,N_Index), OPTIONAL :: D_Vector_Odd
  
  TYPE (PS_Property)   :: Property  , D_Property
  TYPE (PS_Grid_Index) :: Grid_Index, D_Grid_Index
  TYPE (PS_Domain)     :: Domain    , D_Domain
  TYPE (PS_Mapping)    :: Mapping   , D_Mapping
  
  L_D      = PRESENT(D_Vector     )
  L_D_Even = PRESENT(D_Vector_Even)
  L_D_Odd  = PRESENT(D_Vector_Odd )
  
  N = N_Pnts-1 ; N2 = N/2 
  
#if defined (DEBUG)
  call Data_Integrity_2 (L_D, L_D_Even, L_D_Odd, N, N_Index, Vector_Index)
#endif
  
  D_Property   = Property 
  D_Grid_Index = Grid_Index
  D_Domain     = Domain
  D_Mapping    = Mapping
  
  if (L_D) then
    D_Property%Max_Order_Of_D = Order
    D_Property%Algorithm      = 0 
    D_Property%Symmetry       = 0
  
    call PS_Get_Operator_Size ('D', N_Pnts, M_D, D_Property)
  
      Allocate (D_Operator(M_D), x(0:N))
  
    call PS_Setup_Operator    (N_Pnts, D_Operator, x,                      &
                                       D_Property, D_Grid_Index, D_Domain, &
                                       D_Mapping)
  
      Allocate (D_Matrix(0:N,0:N))
  
      call PS_Get_Matrix (D_Operator, D     =D_Matrix, Order=Order)
  
      do i = 1,N_Index
        D_Vector(:,i) = D_Matrix(Vector_Index(i),:)
      enddo
  
      DeAllocate (D_Matrix)
  
    DeAllocate (x, D_Operator) 
  endif
  
  if (L_D_Even .OR. L_D_Odd) then
    D_Property%Max_Order_Of_D = Order
    D_Property%Algorithm      = 1 
    D_Property%Symmetry       = 1
  
    call PS_Get_Operator_Size ('D', N_Pnts, M_D, D_Property)
  
      Allocate (D_Operator(M_D), x(0:N))
  
    call PS_Setup_Operator    (N_Pnts, D_Operator, x,                      &
                                       D_Property, D_Grid_Index, D_Domain, &
                                       D_Mapping)
  
    if (L_D_Even) then
      Allocate (D_Matrix(0:N2,0:N2))
  
      call PS_Get_Matrix (D_Operator, D_Even=D_Matrix, Order=Order)
  
      do i = 1,N_Index
        D_Vector_Even(:,i) = D_Matrix(Vector_Index(i),:)
      enddo
  
      DeAllocate (D_Matrix)
    endif
  
    if (L_D_Odd ) then
      Allocate (D_Matrix(0:N2,0:N2))
  
     call PS_Get_Matrix (D_Operator, D_Odd =D_Matrix, Order=Order)
  
      do i = 1,N_Index
        D_Vector_Odd (:,i) = D_Matrix(Vector_Index(i),:)
      enddo
  
      DeAllocate (D_Matrix)
    endif
  
    DeAllocate (x, D_Operator) 
  endif
  
  END Subroutine PS_Get_Differentiation_Vector_4
!
!------------------------------------------------------------------------------
!
#if defined (DEBUG)
  Subroutine Data_Integrity_0 (N, Vector_Index)
  
  integer                               :: N
  logical                               :: Error
  
  integer                               :: Vector_Index
  
  Error = .FALSE.
  
    if ((Vector_Index    <  0) .OR. (Vector_Index    >  N)) then
      write (6 ,100) 1, Vector_Index, N
      write (99,100) 1, Vector_Index, N
      Error = .TRUE.
    endif
  
  if (Error) STOP
  
  100 format (1x,72('*')/ &
              1x,'Input Error in PS_Get_Differentiation_Vector! ' / &
              1x,'        Index = ',i5,                           / &
              1x,' Vector_Index = ',i5, ' is < 0 or > ',i5        / &
              1x,72('*')//)
  
  END Subroutine Data_Integrity_0
!
!------------------------------------------------------------------------------
!
  Subroutine Data_Integrity_1 (N, N_Index, Vector_Index)
  
  integer                               :: N_Index, N, i
  logical                               :: Error
  
  integer , dimension(N_Index)          :: Vector_Index
  
  Error = .FALSE.
  
  do i = 1,N_Index
    if ((Vector_Index(i) <  0) .OR. (Vector_Index(i) >  N)) then
      write (6 ,100) i, Vector_Index(i), N
      write (99,100) i, Vector_Index(i), N
      Error = .TRUE.
    endif
  enddo
  
  if (Error) STOP
  
  100 format (1x,72('*')/ &
              1x,'Input Error in PS_Get_Differentiation_Vector! ' / &
              1x,'        Index = ',i5,                           / &
              1x,' Vector_Index = ',i5, ' is < 0 or > ',i5        / &
              1x,72('*')//)
  
  END Subroutine Data_Integrity_1
!
!------------------------------------------------------------------------------
!
  Subroutine Data_Integrity_2 (L_D, L_D_Even, L_D_Odd, N, N_Index, Vector_Index)
  
  integer                               :: N_Index, N, i
  logical                               :: L_D, L_D_Even, L_D_Odd
  logical                               :: Error
  
  integer , dimension(N_Index)          :: Vector_Index
  
  Error = .FALSE.
  
  do i = 1,N_Index
    if ((Vector_Index(i) <  0) .OR. (Vector_Index(i) >  N)) then
      write (6 ,100) i, Vector_Index(i), N
      write (99,100) i, Vector_Index(i), N
      Error = .TRUE.
    endif
  enddo
  
  if (.NOT. (L_D .OR. L_D_Even .OR. L_D_Odd)) then
    write (6 ,101)
    write (99,101)
    Error = .TRUE.
  endif
  
  if (Error) STOP
  
  100 format (1x,72('*')/ &
              1x,'Input Error in PS_Get_Differentiation_Vector! ' / &
              1x,'        Index = ',i5,                           / &
              1x,' Vector_Index = ',i5, ' is < 0 or > ',i5        / &
              1x,72('*')//)
  101 format (1x,72('*')/ &
              1x,'Input Error in PS_Get_Differentiation_Vector! '      / &
              1x,'Please specify either D_Vector                     ' / &
              1x,'  or at least one of D_Vector_Even and D_Vector_Odd' / &
              1x,72('*')//)
  
  END Subroutine Data_Integrity_2
#endif
  
END MODULE Differentiation_Vector
