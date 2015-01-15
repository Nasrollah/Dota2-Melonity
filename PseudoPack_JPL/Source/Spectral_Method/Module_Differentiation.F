#include "pseudopack.h"


! ----------------------------------------------------------------------
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!  Index      = 1 if Differentiating along the first index,
!               2 if Differentiating along the second index.
!
!  Method     = 0-Fourier, 1-Chebyshev, 2-Legendre
!
!  Algorithm  = 0-MXM, 1-EOD, 2-CFT
!
!  Max_Order_Of_D = Maximum order of derivative, must be < MAX_DERIVATIVE
!
!  N          = Total number of points
!
!  M          = total number of vectors to be differentiated.
!
!  x          = Collocated points (Mapped/UnMapped)
!
!  Point_Type = 1 Lobatto Point
!             = 2 Radau   Point
!             = 3 Gauss   Point
!
!  ldy        = the leading dimension of the array to be processed.
!
!  BlockSize  = number of vectors to be processed at one time.
!               if BlockSize=0, BlockSize=M  for (EOD and CFT)
!
!  Symmetry   = 0 if the data given does not have any symmetry properties
!             = 1 if the data is either symmetric/anti-symmetry about the
!                    mid-points int((N+1)/2) and want to make use of
!                    special algorithm for computations.
!
!  x_Left     = Left /Bottom  Physical Domain Point
!  x_Right    = Right/Top     Physical Domain Point
!                (Note : For the Fourier Method,
!                        see parameters Map_F for additional info.)
!
!  For (Fourier Method)
!  Map_F      = -1 if Map [0,2 Pi] -> [x_Left, x_Right] in Compuational Space.
!
!                   In Computatonal Space, (x_Left-x_Right = 360 degrees).
!                   It is a simple translation.
!                   Both x_Left and x_Right are input in DEGREES.
!                   Upon return from the subroutine call,
!                   they will be return in RADIAN.
!
!             =  0 then [x_Left, x_Right] is defaulted to [0, 2*PI]
!                  x_Left and x_Right are not used.
!
!             =  1 if Map [0,2 Pi] -> [x_Left, x_Right] in Physical Space.
!
!                  In Physical     Space, [0, 2 Pi] <-> [x_Left, x_Right].
!                    x_Left and x_Right can be any real number.
!
!  if Map  <  0  use user defined Mapping functions. (See C_Mapping.F)
!  if Map  =  0  No mapping is used
!  if Map  >  0  use built in Mapping functions.
!
!  For Fourier Collocation methods (Method=0)
!          =  1  Angular Map (alpha, beta)
!                  if Map_F = -1,0, then
!                    As  Input, beta is in Degrees.
!                    As Output, beta is in Radians.
!                  if Map_F = 1   , then
!                    beta is a real number in [x_Left, x_Right].
!
!  For Chebyshev/Legendre Collocation methods (Method=1,2)
!
!  if  1 < abs(Map) <= 10, It is mapping for FINITE        domain mapping
!          =  1  Tal-Ezer Map (Manual, alpha)
!          =  2  Tangent Map  (alpha, beta)
!
!  beta is the point in the physical domain
!          in which points are clustered around
!
!  if 10 < abs(Map) <= 20, It is mapping for SEMI-INFINITE domain mapping
!          = 11  Algebrai! Map (alpha)
!          = 12  Exponential Map (alpha)
!          = 13  Logarithmi! Map (alpha, beta)
!
!  if 20 < abs(Map) <= 30, It is mapping for INFINITE      domain mapping
!          = 21  Tangent Map (alpha)
!          = 22  ArcTanH Map (alpha)
!          = 23  Algebrai! Map (alpha)
!
!  For Chebyshev/Legendre Method (Method=1,2) and Tal-Ezer Map (Map=1) ONLY
!  Manual     = 0 if the value of alpha is computed from machine zero.
!               1 if the value of alpha is input.
!
!  alpha      = Mapping parameter (see manual for details)
!
!  beta       = Mapping parameter (see manual for details)
!
!  Aux_Diff   = Storage for Differentiation Matrices, Trigonometric Table,
!                 essential parameters, Filtering and Transformation Metrices.
!
! OPTIONAL Argument :
!  Update     = 0 Keep the original Input data
!             = 1 Updated those input data being modified
!
!  EOD_Code   = 0 Add mapping to the EO matrix  (EOD algorithm only)
!             = 1 Add mapping only after computed the derivative. (EOD only)
!                 (OUTPUT, Internal use only)
!
! For Smoothing the Derivative...
!
!     call Filter (Filter_Choice, N, Mode_CutOff, Omega, Order, Sigma)
!
!  Filter_Choice = Choice of Filtering function
!                    1. Exponential Filter (Mode_CutOff, Omega, Order)
!                    2. Lanczos Filter (2nd order)
!                    3. Raised Cosine Filter (2nd order)
!                    4. 8th order Sharpened Rasied Cosine Filter
!  Mode_CutOff   = Cut-Off Frequency for the Exponential Filter
!  Omega         = constant in which the value of the Nth mode = Exp(-Omega)
!  Order         = Order of exponential filter = Exp(-a(k/N)**Order)
!  Sigma         = Filter function
!  Smooth_1      = Smooth the Flux BEFORE Recursion
!                    0 = NO
!                    1 = Yes
!  Smooth_2      = Smooth the Flux AFTER  Recursion
!                    0 = NO
!                    1 = Yes
!                 (NOT USED FOR Fourier Method)
!
!  Note :  Here is the way the quarter wave's index corresponding to the CFT
!          Mod(N,4) = 0 and Mod(N,2) = 0
!
!           CFT    index k = 0, 1, 2, 3, 4, 5,., N/2,....,N-1  , N
!           Q_Even index k = 0,    1,  , 2,  ,., N/4,....,     , N/2
!           Q_Odd  index k =    0,  , 1,  , 2,...........,N/2-1,     ,N/2
!
!            CFT(2k  ) = Q_Even(k)    k = 0,N/2
!            CFT(2k+1) = Q_Odd(k)     k = 0,N/2-1
!
!         Here the definition of N is not the same as above..
!
! ----------------------------------------------------------------------


MODULE PseudoPack_Differentiation

  USE MPI_Data_Distribution

  USE   Fourier_Differentiation_MXM
  USE   Fourier_Differentiation_EOD
  USE   Fourier_Differentiation_FFT
  USE Chebyshev_Differentiation_MXM
  USE Chebyshev_Differentiation_EOD
  USE Chebyshev_Differentiation_CFT
  USE  Legendre_Differentiation_MXM
  USE  Legendre_Differentiation_EOD

implicit NONE

PRIVATE

PUBLIC  :: Setup_Aux_Diff
PUBLIC  :: Setup_Aux_Diff_Base
PUBLIC  :: Setup_Aux_Diff_Map
PUBLIC  :: Setup_Aux_Diff_Filter

PUBLIC  :: Diff

CONTAINS

! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Diff_Base
! DESCRIPTION  : sets up the aux Diff structure for the Derivative.
!                Also computes the x that functions should be
!                sampled at WITHOUT Mapping and Filtering.
!
!                For Fourier Method, x_Left=0, x_Right=360
!                For Chebyshev and Legendre Method, x_Left=-1, x_Right=1
! ----------------------------------------------------------------------
  Subroutine Setup_Aux_Diff_Base   (Index, Method, Point_Type,            &
                                    Max_Order_Of_D, Algorithm, Symmetry,  &
                                    N, M, ldy, BlockSize,                 &
                                    x,                                    &
                                    Aux_Diff, Update)

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry
  integer  :: N, M, ldy, BlockSize, Map, Manual, Map_F
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: alpha, beta, Omega, Order, x_Left, x_Right

  integer, OPTIONAL  :: Update

  REALTYPE, dimension(:) :: x
  REALTYPE, dimension(*) :: Aux_Diff

  Map_F = 0 ; 
  Map   = 0 ;  Manual = 0 ; alpha  = ZERO ; beta   = ZERO

  if ((Method == 1) .OR. (Method == 2)) then
    x_Left  = -ONE ; x_Right =  ONE
  endif

  Filter_Choice = 0
  Mode_CutOff   = 0 ; Omega = ZERO ; Order = ZERO ; Smooth_1 = 0 ; Smooth_2 = 0

  call Setup_Aux_Diff        (Index, Method, Point_Type,            &
                              Max_Order_Of_D, Algorithm, Symmetry,  &
                              N, M, ldy, BlockSize,                 &
                              x, x_Left, x_Right, Map_F,            &
                              Map, Manual, alpha, beta,             &
                              Filter_Choice, Mode_CutOff,           &
                              Omega, Order, Smooth_1, Smooth_2,     &
                              Aux_Diff, Update)

  END Subroutine Setup_Aux_Diff_Base

! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Diff_Map
! DESCRIPTION  : sets up the aux Diff structure for the Derivative.
!                Also computes the x that functions should be
!                sampled at WITHOUT Filtering.
! ----------------------------------------------------------------------
  Subroutine Setup_Aux_Diff_Map    (Index, Method, Point_Type,            &
                                    Max_Order_Of_D, Algorithm, Symmetry,  &
                                    N, M, ldy, BlockSize,                 &
                                    x, x_Left, x_Right, Map_F,            &
                                    Map, Manual, alpha, beta,             &
                                    Aux_Diff, Update)

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry
  integer  :: N, M, ldy, BlockSize, Map, Manual, Map_F
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: alpha, beta, Omega, Order, x_Left, x_Right

  integer, OPTIONAL  :: Update

  REALTYPE, dimension(:) :: x
  REALTYPE, dimension(*) :: Aux_Diff

  Filter_Choice = 0
  Mode_CutOff   = 0 ; Omega = ZERO ; Order = ZERO ; Smooth_1 = 0 ; Smooth_2 = 0

  call Setup_Aux_Diff        (Index, Method, Point_Type,            &
                              Max_Order_Of_D, Algorithm, Symmetry,  &
                              N, M, ldy, BlockSize,                 &
                              x, x_Left, x_Right, Map_F,            &
                              Map, Manual, alpha, beta,             &
                              Filter_Choice, Mode_CutOff,           &
                              Omega, Order, Smooth_1, Smooth_2,     &
                              Aux_Diff, Update)

  END Subroutine Setup_Aux_Diff_Map    

! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Diff_Filter
! DESCRIPTION  : sets up the aux Diff structure for the Derivative.
!                Also computes the x that functions should be
!                sampled at WITHOUT Mapping.
!
!                For Fourier Method, x_Left=0, x_Right=360
!                For Chebyshev and Legendre Method, x_Left=-1, x_Right=1
! ----------------------------------------------------------------------
  Subroutine Setup_Aux_Diff_Filter (Index, Method, Point_Type,            &
                                    Max_Order_Of_D, Algorithm, Symmetry,  &
                                    N, M, ldy, BlockSize,                 &
                                    x,                                    &
                                    Filter_Choice, Mode_CutOff,           &
                                    Omega, Order, Smooth_1, Smooth_2,     &
                                    Aux_Diff, Update)

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry
  integer  :: N, M, ldy, BlockSize, Map, Manual, Map_F
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: alpha, beta, Omega, Order, x_Left, x_Right

  integer, OPTIONAL  :: Update

  REALTYPE, dimension(:) :: x
  REALTYPE, dimension(*) :: Aux_Diff

  Map_F = 0 ; 
  Map   = 0 ;  Manual = 0 ; alpha  = ZERO ; beta   = ZERO

  if ((Method == 1) .OR. (Method == 2)) then
    x_Left  = -ONE ; x_Right =  ONE
  endif

  call Setup_Aux_Diff        (Index, Method, Point_Type,            &
                              Max_Order_Of_D, Algorithm, Symmetry,  &
                              N, M, ldy, BlockSize,                 &
                              x, x_Left, x_Right, Map_F,            &
                              Map, Manual, alpha, beta,             &
                              Filter_Choice, Mode_CutOff,           &
                              Omega, Order, Smooth_1, Smooth_2,     &
                              Aux_Diff, Update)

  END Subroutine Setup_Aux_Diff_Filter 

! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Diff
! DESCRIPTION  : sets up the aux Diff structure for the Derivative.
!                Also computes the x that functions should be
!                sampled at.
! ----------------------------------------------------------------------
  Subroutine Setup_Aux_Diff        (Index, Method, Point_Type,            &
                                    Max_Order_Of_D, Algorithm, Symmetry,  &
                                    N, M, ldy, BlockSize,                 &
                                    x, x_Left, x_Right, Map_F,            &
                                    Map, Manual, alpha, beta,             &
                                    Filter_Choice, Mode_CutOff,           &
                                    Omega, Order, Smooth_1, Smooth_2,     &
                                    Aux_Diff, Update)

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry
  integer  :: N, M, ldy, BlockSize, Map, Manual, Map_F
  integer  :: Filter_Choice, Mode_CutOff, Smooth_1, Smooth_2
  REALTYPE :: alpha, beta, Omega, Order, x_Left, x_Right

  integer, OPTIONAL  :: Update
  logical            :: L_Update

  REALTYPE, dimension(:) :: x
  REALTYPE, dimension(*) :: Aux_Diff

  integer  :: MPI_Distributed, n_0, n_Remain, n_Star, N_Block, NN
  integer  :: Storage, EOD_Code, Rank

  integer  :: A_Algorithm, A_BlockSize
  integer  :: A_Map, A_Filter_Choice, A_Diff_1, A_Diff_2
  REALTYPE :: A_alpha, A_beta, A_Omega 

#if defined (DEBUG)
  integer  :: IO_Error
#endif

  integer  :: i, N_Threads = 1
#if defined (PARALLEL_OPENMP)
  integer  :: OMP_GET_NUM_THREADS

!$OMP PARALLEL DO
  do i = 1,2
    N_Threads = OMP_GET_NUM_THREADS()
  enddo
!$OMP END PARALLEL DO
#endif

#if 0
  if ((Index <  1) .OR. (Index > 2)) then
    write (6 ,*)  'Incorrect Index (1, 2 ONLY) : ', Index 
    write (99,*)  'Incorrect Index (1, 2 ONLY) : ', Index 
    STOP
  endif
#endif

#if defined (PARALLEL_MPI)
                  MPI_Distributed = 0
  if (Index == 2) MPI_Distributed = 1
#endif

  Rank = Index
  if (Index > 2) Index = 2

#if defined (DEBUG)
  IO_Error = 0

  call Data_Integrity (6) ; call Data_Integrity (99)

  if (IO_Error == 1) STOP
#endif

  L_Update = PRESENT(Update)

  if (L_Update) then 
    if (Update == 0) then
      A_Algorithm     = Algorithm
      A_BlockSize     = BlockSize
      A_Map           = Map
      A_alpha         = alpha
      A_beta          = beta
      A_Filter_Choice = Filter_Choice
      A_Omega         = Omega
      A_Diff_1        = Smooth_1
      A_Diff_2        = Smooth_2
    endif
  endif

  if ((Symmetry == 1) .AND.  &
      ((ABS(Map) >  10) .AND. (ABS(Map) <= 20))) Symmetry = 0

  if (BlockSize == 0) BlockSize = M
                      BlockSize = MIN(BlockSize, M)

#if defined(ESSL)
#if defined(PARALLEL_OPENMP)
  if (Algorithm == 2) BlockSize = M/N_Threads
#endif
#endif

  if ((Symmetry /= 0) .AND. (Algorithm == 0)) Algorithm = 1

  if (Method == 0) then
    if  (Smooth_1 == 0)                        Filter_Choice = 0
  else
    if ((Smooth_1 == 0) .AND. (Smooth_2 == 0)) Filter_Choice = 0
  endif

   if ( (Algorithm == 0) .OR. ((Algorithm == 1) .AND. (Symmetry == 1)) ) then
     Storage = 0
   else
     Storage = 1
   endif

  if (Method == 0) then
    if ((Point_Type /= 2) .AND. (MOD(N,2) /= 0)) N = N-1
    if ((Point_Type == 2) .AND. (MOD(N,2) /= 1)) N = N-1
  endif

  call Setup_Aux_Diff_Parameters 

  SELECT CASE (Method)
    CASE (0)                        !  Fourier
      SELECT CASE (Algorithm)
        CASE (0)                    !    MXM
          call F_Setup_Aux_Diff_MXM (Index, Method, Point_Type,           &
                                     Max_Order_Of_D, Algorithm, Symmetry, &
                                     N, M, ldy, BlockSize,                &
                                     x,  x_Left, x_Right,                 &
                                     Map, Map_F , alpha, beta,            &
                                     Filter_Choice, Mode_CutOff,          &
                                     Omega, Order, Smooth_1, Smooth_2,    &
                                     Aux_Diff, EOD_Code)

        CASE (1)                    !    EOD
          call F_Setup_Aux_Diff_EOD (Index, Method, Point_Type,           &
                                     Max_Order_Of_D, Algorithm, Symmetry, &
                                     N, M, ldy, BlockSize,                &
                                     x,  x_Left, x_Right,                 &
                                     Map, Map_F , alpha, beta,            &
                                     Filter_Choice, Mode_CutOff,          &
                                     Omega, Order, Smooth_1, Smooth_2,    &
                                     Aux_Diff, EOD_Code)

        CASE (2)                    !    FFT
          call F_Setup_Aux_Diff_FFT (Index, Method, Point_Type,           &
                                     Max_Order_Of_D, Algorithm, Symmetry, &
                                     N, M, ldy, BlockSize,                &
                                     x,  x_Left, x_Right,                 &
                                     Map, Map_F , alpha, beta,            &
                                     Filter_Choice, Mode_CutOff,          &
                                     Omega, Order, Smooth_1, Smooth_2,    &
                                     Aux_Diff, EOD_Code)
      END SELECT

    CASE (1)                        !  Chebyshev
      SELECT CASE (Algorithm)
        CASE (0)                    !    MXM
          call C_Setup_Aux_Diff_MXM (Index, Method, Point_Type,           &
                                     Max_Order_Of_D, Algorithm, Symmetry, &
                                     N, M, ldy, BlockSize,                &
                                     x,  x_Left, x_Right,                 &
                                     Map, Manual, alpha, beta,            &
                                     Filter_Choice, Mode_CutOff,          &
                                     Omega, Order, Smooth_1, Smooth_2,    &
                                     Aux_Diff, EOD_Code)

        CASE (1)                    !    EOD
          call C_Setup_Aux_Diff_EOD (Index, Method, Point_Type,           &
                                     Max_Order_Of_D, Algorithm, Symmetry, &
                                     N, M, ldy, BlockSize,                &
                                     x,  x_Left, x_Right,                 &
                                     Map, Manual, alpha, beta,            &
                                     Filter_Choice, Mode_CutOff,          &
                                     Omega, Order, Smooth_1, Smooth_2,    &
                                     Aux_Diff, EOD_Code)

        CASE (2)                    !    CFT
          call C_Setup_Aux_Diff_CFT (Index, Method, Point_Type,           &
                                     Max_Order_Of_D, Algorithm, Symmetry, &
                                     N, M, ldy, BlockSize,                &
                                     x,  x_Left, x_Right,                 &
                                     Map, Manual, alpha, beta,            &
                                     Filter_Choice, Mode_CutOff,          &
                                     Omega, Order, Smooth_1, Smooth_2,    &
                                     Aux_Diff, EOD_Code)
      END SELECT

    CASE (2)                        !  Legendre
      SELECT CASE (Algorithm)
        CASE (0)                    !    MXM
          call L_Setup_Aux_Diff_MXM (Index, Method, Point_Type,           &
                                     Max_Order_Of_D, Algorithm, Symmetry, &
                                     N, M, ldy, BlockSize,                &
                                     x,  x_Left, x_Right,                 &
                                     Map, Manual, alpha, beta,            &
                                     Filter_Choice, Mode_CutOff,          &
                                     Omega, Order, Smooth_1, Smooth_2,    &
                                     Aux_Diff, EOD_Code)

        CASE (1)                    !    EOD
          call L_Setup_Aux_Diff_EOD (Index, Method, Point_Type,           &
                                     Max_Order_Of_D, Algorithm, Symmetry, &
                                     N, M, ldy, BlockSize,                &
                                     x,  x_Left, x_Right,                 &
                                     Map, Manual, alpha, beta,            &
                                     Filter_Choice, Mode_CutOff,          &
                                     Omega, Order, Smooth_1, Smooth_2,    &
                                     Aux_Diff, EOD_Code)

        CASE (2)                    !    NOT SUPPORTED
          write (6 ,*) 'L_Setup_Aux_Diff: Invalid Legendre Algorithm=',Algorithm
          write (99,*) 'L_Setup_Aux_Diff: Invalid Legendre Algorithm=',Algorithm
          STOP

      END SELECT

  END SELECT

  call PS_MPI_Block_Distribution (N, n_Remain, n_0, n_Star, N_Block, NN)

  call Setup_Aux_Diff_Parameters 

! if (Map == 99) Map = 0

  if (L_Update) then 
    if (Update == 0) then
      Algorithm     = A_Algorithm
      BlockSize     = A_BlockSize
      Map           = A_Map
      alpha         = A_alpha
      beta          = A_beta
      Filter_Choice = A_Filter_Choice
      Omega         = A_Omega
      Smooth_1      = A_Diff_1
      Smooth_2      = A_Diff_2
    endif
  endif

  CONTAINS

    Subroutine Setup_Aux_Diff_Parameters 

#if defined (ESSL)
    Aux_Diff(_I_LIBRARY)        = 1
#endif
#if defined (VFFT)
    Aux_Diff(_I_LIBRARY)        = 2
#endif
#if defined (CRAY)
    Aux_Diff(_I_LIBRARY)        = 3
#endif

#if defined (CRAY)
    Aux_Diff(_I_PRECISION)      = 0
#else
# if defined (SINGLE_PRECISION)
    Aux_Diff(_I_PRECISION)      = 0
# else
    Aux_Diff(_I_PRECISION)      = 1
# endif
#endif

    Aux_Diff(_I_CASE)            = 1       ! 0=Smoothing, 1=Differentiation
    Aux_Diff(_I_RANK)            = Rank
    Aux_Diff(_I_INDEX)           = Index
    Aux_Diff(_I_METHOD)          = Method
    Aux_Diff(_I_POINT_TYPE)      = Point_Type
    Aux_Diff(_I_MAX_ORDER_OF_D)  = Max_Order_Of_D
    Aux_Diff(_I_ALGORITHM)       = Algorithm
    Aux_Diff(_I_SYMMETRY)        = Symmetry
    Aux_Diff(_I_N_PNTS)          = N
    Aux_Diff(_I_N_VECTORS)       = M
    Aux_Diff(_I_BLOCKSIZE)       = BlockSize
    Aux_Diff(_I_X_LEFT)          = x_Left
    Aux_Diff(_I_X_RIGHT)         = x_Right
    Aux_Diff(_I_MAP_F)           = Map_F
    Aux_Diff(_I_MAP)             = Map
    Aux_Diff(_I_MANUAL)          = Manual
    Aux_Diff(_I_ALPHA)           = alpha
    Aux_Diff(_I_BETA)            = beta
    Aux_Diff(_I_SMOOTH_1)        = Smooth_1
    Aux_Diff(_I_SMOOTH_2)        = Smooth_2
    Aux_Diff(_I_EOD_CODE)        = EOD_Code
    Aux_Diff(_I_EPSILON)         = EPSILON(ONE)
    Aux_Diff(_I_STORAGE)         = Storage
    Aux_Diff(_I_MPI_DISTRIBUTED) = MPI_Distributed
    Aux_Diff(_I_N_REMAIN)        = n_Remain
    Aux_Diff(_I_N_0)             = n_0
    Aux_Diff(_I_N_STAR)          = n_Star
    Aux_Diff(_I_N_BLOCK)         = n_Block
    Aux_Diff(_I_NN)              = NN
    Aux_Diff(_I_NUMBER_THREADS)  = N_Threads

    END Subroutine Setup_Aux_Diff_Parameters 

! ----------------------------------------------------------------------

#if defined (DEBUG)
    Subroutine Data_Integrity (lid)

    integer   :: lid, NN
    character(LEN=20) :: xy

    if (Max_Order_Of_D >  MAX_DERIVATIVE) then
      write (lid,*) 'Error : Max_Order_Of_D =', Max_Order_Of_D,   &
                                            ' > MAX_DERIVATIVE'
      IO_Error = 1
    endif

    if (N <= 0) then
      write (lid,*) 'Error : N = ', N,' <= 0'
      IO_Error = 1
    endif

    if (M <  0) then
      write (lid,*) 'Error : M = ', M,' < 0'
      IO_Error = 1
    endif
 
    if (Symmetry == 0) then
      if (Index == 1) NN = N
      if (Index == 2) NN = M
    else
      if (Method == 0) then
        if (Index == 1) NN = N/2+1
        if (Index == 2) NN = M/2+1
      else
        if (Index == 1) NN = (N+1)/2
        if (Index == 2) NN = (M+1)/2
      endif
    endif

    if (Index == 1) xy = 'First  Index'
    if (Index == 2) xy = 'Second Index'
  
    if (ldy <  NN) then
      write (lid,*) 'Error : ',TRIM(xy),' - ldy < NN => ',ldy,'<',NN
      write (lid,*) '----------------------------------------'
      IO_Error = 1
    endif

    if ((ABS(Map) <= 10) .AND. (x_Right == x_Left)) then
      write (lid,*) 'Error :          Index = ', Index
      write (lid,*) '                   Map = ', Map,' < 10'
      write (lid,*) '        x_Left=x_Right = ', x_Right
      write (lid,*) '----------------------------------------'
      IO_Error = 1 
    endif

    if ((Symmetry ==  1) .AND. ((ABS(Map) >  10) .AND. (ABS(Map) <= 20))) then
      write (lid,*) 'Error : Semi-Infinite Map = ',Map
      write (lid,*) '                 Symmetry = ',Symmetry
      write (lid,*) '        No Symmetry is allowed for Semi-Infinte Map'
      write (lid,*) '----------------------------------------'
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

    if ((Algorithm <  0) .OR. (Algorithm >  2)) then
      write (lid,*) 'Error : Algorithm = ', Algorithm
      write (lid,*) '        Algorithm = 0,  MXM       Algorithm'
      write (lid,*) '                  = 1,  EOD       Algorithm'
      write (lid,*) '                  = 2,  Transform Algorithm'
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
        IO_Error = 1
      endif

      if ((Point_Type == 2) .AND. (MOD(N,2) /= 1)) then
        write (lid,*) 'Error : Method     = ', Method,' Fourier'
        write (lid,*) '        Point_Type = ', Point_Type,' Radau'
        write (lid,*) '                 N = ', N,' must be an ODD  number'
        write (lid,*) '----------------------------------------'
        IO_Error = 1
      endif
    endif
   
    if ((Symmetry <  0) .OR. (Symmetry >  1)) then
      write (lid,*) 'Error : Symmetry = ', Symmetry
      write (lid,*) '        Symmetry = 0, No Symmetry'
      write (lid,*) '                 = 1, Even/Odd Symmetric Data'
      write (lid,*) '----------------------------------------'
      IO_Error = 1
    endif

    if ((Method /= 2) .AND. (Algorithm == 2)) then
      if (Point_Type == 2) then
        if (Method == 0) write (lid,*) 'FFT for Fourier Gauss Radau point ', &
                                       'is not supported!'
        if (Method == 1) write (lid,*) 'CFT for Fourier Gauss Radau point ', &
                                       'is not supported!'
        IO_Error = 1
      endif

      if (Point_Type == 3) then
        if (Method == 1) write (lid,*) 'CFT for Chebyshev Gauss point ',  &
                                       'is not supported!'
        IO_Error = 1
      endif
    endif

    if ((Method == 0) .AND. ((Map_F <  -1) .OR. (Map_F >  1))) then
      write (lid,*) 'Error : Map_F = ', Map_F
      write (lid,*) '        Map_F =-1, Fourier Linear Translation'
      write (lid,*) '              = 0, Defaulted [0,2 Pi]'
      write (lid,*) '              = 1, Fourier Linear Mapping'
      write (lid,*) '----------------------------------------'
      IO_Error = 1
    endif

    if ((Method /= 0) .AND. ((Manual <  0) .OR. (Manual >  1))) then
      write (lid,*) 'Error : Manual = ', Manual
      write (lid,*) '        Manual = 0, machine computed alpha'
      write (lid,*) '               = 1, manually specify alpha'
      write (lid,*) '----------------------------------------'
      IO_Error = 1
    endif

    if ((Method == 2) .AND. (Algorithm == 2)) then
      write (lid,*) 'Error : Method = ',Method,' Algorithm =', Algorithm
      write (lid,*) '        Legendre Method (=2) does not have '
      write (lid,*) '        Transform Algorithm (=2)'
      write (lid,*) '----------------------------------------'
      IO_Error = 1
    endif

    if ((Filter_Choice <  0) .OR. (Filter_Choice >  4)) then
      write (lid,*) 'Error : Filter_Choice = ', Filter_Choice
      write (lid,*) '        Filter_Choice = 0, No Filtering'
      write (lid,*) '                      = 1, Exponential filter'
      write (lid,*) '                      = 2, Lanczos filter'
      write (lid,*) '                      = 3, Raised Cosine'
      write (lid,*) '                      = 4, O(8) Raised Cosine'
      write (lid,*) '----------------------------------------'
      IO_Error = 1
    else
      if ((Smooth_1 <  0) .OR. (Smooth_1 >  1)) then
        write (lid,*) 'Error : Smooth_1 = ', Smooth_1
        write (lid,*) '        Smooth_1 = 0, No'
        write (lid,*) '                 = 1, Yes'
        write (lid,*) '----------------------------------------'
        IO_Error = 1
      endif

      if ((Method /= 0) .AND. ((Smooth_2 <  0) .OR. (Smooth_2 >  1))) then
        write (lid,*) 'Error : Smooth_2 = ', Smooth_2
        write (lid,*) '        Smooth_2 = 0, No'
        write (lid,*) '                 = 1, Yes'
        write (lid,*) '----------------------------------------'
        IO_Error = 1
      endif
    endif

    if (IO_Error == 1) then
      write (lid,*) 
      write (lid,*) 'Arguments Error in Setup_Aux_Diff Subroutine'
      write (lid,*) '==============================================='
    endif

    END Subroutine Data_Integrity
#endif

  END Subroutine Setup_Aux_Diff        



! ----------------------------------------------------------------------
! FUNCTION NAME: Diff
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : Branch to proper differentiation routines.
!
!  M  = Total number of Vectors  to be differentiated.
!
! Optional Argument :
!
!  Even_Or_Odd= -1 if Anti-Symmetric function
!                0 if None
!                1 if Symmetric function
!
!                Ignored if Symmetry=0
! ----------------------------------------------------------------------

  Subroutine Diff (M, O_Of_D, y, dy, ldy, Aux_Diff, Even_Or_Odd)

  integer           :: M, O_Of_D, ldy

  integer, OPTIONAL :: Even_Or_Odd
  logical           :: L_Even_Or_Odd

  REALTYPE, dimension(ldy,*) :: y, dy
  REALTYPE, dimension(*)     :: Aux_Diff

  integer  :: Index, Method, Algorithm, Symmetry, N, BlockSize, IEO

  L_Even_Or_Odd = PRESENT(Even_Or_Odd)

  Index      = OPERATOR_INDEX(Aux_Diff)
  Method     = OPERATOR_METHOD(Aux_Diff)
  Algorithm  = OPERATOR_ALGORITHM(Aux_Diff)
  Symmetry   = OPERATOR_SYMMETRY(Aux_Diff)
  N          = OPERATOR_N_PNTS(Aux_Diff)
  BlockSize  = OPERATOR_BLOCKSIZE(Aux_Diff)

  if (BlockSize == 0) BlockSize = M
                      BlockSize = MIN(BlockSize, M)

#if defined (DEBUG)
  call Data_Integrity
#endif

  if (M == 0) RETURN
  if (N == 1) then ; dy(1:N,1:M) = ZERO ; RETURN ; endif 

  if (L_Even_Or_Odd .AND. (Symmetry /= 0)) then
    IEO = Even_Or_Odd
  else
    IEO = 0
  endif

  SELECT CASE (Index)
    CASE (1)                       ! Index 1
      SELECT CASE (Method)
        CASE (0)                   !   Fourier Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call F_Diff_MXM_1 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (1)               !     EOD
              call F_Diff_EOD_1 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (2)               !     FFT
              call F_Diff_FFT_1 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

          END SELECT

        CASE (1)                   !   Chebyshev Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call C_Diff_MXM_1 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (1)               !     EOD
              call C_Diff_EOD_1 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (2)               !     CFT
              call C_Diff_CFT_1 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

          END SELECT

        CASE (2)                   !   Legendre  Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call L_Diff_MXM_1 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (1)               !     EOD
              call L_Diff_EOD_1 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (2)               !     NOT SUPPORTED
              write (6 ,*) 'Diff: Invalid Legendre Algorithm=',Algorithm
              write (99,*) 'Diff: Invalid Legendre Algorithm=',Algorithm
              STOP

          END SELECT

      END SELECT

    CASE (2)                       ! Index 2
      SELECT CASE (Method)
        CASE (0)                   !   Fourier Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call F_Diff_MXM_2 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (1)               !     EOD
              call F_Diff_EOD_2 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (2)               !     FFT
              call F_Diff_FFT_2 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

          END SELECT

        CASE (1)                   !   Chebyshev Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call C_Diff_MXM_2 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (1)               !     EOD
              call C_Diff_EOD_2 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (2)               !     CFT
              call C_Diff_CFT_2 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

          END SELECT

        CASE (2)                   !   Legendre  Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call L_Diff_MXM_2 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (1)               !     EOD
              call L_Diff_EOD_2 (N, M, O_Of_D, y, dy, ldy, Aux_Diff, IEO)

            CASE (2)               !     NOT SUPPORTED
              write (6 ,*) 'Diff: Invalid Legendre Algorithm=',Algorithm
              write (99,*) 'Diff: Invalid Legendre Algorithm=',Algorithm
              STOP

          END SELECT

      END SELECT

  END SELECT

#if defined(DEBUG)
  CONTAINS

    Subroutine Data_Integrity

    ASSERT(N> 0,"Diff:N<=0")
    ASSERT(M> 0,"Diff:M<=0")

    ASSERT(Index>=1,"Diff:Index<1")
    ASSERT(Index<=2,"Diff:Index>2")

    ASSERT(Method>=0,"Diff:Method<0")
    ASSERT(Method<=2,"Diff:Method>2")

    ASSERT(Algorithm>=0,"Diff:Algorithm<0")
    ASSERT(Algorithm<=2,"Diff:Algorithm>2")

    if (L_Even_Or_Odd) then
      ASSERT(ABS(Even_Or_Odd)>=0,"Diff:ABS(Even_Or_Odd)<0")
      ASSERT(ABS(Even_Or_Odd)<=1,"Diff:ABS(Even_Or_Odd)>1")

      if ((Symmetry == 1) .AND. (Even_Or_Odd == 0)) then
        write (6 ,*) 'FATAL ERROR : Diff: Symmetry=1 and Even_Or_Odd=0'
        write (99,*) 'FATAL ERROR : Diff: Symmetry=1 and Even_Or_Odd=0'
        STOP
      endif
    endif

# if defined (ESSL)
    if (M <  BlockSize) then
      write (6 ,*) 'FATAL ERROR : Diff: ESSL Library :'
      write (6 ,*) '              M < BlockSize', M, BlockSize
      write (99,*) 'FATAL ERROR : Diff: ESSL Library :'
      write (99,*) '              M < BlockSize', M, BlockSize
      STOP
    endif
# endif

    END Subroutine Data_Integrity
#endif

  END Subroutine Diff 

END MODULE PseudoPack_Differentiation
