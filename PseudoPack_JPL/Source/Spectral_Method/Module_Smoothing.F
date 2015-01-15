#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: Setup_Aux_Smth
! AUTHOR       : Wai-Sun Don and Alex Solomonoff
! AUTHOR       : Wai-Sun Don and Bruno Costa      (Rio)
! DESCRIPTION  : sets up the aux Smooth structure for the (possibly Mapped)
!                Chebyshev Solution Smoothing.
!
!  Index      = 1 if Smoothing       along the first index,
!               2 if Smoothing       along the second index.
!
!  Method     = 0-Fourier, 1-Chebyshev, 2-Legendre
!
!  Point_Type = 1-Gauss-Lobatto, 2-Gauss-Radau, 3-Gauss
!
!  Algorithm  = 0-MXM, 1-EOD, 2-CFT
!
!  N          = Total number of points
!
!  M          = total number of vectors to be differentiated.
!
!  ldy        = the leading dimension of the array to be processed.
!
!  Aux_Smth   = Storage for Smoothing Matrices, Trigonometric Table,
!                 essential parameters and Filtering.
!
!  BlockSize  = number of vectors to be processed at one time.
!               if BlockSize=0, BlockSize=M  for (EOD and CFT)
!
!  Symmetry   = 0 if the data given does not have any symmetry properties
!             = 1 if the data is either symmetric/anti-symmetry about the
!                    mid-points int((N+1)/2) and want to make use of
!                    special algorithm for computations.
!
!  Update     = 0 Return the original input data
!             = 1 update any input data that has been modified
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
!
!         Here the definition of N is not the same as above..
!
! ----------------------------------------------------------------------

MODULE PseudoPack_Smoothing

  USE MPI_Data_Distribution

  USE   Fourier_Smoothing_MXM
  USE   Fourier_Smoothing_EOD
  USE   Fourier_Smoothing_FFT
  USE Chebyshev_Smoothing_MXM
  USE Chebyshev_Smoothing_EOD
  USE Chebyshev_Smoothing_CFT
  USE  Legendre_Smoothing

implicit NONE

PRIVATE

PUBLIC  :: Setup_Aux_Smth

PUBLIC  :: Smooth

CONTAINS

  Subroutine Setup_Aux_Smth     (Index, Method, Point_Type,     &
                                 Algorithm, Symmetry,           &
                                 N, M, ldy, BlockSize,          &
                                 Filter_Choice, Mode_CutOff,    &
                                 Omega, Order,                  &
                                 Aux_Smth, Update)

  integer  :: Index, Method, Point_Type, Algorithm, Symmetry 
  integer  :: N, M, ldy, BlockSize
  integer  :: Filter_Choice, Mode_CutOff
  REALTYPE :: Omega, Order

  integer           :: Rank
  integer, OPTIONAL :: Update
  logical           :: L_Update

  REALTYPE, dimension(*) :: Aux_Smth

  integer  :: MPI_Distributed, n_0, n_Remain, n_Star, N_Block, NN
  integer  :: A_Algorithm, A_BlockSize, A_Filter_Choice
  REALTYPE :: A_Omega

#if defined(DEBUG)
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
    write (6 ,*)  'Incorrect Index (1,2 Only) : ', Index 
    write (99,*)  'Incorrect Index (1,2 Only) : ', Index 
    STOP
  endif
#endif

#if defined (PARALLEL_MPI)
                  MPI_Distributed = 0
  if (Index == 2) MPI_Distributed = 1
#endif

  Rank = Index
  if (Index >  2) Index = 2

#if defined(DEBUG)
  IO_Error = 0

  call Data_Integrity (6) ; call Data_Integrity (99)

  if (IO_Error /= 0) STOP
#endif

  L_Update = PRESENT(Update)

  if (L_Update) then 
    if (Update == 0) then
      A_Algorithm     = Algorithm
      A_BlockSize     = BlockSize
      A_Filter_Choice = Filter_Choice
      A_Omega         = Omega
    endif
  endif

  if (BlockSize == 0) BlockSize = M
                      BlockSize = MIN(BlockSize, M)

#if defined(ESSL)
#if defined(PARALLEL_OPENMP)
  if (Algorithm == 2) BlockSize = M/N_Threads
#endif
#endif

  if ((Symmetry /= 0) .AND. (Algorithm == 0)) Algorithm = 1

  if (Method == 0) then
    if ((Point_Type /= 2) .AND. (MOD(N,2) /= 0)) N = N-1
    if ((Point_Type == 2) .AND. (MOD(N,2) /= 1)) N = N-1
  endif

  call Set_Aux_Smth_Parameters 

  if (Filter_Choice == 0) RETURN

  SELECT CASE (Method)
    CASE (0)                        !  Fourier   Method
      SELECT CASE (Algorithm)
        CASE (0)                    !    MXM
          call F_Setup_Aux_Smth_MXM (Index, Method, Point_Type,   &
                                     Algorithm, Symmetry,         &
                                     N, M, ldy, BlockSize,        &
                                     Filter_Choice, Mode_CutOff,  &
                                     Omega, Order, Aux_Smth)

        CASE (1)                    !    EOD
          call F_Setup_Aux_Smth_EOD (Index, Method, Point_Type,   &
                                     Algorithm, Symmetry,         &
                                     N, M, ldy, BlockSize,        &
                                     Filter_Choice, Mode_CutOff,  &
                                     Omega, Order, Aux_Smth)
     
        CASE DEFAULT                !    FFT
          call F_Setup_Aux_Smth_FFT (Index, Method, Point_Type,   &
                                     Algorithm, Symmetry,         &
                                     N, M, ldy, BlockSize,        &
                                     Filter_Choice, Mode_CutOff,  &
                                     Omega, Order, Aux_Smth)
     
      END SELECT

    CASE (1)                        !  Chebyshev Method
      SELECT CASE (Algorithm)
        CASE (0)                    !    MXM
          call C_Setup_Aux_Smth_MXM (Index, Method, Point_Type,   &
                                     Algorithm, Symmetry,         &
                                     N, M, ldy, BlockSize,        &
                                     Filter_Choice, Mode_CutOff,  &
                                     Omega, Order, Aux_Smth)

        CASE (1)                    !    EOD
          call C_Setup_Aux_Smth_EOD (Index, Method, Point_Type,   &
                                     Algorithm, Symmetry,         &
                                     N, M, ldy, BlockSize,        &
                                     Filter_Choice, Mode_CutOff,  &
                                     Omega, Order, Aux_Smth)

        CASE DEFAULT                !    CFT
          call C_Setup_Aux_Smth_CFT (Index, Method, Point_Type,   &
                                     Algorithm, Symmetry,         &
                                     N, M, ldy, BlockSize,        &
                                     Filter_Choice, Mode_CutOff,  &
                                     Omega, Order, Aux_Smth)

      END SELECT

    CASE (2)                        !  Legendre  Method
      SELECT CASE (Algorithm)
        CASE (0)                    !    MXM
          call L_Setup_Aux_Smth_MXM (Index, Method, Point_Type,   &
                                     Algorithm, Symmetry,         &
                                     N, M, ldy, BlockSize,        &
                                     Filter_Choice, Mode_CutOff,  &
                                     Omega, Order, Aux_Smth)

        CASE DEFAULT                !    EOD
          call L_Setup_Aux_Smth_EOD (Index, Method, Point_Type,   &
                                     Algorithm, Symmetry,         &
                                     N, M, ldy, BlockSize,        &
                                     Filter_Choice, Mode_CutOff,  &
                                     Omega, Order, Aux_Smth)

        CASE (2)                    !    NOT SUPPORTED
          write (6 ,*) 'L_Setup_Aux_Smth: Invalid Legendre Algorithm=',Algorithm
          write (99,*) 'L_Setup_Aux_Smth: Invalid Legendre Algorithm=',Algorithm
          STOP

      END SELECT

  END SELECT

  call PS_MPI_Block_Distribution (N, n_Remain, n_0, n_Star, N_Block, NN)

  call Set_Aux_Smth_Parameters 

  if (L_Update) then
    if (Update == 0) then
      Algorithm     = A_Algorithm
      BlockSize     = A_BlockSize
      Filter_Choice = A_Filter_Choice
      Omega         = A_Omega
    endif
  endif

  CONTAINS

    Subroutine Set_Aux_Smth_Parameters 

#if defined (ESSL)
    Aux_Smth(_I_LIBRARY)    = 1
#endif
#if defined (VFFT)
    Aux_Smth(_I_LIBRARY)    = 2
#endif
#if defined (CRAY)
    Aux_Smth(_I_LIBRARY)    = 3
#endif

#if defined (CRAY)
    Aux_Smth(_I_PRECISION)  = 0
#else
# if defined (SINGLE_PRECISION)
    Aux_Smth(_I_PRECISION)  = 0
# else
    Aux_Smth(_I_PRECISION)  = 1
# endif
#endif

    Aux_Smth(_I_CASE)            = 0           ! 0=Smoothing, 1=Differentiation

    Aux_Smth(_I_RANK)            = Rank
    Aux_Smth(_I_INDEX)           = Index
    Aux_Smth(_I_METHOD)          = Method
    Aux_Smth(_I_POINT_TYPE)      = Point_Type
    Aux_Smth(_I_ALGORITHM)       = Algorithm
    Aux_Smth(_I_SYMMETRY)        = Symmetry
    Aux_Smth(_I_N_PNTS)          = N
    Aux_Smth(_I_N_VECTORS)       = M
    Aux_Smth(_I_BLOCKSIZE)       = BlockSize
    Aux_Smth(_I_SMOOTH_1)        = Filter_Choice
    Aux_Smth(_I_SMOOTH_2)        = Filter_Choice
    Aux_Smth(_I_MPI_DISTRIBUTED) = MPI_Distributed
    Aux_Smth(_I_N_REMAIN)        = n_Remain
    Aux_Smth(_I_N_0)             = n_0
    Aux_Smth(_I_N_STAR)          = n_Star
    Aux_Smth(_I_N_BLOCK)         = N_Block
    Aux_Smth(_I_NN)              = NN
    Aux_Smth(_I_NUMBER_THREADS)  = N_Threads

    END Subroutine Set_Aux_Smth_Parameters 

! ----------------------------------------------------------------------
#if defined (DEBUG)
    Subroutine Data_Integrity (lid)

    integer  :: lid, NN
    character(LEN=20) :: xy

    if (N <  0) then
      write (lid,*) 'Error : N = ', N,' < 0'
      write (lid,*) '----------------------------------------'
      IO_Error = 1
    endif

    if (M <  0) then
      write (lid,*) 'Error : M = ', M,' < 0'
      write (lid,*) '----------------------------------------'
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

    if ((Index <  1) .OR. (Index >  2)) then
      write (lid,*) 'Error : Index = ', Index
      write (lid,*) '        Index = 1, do the First  index'
      write (lid,*) '                   = 2, do the Second index'
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

    if ((Point_Type <  1) .OR. (Point_Type >  3)) then
      write (lid,*) 'Error : Point_Type = ', Point_Type
      write (lid,*) '        Point_Type = 1, Gauss Lobatto Point '
      write (lid,*) '                   = 2, Gauss Radau Point'
      write (lid,*) '                   = 3, Gauss Point'
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

    if ((Symmetry <  0) .OR. (Symmetry >  1)) then
      write (lid,*) 'Error : Symmetry = ', Symmetry
      write (lid,*) '        Symmetry = 0, No Symmetry'
      write (lid,*) '                 = 1, Even/Odd Symmetric Data'
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
    endif

    if (IO_Error == 1) then
      write (lid,*) 
      write (lid,*) 'Arguments Error in Setup_Aux_Smth Subroutine'
      write (lid,*) '==============================================='
    endif

    END Subroutine Data_Integrity 
#endif

  END Subroutine Setup_Aux_Smth     


! ----------------------------------------------------------------------
! FUNCTION NAME: Smooth
! DESCRIPTION  : Perform Solution Smoothing.
! 
!  M          = number of vectors to be smoothed.
!
! Optional Argument :
!
!  Even_Or_Odd= -1 if Anti-Symmetric function
!                0 if None
!                1 if Symmetric function
!                ignored if symmetry=0
! ----------------------------------------------------------------------

  Subroutine Smooth (M, y, ldy, Aux_Smth, Even_Or_Odd)

  integer  :: M, ldy

  integer, OPTIONAL :: Even_Or_Odd
  logical           :: L_Even_Or_Odd

  REALTYPE, dimension(ldy,*) :: y
  REALTYPE, dimension(*)     :: Aux_Smth

  integer  :: Index, Method, Algorithm, Symmetry, N, BlockSize, IEO

  if (OPERATOR_SMOOTH_1(Aux_Smth) == 0) RETURN

  L_Even_Or_Odd = PRESENT(Even_Or_Odd)

  Index      = OPERATOR_INDEX(Aux_Smth)
  Method     = OPERATOR_METHOD(Aux_Smth)
  Algorithm  = OPERATOR_ALGORITHM(Aux_Smth)
  Symmetry   = OPERATOR_SYMMETRY(Aux_Smth)
  N          = OPERATOR_N_PNTS(Aux_Smth)
  BlockSize  = OPERATOR_BLOCKSIZE(Aux_Smth)

  if (BlockSize == 0) BlockSize = M
                      BlockSize = MIN(BlockSize, M)

#if defined(DEBUG)
  call Data_Integrity
#endif

  if ((M == 0) .OR. (N == 1)) RETURN
    
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
              call F_Smooth_MXM_1 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (1)               !     EOD
              call F_Smooth_EOD_1 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (2)               !     FFT
              call F_Smooth_FFT_1 (N, M, y, ldy, Aux_Smth, IEO)

          END SELECT

        CASE (1)                   !   Chebyshev Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call C_Smooth_MXM_1 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (1)               !     EOD
              call C_Smooth_EOD_1 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (2)               !     CFT
              call C_Smooth_CFT_1 (N, M, y, ldy, Aux_Smth, IEO)

          END SELECT

        CASE (2)                   !   Legendre  Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call L_Smooth_MXM_1 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (1)               !     EOD
              call L_Smooth_EOD_1 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (2)               !     NOT SUPPORTED
              write (6 ,*) 'Smooth: Invalid Legendre Algorithm=',Algorithm
              write (99,*) 'Smooth: Invalid Legendre Algorithm=',Algorithm
              STOP

          END SELECT

      END SELECT

    CASE (2)                       ! Index 2
      SELECT CASE (Method)
        CASE (0)                   !   Fourier Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call F_Smooth_MXM_2 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (1)               !     EOD
              call F_Smooth_EOD_2 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (2)               !     FFT
              call F_Smooth_FFT_2 (N, M, y, ldy, Aux_Smth, IEO)

          END SELECT

        CASE (1)                   !   Chebyshev Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call C_Smooth_MXM_2 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (1)               !     EOD
              call C_Smooth_EOD_2 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (2)               !     CFT
              call C_Smooth_CFT_2 (N, M, y, ldy, Aux_Smth, IEO)

          END SELECT

        CASE (2)                   !   Legendre  Method
          SELECT CASE (Algorithm)
            CASE (0)               !     MXM
              call L_Smooth_MXM_2 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (1)               !     EOD
              call L_Smooth_EOD_2 (N, M, y, ldy, Aux_Smth, IEO)

            CASE (2)               !     NOT SUPPORTED
              write (6 ,*) 'Smooth: Invalid Legendre Algorithm=',Algorithm
              write (99,*) 'Smooth: Invalid Legendre Algorithm=',Algorithm
              STOP

          END SELECT

      END SELECT

  END SELECT

#if defined(DEBUG)
  CONTAINS

    Subroutine Data_Integrity

    ASSERT(N> 0,"Smooth:N<=0")
    ASSERT(M> 0,"Smooth:M<=0")

    ASSERT(Index>=1,"Smooth:Index<1")
    ASSERT(Index<=2,"Smooth:Index>2")

    ASSERT(Method>=0,"Smooth:Method<0")
    ASSERT(Method<=2,"Smooth:Method>2")

    ASSERT(Algorithm>=0,"Smooth:Algorithm<0")
    ASSERT(Algorithm<=2,"Smooth:Algorithm>2")

    if (L_Even_Or_Odd) then
      ASSERT(ABS(Even_Or_Odd)>=0,"Smooth:ABS(Even_Or_Odd)<0")
      ASSERT(ABS(Even_Or_Odd)<=1,"Smooth:ABS(Even_Or_Odd)>1")

      if ((Symmetry == 1) .AND. (Even_Or_Odd == 0)) then
        write (6 ,*) 'FATAL ERROR : Smooth: Symmetry=1 and Even_Or_Odd=0'
        write (99,*) 'FATAL ERROR : Smooth: Symmetry=1 and Even_Or_Odd=0'
        STOP
      endif
    endif

# if defined (ESSL)
    if (M <  BlockSize) then
      write (6 ,*) 'FATAL ERROR : Smooth: ESSL Library :'
      write (6 ,*) '              M < BlockSize', M, BlockSize
      write (99,*) 'FATAL ERROR : Smooth: ESSL Library :'
      write (99,*) '              M < BlockSize', M, BlockSize
      STOP
    endif
# endif

    END Subroutine Data_Integrity
#endif

  END Subroutine Smooth 

END MODULE PseudoPack_Smoothing
