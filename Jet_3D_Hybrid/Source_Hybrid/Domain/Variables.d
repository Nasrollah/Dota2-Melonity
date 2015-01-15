  integer , parameter :: Dimensions    = 3

  REALTYPE, dimension(:,:,:,:), TARGET, ALLOCATABLE :: Q
  REALTYPE, dimension(:,:,:,:), TARGET, ALLOCATABLE :: S1
  REALTYPE, dimension(:,:,:)  , TARGET, ALLOCATABLE :: Specific_Heat

  REALTYPE, dimension(:,:)    ,         ALLOCATABLE :: Q_0

#if defined (MULTI_RESOLUTION)
  integer , dimension(:,:,:,:)        , ALLOCATABLE :: Flag_WENO
#endif

