integer  ::  Index_x, Method_x, Point_Type_x, Max_Order_Of_D_x, &
               Algorithm_x, Symmetry_x,                         &
             N_x, BlockSize_x, Map_x, Manual_x, Map_F_x,        &
             Filter_Choice_D_x, Mode_CutOff_D_x,                &
               Smooth_1_x, Smooth_2_x,                          &
             Filter_Choice_S_x, Mode_CutOff_S_x
REALTYPE ::  Left_x, Right_x, alpha_x, Angle_x,                 &
             Omega_D_x, Order_D_x, Omega_S_x, Order_S_x

logical  ::   Distributed_x
integer  :: Spatial_Order_x, M_s_x, O_s_x, Degree_x
integer  :: NN_g, N_g, N0_g, N1_g, N2_g, N3_g, N4_g, N5_g
integer  :: NN, N0, N1, N2, N3, N4, N5, n_Remain, n_0, n_Star

logical  :: x_Symmetry
REALTYPE :: x0, x1
REALTYPE :: dx = ZERO

integer  :: WENO_Order_x
integer  :: WENO_Flux_Type_x
integer  :: WENO_Weight_Type_x
integer  :: WENO_Unroll_Loop_x
REALTYPE :: WENO_Epsilon_x
REALTYPE :: WENO_Power_x
REALTYPE :: WENO_Lambda_Factor_x
REALTYPE :: WENO_Sound_Speed_Min_x

  REALTYPE, dimension(:)    , ALLOCATABLE         :: xi, x_1D
  REALTYPE, dimension(:)    , ALLOCATABLE         :: I_x, D_x, S_x, R_x

  REALTYPE, dimension(:,:,:)  , ALLOCATABLE, TARGET ::  x
  REALTYPE, dimension(:,:,:,:), ALLOCATABLE, TARGET :: Q0_N0, Q0_N5
  REALTYPE, dimension(:,:,:,:), ALLOCATABLE, TARGET :: Normal_N0, Normal_N5
