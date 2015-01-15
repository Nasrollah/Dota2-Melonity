integer  ::  Index_z, Method_z, Point_Type_z, Max_Order_Of_D_z, &
               Algorithm_z, Symmetry_z,                         &
             N_z, BlockSize_z, Map_z, Manual_z, Map_F_z,        &
             Filter_Choice_D_z, Mode_CutOff_D_z,                &
               Smooth_1_z, Smooth_2_z,                          &
             Filter_Choice_S_z, Mode_CutOff_S_z
REALTYPE ::  Left_z, Right_z, alpha_z, Angle_z,                 &
             Omega_D_z, Order_D_z, Omega_S_z, Order_S_z

logical  ::   Distributed_z
integer  :: Spatial_Order_z, M_s_z, O_s_z, Degree_z
integer  :: KK_g, K_g, K0_g, K1_g, K2_g, K3_g, K4_g, K5_g
integer  :: KK, K0, K1, K2, K3, K4, K5, k_Remain, k_0, k_Star

logical  :: z_Symmetry
REALTYPE :: z0, z1
REALTYPE :: dz = ZERO

integer  :: WENO_Order_z
integer  :: WENO_Flux_Type_z
integer  :: WENO_Weight_Type_z
integer  :: WENO_Unroll_Loop_z
REALTYPE :: WENO_Epsilon_z
REALTYPE :: WENO_Power_z
REALTYPE :: WENO_Lambda_Factor_z
REALTYPE :: WENO_Sound_Speed_Min_z

  REALTYPE, dimension(:)      , ALLOCATABLE         :: fi, z_1D
  REALTYPE, dimension(:)      , ALLOCATABLE         :: I_z, D_z, S_z, R_z

  REALTYPE, dimension(:,:,:)  , ALLOCATABLE, TARGET ::  z
  REALTYPE, dimension(:,:,:,:), ALLOCATABLE, TARGET :: Q0_K0, Q0_K5
  REALTYPE, dimension(:,:,:,:), ALLOCATABLE, TARGET :: Normal_K0, Normal_K5
