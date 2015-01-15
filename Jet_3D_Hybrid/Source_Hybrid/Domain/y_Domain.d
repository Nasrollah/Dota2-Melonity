integer  ::  Index_y, Method_y, Point_Type_y, Max_Order_Of_D_y, &
               Algorithm_y, Symmetry_y,                         &
             N_y, BlockSize_y, Map_y, Manual_y, Map_F_y,        &
             Filter_Choice_D_y, Mode_CutOff_D_y,                &
               Smooth_1_y, Smooth_2_y,                          &
             Filter_Choice_S_y, Mode_CutOff_S_y
REALTYPE ::  Left_y, Right_y, alpha_y, Angle_y,                 &
             Omega_D_y, Order_D_y, Omega_S_y, Order_S_y

logical  ::   Distributed_y
integer  :: Spatial_Order_y, M_s_y, O_s_y, Degree_y
integer  :: MM_g, M_g, M0_g, M1_g, M2_g, M3_g, M4_g, M5_g
integer  :: MM, M0, M1, M2, M3, M4, M5, m_Remain, m_0, m_Star

logical  :: y_Symmetry
REALTYPE :: y0, y1
REALTYPE :: dy = ZERO

integer  :: WENO_Order_y
integer  :: WENO_Flux_Type_y
integer  :: WENO_Weight_Type_y
integer  :: WENO_Unroll_Loop_y
REALTYPE :: WENO_Epsilon_y
REALTYPE :: WENO_Power_y
REALTYPE :: WENO_Lambda_Factor_y
REALTYPE :: WENO_Sound_Speed_Min_y

  REALTYPE, dimension(:)      , ALLOCATABLE         :: eta, y_1D
  REALTYPE, dimension(:)      , ALLOCATABLE         :: I_y, D_y, S_y, R_y

  REALTYPE, dimension(:,:,:)  , ALLOCATABLE, TARGET :: y
  REALTYPE, dimension(:,:,:,:), ALLOCATABLE, TARGET :: Q0_M0, Q0_M5
  REALTYPE, dimension(:,:,:,:), ALLOCATABLE, TARGET :: Normal_M0, Normal_M5

