#if defined (THREE_DIMENSION)
#define _GRID_G			N0_g:N5_g,M0_g:M5_g,K0_g:K5_g
#define _GRID_G_NV		N0_g:N5_g,M0_g:M5_g,K0_g:K5_g,1:NV
#define _SUB_GRID_G		N2_g:N3_g,M2_g:M3_g,K2_g:K3_g
#define _SUB_GRID_G_NV		N2_g:N3_g,M2_g:M3_g,K2_g:K3_g,1:NV

#define _GRID			N0:N5,M0:M5,K0:K5
#define _GRID_NV		N0:N5,M0:M5,K0:K5,1:NV
#define _SUB_GRID		N2:N3,M2:M3,K2:K3
#define _SUB_GRID_NV		N2:N3,M2:M3,K2:K3,1:NV
#define _GRID_ND		N0:N5,M0:M5,K0:K5,1:ND

#define _REGION_G		N0_g,N5_g,N2_g,N3_g, M0_g,M5_g,M2_g,M3_g, K0_g,K5_g,K2_g,K3_g
#define _REGION			N0,N5,N2,N3, M0,M5,M2,M3, K0,K5,K2,K3
#define _REGION_DIFF		N0,N5,N2,N3, WENO_Order_x,dx, M0,M5,M2,M3, WENO_Order_y,dy, K0,K5,K2,K3, WENO_Order_z,dz

#define _REGION_DIFF_WENO_X	N0,N5,N2,N3, WENO_Order_x,dx
#define _REGION_DIFF_WENO_Y	M0,M5,M2,M3, WENO_Order_y,dy
#define _REGION_DIFF_WENO_Z	K0,K5,K2,K3, WENO_Order_z,dz
#define _REGION_DIFF_WENO	_REGION_DIFF_WENO_X, _REGION_DIFF_WENO_Y, _REGION_DIFF_WENO_Z

#define _REGION_DIFF_CENTRAL_X	_REGION_DIFF_WENO_X, D_Order_x,Central_Degree_x
#define _REGION_DIFF_CENTRAL_Y	_REGION_DIFF_WENO_Y, D_Order_y,Central_Degree_y
#define _REGION_DIFF_CENTRAL_Z	_REGION_DIFF_WENO_Z, D_Order_z,Central_Degree_z
#define _REGION_DIFF_CENTRAL  	_REGION_DIFF_CENTRAL_X, _REGION_DIFF_CENTRAL_Y, _REGION_DIFF_CENTEAL_Z

#define _GHOST_CELL		Index_x,Distributed_x,M_s_x, Index_y,Distributed_y,M_s_y, Index_z,Distributed_z,M_s_z
#endif

#if defined (TWO_DIMENSION)
#define _GRID_G			N0_g:N5_g,M0_g:M5_g
#define _GRID_G_NV		N0_g:N5_g,M0_g:M5_g,1:NV
#define _SUB_GRID_G		N2_g:N3_g,M2_g:M3_g
#define _SUB_GRID_G_NV		N2_g:N3_g,M2_g:M3_g,1:NV

#define _GRID			N0:N5,M0:M5
#define _GRID_NV		N0:N5,M0:M5,1:NV
#define _SUB_GRID		N2:N3,M2:M3
#define _SUB_GRID_NV		N2:N3,M2:M3,1:NV
#define _GRID_ND		N0:N5,M0:M5,1:ND

#define _REGION			N0,N5,N2,N3, M0,M5,M2,M3
#define _REGION_DIFF		N0,N5,N2,N3, WENO_Order_x,dx, M0,M5,M2,M3, WENO_Order_y,dy

#define _REGION_DIFF_WENO_X	N0,N5,N2,N3, WENO_Order_x,dx
#define _REGION_DIFF_WENO_Y	M0,M5,M2,M3, WENO_Order_y,dy
#define _REGION_DIFF_WENO	_REGION_DIFF_WENO_X, _REGION_DIFF_WENO_Y

#define _REGION_DIFF_CENTRAL_X	_REGION_DIFF_WENO_X, D_Order_x,Central_Degree_x
#define _REGION_DIFF_CENTRAL_Y	_REGION_DIFF_WENO_Y, D_Order_y,Central_Degree_y
#define _REGION_DIFF_CENTRAL  	_REGION_DIFF_CENTRAL_X, _REGION_DIFF_CENTRAL_Y

#define _GHOST_CELL		Index_x,Distributed_x,M_s_x, Index_y,Distributed_y,M_s_y
#endif

#if defined (ONE_DIMENSION)
#define _GRID_G			N0_g:N5_g
#define _GRID_G_NV		N0_g:N5_g,1:NV
#define _SUB_GRID_G		N2_g:N3_g
#define _SUB_GRID_G_NV		N2_g:N3_g,1:NV

#define _GRID			N0:N5
#define _GRID_NV		N0:N5,1:NV
#define _SUB_GRID		N2:N3
#define _SUB_GRID_NV		N2:N3,1:NV
#define _GRID_ND		N0:N5,1:ND

#define _REGION			N0,N5,N2,N3
#define _REGION_DIFF		N0,N5,N2,N3, WENO_Order_x,dx

#define _REGION_DIFF_WENO_X	N0,N5,N2,N3, WENO_Order_x,dx
#define _REGION_DIFF_WENO	_REGION_DIFF_WENO_X

#define _REGION_DIFF_CENTRAL_X	_REGION_DIFF_WENO_X, D_Order_x,Central_Degree_x
#define _REGION_DIFF_CENTRAL  	_REGION_DIFF_CENTRAL_X

#define _GHOST_CELL		Index_x,Distributed_x,M_s_x
#endif

#define _PS_INPUT_X_1	Index_x, Method_x, Point_Type_x, Max_Order_Of_D_x, Algorithm_x, Symmetry_x
#define _PS_INPUT_X_2	N_x, BlockSize_x, Left_x, Right_x, Map_F_x, Map_x, Manual_x, alpha_x, beta_x
#define _PS_INPUT_X_3	Filter_Choice_D_x, Mode_CutOff_D_x, Omega_D_x, Order_D_x, Smooth_1_x, Smooth_2_x
#define _PS_INPUT_X_4	Filter_Choice_S_x, Mode_CutOff_S_x, Omega_S_x, Order_S_x
#define _PS_INPUT_X_5	Filter_Choice_R_x, Mode_CutOff_R_x, Omega_R_x, Order_R_x

#define _PS_INPUT_Y_1	Index_y, Method_y, Point_Type_y, Max_Order_Of_D_y, Algorithm_y, Symmetry_y
#define _PS_INPUT_Y_2	N_y, BlockSize_y, Left_y, Right_y, Map_F_y, Map_y, Manual_y, alpha_y, beta_y
#define _PS_INPUT_Y_3	Filter_Choice_D_y, Mode_CutOff_D_y, Omega_D_y, Order_D_y, Smooth_1_y, Smooth_2_y
#define _PS_INPUT_Y_4	Filter_Choice_S_y, Mode_CutOff_S_y, Omega_S_y, Order_S_y
#define _PS_INPUT_Y_5	Filter_Choice_R_y, Mode_CutOff_R_y, Omega_R_y, Order_R_y

#define _PS_INPUT_Z_1	Index_z, Method_z, Point_Type_z, Max_Order_Of_D_z, Algorithm_z, Symmetry_z
#define _PS_INPUT_Z_2	N_z, BlockSize_z, Left_z, Right_z, Map_F_z, Map_z, Manual_z, alpha_z, beta_z
#define _PS_INPUT_Z_3	Filter_Choice_D_z, Mode_CutOff_D_z, Omega_D_z, Order_D_z, Smooth_1_z, Smooth_2_z
#define _PS_INPUT_Z_4	Filter_Choice_S_z, Mode_CutOff_S_z, Omega_S_z, Order_S_z
#define _PS_INPUT_Z_5	Filter_Choice_R_z, Mode_CutOff_R_z, Omega_R_z, Order_R_z

#define _PS_INPUT_X_1_I	Index_x,Method_x,Point_Type_x,Max_Order_Of_D_x, Algorithm_x,Symmetry_x
#define _PS_INPUT_X_2_I	N_x,BlockSize_x,Map_x,Manual_x,Map_F_x
#define _PS_INPUT_X_2_R	Left_x, Right_x, alpha_x, beta_x, Angle_x
#define _PS_INPUT_X_3_I	Filter_Choice_D_x,Mode_CutOff_D_x,Smooth_1_x,Smooth_2_x
#define _PS_INPUT_X_3_R	Omega_D_x, Order_D_x
#define _PS_INPUT_X_4_I	Filter_Choice_S_x,Mode_CutOff_S_x
#define _PS_INPUT_X_4_R	Omega_S_x, Order_S_x
#define _PS_INPUT_X_5_I	Filter_Choice_R_x,Mode_CutOff_R_x
#define _PS_INPUT_X_5_R	Omega_R_x, Order_R_x

#define _PS_INPUT_Y_1_I	Index_y,Method_y,Point_Type_y,Max_Order_Of_D_y, Algorithm_y,Symmetry_y
#define _PS_INPUT_Y_2_I	N_y,BlockSize_y,Map_y,Manual_y,Map_F_y
#define _PS_INPUT_Y_2_R	Left_y, Right_y, alpha_y, beta_y, Angle_y
#define _PS_INPUT_Y_3_I	Filter_Choice_D_y,Mode_CutOff_D_y,Smooth_1_y,Smooth_2_y
#define _PS_INPUT_Y_3_R	Omega_D_y, Order_D_y
#define _PS_INPUT_Y_4_I	Filter_Choice_S_y,Mode_CutOff_S_y
#define _PS_INPUT_Y_4_R	Omega_S_y, Order_S_y
#define _PS_INPUT_Y_5_I	Filter_Choice_R_y,Mode_CutOff_R_y
#define _PS_INPUT_Y_5_R	Omega_R_y, Order_R_y

#define _PS_INPUT_Z_1_I	Index_z,Method_z,Point_Type_z,Max_Order_Of_D_z, Algorithm_z,Symmetry_z
#define _PS_INPUT_Z_2_I	N_z,BlockSize_z,Map_z,Manual_z,Map_F_z
#define _PS_INPUT_Z_2_R	Left_z, Right_z, alpha_z, beta_z, Angle_z
#define _PS_INPUT_Z_3_I	Filter_Choice_D_z,Mode_CutOff_D_z,Smooth_1_z,Smooth_2_z
#define _PS_INPUT_Z_3_R	Omega_D_z, Order_D_z
#define _PS_INPUT_Z_4_I	Filter_Choice_S_z,Mode_CutOff_S_z
#define _PS_INPUT_Z_4_R	Omega_S_z, Order_S_z
#define _PS_INPUT_Z_5_I	Filter_Choice_R_z,Mode_CutOff_R_z
#define _PS_INPUT_Z_5_R	Omega_R_z, Order_R_z

#define _PS_WENO_GRID_X_0	Index_x,              Method_x,Symmetry_x, WENO_Order_x
#define _PS_WENO_GRID_X_1	Index_x,Distributed_x,Method_x,Symmetry_x, WENO_Order_x
#define _PS_WENO_GRID_X_2       NN_g,N_g,N0_g,N1_g,N2_g,N3_g,N4_g,N5_g
#define _PS_WENO_GRID_X_3	NN  ,N  ,N0  ,N1  ,N2  ,N3  ,N4  ,N5  
#define _PS_WENO_GRID_X_4	n_Remain, n_0, n_Star 
#define _PS_WENO_GRID_X_X	_PS_WENO_GRID_X_2, _PS_WENO_GRID_X_3

#define _PS_WENO_GRID_Y_0	Index_y,              Method_y,Symmetry_y, WENO_Order_y
#define _PS_WENO_GRID_Y_1	Index_y,Distributed_y,Method_y,Symmetry_y, WENO_Order_y
#define _PS_WENO_GRID_Y_2       MM_g,M_g,M0_g,M1_g,M2_g,M3_g,M4_g,M5_g
#define _PS_WENO_GRID_Y_3	MM  ,M  ,M0  ,M1  ,M2  ,M3  ,M4  ,M5  
#define _PS_WENO_GRID_Y_4	m_Remain, m_0, m_Star
#define _PS_WENO_GRID_Y_Y	_PS_WENO_GRID_Y_2, _PS_WENO_GRID_Y_3

#define _PS_WENO_GRID_Z_0	Index_z,              Method_z,Symmetry_z, WENO_Order_z
#define _PS_WENO_GRID_Z_1	Index_z,Distributed_z,Method_z,Symmetry_z, WENO_Order_z
#define _PS_WENO_GRID_Z_2       KK_g,K_g,K0_g,K1_g,K2_g,K3_g,K4_g,K5_g
#define _PS_WENO_GRID_Z_3	KK  ,K  ,K0  ,K1  ,K2  ,K3  ,K4  ,K5  
#define _PS_WENO_GRID_Z_4	k_Remain, k_0, k_Star
#define _PS_WENO_GRID_Z_Z	_PS_WENO_GRID_Z_2, _PS_WENO_GRID_Z_3

#define _PS_SPECTRAL_GRID_X_0	Index_x,              Method_x,Symmetry_x, Point_Type_x
#define _PS_SPECTRAL_GRID_X_1	Index_x,Distributed_x,Method_x,Symmetry_x, Point_Type_x
#define _PS_SPECTRAL_GRID_X_2	NN_g,N_g,N0_g,N1_g,N2_g,N3_g,N4_g,N5_g
#define _PS_SPECTRAL_GRID_X_3	NN  ,N  ,N0  ,N1  ,N2  ,N3  ,N4  ,N5   
#define _PS_SPECTRAL_GRID_X_4	n_Remain,n_0,n_Star
#define _PS_SPECTRAL_GRID_X_X	_PS_SPECTRAL_GRID_X_2,_PS_SPECTRAL_GRID_X_3, _PS_SPECTRAL_GRID_X_4

#define _PS_SPECTRAL_GRID_Y_0	Index_y,              Method_y,Symmetry_y, Point_Type_y
#define _PS_SPECTRAL_GRID_Y_1	Index_y,Distributed_y,Method_y,Symmetry_y, Point_Type_y
#define _PS_SPECTRAL_GRID_Y_2	MM_g,M_g,M0_g,M1_g,M2_g,M3_g,M4_g,M5_g
#define _PS_SPECTRAL_GRID_Y_3	MM  ,M  ,M0  ,M1  ,M2  ,M3  ,M4  ,M5
#define _PS_SPECTRAL_GRID_Y_4	m_Remain,m_0,m_Star
#define _PS_SPECTRAL_GRID_Y_Y	_PS_SPECTRAL_GRID_Y_2,_PS_SPECTRAL_GRID_Y_3, _PS_SPECTRAL_GRID_Y_4

#define _PS_SPECTRAL_GRID_Z_0	Index_z,              Method_z,Symmetry_z, Point_Type_z
#define _PS_SPECTRAL_GRID_Z_1	Index_z,Distributed_z,Method_z,Symmetry_z, Point_Type_z
#define _PS_SPECTRAL_GRID_Z_2	KK_g,K_g,K0_g,K1_g,K2_g,K3_g,K4_g,K5_g
#define _PS_SPECTRAL_GRID_Z_3	KK  ,K  ,K0  ,K1  ,K2  ,K3  ,K4  ,K5
#define _PS_SPECTRAL_GRID_Z_4	k_Remain,k_0,k_Star
#define _PS_SPECTRAL_GRID_Z_Z	_PS_SPECTRAL_GRID_Z_2,_PS_SPECTRAL_GRID_Z_3, _PS_SPECTRAL_GRID_Z_4
