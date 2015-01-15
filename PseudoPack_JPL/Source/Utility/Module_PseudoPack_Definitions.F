#include "pseudopack.h"


Module PseudoPack_Definition

  USE PseudoPack_Input
  USE MPI_Processor_Topology

IMPLICIT NONE

#if 0
integer  ::  Index_x, Method_x, Point_Type_x, Max_Order_Of_D_x, &
               Algorithm_x, Symmetry_x,                         &
             N_x, BlockSize_x, Map_x, Manual_x, Map_F_x,        &
             Filter_Choice_D_x, Mode_CutOff_D_x,                &
               Smooth_1_x, Smooth_2_x,                          &
             Filter_Choice_S_x, Mode_CutOff_S_x,                &
             Filter_Choice_R_x, Mode_CutOff_R_x

REALTYPE ::  Left_x, Right_x, alpha_x, beta_x,                  &
             Omega_D_x, Order_D_x, Omega_S_x, Order_S_x,        &
             Omega_R_x, Order_R_x

REALTYPE :: x0, x1

integer  :: NN_g, N_g, N0_g, N1_g, N2_g, N3_g, N4_g, N5_g
integer  :: NN  , N  , N0  , N1  , N2  , N3  , N4  , N5   
integer  :: n_Remain, n_0, n_Star

integer  ::     WENO_Order_x, N_Ghost_Cell_x       , M_s_x
integer  :: Central_Degree_x, N_Collocation_Point_x, O_s_x

logical  :: Distributed_x, Periodicity_x
integer  :: Topology_Grid_x

#else

integer  ::           Index_x = 1,       &
                     Method_x = 0,       &
	         Point_Type_x = 1,       &
	     Max_Order_Of_D_x = 1,       &
                  Algorithm_x = 2,       &
	           Symmetry_x = 0

integer  ::               N_x = 32,      &  
	          BlockSize_x = 0 

REALTYPE ::            Left_x = ZERO  ,  &
                      Right_x = TWO*PI
integer  ::           Map_F_x = 0
		      
integer  ::             Map_x = 0,       & 
	             Manual_x = 0
REALTYPE ::           alpha_x = 0.99999, &
                       beta_x = ZERO

integer  :: Filter_Choice_D_x = 0,       &
	      Mode_CutOff_D_x = 0,       &
                   Smooth_1_x = 1,       &
	           Smooth_2_x = 1
REALTYPE ::         Omega_D_x = 0,       &
                    Order_D_x = FOUR*FOUR

integer  :: Filter_Choice_S_x = 0,       &
	      Mode_CutOff_S_x = 0
REALTYPE ::         Omega_S_x = 0,       &
                    Order_S_x = FOUR*FOUR

integer  :: Filter_Choice_R_x = 0,       & 
	      Mode_CutOff_R_x = 0 
REALTYPE ::         Omega_R_x = 0,       &
                    Order_R_x = FOUR*FOUR

integer  ::          WENO_Order_x = 5,   &  
                   N_Ghost_Cell_x = 3,   &
	                    M_s_x = 3

integer  ::      Central_Degree_x = 5,   & 
            N_Collocation_Point_x = 7,   &
	                    O_s_x = 7

logical  ::   Distributed_x = .TRUE.,    &
              Periodicity_x = .TRUE.  
integer  :: Topology_Grid_x = 0

REALTYPE :: x0 = ZERO,                   &
            x1 = TWO*PI

integer  :: NN_g=32, N_g=32, N0_g=0, N1_g=0, N2_g=0, N3_g=31, N4_g=31, N5_g=31
integer  :: NN  =32, N  =32, N0  =0, N1  =0, N2  =0, N3  =31, N4  =31, N5  =31
integer  :: n_Remain=0, n_0=0, n_Star=32
#endif

REALTYPE, dimension(:), ALLOCATABLE, TARGET :: xi, D_x, S_x, R_x

TYPE (PS_Property)  , TARGET ::   Property_x
TYPE (PS_Grid_Index), TARGET :: Grid_Index_x
TYPE (PS_Domain)    , TARGET ::     Domain_x
TYPE (PS_Mapping)   , TARGET ::    Mapping_x
TYPE (PS_Filtering) , TARGET ::   Filter_D_x, Filter_S_x, Filter_R_x
!
! =====================================================================
!
#if 0
integer  ::  Index_y, Method_y, Point_Type_y, Max_Order_Of_D_y, &
               Algorithm_y, Symmetry_y,                         &
             N_y, BlockSize_y, Map_y, Manual_y, Map_F_y,        &
             Filter_Choice_D_y, Mode_CutOff_D_y,                &
               Smooth_1_y, Smooth_2_y,                          &
             Filter_Choice_S_y, Mode_CutOff_S_y,                &
             Filter_Choice_R_y, Mode_CutOff_R_y

REALTYPE ::  Left_y, Right_y, alpha_y, beta_y,                  &
             Omega_D_y, Order_D_y, Omega_S_y, Order_S_y,        &
             Omega_R_y, Order_R_y

REALTYPE :: y0, y1

integer  :: MM_g, M_g, M0_g, M1_g, M2_g, M3_g, M4_g, M5_g
integer  :: MM  , M  , M0  , M1  , M2  , M3  , M4  , M5  
integer  :: m_Remain, m_0, m_Star

integer  ::     WENO_Order_y, N_Ghost_Cell_y       , M_s_y
integer  :: Central_Degree_y, N_Collocation_Point_y, O_s_y

logical  :: Distributed_y, Periodicity_y
integer  :: Topology_Grid_y

#else

integer  ::           Index_y = 2,       &
                     Method_y = 0,       &
	         Point_Type_y = 1,       &
	     Max_Order_Of_D_y = 1,       &
                  Algorithm_y = 2,       &
	           Symmetry_y = 0

integer  ::               N_y = 32,      &  
	          BlockSize_y = 0 

REALTYPE ::            Left_y = ZERO  ,  &
                      Right_y = TWO*PI
integer  ::           Map_F_y = 0
		      
integer  ::             Map_y = 0,       & 
	             Manual_y = 0
REALTYPE ::           alpha_y = 0.99999, &
                       beta_y = ZERO

integer  :: Filter_Choice_D_y = 0,       &
	      Mode_CutOff_D_y = 0,       &
                   Smooth_1_y = 1,       &
	           Smooth_2_y = 1
REALTYPE ::         Omega_D_y = 0,       &
                    Order_D_y = FOUR*FOUR

integer  :: Filter_Choice_S_y = 0,       &
	      Mode_CutOff_S_y = 0
REALTYPE ::         Omega_S_y = 0,       &
                    Order_S_y = FOUR*FOUR

integer  :: Filter_Choice_R_y = 0,       & 
	      Mode_CutOff_R_y = 0 
REALTYPE ::         Omega_R_y = 0,       &
                    Order_R_y = FOUR*FOUR

integer  ::          WENO_Order_y = 5,   &  
                   N_Ghost_Cell_y = 3,   &
	                    M_s_y = 3

integer  ::      Central_Degree_y = 5,   & 
            N_Collocation_Point_y = 7,   &
	                    O_s_y = 7

logical  ::   Distributed_y = .TRUE.,    &
              Periodicity_y = .TRUE.  
integer  :: Topology_Grid_y = 0

REALTYPE :: y0 = ZERO,                   &
            y1 = TWO*PI

integer  :: MM_g=32, M_g=32, M0_g=0, M1_g=0, M2_g=0, M3_g=31, M4_g=31, M5_g=31
integer  :: MM  =32, M  =32, M0  =0, M1  =0, M2  =0, M3  =31, M4  =31, M5  =31
integer  :: m_Remain=0, m_0=0, m_Star=32
#endif

REALTYPE, dimension(:), ALLOCATABLE, TARGET :: eta, D_y, S_y, R_y

TYPE (PS_Property)  , TARGET ::   Property_y
TYPE (PS_Grid_Index), TARGET :: Grid_Index_y
TYPE (PS_Domain)    , TARGET ::     Domain_y
TYPE (PS_Mapping)   , TARGET ::    Mapping_y
TYPE (PS_Filtering) , TARGET ::   Filter_D_y, Filter_S_y, Filter_R_y
!
! =====================================================================
!
#if 0
integer  ::  Index_z, Method_z, Point_Type_z, Max_Order_Of_D_z, &
               Algorithm_z, Symmetry_z,                         &
             N_z, BlockSize_z, Map_z, Manual_z, Map_F_z,        &
             Filter_Choice_D_z, Mode_CutOff_D_z,                &
               Smooth_1_z, Smooth_2_z,                          &
             Filter_Choice_S_z, Mode_CutOff_S_z,                &
             Filter_Choice_R_z, Mode_CutOff_R_z

REALTYPE ::  Left_z, Right_z, alpha_z, beta_z,                  &
             Omega_D_z, Order_D_z, Omega_S_z, Order_S_z,        &
             Omega_R_z, Order_R_z

REALTYPE :: z0, z1

integer  :: KK_g, K_g, K0_g, K1_g, K2_g, K3_g, K4_g, K5_g
integer  :: KK  , K  , K0  , K1  , K2  , K3  , K4  , K5  
integer  :: k_Remain, k_0, k_Star

integer  ::     WENO_Order_z, N_Ghost_Cell_z       , M_s_z
integer  :: Central_Degree_z, N_Collocation_Point_z, O_s_z

logical  :: Distributed_z, Periodicity_z
integer  :: Topology_Grid_z

#else

integer  ::           Index_z = 3,       &
                     Method_z = 0,       &
	         Point_Type_z = 1,       &
	     Max_Order_Of_D_z = 1,       &
                  Algorithm_z = 2,       &
	           Symmetry_z = 0

integer  ::               N_z = 32,      &  
	          BlockSize_z = 0 

REALTYPE ::            Left_z = ZERO  ,  &
                      Right_z = TWO*PI
integer  ::           Map_F_z = 0
		      
integer  ::             Map_z = 0,       & 
	             Manual_z = 0
REALTYPE ::           alpha_z = 0.99999, &
                       beta_z = ZERO

integer  :: Filter_Choice_D_z = 0,       &
	      Mode_CutOff_D_z = 0,       &
                   Smooth_1_z = 1,       &
	           Smooth_2_z = 1
REALTYPE ::         Omega_D_z = 0,       &
                    Order_D_z = FOUR*FOUR

integer  :: Filter_Choice_S_z = 0,       &
	      Mode_CutOff_S_z = 0
REALTYPE ::         Omega_S_z = 0,       &
                    Order_S_z = FOUR*FOUR

integer  :: Filter_Choice_R_z = 0,       & 
	      Mode_CutOff_R_z = 0 
REALTYPE ::         Omega_R_z = 0,       &
                    Order_R_z = FOUR*FOUR

integer  ::          WENO_Order_z = 5,   &  
                   N_Ghost_Cell_z = 3,   &
	                    M_s_z = 3

integer  ::      Central_Degree_z = 5,   & 
            N_Collocation_Point_z = 7,   &
	                    O_s_z = 7

logical  ::   Distributed_z = .TRUE.,    &
              Periodicity_z = .TRUE.  
integer  :: Topology_Grid_z = 0

REALTYPE :: z0 = ZERO,                   &
            z1 = TWO*PI

integer  :: KK_g=32, K_g=32, K0_g=0, K1_g=0, K2_g=0, K3_g=31, K4_g=31, K5_g=31
integer  :: KK  =32, K  =32, K0  =0, K1  =0, K2  =0, K3  =31, K4  =31, K5  =31
integer  :: k_Remain=0, k_0=0, k_Star=32
#endif

REALTYPE, dimension(:), ALLOCATABLE, TARGET ::  fi, D_z, S_z, R_z

TYPE (PS_Property)  , TARGET ::   Property_z
TYPE (PS_Grid_Index), TARGET :: Grid_Index_z
TYPE (PS_Domain)    , TARGET ::     Domain_z
TYPE (PS_Mapping)   , TARGET ::    Mapping_z
TYPE (PS_Filtering) , TARGET ::   Filter_D_z, Filter_S_z, Filter_R_z

!
! =====================================================================
!
TYPE (Processor_Topology) :: Topology

!
! =====================================================================
!
TYPE PS_Description

  integer  ::  Index, Method, Point_Type, Max_Order_Of_D,     &
                 Algorithm, Symmetry,                         &
               N, BlockSize, Map, Manual, Map_F,              &
               Filter_Choice_D, Mode_CutOff_D,                &
                 Smooth_1, Smooth_2,                          &
               Filter_Choice_S, Mode_CutOff_S,                &
               Filter_Choice_R, Mode_CutOff_R

  REALTYPE ::  Left, Right, alpha, beta,                      &
               Omega_D, Order_D, Omega_S, Order_S,            &
               Omega_R, Order_R

  integer  ::    WENO_Order , N_Ghost_Cell       , M_s
  integer  :: Central_Degree, N_Collocation_Point, O_s

  integer  :: NN_g, N_g, N0_g, N1_g, N2_g, N3_g, N4_g, N5_g
  integer  :: NN       , N0  , N1  , N2  , N3  , N4  , N5
  integer  :: n_Remain, n_0, n_Star

  logical  :: Distributed, Periodicity
  integer  :: Topology_Grid

END TYPE PS_Description

TYPE (PS_Description), dimension(:), ALLOCATABLE :: PS_Descriptor

PUBLIC 

PUBLIC  :: PS_Description
PUBLIC  :: PS_Descriptor

PUBLIC  :: PS_Input_x
PUBLIC  :: PS_Input_y
PUBLIC  :: PS_Input_z
PUBLIC  :: PS_Input_MPI_Topology

CONTAINS

  Subroutine PS_Input_x (lid)

  integer  :: lid

  call PS_Input (Index_x, Method_x, Point_Type_x, Max_Order_Of_D_x,  &
                   Algorithm_x, Symmetry_x,                          &
                 N_x, BlockSize_x, Left_x, Right_x, Map_F_x,         &
                   Map_x, Manual_x, alpha_x, beta_x,                 &
                 Filter_Choice_D_x, Mode_CutOff_D_x,                 &
                   Omega_D_x, Order_D_x, Smooth_1_x, Smooth_2_x,     &
                 Filter_Choice_S_x, Mode_CutOff_S_x,                 &
                   Omega_S_x, Order_S_x,                             &
                 Filter_Choice_R_x, Mode_CutOff_R_x,                 &
                   Omega_R_x, Order_R_x,                             &
                 lid)

  read (lid,'(//)')
  read (lid,*) WENO_Order_x, Central_Degree_x
  read (lid,'(//)')
  read (lid,*) Distributed_x, Periodicity_x, Topology_Grid_x
  read (lid,*)

  END Subroutine PS_Input_x 
!
! =========================================================================
!
  Subroutine PS_Input_y (lid)

  integer  :: lid

  call PS_Input (Index_y, Method_y, Point_Type_y, Max_Order_Of_D_y,  &
                   Algorithm_y, Symmetry_y,                          &
                 N_y, BlockSize_y, Left_y, Right_y, Map_F_y,         &
                   Map_y, Manual_y, alpha_y, beta_y,                 &
                 Filter_Choice_D_y, Mode_CutOff_D_y,                 &
                   Omega_D_y, Order_D_y, Smooth_1_y, Smooth_2_y,     &
                 Filter_Choice_S_y, Mode_CutOff_S_y,                 &
                   Omega_S_y, Order_S_y,                             &
                 Filter_Choice_R_y, Mode_CutOff_R_y,                 &
                   Omega_R_y, Order_R_y,                             &
                 lid)

  read (lid,'(//)')
  read (lid,*) WENO_Order_y, Central_Degree_y
  read (lid,'(//)')
  read (lid,*) Distributed_y, Periodicity_y, Topology_Grid_y
  read (lid,*)

  END Subroutine PS_Input_y 
!
! =========================================================================
!
  Subroutine PS_Input_z (lid)

  integer  :: lid

  call PS_Input (Index_z, Method_z, Point_Type_z, Max_Order_Of_D_z,  &
                   Algorithm_z, Symmetry_z,                          &
                 N_z, BlockSize_z, Left_z, Right_z, Map_F_z,         &
                   Map_z, Manual_z, alpha_z, beta_z,                 &
                 Filter_Choice_D_z, Mode_CutOff_D_z,                 &
                   Omega_D_z, Order_D_z, Smooth_1_z, Smooth_2_z,     &
                 Filter_Choice_S_z, Mode_CutOff_S_z,                 &
                   Omega_S_z, Order_S_z,                             &
                 Filter_Choice_R_z, Mode_CutOff_R_z,                 &
                   Omega_R_z, Order_R_z,                             &
                 lid)

  read (lid,'(//)')
  read (lid,*) WENO_Order_z, Central_Degree_z
  read (lid,'(//)')
  read (lid,*) Distributed_z, Periodicity_z, Topology_Grid_z
  read (lid,*)

  END Subroutine PS_Input_z 
!
! =========================================================================
!
  Subroutine PS_Input_MPI_Topology (lid)

  integer  :: lid

  read (lid,'(//)')
  read (lid,*) Distributed_x, Periodicity_x, Topology_Grid_x
  read (lid,*) Distributed_y, Periodicity_y, Topology_Grid_y
  read (lid,*) Distributed_z, Periodicity_z, Topology_Grid_z

  END Subroutine PS_Input_MPI_Topology 
!
! =========================================================================
!
  Subroutine PS_Input_2 (PS, lid)

  TYPE (PS_Description) :: PS

  integer  :: lid

  call PS_Input (PS%Index, PS%Method, PS%Point_Type, PS%Max_Order_Of_D, &
                   PS%Algorithm, PS%Symmetry,                           &
                 PS%N, PS%BlockSize, PS%Left, PS%Right, PS%Map_F,       &
                   PS%Map, PS%Manual, PS%alpha, PS%beta,                &
                 PS%Filter_Choice_D, PS%Mode_CutOff_D,                  &
                   PS%Omega_D, PS%Order_D, PS%Smooth_1, PS%Smooth_2,    &
                 PS%Filter_Choice_S, PS%Mode_CutOff_S,                  &
                   PS%Omega_S, PS%Order_S,                              &
                 PS%Filter_Choice_R, PS%Mode_CutOff_R,                  &
                   PS%Omega_R, PS%Order_R,                              &
                 lid)

  read (lid,'(//)')
  read (lid,*) PS%WENO_Order, PS%Central_Degree
  read (lid,'(//)')
  read (lid,*) PS%Distributed, PS%Periodicity, PS%Topology_Grid
  read (lid,*)

  END Subroutine PS_Input_2 

END Module PseudoPack_Definition
