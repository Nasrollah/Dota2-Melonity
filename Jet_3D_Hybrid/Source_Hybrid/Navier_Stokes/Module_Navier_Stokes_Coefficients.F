#include "pseudopack.h"


MODULE Navier_Stokes_Coefficients

implicit NONE

INTERFACE PS_Grid_Index_Update_NS
  MODULE PROCEDURE Grid_Index_Update_NS_1D 
  MODULE PROCEDURE Grid_Index_Update_NS_2D 
  MODULE PROCEDURE Grid_Index_Update_NS_3D 
END INTERFACE

Logical  :: Navier_Stokes         = .TRUE.

Logical  ::   Thermal_Conduction  = .TRUE.
Logical  ::   Viscous_Dissipation = .TRUE.
Logical  :: Molecular_Diffusion   = .TRUE.

integer  :: NS_TC_Order_x, NS_VD_Order_x, NS_MD_Order_x, NS_Spatial_Order_x
integer  :: NS_TC_Order_y, NS_VD_Order_y, NS_MD_Order_y, NS_Spatial_Order_y
integer  :: NS_TC_Order_z, NS_VD_Order_z, NS_MD_Order_z, NS_Spatial_Order_z

integer  :: NS_N0, NS_N2, NS_N3, NS_N5, NS_M_s_x, N_Ghost_x
integer  :: NS_M0, NS_M2, NS_M3, NS_M5, NS_M_s_y, N_Ghost_y
integer  :: NS_K0, NS_K2, NS_K3, NS_K5, NS_M_s_z, N_Ghost_z

REALTYPE, dimension(-20:20,-1:4) :: CD_TC_x, CD_TC_y, CD_TC_z
REALTYPE, dimension(-20:20,-1:4) :: CD_VD_x, CD_VD_y, CD_VD_z
REALTYPE, dimension(-20:20,-1:4) :: CD_MD_x, CD_MD_y, CD_MD_z

PUBLIC

PUBLIC  :: Navier_Stokes

PUBLIC  ::   Thermal_Conduction 
PUBLIC  ::   Viscous_Dissipation
PUBLIC  :: Molecular_Diffusion

PUBLIC  :: PS_Grid_Index_Update_NS

CONTAINS
!
!===================================================================
!
  Subroutine Grid_Index_Update_NS_1D (N0, N2, N3, N5)

  integer  :: N0, N2, N3, N5

  NS_M_s_x = (MAX(NS_TC_Order_x, NS_VD_Order_x, NS_MD_Order_x)/2)*2

  N_Ghost_x = MAX(N2-N0, N5-N3, NS_M_s_x) 

  END Subroutine Grid_Index_Update_NS_1D 
!
!===================================================================
!
  Subroutine Grid_Index_Update_NS_2D (N0, N2, N3, N5,   &
                                      M0, M2, M3, M5)

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5

  NS_M_s_x = (MAX(NS_TC_Order_x, NS_VD_Order_x, NS_MD_Order_x)/2)*2
  NS_M_s_y = (MAX(NS_TC_Order_y, NS_VD_Order_y, NS_MD_Order_y)/2)*2

  N_Ghost_x = MAX(N2-N0, N5-N3, NS_M_s_x) 
  N_Ghost_y = MAX(M2-M0, M5-M3, NS_M_s_y) 

  END Subroutine Grid_Index_Update_NS_2D 
!
!===================================================================
!
  Subroutine Grid_Index_Update_NS_3D (N0, N2, N3, N5,   &
                                      M0, M2, M3, M5,   &
                                      K0, K2, K3, K5)

  integer  :: N0, N2, N3, N5
  integer  :: M0, M2, M3, M5
  integer  :: K0, K2, K3, K5

  NS_M_s_x = (MAX(NS_TC_Order_x, NS_VD_Order_x, NS_MD_Order_x)/2)*2
  NS_M_s_y = (MAX(NS_TC_Order_y, NS_VD_Order_y, NS_MD_Order_y)/2)*2
  NS_M_s_z = (MAX(NS_TC_Order_z, NS_VD_Order_z, NS_MD_Order_z)/2)*2

  N_Ghost_x = MAX(N2-N0, N5-N3, NS_M_s_x) 
  N_Ghost_y = MAX(M2-M0, M5-M3, NS_M_s_y) 
  N_Ghost_z = MAX(K2-K0, K5-K3, NS_M_s_z) 

  END Subroutine Grid_Index_Update_NS_3D 
  
END MODULE Navier_Stokes_Coefficients
