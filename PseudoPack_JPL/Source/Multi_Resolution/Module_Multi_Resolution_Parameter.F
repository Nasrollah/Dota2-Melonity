#include "pseudopack.h"


Module Multi_Resolution_Parameter

IMPLICIT NONE  

integer  :: MR_Order_x, MR_Index_x, MR_N0, MR_N2, MR_N3, MR_N5
integer  :: MR_Order_y, MR_Index_y, MR_M0, MR_M2, MR_M3, MR_M5
integer  :: MR_Order_z, MR_Index_z, MR_K0, MR_K2, MR_K3, MR_K5

REALTYPE, dimension(:,:,:,:), ALLOCATABLE :: MR_Q
REALTYPE, dimension(:,:,:,:), ALLOCATABLE :: MR_Scale
integer , dimension(:,:,:,:), ALLOCATABLE :: MR_Flag

PRIVATE

PUBLIC  :: MR_Order_x, MR_Index_x, MR_N0, MR_N2, MR_N3, MR_N5
PUBLIC  :: MR_Order_y, MR_Index_y, MR_M0, MR_M2, MR_M3, MR_M5
PUBLIC  :: MR_Order_z, MR_Index_z, MR_K0, MR_K2, MR_K3, MR_K5

PUBLIC  :: MR_Q, MR_Scale, MR_Flag

PUBLIC  :: PS_MR_Allocate_Flag_And_Scale
PUBLIC  :: PS_MR_DeAllocate_Flag_And_Scale

CONTAINS

  Subroutine PS_MR_Allocate_Flag_And_Scale (Dimensions, Flag, Scale, Q)

  integer            :: Dimensions
  logical , OPTIONAL ::  Flag,  Scale,  Q
  logical            :: IFlag, IScale, IQ

  IQ     = .TRUE. ; IF (PRESENT(Q    )) IQ     = Q
  IFlag  = .TRUE. ; IF (PRESENT(Flag )) IFlag  = Flag
  IScale = .TRUE. ; IF (PRESENT(Scale)) IScale = Scale

  call PS_MR_DeAllocate_Flag_And_Scale (Dimensions, Flag, Scale, Q)

  SELECT CASE (Dimensions)
    CASE DEFAULT
      if (IQ    ) Allocate (MR_Q    (MR_N0:MR_N5,    1:1   ,     1:1    ,1:1))
      if (IFlag ) Allocate (MR_Flag (MR_N0:MR_N5,    1:1   ,     1:1    ,1:1))
      if (IScale) Allocate (MR_Scale(MR_N0:MR_N5,    1:1   ,     1:1    ,1:1))

    CASE (2)
      if (IQ    ) Allocate (MR_Q    (MR_N0:MR_N5,MR_M0:MR_M5,    1:1    ,1:1))
      if (IFlag ) Allocate (MR_Flag (MR_N0:MR_N5,MR_M0:MR_M5,    1:1    ,1:2))
      if (IScale) Allocate (MR_Scale(MR_N0:MR_N5,MR_M0:MR_M5,    1:1    ,1:2))

    CASE (3)
      if (IQ    ) Allocate (MR_Q    (MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,1:1))
      if (IFlag ) Allocate (MR_Flag (MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,1:3))
      if (IScale) Allocate (MR_Scale(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,1:3))

  END SELECT

      if (IQ    ) MR_Q     = ZERO
      if (IFlag ) MR_Flag  = 1
      if (IScale) MR_Scale = ZERO
  
  END Subroutine PS_MR_Allocate_Flag_And_Scale 
! 
! ==============================================================
!
  Subroutine PS_MR_DeAllocate_Flag_And_Scale (Dimensions, Flag, Scale, Q)

  integer , OPTIONAL :: Dimensions
  logical , OPTIONAL ::  Flag,  Scale,  Q
  logical            :: IFlag, IScale, IQ

  IQ     = .TRUE. ; IF (PRESENT(Q    )) IQ     = Q
  IFlag  = .TRUE. ; IF (PRESENT(Flag )) IFlag  = Flag
  IScale = .TRUE. ; IF (PRESENT(Scale)) IScale = Scale

  if (IQ     .AND. ALLOCATED(MR_Q    )) DeAllocate (MR_Q    )
  if (IFlag  .AND. ALLOCATED(MR_Flag )) DeAllocate (MR_Flag )
  if (ISCale .AND. ALLOCATED(MR_Scale)) DeAllocate (MR_Scale)

  END Subroutine PS_MR_DeAllocate_Flag_And_Scale 

END Module Multi_Resolution_Parameter
