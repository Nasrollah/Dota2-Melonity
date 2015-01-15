#include "pseudopack.h"


MODULE WENO_Option

implicit NONE

REALTYPE, parameter ::         Epsilon_Default = 1.0d-10
REALTYPE, parameter ::   Lambda_Factor_Default = ELEVEN/TEN
REALTYPE, parameter :: Sound_Speed_Min_Default = TEN*EPSILON(ONE)
REALTYPE, parameter ::           Power_Default = ONE
integer , parameter ::     Unroll_Loop_Default = 2                 ! 2
integer , parameter ::       Flux_Type_Default = 1
integer , parameter ::     Weight_Type_Default = 0

#if 0
REALTYPE, dimension(:)    , ALLOCATABLE :: d, W
REALTYPE, dimension(:,:)  , ALLOCATABLE :: C
REALTYPE, dimension(:,:,:), ALLOCATABLE :: S
#endif

TYPE PS_WENO_Option_Type
  integer  :: Order
  REALTYPE :: Epsilon
  REALTYPE :: Power
  REALTYPE :: Lambda_Factor
  REALTYPE :: Sound_Speed_Min
  integer  :: Unroll_Loop
  integer  :: Flux_Type
  integer  :: Weight_Type
END TYPE

TYPE (PS_WENO_Option_Type) :: PS_WENO_Option

PRIVATE
PUBLIC :: PS_WENO_Option_Type
PUBLIC :: PS_WENO_Option_Setup
PUBLIC :: PS_WENO_Option
PUBLIC :: PS_WENO_Mapped_Weight

CONTAINS

  Subroutine PS_WENO_Option_Setup (Order, Epsilon, Power, Lambda_Factor,           &
                                   Sound_Speed_Min, Unroll_Loop, Flux_Type, &
                                   Weight_Type)

  integer            :: Order, M
  REALTYPE, OPTIONAL :: Epsilon, Power, Lambda_Factor, Sound_Speed_Min
  integer , OPTIONAL :: Unroll_Loop, Flux_Type, Weight_Type

    PS_WENO_Option%Epsilon             = Epsilon_Default 
    PS_WENO_Option%Power               = Power_Default
    PS_WENO_Option%Lambda_Factor       = Lambda_Factor_Default 
    PS_WENO_Option%Sound_Speed_Min     = Sound_Speed_Min_Default
    PS_WENO_Option%Unroll_Loop         = Unroll_Loop_Default
    PS_WENO_Option%Flux_Type           = Flux_Type_Default
    PS_WENO_Option%Weight_Type         = Weight_Type_Default

  if (PRESENT(Epsilon))       PS_WENO_Option%Epsilon       = Epsilon
  if (PRESENT(Power))         PS_WENO_Option%Power         = Power
  if (PRESENT(Lambda_Factor)) PS_WENO_Option%Lambda_Factor = Lambda_Factor
  if (PRESENT(Unroll_Loop))   PS_WENO_Option%Unroll_Loop   = Unroll_Loop
  if (PRESENT(Flux_Type))     PS_WENO_Option%Flux_Type     = Flux_Type
  if (PRESENT(Weight_Type))   PS_WENO_Option%Weight_Type   = Weight_type

  if (PRESENT(Sound_Speed_Min))                              &
    PS_WENO_Option%Sound_Speed_Min = Sound_Speed_Min

#if 0
  if (ALLOCATED(d)) DEALLOCATE (d) 
  if (ALLOCATED(C)) DEALLOCATE (C) 
  if (ALLOCATED(S)) DEALLOCATE (S)
  if (ALLOCATED(W)) DEALLOCATE (W)

  M = (Order-1)/2

  ALLOCATE (d(0:M), C(0:M,0:M), S(0:M,0:M,0:M), W(-M:M))

  call PS_WENO_Coefficient (Order, d, C, S, W)
#endif

  END Subroutine PS_WENO_Option_Setup 
!
! ========================================================
!
  FUNCTION PS_WENO_Mapped_Weight (d, W_k) RESULT (W)
 
  REALTYPE                     :: S_W_k
  REALTYPE, dimension(:)       :: d, W_k
  REALTYPE, dimension(SIZE(d)) :: W, Omega

  S_W_k = SUM(W_k) 
  Omega = W_k/S_W_k
  W     = Omega*(d*(d+ONE)+Omega*(Omega-THREE*d))/                 &
                (d*d      +Omega*(ONE  -  TWO*d))

  END FUNCTION PS_WENO_Mapped_Weight 
 
END MODULE WENO_Option
