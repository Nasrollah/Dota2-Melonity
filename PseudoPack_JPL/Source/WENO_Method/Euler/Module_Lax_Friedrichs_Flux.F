#include "pseudopack.h"


MODULE Lax_Friedrichs_Flux_Splitting

implicit NONE

PRIVATE
PUBLIC  :: PS_Global_Lax_Friedrichs
PUBLIC  :: PS_Local_Lax_Friedrichs

CONTAINS
  
  Subroutine PS_Global_Lax_Friedrichs (N0, N5, NV, Lambda, F, Q, g)

  integer  :: N0, N5, NV

  REALTYPE                        :: Lambda
  REALTYPE, dimension(N0:N5,NV)   :: F, Q
  REALTYPE, dimension(N0:N5,NV,2) :: g

  integer :: i

  do i = N0,N5
    g(i,:,1) = HALF*(F(i,:) + Lambda*Q(i,:))
    g(i,:,2) = HALF*(F(i,:) - Lambda*Q(i,:))
  enddo

  END Subroutine PS_Global_Lax_Friedrichs 
!
! ==================================================================
!
  Subroutine PS_Local_Lax_Friedrichs (N0, N5, NV, Lambda, F, Q, g)

  integer  :: N0, N5, NV

  REALTYPE, dimension(N0:N5)      :: Lambda
  REALTYPE, dimension(N0:N5,NV)   :: F, Q
  REALTYPE, dimension(N0:N5,NV,2) :: g

  integer :: i

  do i = N0,N5
    g(i,:,1) = HALF*(F(i,:) + Lambda(i)*Q(i,:))
    g(i,:,2) = HALF*(F(i,:) - Lambda(i)*Q(i,:))
  enddo

  END Subroutine PS_Local_Lax_Friedrichs

END MODULE Lax_Friedrichs_Flux_Splitting
