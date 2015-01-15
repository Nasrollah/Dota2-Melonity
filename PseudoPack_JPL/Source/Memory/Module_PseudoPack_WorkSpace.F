#include "pseudopack.h"


MODULE PS_WorkSpace_Allocation
  USE Memory_Allocation_Status

implicit NONE

REALTYPE, dimension(:), ALLOCATABLE :: WorkSpace

PRIVATE

PUBLIC  :: WorkSpace
PUBLIC  :: Allocate_WorkSpace, DeAllocate_WorkSpace

CONTAINS

  Subroutine Allocate_WorkSpace (Size)

  integer  :: Size

  if (.NOT. ALLOCATED(WorkSpace)) then
    ALLOCATE (WorkSpace(Size), STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = Size
    Subroutine_Name  = '  Allocate_WorkSpace'

    call Memory_Error_Check (0)
#endif
  endif

  END Subroutine Allocate_WorkSpace
!
!=======================================================================
!
  Subroutine DeAllocate_WorkSpace (Size)

  integer  :: Size

  if (ALLOCATED(WorkSpace)) then
    DEALLOCATE (WorkSpace, STAT=M_Error)

#if defined (DEBUG)
    Domain_Number    = 0
    Memory_Requested = Size
    Subroutine_Name  = 'DeAllocate_WorkSpace'

    call Memory_Error_Check (1)
#endif
  endif

  END Subroutine DeAllocate_WorkSpace

END MODULE PS_WorkSpace_Allocation
