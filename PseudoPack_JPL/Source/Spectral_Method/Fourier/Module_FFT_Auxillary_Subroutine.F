#include "pseudopack.h"


! ----------------------------------------------------------------------
! FUNCTION NAME: FFT_Filtering
! DESCRIPTION  : Compute the Fourier Differentiation coefficients with
!                or without smoothing for the mapping and no mapping case.
! ----------------------------------------------------------------------

MODULE FFT_Auxillary_Subroutine

implicit NONE

PRIVATE

PUBLIC   :: FFT_D_Coefficients

CONTAINS

  Subroutine FFT_D_Coefficients (Map, Smooth, Even_Or_Odd, O_Of_D, N, Sigma, b)

  integer  :: Map, Smooth, Even_Or_Odd, O_Of_D, N, i, k

  REALTYPE, dimension(0:N)   :: Sigma
  REALTYPE, dimension(0:N,4) :: b

  if (Map == 0) then
    call No_Mapping_Coefficients
  else
    call    Mapping_Coefficients
  endif

  CONTAINS

    Subroutine No_Mapping_Coefficients

    SELECT CASE (Even_Or_Odd)
      CASE DEFAULT
        if (O_Of_D == 1) b(:,1) = (/ ( i   , i=0,N) /)
        if (O_Of_D == 2) b(:,2) = (/ (-i**2, i=0,N) /)
        if (O_Of_D == 3) b(:,3) = (/ (-i**3, i=0,N) /)
        if (O_Of_D == 4) b(:,4) = (/ ( i**4, i=0,N) /)

      CASE (1)
        if (O_Of_D == 1) b(:,1) = (/ (-i   , i=0,N) /)
        if (O_Of_D == 2) b(:,2) = (/ (-i**2, i=0,N) /)
        if (O_Of_D == 3) b(:,3) = (/ ( i**3, i=0,N) /)
        if (O_Of_D == 4) b(:,4) = (/ ( i**4, i=0,N) /)

    END SELECT

    if (Smooth == 1) b(:,O_Of_D) = b(:,O_Of_D)*Sigma

    END Subroutine No_Mapping_Coefficients

! ----------------------------------------------------------------------

    Subroutine    Mapping_Coefficients

    SELECT CASE (Even_Or_Odd)
      CASE DEFAULT
        if (O_Of_D >= 1) b(:,1) = (/ ( i   , i=0,N) /)
        if (O_Of_D >= 2) b(:,2) = (/ (-i**2, i=0,N) /)
        if (O_Of_D >= 3) b(:,3) = (/ (-i**3, i=0,N) /)
        if (O_Of_D >= 4) b(:,4) = (/ ( i**4, i=0,N) /)

      CASE (1)
        if (O_Of_D >= 1) b(:,1) = (/ (-i   , i=0,N) /)
        if (O_Of_D >= 2) b(:,2) = (/ (-i**2, i=0,N) /)
        if (O_Of_D >= 3) b(:,3) = (/ ( i**3, i=0,N) /)
        if (O_Of_D >= 4) b(:,4) = (/ ( i**4, i=0,N) /)

    END SELECT

    if (Smooth == 1) then
      do k = 1,O_Of_D
        b(:,k) = b(:,k)*Sigma
      enddo
    endif

    END Subroutine    Mapping_Coefficients

  END Subroutine FFT_D_Coefficients 

END MODULE FFT_Auxillary_Subroutine
