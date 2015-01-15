#include "pseudopack.h"


!-------------------------------------------------------------------------
! N = Total number of grid points
! 
! Fourier            : F = N+1,  H = (F+1)/2 = N/2+1  , Q = (H+1)/2
! Chebyshev/Legendre : F = N  ,  H = (F+1)/2 = (N+1)/2, Q = (H+1)/2
!
!   Aux_Size = HEADER_SIZE
!
!   if (Transform Algorithm) then
!     if (Fourier   Method) Filter_Size = H
!     if (Chebyshev Method) Filter_Size = N
!
!     Aux_Size = Aux_Size+Filter_Size
!   endif
!
!   if (Differentiation) then
!     if (Full_Domain) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*N
!     if (Half_Domain) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*H
!
!     Aux_Size = Aux_Size+Metric_Size
!   endif
!
!   For Non-Transform Algorithm, 
!     Aux_Size = Aux_Size+D_Size*MAX_DERIVATIVE
!   else
!     Aux_Size = Aux_Size+Trig_Table_Size
!   endif
!
! ----------------------------------------------------------------
!                 Differentiation                  Smooth
! ----------------------------------------------------------------
!                            D_Size
!                            ------
! EOD              2*H*H                            Same
! ----------------------------------------------------------------
!
!                            Fourier Method
! ----------------------------------------------------------------
!                            D_Size
!                            ------
! MXM              N*N                             Same
! ----------------------------------------------------------------
!
!                            Trig_Table_Size 
!                            ---------------
#if defined (ESSL)
!      Full      [ 50000,44000]                   Same
!      Half      [100000,85000]                   [60000,50000]
! ----------------------------------------------------------------
#endif
#if defined (VFFT)
!      Full      N+15                             Same
!      Half      5*H+Q+30                         Same
! ----------------------------------------------------------------
#endif
#if defined (CRAY)
!      Full      2*N+19                           Same
!      Half      2*(Q+2*H+19)                     Same
! ----------------------------------------------------------------
#endif
!
!                            Chebyshev or Legendre
! ----------------------------------------------------------------
!                            D_Size
!                            ------
! MXM              N*N                            Same
! ----------------------------------------------------------------
!
!                            Chebyshev Only
! ----------------------------------------------------------------
!
!                            Trig_Table_Size
!                            ---------------
#if defined (ESSL)
!      Full      [40000,35000]                    Same
!      Half      [65000,57000]+H                  Same
! ----------------------------------------------------------------
#endif
#if defined (VFFT)
!      Full         3*N+15                        Same
!      Half      2*(3*H+15)                       Same
! ----------------------------------------------------------------
#endif
#if defined (CRAY)
!      Full      H+2*N+19                         Same
!      Half      Q+3*H+19                         Same
! ----------------------------------------------------------------
#endif

MODULE PseudoPack_Aux_Storage

Implicit NONE

PRIVATE

PUBLIC  :: Get_Aux_XXXX_Size

PUBLIC  ::  MXM_Aux
PUBLIC  ::  EOD_Aux
#if defined (ESSL)
PUBLIC  :: ESSL_Aux
#endif
#if defined (VFFT)
PUBLIC  :: VFFT_Aux
#endif
#if defined (CRAY)
PUBLIC  :: CRAY_Aux
#endif

CONTAINS

  Subroutine Get_Aux_XXXX_Size (D_Or_S, Method, Algorithm, Symmetry, &
                                        S_Smooth, N, M_Required)

  character :: D_Or_S
  integer   :: Method, Algorithm, Symmetry, S_Smooth, N, M_Required

  if (D_Or_S == 'd') D_Or_S = 'D'
  if (D_Or_S == 's') D_Or_S = 'S'

  SELECT CASE (Algorithm)
    CASE (0)
      M_Required =  MXM_Aux (D_Or_S, Method, Symmetry, S_Smooth, N)

    CASE (1)
      M_Required =  EOD_Aux (D_Or_S, Method, Symmetry, S_Smooth, N)

    CASE (2)
#if defined (CRAY)
      M_Required = CRAY_Aux (D_Or_S, Method, Symmetry, S_Smooth, N)
#endif
#if defined (VFFT)
      M_Required = VFFT_Aux (D_Or_S, Method, Symmetry, S_Smooth, N)
#endif
#if defined (ESSL)
      M_Required = ESSL_Aux (D_Or_S, Method, Symmetry, S_Smooth, N)
#endif

  END SELECT

  END Subroutine Get_Aux_XXXX_Size 
!
!========================================================================
!
  Integer FUNCTION  MXM_Aux  (D_Or_S, Method, Symmetry, S_Smooth, N)

  character :: D_Or_S
  integer   :: Method, Symmetry, S_Smooth, N
  integer   :: Aux_Size, Metric_Size

  integer   :: MXM_Full_Size

  MXM_Full_Size(n) = n*n

  Aux_Size = HEADER_SIZE

  if (D_Or_S == 'S') then
    if (S_Smooth == 0) then ; MXM_Aux = Aux_Size ; RETURN ; endif

    Aux_Size = Aux_Size +                   MXM_Full_Size(N)
  endif

  if (D_Or_S == 'D') then
       Aux_Size = Aux_Size + MAX_DERIVATIVE*MXM_Full_Size(N)

    Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*N

    Aux_Size = Aux_Size + Metric_Size
  endif

  MXM_Aux = Aux_Size 

  END FUNCTION  MXM_Aux 
!
!========================================================================
!
  Integer FUNCTION  EOD_Aux  (D_Or_S, Method, Symmetry, S_Smooth, N)

  character :: D_Or_S
  integer   :: Method, Symmetry, S_Smooth, N
  integer   :: Aux_Size, Metric_Size, F, H

  integer   :: EOD_Diff_Size

  EOD_Diff_Size(n) = 2*n*n

  if (Method == 0) then ; F = N+1 ; else ; F = N       ; endif

    H = (F+1)/2

  Aux_Size = HEADER_SIZE

  if (D_Or_S == 'S') then
    if (S_Smooth == 0) then ; EOD_Aux = Aux_Size ; RETURN ; endif

    Aux_Size = Aux_Size + EOD_Diff_Size(H)
  endif

  if (D_Or_S == 'D') then
    Aux_Size = Aux_Size + MAX_DERIVATIVE*EOD_Diff_Size(H)

    if (Symmetry == 0) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*N
    if (Symmetry == 1) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*H

    Aux_Size = Aux_Size + Metric_Size
  endif

  EOD_Aux = Aux_Size 

  END FUNCTION  EOD_Aux 
!
!========================================================================
!
#if defined (ESSL)
  Integer FUNCTION ESSL_Aux  (D_Or_S, Method, Symmetry, S_Smooth, N)

  character :: D_Or_S
  integer   :: Method, Symmetry, S_Smooth, N
  integer   :: Aux_Size, Metric_Size, Filter_Size, F, H, Q

  integer   :: Fourier_Full_Size
  integer   :: Fourier_Diff_Half_Size
  integer   :: Fourier_Smth_Half_Size
  integer   :: Chebyshev_Full_Size
  integer   :: Chebyshev_Half_Size

#if defined (SINGLE_PRECISION)
         Fourier_Full_Size(n,h) =  25000+25000
    Fourier_Diff_Half_Size(h,q) =  40000+60000
    Fourier_Smth_Half_Size(h,q) =  60000

       Chebyshev_Full_Size(n,h) =  40000
       Chebyshev_Half_Size(h,q) =  40000+25000+h
#else
         Fourier_Full_Size(n,h) =  22000+22000
    Fourier_Diff_Half_Size(h,q) =  35000+50000
    Fourier_Smth_Half_Size(h,q) =  50000

       Chebyshev_Full_Size(n,h) =  35000
       Chebyshev_Half_Size(h,q) =  35000+22000+h
#endif

  if (Method == 0) then ; F = N+1 ; else ; F = N       ; endif

    H = (F+1)/2
    Q = (H+1)/2

  Aux_Size = HEADER_SIZE

  if ((D_Or_S == 'S') .AND. (S_Smooth == 0)) then
    ESSL_Aux = Aux_Size ; RETURN
  endif

  if (Method == 0) then ; Filter_Size = (F+1)/2 ; else ; Filter_Size = F ; endif

  Aux_Size = Aux_Size+Filter_Size

  if (D_Or_S == 'D') then
    if (Symmetry == 0) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*N
    if (Symmetry == 1) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*H

    Aux_Size = Aux_Size+Metric_Size
  endif

  if (Method == 0) then
    if (Symmetry == 0) then
      Aux_Size = Aux_Size+Fourier_Full_Size(N,H)
    else
      if (D_Or_S == 'D') then
        Aux_Size = Aux_Size+Fourier_Diff_Half_Size(H,Q)
      else
        Aux_Size = Aux_Size+Fourier_Smth_Half_Size(H,Q)
      endif
    endif
  endif

  if (Method == 1) then
    if (Symmetry == 0) then
      Aux_Size = Aux_Size+Chebyshev_Full_Size(N,H)
    else
      Aux_Size = Aux_Size+Chebyshev_Half_Size(H,Q)
    endif
  endif

  ESSL_Aux = Aux_Size 

  END FUNCTION ESSL_Aux
#endif



#if defined (VFFT)
  Integer FUNCTION VFFT_Aux  (D_Or_S, Method, Symmetry, S_Smooth, N)

  character :: D_Or_S
  integer   :: Method, Symmetry, S_Smooth, N
  integer   :: Aux_Size, Metric_Size, Filter_Size, F, H, Q

  integer   ::   Fourier_Full_Size
  integer   ::   Fourier_Half_Size
  integer   :: Chebyshev_Full_Size
  integer   :: Chebyshev_Half_Size

    Fourier_Full_Size(n,h) = n+15
    Fourier_Half_Size(h,q) = q+5*h+30
  Chebyshev_Full_Size(n,h) = 3*n+15
  Chebyshev_Half_Size(h,q) = 2*(3*h+15)

  if (Method == 0) then ; F = N+1 ; else ; F = N       ; endif

    H = (F+1)/2
    Q = (H+1)/2

  Aux_Size = HEADER_SIZE

  if ((D_Or_S == 'S') .AND. (S_Smooth == 0)) then
    VFFT_Aux = Aux_Size ; RETURN
  endif

  if (Method == 0) then ; Filter_Size = (F+1)/2 ; else ; Filter_Size = F ; endif

  Aux_Size = Aux_Size+Filter_Size

  if (D_Or_S == 'D') then
    if (Symmetry == 0) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*N
    if (Symmetry == 1) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*H

    Aux_Size = Aux_Size+Metric_Size
  endif

  if (Method == 0) then
    if (Symmetry == 0) then
      Aux_Size = Aux_Size+  Fourier_Full_Size(N,H)
    else
      Aux_Size = Aux_Size+  Fourier_Half_Size(H,Q)
    endif
  endif

  if (Method == 1) then
    if (Symmetry == 0) then
      Aux_Size = Aux_Size+Chebyshev_Full_Size(N,H)
    else
      Aux_Size = Aux_Size+Chebyshev_Half_Size(H,Q)
    endif
  endif

  VFFT_Aux = Aux_Size 

  END FUNCTION VFFT_Aux
#endif



#if defined (CRAY)
  Integer FUNCTION CRAY_Aux  (D_Or_S, Method, Symmetry, S_Smooth, N)

  character :: D_Or_S
  integer   :: Method, Symmetry, S_Smooth, N
  integer   :: Aux_Size, Metric_Size, Filter_Size, F, H, Q

  integer   ::   Fourier_Full_Size
  integer   ::   Fourier_Half_Size
  integer   :: Chebyshev_Full_Size
  integer   :: Chebyshev_Half_Size

    Fourier_Full_Size(n,h) = 2*n+19
    Fourier_Half_Size(h,q) = 2*(q+2*h+19)
  Chebyshev_Full_Size(n,h) = h+2*n+19
  Chebyshev_Half_Size(h,q) = q+3*h+19

  if (Method == 0) then ; F = N+1 ; else ; F = N       ; endif

    H = (F+1)/2
    Q = (H+1)/2

  Aux_Size = HEADER_SIZE

  if ((D_Or_S == 'S') .AND. (S_Smooth == 0)) then
    CRAY_Aux = Aux_Size ; RETURN
  endif

  if (Method == 0) then ; Filter_Size = (F+1)/2 ; else ; Filter_Size = F ; endif

  Aux_Size = Aux_Size+Filter_Size

  if (D_Or_S == 'D') then
    if (Symmetry == 0) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*N
    if (Symmetry == 1) Metric_Size = METRIC_SIZE(MAX_DERIVATIVE)*H

    Aux_Size = Aux_Size+Metric_Size
  endif

  if (Method == 0) then
    if (Symmetry == 0) then
      Aux_Size = Aux_Size+  Fourier_Full_Size(N,H)
    else
      Aux_Size = Aux_Size+  Fourier_Half_Size(H,Q)
    endif
  endif

  if (Method == 1) then
    if (Symmetry == 0) then
      Aux_Size = Aux_Size+Chebyshev_Full_Size(N,H)
    else
      Aux_Size = Aux_Size+Chebyshev_Half_Size(H,Q)
    endif
  endif

  CRAY_Aux = Aux_Size 

  END FUNCTION CRAY_Aux
#endif

END MODULE PseudoPack_Aux_Storage
