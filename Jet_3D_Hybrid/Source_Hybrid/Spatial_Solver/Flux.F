!--------------------------------------
! KAUSHIK COMMENTED THE LINES BELOW
!--------------------------------------


!#if defined (HYBRID)

Subroutine Flux (Q, D_Flux, Flag, Flag_WENO)

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, D_Flux
  integer , dimension(N0:N5,M0:M5,K0:K5,3) , OPTIONAL :: Flag_WENO

  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) , OPTIONAL :: Flag

!  if (MR_Activation == 1) then
    call Hybrid_Flux (Q, D_Flux, Flag, Flag_WENO)
!  else
!    call   WENO_Flux (Q, D_Flux) 
!  endif

END Subroutine Flux

!#else

!Subroutine Flux (Q, D_Flux)

!  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, D_Flux

!    call   WENO_Flux (Q, D_Flux) 

!END Subroutine Flux

!#endif
