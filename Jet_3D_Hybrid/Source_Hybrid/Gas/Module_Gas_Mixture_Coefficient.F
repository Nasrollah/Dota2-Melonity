!------------
#if 1
!------------

#include "pseudopack.h"


MODULE Gas_Mixture_Coefficient

  USE PseudoPack

  USE Gas_Properties
  USE Gas_Mixture_Formulation

implicit NONE

INTERFACE PS_Gas_Mixture_Coefficient
  MODULE PROCEDURE Gas_Mixture_Coefficient_1D
  MODULE PROCEDURE Gas_Mixture_Coefficient_2D
  MODULE PROCEDURE Gas_Mixture_Coefficient_3D
END INTERFACE

REALTYPE :: S_Power = 0.7860d0
REALTYPE :: T1                 

PRIVATE
PUBLIC  :: PS_Gas_Mixture_Coefficient

CONTAINS
!
!===================================================================
!
  Subroutine Gas_Mixture_Coefficient_1D (Q, T,                                 &
                                           Viscous_Dissipation,   mu_0, mu    ,&
                                           Thermal_Conduction , kappa_0, kappa,&
                                         Molecular_Diffusion  ,  zeta_0, zeta  )

  logical  ::   Viscous_Dissipation 
  logical  ::   Thermal_Conduction
  logical  :: Molecular_Diffusion 

  REALTYPE :: mu_0, kappa_0, zeta_0

  REALTYPE, dimension(:,:)      :: Q
  REALTYPE, dimension(:)        :: T, mu, kappa, zeta

  REALTYPE, dimension(SIZE(Q,DIM=1),2)  :: f, g
  REALTYPE, dimension(SIZE(Q,DIM=1))    :: T_S_Power

  if ((mu_0 <  ZERO) .OR. (kappa_0 <  ZERO) .OR. ( zeta_0 <  ZERO)) then 
#if defined (DEBUG)
    call Debug_1D
#endif

    f(:,1) = Q(:,4)/Q(:,1) ; f(:,2) = ONE-f(:,1)

    T_S_Power = T**S_Power
  endif

  call   Viscous_Dissipation_1D
  call   Thermal_Conduction_1D
  call Molecular_Diffusion_1D

  CONTAINS
!
! -----------------------------------------------------
!
  SUBROUTINE Viscous_Dissipation_1D

  if (.NOT. Viscous_Dissipation) RETURN

  if (mu_0 >  ZERO) then ; mu = mu_0 ; RETURN ; endif

#if 0
  call PS_Viscous_Dissipativity (Gas(1)%T_mu   , Gas(1)%mu   , T, g(:,1))
  call PS_Viscous_Dissipativity (Gas(2)%T_mu   , Gas(2)%mu   , T, g(:,2))
#else
  T1 = Gas(1)%mu/Gas(1)%T_mu**S_Power ; g(:,1) = T1*T_S_Power
  T1 = Gas(2)%mu/Gas(2)%T_mu**S_Power ; g(:,2) = T1*T_S_Power
#endif

  call PS_Mixture_Coefficient (g(:,1), f(:,1), Gas(1)%MW, &
                               g(:,2), f(:,2), Gas(2)%MW, mu   )

  END SUBROUTINE Viscous_Dissipation_1D
!
! -----------------------------------------------------
!
  SUBROUTINE Thermal_Conduction_1D

  if (.NOT. Thermal_Conduction) RETURN

  if (kappa_0 >  ZERO) then ; kappa = kappa_0 ; RETURN ; endif

#if 0
  call PS_Thermal_Conductivity  (Gas(1)%T_kappa, Gas(1)%kappa, T, g(:,1))
  call PS_Thermal_Conductivity  (Gas(2)%T_kappa, Gas(2)%kappa, T, g(:,2))
#else
  T1 = Gas(1)%kappa/Gas(1)%T_kappa**S_Power ; g(:,1) = T1*T_S_Power
  T1 = Gas(2)%kappa/Gas(2)%T_kappa**S_Power ; g(:,2) = T1*T_S_Power
#endif

  call PS_Mixture_Coefficient (g(:,1), f(:,1), Gas(1)%MW, &
                               g(:,2), f(:,2), Gas(2)%MW, kappa)

  END SUBROUTINE Thermal_Conduction_1D
!
! -----------------------------------------------------
!
  SUBROUTINE Molecular_Diffusion_1D 

  if (.NOT. Molecular_Diffusion) RETURN

  if (zeta_0 >  ZERO) then ; zeta = zeta_0 ; RETURN ; endif

#if 0
  call PS_Molecular_Diffusivity (Gas(1)%T_zeta , Gas(1)%zeta , T, g(:,1))
  call PS_Molecular_Diffusivity (Gas(2)%T_zeta , Gas(2)%zeta , T, g(:,2))
#else
  T1 = Gas(1)%zeta/Gas(1)%T_zeta**S_Power ; g(:,1) = T1*T_S_Power
  T1 = Gas(2)%zeta/Gas(2)%T_zeta**S_Power ; g(:,2) = T1*T_S_Power
#endif

  call PS_Mixture_Coefficient (g(:,1), f(:,1), Gas(1)%MW, &
                               g(:,2), f(:,2), Gas(2)%MW, zeta )

  END SUBROUTINE Molecular_Diffusion_1D 
#if defined (DEBUG)
!
! -----------------------------------------------------
!
  SUBROUTINE Debug_1D

  integer :: N1, NS

  if (I_Am == First_Processor) then
    N1 = SIZE(Q,DIM=1) 
    NS = SIZE(Q,DIM=2)-3+1

    if (NS /= 2) write (6,*) 'Gas_Mixture_Coef : NS ',NS,' /= 2'
    if (NS /= 2) STOP 
  endif

  END SUBROUTINE Debug_1D
#endif

  END Subroutine Gas_Mixture_Coefficient_1D 
!
!===================================================================
!
  Subroutine Gas_Mixture_Coefficient_2D (Q, T,                                 &
                                           Viscous_Dissipation,   mu_0, mu    ,&
                                           Thermal_Conduction , kappa_0, kappa,&
                                         Molecular_Diffusion  ,  zeta_0, zeta  )

  logical  ::   Viscous_Dissipation 
  logical  ::   Thermal_Conduction
  logical  :: Molecular_Diffusion 

  REALTYPE :: mu_0, kappa_0, zeta_0

  REALTYPE, dimension(:,:,:)  :: Q
  REALTYPE, dimension(:,:)    :: T, mu, kappa, zeta

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),2) :: f, g
  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2))   :: T_S_Power

  if ((mu_0 <  ZERO) .OR. (kappa_0 <  ZERO) .OR. ( zeta_0 <  ZERO)) then 
#if defined (DEBUG)
    call Debug_2D
#endif

    f(:,:,1) = Q(:,:,5)/Q(:,:,1) ; f(:,:,2) = ONE-f(:,:,1)

    T_S_Power = T**S_Power
  endif

  call   Viscous_Dissipation_2D
  call   Thermal_Conduction_2D
  call Molecular_Diffusion_2D

  CONTAINS
!
! -----------------------------------------------------
!
  SUBROUTINE Viscous_Dissipation_2D

  if (.NOT. Viscous_Dissipation) RETURN

  if (mu_0 >  ZERO) then ; mu = mu_0 ; RETURN ; endif

#if 0
  call PS_Viscous_Dissipativity (Gas(1)%T_mu   , Gas(1)%mu   , T, g(:,:,1))
  call PS_Viscous_Dissipativity (Gas(2)%T_mu   , Gas(2)%mu   , T, g(:,:,2))
#else
  T1 = Gas(1)%mu/Gas(1)%T_mu**S_Power ; g(:,:,1) = T1*T_S_Power
  T1 = Gas(2)%mu/Gas(2)%T_mu**S_Power ; g(:,:,2) = T1*T_S_Power
#endif

  call PS_Mixture_Coefficient (g(:,:,1), f(:,:,1), Gas(1)%MW, &
                               g(:,:,2), f(:,:,2), Gas(2)%MW, mu   )

  END SUBROUTINE Viscous_Dissipation_2D
!
! -----------------------------------------------------
!
  SUBROUTINE Thermal_Conduction_2D

  if (.NOT. Thermal_Conduction) RETURN

  if (kappa_0 >  ZERO) then ; kappa = kappa_0 ; RETURN ; endif

#if 0
  call PS_Thermal_Conductivity  (Gas(1)%T_kappa, Gas(1)%kappa, T, g(:,:,1))
  call PS_Thermal_Conductivity  (Gas(2)%T_kappa, Gas(2)%kappa, T, g(:,:,2))
#else
  T1 = Gas(1)%kappa/Gas(1)%T_kappa**S_Power ; g(:,:,1) = T1*T_S_Power
  T1 = Gas(2)%kappa/Gas(2)%T_kappa**S_Power ; g(:,:,2) = T1*T_S_Power
#endif

  call PS_Mixture_Coefficient (g(:,:,1), f(:,:,1), Gas(1)%MW, &
                               g(:,:,2), f(:,:,2), Gas(2)%MW, kappa)

  END SUBROUTINE Thermal_Conduction_2D
!
! -----------------------------------------------------
!
  SUBROUTINE Molecular_Diffusion_2D 

  if (.NOT. Molecular_Diffusion) RETURN

  if (zeta_0 >  ZERO) then ; zeta = zeta_0 ; RETURN ; endif

#if 0
  call PS_Molecular_Diffusivity (Gas(1)%T_zeta , Gas(1)%zeta , T, g(:,:,1))
  call PS_Molecular_Diffusivity (Gas(2)%T_zeta , Gas(2)%zeta , T, g(:,:,2))
#else
  T1 = Gas(1)%zeta/Gas(1)%T_zeta**S_Power ; g(:,:,1) = T1*T_S_Power
  T1 = Gas(2)%zeta/Gas(2)%T_zeta**S_Power ; g(:,:,2) = T1*T_S_Power
#endif

  call PS_Mixture_Coefficient (g(:,:,1), f(:,:,1), Gas(1)%MW, &
                               g(:,:,2), f(:,:,2), Gas(2)%MW, zeta )

  END SUBROUTINE Molecular_Diffusion_2D 
#if defined (DEBUG)
!
! -----------------------------------------------------
!
  SUBROUTINE Debug_2D

  integer  :: N1, N2, NS

  if (I_Am == First_Processor) then
    N1 = SIZE(Q,DIM=1) ; N2 = SIZE(Q,DIM=2) 
    NS = SIZE(Q,DIM=3)-4+1

    if (NS /= 2) write (6,*) 'Gas_Mixture_Coef : NS ',NS,' /= 2'
    if (NS /= 2) STOP
  endif

  END SUBROUTINE Debug_2D
#endif

  END Subroutine Gas_Mixture_Coefficient_2D 
!
!===================================================================
!
  Subroutine Gas_Mixture_Coefficient_3D (Q, T,                                 &
                                           Viscous_Dissipation,   mu_0, mu    ,&
                                           Thermal_Conduction , kappa_0, kappa,&
                                         Molecular_Diffusion  ,  zeta_0, zeta  )

  logical  ::   Viscous_Dissipation 
  logical  ::   Thermal_Conduction
  logical  :: Molecular_Diffusion 

  REALTYPE :: mu_0, kappa_0, zeta_0

  REALTYPE, dimension(:,:,:,:)  :: Q
  REALTYPE, dimension(:,:,:)    :: T, mu, kappa, zeta

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3),2) :: f, g
  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3))   :: T_S_Power

  if ((mu_0 <  ZERO) .OR. (kappa_0 <  ZERO) .OR. ( zeta_0 <  ZERO)) then 
#if defined (DEBUG)
    call Debug_3D
#endif

    f(:,:,:,1) = Q(:,:,:,6)/Q(:,:,:,1) ; f(:,:,:,2) = ONE-f(:,:,:,1)

    T_S_Power = T**S_Power
  endif

  call   Viscous_Dissipation_3D
  call   Thermal_Conduction_3D
  call Molecular_Diffusion_3D

  CONTAINS
!
! -----------------------------------------------------
!
  SUBROUTINE Viscous_Dissipation_3D

  if (.NOT. Viscous_Dissipation) RETURN


  if (mu_0 >  ZERO) then ; mu = mu_0 ; RETURN ; endif

          ! KAUSHIK ADDED THIS
             write(*,*) 'should not be here '
               stop   


#if 0
  call PS_Viscous_Dissipativity (Gas(1)%T_mu   , Gas(1)%mu   , T, g(:,:,:,1))
  call PS_Viscous_Dissipativity (Gas(2)%T_mu   , Gas(2)%mu   , T, g(:,:,:,2))
#else
  T1 = Gas(1)%mu/Gas(1)%T_mu**S_Power ; g(:,:,:,1) = T1*T_S_Power
  T1 = Gas(2)%mu/Gas(2)%T_mu**S_Power ; g(:,:,:,2) = T1*T_S_Power
#endif

  call PS_Mixture_Coefficient (g(:,:,:,1), f(:,:,:,1), Gas(1)%MW, &
                               g(:,:,:,2), f(:,:,:,2), Gas(2)%MW, mu   )

  END SUBROUTINE Viscous_Dissipation_3D
!
! -----------------------------------------------------
!
  SUBROUTINE Thermal_Conduction_3D

  if (.NOT. Thermal_Conduction) RETURN


  if (kappa_0 >  ZERO) then ; kappa = kappa_0 ; RETURN ; endif

          ! KAUSHIK ADDED THIS
             write(*,*) 'should not be here '
               stop   

#if 0
  call PS_Thermal_Conductivity  (Gas(1)%T_kappa, Gas(1)%kappa, T, g(:,:,:,1))
  call PS_Thermal_Conductivity  (Gas(2)%T_kappa, Gas(2)%kappa, T, g(:,:,:,2))
#else
  T1 = Gas(1)%kappa/Gas(1)%T_kappa**S_Power ; g(:,:,:,1) = T1*T_S_Power
  T1 = Gas(2)%kappa/Gas(2)%T_kappa**S_Power ; g(:,:,:,2) = T1*T_S_Power
#endif

  call PS_Mixture_Coefficient (g(:,:,:,1), f(:,:,:,1), Gas(1)%MW, &
                               g(:,:,:,2), f(:,:,:,2), Gas(2)%MW, kappa)

  END SUBROUTINE Thermal_Conduction_3D
!
! -----------------------------------------------------
!
  SUBROUTINE Molecular_Diffusion_3D 

  if (.NOT. Molecular_Diffusion) RETURN

  if (zeta_0 >  ZERO) then ; zeta = zeta_0 ; RETURN ; endif

#if 0
  call PS_Molecular_Diffusivity (Gas(1)%T_zeta , Gas(1)%zeta , T, g(:,:,:,1))
  call PS_Molecular_Diffusivity (Gas(2)%T_zeta , Gas(2)%zeta , T, g(:,:,:,2))
#else
  T1 = Gas(1)%zeta/Gas(1)%T_zeta**S_Power ; g(:,:,:,1) = T1*T_S_Power
  T1 = Gas(2)%zeta/Gas(2)%T_zeta**S_Power ; g(:,:,:,2) = T1*T_S_Power
#endif

  call PS_Mixture_Coefficient (g(:,:,:,1), f(:,:,:,1), Gas(1)%MW, &
                               g(:,:,:,2), f(:,:,:,2), Gas(2)%MW, zeta )

  END SUBROUTINE Molecular_Diffusion_3D 
#if defined (DEBUG)
!
! -----------------------------------------------------
!
  SUBROUTINE Debug_3D

  integer  :: N1, N2, N3, NS

  if (I_Am == First_Processor) then
    N1 = SIZE(Q,DIM=1) ; N2 = SIZE(Q,DIM=2) ; N3 = SIZE(Q,DIM=3)
    NS = SIZE(Q,DIM=4)-5+1

    if (NS /= 2) write (6,*) 'Gas_Mixture_Coef : NS ',NS,' /= 2'
    if (NS /= 2) STOP
  endif

  END SUBROUTINE Debug_3D
#endif

  END Subroutine Gas_Mixture_Coefficient_3D 

END MODULE Gas_Mixture_Coefficient

!------------
#else
!------------

#include "pseudopack.h"


MODULE Gas_Mixture_Coefficient

  USE PseudoPack
  USE Gas_Properties
  USE Gas_Mixture_Formulation

implicit NONE

INTERFACE PS_Gas_Mixture_Coefficient
  MODULE PROCEDURE Gas_Mixture_Coefficient_1D
  MODULE PROCEDURE Gas_Mixture_Coefficient_2D
  MODULE PROCEDURE Gas_Mixture_Coefficient_3D
END INTERFACE

REALTYPE :: S_Power = 0.7860d0
REALTYPE :: T1                 

PRIVATE
PUBLIC  :: PS_Gas_Mixture_Coefficient

CONTAINS
!
!===================================================================
!
  Subroutine Gas_Mixture_Coefficient_1D (Q, T,                          &
                                           Viscous_Dissipation, mu   ,  &
                                           Thermal_Conduction , kappa,  & 
                                         Molecular_Diffusion  , zeta  )

  logical  ::   Viscous_Dissipation 
  logical  ::   Thermal_Conduction
  logical  :: Molecular_Diffusion 

  REALTYPE, dimension(:,:)      :: Q
  REALTYPE, dimension(:)        :: T, mu, kappa, zeta

  integer                               :: N1, NS
  REALTYPE, dimension(SIZE(Q,DIM=1),2)  :: f, g
  REALTYPE, dimension(SIZE(Q,DIM=1))    :: T_S_Power

  N1 = SIZE(Q,DIM=1) 
  NS = SIZE(Q,DIM=2)-3+1

#if defined (DEBUG)
  if (I_Am == First_Processor) then
    if (NS /= 2) write (6,*) 'Gas_Mixture_Coef : NS ',NS,' /= 2'
    if (NS /= 2) STOP
  endif
#endif

  f(:,1) = Q(:,4)/Q(:,1) ; f(:,2) = ONE-f(:,1)

  T_S_Power = T**S_Power

  if (Viscous_Dissipation) then
#if 0
    call PS_Viscous_Dissipativity (Gas(1)%T_mu   , Gas(1)%mu   , T, g(:,1))
    call PS_Viscous_Dissipativity (Gas(2)%T_mu   , Gas(2)%mu   , T, g(:,2))
#else
    T1 = Gas(1)%mu/Gas(1)%T_mu**S_Power ; g(:,1) = T1*T_S_Power
    T1 = Gas(2)%mu/Gas(2)%T_mu**S_Power ; g(:,2) = T1*T_S_Power
#endif

    call PS_Mixture_Coefficient (g(:,1), f(:,1), Gas(1)%MW, &
                                 g(:,2), f(:,2), Gas(2)%MW, mu   )
  endif

  if (Thermal_Conduction) then
#if 0
    call PS_Thermal_Conductivity  (Gas(1)%T_kappa, Gas(1)%kappa, T, g(:,1))
    call PS_Thermal_Conductivity  (Gas(2)%T_kappa, Gas(2)%kappa, T, g(:,2))
#else
    T1 = Gas(1)%kappa/Gas(1)%T_kappa**S_Power ; g(:,1) = T1*T_S_Power
    T1 = Gas(2)%kappa/Gas(2)%T_kappa**S_Power ; g(:,2) = T1*T_S_Power
#endif

    call PS_Mixture_Coefficient (g(:,1), f(:,1), Gas(1)%MW, &
                                 g(:,2), f(:,2), Gas(2)%MW, kappa)
  endif

  if (Molecular_Diffusion) then 
#if 0
    call PS_Molecular_Diffusivity (Gas(1)%T_zeta , Gas(1)%zeta , T, g(:,1))
    call PS_Molecular_Diffusivity (Gas(2)%T_zeta , Gas(2)%zeta , T, g(:,2))
#else
    T1 = Gas(1)%zeta/Gas(1)%T_zeta**S_Power ; g(:,1) = T1*T_S_Power
    T1 = Gas(2)%zeta/Gas(2)%T_zeta**S_Power ; g(:,2) = T1*T_S_Power
#endif

    call PS_Mixture_Coefficient (g(:,1), f(:,1), Gas(1)%MW, &
                                 g(:,2), f(:,2), Gas(2)%MW, zeta )
  endif

  END Subroutine Gas_Mixture_Coefficient_1D 
!
!===================================================================
!
  Subroutine Gas_Mixture_Coefficient_2D (Q, T,                          &
                                           Viscous_Dissipation, mu   ,  &
                                           Thermal_Conduction , kappa,  & 
                                         Molecular_Diffusion  , zeta  )

  logical  ::   Viscous_Dissipation 
  logical  ::   Thermal_Conduction
  logical  :: Molecular_Diffusion 

  REALTYPE, dimension(:,:,:)      :: Q
  REALTYPE, dimension(:,:)        :: T, mu, kappa, zeta

  integer                                   :: N1, N2, NS

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),2) :: f, g
  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2))   :: T_S_Power

  N1 = SIZE(Q,DIM=1) ; N2 = SIZE(Q,DIM=2) 
  NS = SIZE(Q,DIM=3)-4+1

#if defined (DEBUG)
  if (I_Am == First_Processor) then
    if (NS /= 2) write (6,*) 'Gas_Mixture_Coef : NS ',NS,' /= 2'
    if (NS /= 2) STOP
  endif
#endif

  f(:,:,1) = Q(:,:,5)/Q(:,:,1) ; f(:,:,2) = ONE-f(:,:,1)

  T_S_Power = T**S_Power

  if (Viscous_Dissipation) then
#if 0
    call PS_Viscous_Dissipativity (Gas(1)%T_mu   , Gas(1)%mu   , T, g(:,:,1))
    call PS_Viscous_Dissipativity (Gas(2)%T_mu   , Gas(2)%mu   , T, g(:,:,2))
#else
    T1 = Gas(1)%mu/Gas(1)%T_mu**S_Power ; g(:,:,1) = T1*T_S_Power
    T1 = Gas(2)%mu/Gas(2)%T_mu**S_Power ; g(:,:,2) = T1*T_S_Power
#endif

    call PS_Mixture_Coefficient (g(:,:,1), f(:,:,1), Gas(1)%MW, &
                                 g(:,:,2), f(:,:,2), Gas(2)%MW, mu   )
  endif

  if (Thermal_Conduction) then
#if 0
    call PS_Thermal_Conductivity  (Gas(1)%T_kappa, Gas(1)%kappa, T, g(:,:,1))
    call PS_Thermal_Conductivity  (Gas(2)%T_kappa, Gas(2)%kappa, T, g(:,:,2))
#else
    T1 = Gas(1)%kappa/Gas(1)%T_kappa**S_Power ; g(:,:,1) = T1*T_S_Power
    T1 = Gas(2)%kappa/Gas(2)%T_kappa**S_Power ; g(:,:,2) = T1*T_S_Power
#endif

    call PS_Mixture_Coefficient (g(:,:,1), f(:,:,1), Gas(1)%MW, &
                                 g(:,:,2), f(:,:,2), Gas(2)%MW, kappa)
  endif

  if (Molecular_Diffusion) then 
#if 0
    call PS_Molecular_Diffusivity (Gas(1)%T_zeta , Gas(1)%zeta , T, g(:,:,1))
    call PS_Molecular_Diffusivity (Gas(2)%T_zeta , Gas(2)%zeta , T, g(:,:,2))
#else
    T1 = Gas(1)%zeta/Gas(1)%T_zeta**S_Power ; g(:,:,1) = T1*T_S_Power
    T1 = Gas(2)%zeta/Gas(2)%T_zeta**S_Power ; g(:,:,2) = T1*T_S_Power
#endif

    call PS_Mixture_Coefficient (g(:,:,1), f(:,:,1), Gas(1)%MW, &
                                 g(:,:,2), f(:,:,2), Gas(2)%MW, zeta )
  endif

  END Subroutine Gas_Mixture_Coefficient_2D 
!
!===================================================================
!
  Subroutine Gas_Mixture_Coefficient_3D (Q, T,                          &
                                           Viscous_Dissipation, mu   ,  &
                                           Thermal_Conduction , kappa,  & 
                                         Molecular_Diffusion  , zeta  )

  logical  ::   Viscous_Dissipation 
  logical  ::   Thermal_Conduction
  logical  :: Molecular_Diffusion 

  REALTYPE, dimension(:,:,:,:)      :: Q
  REALTYPE, dimension(:,:,:)        :: T, mu, kappa, zeta

  integer                                   :: N1, N2, N3, NS

  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3),2) :: f, g
  REALTYPE, dimension(SIZE(Q,DIM=1),SIZE(Q,DIM=2),SIZE(Q,DIM=3))   :: T_S_Power

  N1 = SIZE(Q,DIM=1) ; N2 = SIZE(Q,DIM=2) ; N3 = SIZE(Q,DIM=3)
  NS = SIZE(Q,DIM=4)-5+1

#if defined (DEBUG)
  if (I_Am == First_Processor) then
    if (NS /= 2) write (6,*) 'Gas_Mixture_Coef : NS ',NS,' /= 2'
    if (NS /= 2) STOP
  endif
#endif

  f(:,:,:,1) = Q(:,:,:,6)/Q(:,:,:,1) ; f(:,:,:,2) = ONE-f(:,:,:,1)

  T_S_Power = T**S_Power

  if (Viscous_Dissipation) then
#if 0
    call PS_Viscous_Dissipativity (Gas(1)%T_mu   , Gas(1)%mu   , T, g(:,:,:,1))
    call PS_Viscous_Dissipativity (Gas(2)%T_mu   , Gas(2)%mu   , T, g(:,:,:,2))
#else
    T1 = Gas(1)%mu/Gas(1)%T_mu**S_Power ; g(:,:,:,1) = T1*T_S_Power
    T1 = Gas(2)%mu/Gas(2)%T_mu**S_Power ; g(:,:,:,2) = T1*T_S_Power
#endif

    call PS_Mixture_Coefficient (g(:,:,:,1), f(:,:,:,1), Gas(1)%MW, &
                                 g(:,:,:,2), f(:,:,:,2), Gas(2)%MW, mu   )
  endif

  if (Thermal_Conduction) then
#if 0
    call PS_Thermal_Conductivity  (Gas(1)%T_kappa, Gas(1)%kappa, T, g(:,:,:,1))
    call PS_Thermal_Conductivity  (Gas(2)%T_kappa, Gas(2)%kappa, T, g(:,:,:,2))
#else
    T1 = Gas(1)%kappa/Gas(1)%T_kappa**S_Power ; g(:,:,:,1) = T1*T_S_Power
    T1 = Gas(2)%kappa/Gas(2)%T_kappa**S_Power ; g(:,:,:,2) = T1*T_S_Power
#endif

    call PS_Mixture_Coefficient (g(:,:,:,1), f(:,:,:,1), Gas(1)%MW, &
                                 g(:,:,:,2), f(:,:,:,2), Gas(2)%MW, kappa)
  endif

  if (Molecular_Diffusion) then 
#if 0
    call PS_Molecular_Diffusivity (Gas(1)%T_zeta , Gas(1)%zeta , T, g(:,:,:,1))
    call PS_Molecular_Diffusivity (Gas(2)%T_zeta , Gas(2)%zeta , T, g(:,:,:,2))
#else
    T1 = Gas(1)%zeta/Gas(1)%T_zeta**S_Power ; g(:,:,:,1) = T1*T_S_Power
    T1 = Gas(2)%zeta/Gas(2)%T_zeta**S_Power ; g(:,:,:,2) = T1*T_S_Power
#endif

    call PS_Mixture_Coefficient (g(:,:,:,1), f(:,:,:,1), Gas(1)%MW, &
                                 g(:,:,:,2), f(:,:,:,2), Gas(2)%MW, zeta )
  endif

  END Subroutine Gas_Mixture_Coefficient_3D 

END MODULE Gas_Mixture_Coefficient

!------------
#endif
!------------
