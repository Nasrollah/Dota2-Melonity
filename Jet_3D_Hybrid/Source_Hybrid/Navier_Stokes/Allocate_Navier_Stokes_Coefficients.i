#if defined (NAVIER_STOKES)

#define SPACE	N0:N5,M0:M5,K0:K5

if (Navier_Stokes) then
                           Allocate (   mu(SPACE))
  if (Molecular_Diffusion) Allocate ( zeta(SPACE))
  if (Thermal_Conduction ) Allocate (kappa(SPACE))

#if 0
       mu =    mu_default
     zeta =  zeta_default
    kappa = kappa_default
#else
       mu = ZERO
     zeta = ZERO
    kappa = ZERO
#endif
endif

#if defined (SPACE)
#undef SPACE
#endif

#endif
