#if defined (NAVIER_STOKES)

if (Navier_Stokes) then
  if (ALLOCATED(mu   )) DEALLOCATE (mu)
  if (ALLOCATED(zeta )) DEALLOCATE (zeta)
  if (ALLOCATED(kappa)) DEALLOCATE (kappa)
endif

#endif
