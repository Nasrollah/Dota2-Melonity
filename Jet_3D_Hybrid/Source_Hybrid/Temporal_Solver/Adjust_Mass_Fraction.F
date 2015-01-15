  Subroutine PS_Adjust_Mass_Fraction (N0, N5, N2, N3,                 &
                                      M0, M5, M2, M3,                 &
                                      K0, K5, K2, K3,                 &
                                      NV, NF, NS, Q,                  &
                                      Adjust, Check, F_Zero) 

  integer  :: N0, N5, N2, N3
  integer  :: M0, M5, M2, M3 
  integer  :: K0, K5, K2, K3 
  integer  :: NV, NF, NS

  integer  :: Adjust, Check
  REALTYPE :: F_Zero

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: f

  REALTYPE :: f_Min, f_Min_g
  REALTYPE :: f_Max, f_Max_g
  INTEGER  :: N_Max, N_Max_g

  integer  :: n

  if (.NOT. (Adjust .OR. Check)) RETURN

  n = NF+1 

  if (Adjust == 1) then
    f     = Q(:,:,:,n)/Q(:,:,:,1)

    WHERE (ABS(f-ONE) <  F_Zero) Q(:,:,:,n) = Q(:,:,:,1)
    WHERE (ABS(f    ) <  F_Zero) Q(:,:,:,n) = ZERO
    WHERE (    f      <  F_Zero) Q(:,:,:,n) = ZERO
    WHERE (    f-ONE  >  F_Zero) Q(:,:,:,n) = Q(:,:,:,1)
  endif
  
  if (Check == 1) then
    f     = Q(:,:,:,n)/Q(:,:,:,1)

    f_Min = MINVAL(f(N2:N3,M2:M3,K2:K3)      )
    f_Max = MAXVAL(f(N2:N3,M2:M3,K2:K3)      )
    N_Max =  COUNT(f(N2:N3,M2:M3,K2:K3) > ONE)

#if defined (PARALLEL_MPI)
  if (N_Processor == 1) then
    f_Min_g = f_Min ; f_Max_g = f_Max ; N_Max_g = N_Max
  else
    call MPI_ALLREDUCE (f_Min, f_Min_g, 1, MPI_REALTYPE, MPI_MIN, &
                                           MPI_Comm_World, MPI_Error_Status)
    call MPI_ALLREDUCE (f_Max, f_Max_g, 1, MPI_REALTYPE, MPI_MAX, &
                                           MPI_Comm_World, MPI_Error_Status)
    call MPI_ALLREDUCE (N_Max, N_Max_g, 1, MPI_INTEGER , MPI_SUM, &
                                           MPI_Comm_World, MPI_Error_Status)
  endif
#else
    f_Min_g = f_Min ; f_Max_g = f_Max ; N_Max_g = N_Max
#endif

    if (I_Am == First_Processor) then
      if ((f_Min_g < ZERO) .OR. (f_Max_g >  ONE)) then
        write (lid6 ,*) 'Massfraction Min/Max/Count : ', &
                        f_Min_g, f_Max_g-ONE, N_Max_g
        write (lid60,*) 'Massfraction Min/Max/Count : ', &
                        f_Min_g, f_Max_g-ONE, N_Max_g
      endif
    endif
  endif

  End Subroutine PS_Adjust_Mass_Fraction
