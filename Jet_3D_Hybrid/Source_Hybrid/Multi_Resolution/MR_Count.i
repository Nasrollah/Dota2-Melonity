  if (MR_Activation == 1) then
      MR_WENO(0) =  SIZE(Flag_WENO(N2:N3,M2:M3,K2:K3,1))
    do k = 1,Dimensions
      MR_WENO(k) = COUNT(Flag_WENO(N2:N3,M2:M3,K2:K3,k) == 1)
    enddo

      MR_Ratio(0           ) = HUNDRED*SUM(MR_WENO(1:Dimensions))/MR_WENO(0)
      MR_Ratio(1:Dimensions) = HUNDRED*    MR_WENO(1:Dimensions) /MR_WENO(0)

    if (MR_Output_Screen == 1) then
      write (6 ,301) I_Am, Step, Time, (MR_WENO(k), MR_Ratio(k), k=0,Dimensions)
    endif
      write (79,302) I_Am, Step, Time, (MR_WENO(k), MR_Ratio(k), k=0,Dimensions)
  endif

 300 format (45x,'         WENO/Ratio '/ &
               2x,'I_Am   Step  ',3x,'Time   ', &
               3x,'Total/Ratio        x            y            z'/  &
               1x,72('-')/)
 301 format ('WENO Coverage :',1x,i5,1x,i5,1x,g12.4, 4(i8,1x,f6.1,:))
 302 format (                  1x,i5,1x,i5,1x,g12.4, 4(i8,1x,f6.1,:))
