  call Time_Step (dt) ; dt_Original = dt

  call Stability_Check

  Step = Step + 1 ; Time = Time + dt

#if 0
  if (Time >  Final_Time) then
     dt_Original = dt ; dt = dt - (Time-Final_Time) ; Time = Final_Time
  endif
#endif
