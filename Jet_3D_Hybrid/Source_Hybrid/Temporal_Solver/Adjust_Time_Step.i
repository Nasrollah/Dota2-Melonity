if (Time <= Final_Time) then
  dt_Original = dt 

  if ( Restart_Resolution >  ZERO) then
    dt_Restart  = Time - PS_Counter(1)* Restart_Resolution
    if (dt_Restart  >  ZERO) dt = MIN(dt_Original-dt_Restart , dt)
  endif

  if (   Movie_Resolution >  ZERO) then
    dt_Movie    = Time - PS_Counter(2)*   Movie_Resolution
    if (dt_Movie    >  ZERO) dt = MIN(dt_Original-dt_Movie   , dt)
  endif

  PS_Switch(1) = dt_Original-dt_Restart  == dt
  PS_Switch(2) = dt_Original-dt_Movie    == dt

  if (ANY(PS_Switch(1:3))) Time = Time - dt_Original + dt
else
  dt_Original = dt ; dt = dt - (Time-Final_Time) ; Time = Final_Time
endif
