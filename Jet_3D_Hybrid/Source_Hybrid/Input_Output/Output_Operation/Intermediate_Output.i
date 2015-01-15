#if defined (MULTI_RESOLUTION) 
#include "MR_Count.i"
#endif

if (PS_Switch(1)) then
  PS_Counter(1) = PS_Counter(1) + 1

#include "Refresh_Output_File.i"

  if (I_Am == First_Processor) write (6,*) 'Restart Output --> Begin ', Step
#include "Write_Restart_File.i"
  if (I_Am == First_Processor) write (6,*) 'Restart Output --> Done  '
endif

if (PS_Switch(2)) then
  PS_Counter(2) = PS_Counter(2) + 1

  if (I_Am == First_Processor) write (6,*) 'Movie Output --> Begin ', Step
  TecPlot_Header = .TRUE. ; call TecPlot (Movie_File)
!                           call MatLab  (Movie_File)
  if (I_Am == First_Processor) write (6,*) 'Movie Output --> Done  '

#if defined (MULTI_RESOLUTION) && defined (MR_OUTPUT)
  if (MR_Activation == 1) then
    if (I_Am == First_Processor) write (6,*) 'MR_Flag Output --> Begin ', Step
    TecPlot_Header = .TRUE. ; call MR_Flag_Plot (MR_Flag_File, MR_Flag, Flag_WENO)
    if (I_Am == First_Processor) write (6,*) 'MR_Flag Output --> Done  '
  endif
#endif
endif

#if defined (ANALYSIS)
if (PS_Switch(3)) then
  PS_Counter(3) = PS_Counter(3) + 1

  if (Analysis) then
    if (I_Am == First_Processor) write (6,*) 'Analysis Output --> Begin ', Step
  
                                call Analysis_Data ()

    if (Analysis_Data_Output == 1) then
      TecPlot_Header = .TRUE. ; call Analysis_Plot (Analysis_Data_File)
    endif

    call Analysis_Resolution_Adjustment (Time, PS_Counter(3))

    if (I_Am == First_Processor) write (6,*) 'Analysis Output --> Done  '
  endif
endif
#endif

  CPU_RK_End   = TIMER

#include "CPU_Timing_Output.i"

