#if defined (MULTI_RESOLUTION) && defined (MR_OUTPUT)
#include "Multi_Resolution.i"

  if (MR_Activation == 1) then
    if (I_Am == First_Processor) write (6,*) 'MR_Flag Output --> Begin ', Step
    TecPlot_Header = .TRUE. ; call MR_Flag_Plot (MR_Flag_File, MR_Flag, Flag_WENO)
    if (I_Am == First_Processor) write (6,*) 'MR_Flag Output --> Done  '
  endif
#endif

  if (I_Am == First_Processor) write (6,*) 'Movie Output --> Begin ', Step
  TecPlot_Header = .TRUE. ; call TecPlot (Movie_File)
!                           call MatLab  (Movie_File)
  if (I_Am == First_Processor) write (6,*) 'Movie Output --> Done  '

if (Write_Final) then
  TecPlot_Header = .TRUE. ; call TecPlot (Plot_Final_File)
endif

#if defined (ANALYSIS)
if (Analysis) then
  if (I_Am == First_Processor) write (6,*) 'Analysis Output --> Begin ', Step
                            call Analysis_Data ()
  TecPlot_Header = .TRUE. ; call Analysis_Plot (Analysis_Data_File)
  if (I_Am == First_Processor) write (6,*) 'Analysis Output --> Done  '
endif
#endif

  if (I_Am == First_Processor) write (6,*) 'Restart  Output --> Begin ', Step
#include "Write_Restart_File.i"
  if (I_Am == First_Processor) write (6,*) 'Restart  Output --> Done  '
