if (Write_IC == 1) then
  TecPlot_Header = .TRUE. ; call TecPlot       (Plot_IC_File) 
endif

  PS_Counter(1:3) = PS_Counter(1:3) + 1

if ((Write_Initial_Movie == 1) .AND. (Restart == 0)) then
  if (I_Am == First_Processor) write (6,*) 'Movie Output --> Begin ', Step

  TecPlot_Header = .TRUE. ; call TecPlot (Movie_File)
!                           call MatLab  (Movie_File)

  if (I_Am == First_Processor) write (6,*) 'Movie Output --> Done  '
endif

#if defined (MULTI_RESOLUTION)
if ((Write_Initial_MR_Flag == 1) .AND. (MR_Activation == 1)) then
  if (I_Am == First_Processor) write (6,*) 'MR_Flag Output --> Begin ', Step

#include "MR_Count_Open_File.i"
#include "Multi_Resolution.i"
#include "Write_MR_Flag_File.i"

  if (I_Am == First_Processor) write (6,*) 'MR_Flag Output --> Done  '
else
#include "MR_Count_Open_File.i"
endif
#endif

#if defined (ANALYSIS)
if (Analysis) then
  if (I_Am == First_Processor) write (6,*) 'Analysis Output --> Begin ', Step

                            call Analysis_Data ()
  TecPlot_Header = .TRUE. ; call Analysis_Plot (Analysis_Data_File)

  if (I_Am == First_Processor) write (6,*) 'Analysis Output --> Done  '
endif
#endif

#include "System_Info.i"
