#if defined (MULTI_RESOLUTION) 
  if (MR_Activation == 1) then
    TecPlot_Header = .TRUE. ; call MR_Flag_Plot (MR_Flag_File, MR_Flag, Flag_WENO)
  endif
#endif
