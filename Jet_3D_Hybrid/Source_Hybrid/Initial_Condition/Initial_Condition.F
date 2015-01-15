!
!=====================================================================
!
  Subroutine Initial_Condition(time,dt3) 

  real(kind=8) :: time, dt3 
  Specific_Heat = Gamma

     call Uniform_Flow(time,dt3)

  call BC (Q,time,dt3)

  END Subroutine Initial_Condition
