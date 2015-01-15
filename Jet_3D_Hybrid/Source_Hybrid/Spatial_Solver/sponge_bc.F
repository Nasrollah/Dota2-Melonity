!----------------------------------------------------------------------
! SPONGE BC BASED ON J. B. FREUND, AIAA J., Vol. 35(4), pp 740-742, 1997
!----------------------------------------------------------------------
        subroutine sponge_zmax(gamma,zz,qq_jpl,dQdz,q_source_jbf)


       implicit none

       real(kind=8) :: gamma, qq_jpl(1:6), q_source_jbf(1:6), zz
       real(kind=8) :: dQdz(1:6)

       real(kind=8), parameter :: rho_target =  0.12d0
       real(kind=8), parameter :: p_target = 0.1d5

       real(kind=8), parameter :: zsponge = 0.25d0   
       real(kind=8), parameter :: zmax = 0.28d0   
       real(kind=8), parameter :: sigma_sponge = 1.0d5 
       real(kind=8), parameter :: beta_sponge = 3.0d0

       real(kind=8) :: term, c_inf, Uconvective 
       real(kind=8) :: sigma_jbf, Usponge



        c_inf = sqrt(gamma*p_target/rho_target)
        Uconvective = 1.3d0*c_inf


        q_source_jbf = 0.0d0 

 
          term = (zz-zsponge)/(zmax-zsponge)
          if(term.lt.0.0d0) then
           write(*,*) 'bug in sponge bc ', term, zz
           stop 
          endif


          sigma_jbf = sigma_sponge*(term**beta_sponge) 

          Usponge = Uconvective*(term**beta_sponge)


          ! DAMPING TERMS
          q_source_jbf(1) = -sigma_jbf*(qq_jpl(1)-rho_target)   
          q_source_jbf(2) = -sigma_jbf*(qq_jpl(2))  
          q_source_jbf(3) = -sigma_jbf*(qq_jpl(3))  
          q_source_jbf(4) = -sigma_jbf*(qq_jpl(4))  
          q_source_jbf(5) = -sigma_jbf*(qq_jpl(5)-p_target/(gamma-1.0d0))  
          q_source_jbf(6) = -sigma_jbf*(qq_jpl(6)-rho_target)  

  
          ! ARTIFICIAL CONVECTION
          q_source_jbf(1:6) = q_source_jbf(1:6) - Usponge*dQdz(1:6) 



             return
            end subroutine 

!------------------------------------------------------------------------

        subroutine sponge_xmin(gamma,xx,qq_jpl,dQdx,q_source_jbf)


       implicit none

       real(kind=8) :: gamma, qq_jpl(1:6), q_source_jbf(1:6), xx
       real(kind=8) :: dQdx(1:6)

       real(kind=8), parameter :: rho_target =  0.12d0
       real(kind=8), parameter :: p_target = 0.1d5

       real(kind=8), parameter :: xsponge = -0.04d0   
       real(kind=8), parameter :: xmin = -0.0525d0   
       real(kind=8), parameter :: sigma_sponge = 1.0d5 
       real(kind=8), parameter :: beta_sponge = 3.0d0

       real(kind=8) :: term, c_inf, Uconvective 
       real(kind=8) :: sigma_jbf, Usponge



        c_inf = sqrt(gamma*p_target/rho_target)
        Uconvective = -1.3d0*c_inf 
        ! Uconvective has a minus sign for xmin sponge bc


        q_source_jbf = 0.0d0 

 
          term = (xx-xsponge)/(xmin-xsponge)
          if(term.lt.0.0d0) then
           write(*,*) 'bug in sponge bc ', term, xx
           stop 
          endif


          sigma_jbf = sigma_sponge*(term**beta_sponge) 

          Usponge = Uconvective*(term**beta_sponge)


          ! DAMPING TERMS
          q_source_jbf(1) = -sigma_jbf*(qq_jpl(1)-rho_target)   
          q_source_jbf(2) = -sigma_jbf*(qq_jpl(2))  
          q_source_jbf(3) = -sigma_jbf*(qq_jpl(3))  
          q_source_jbf(4) = -sigma_jbf*(qq_jpl(4))  
          q_source_jbf(5) = -sigma_jbf*(qq_jpl(5)-p_target/(gamma-1.0d0))  
          q_source_jbf(6) = -sigma_jbf*(qq_jpl(6)-rho_target)  

  
          ! ARTIFICIAL CONVECTION
          q_source_jbf(1:6) = q_source_jbf(1:6) - Usponge*dQdx(1:6) 



             return
            end subroutine 

!------------------------------------------------------------------------

        subroutine sponge_xmax(gamma,xx,qq_jpl,dQdx,q_source_jbf)


       implicit none

       real(kind=8) :: gamma, qq_jpl(1:6), q_source_jbf(1:6), xx
       real(kind=8) :: dQdx(1:6)

       real(kind=8), parameter :: rho_target =  0.12d0
       real(kind=8), parameter :: p_target = 0.1d5

       real(kind=8), parameter :: xsponge = 0.04d0   
       real(kind=8), parameter :: xmax = 0.0525d0   
       real(kind=8), parameter :: sigma_sponge = 1.0d5 
       real(kind=8), parameter :: beta_sponge = 3.0d0

       real(kind=8) :: term, c_inf, Uconvective 
       real(kind=8) :: sigma_jbf, Usponge



        c_inf = sqrt(gamma*p_target/rho_target)
        Uconvective = 1.3d0*c_inf 


        q_source_jbf = 0.0d0 

 
          term = (xx-xsponge)/(xmax-xsponge)
          if(term.lt.0.0d0) then
           write(*,*) 'bug in sponge bc ', term, xx
           stop 
          endif


          sigma_jbf = sigma_sponge*(term**beta_sponge) 

          Usponge = Uconvective*(term**beta_sponge)


          ! DAMPING TERMS
          q_source_jbf(1) = -sigma_jbf*(qq_jpl(1)-rho_target)   
          q_source_jbf(2) = -sigma_jbf*(qq_jpl(2))  
          q_source_jbf(3) = -sigma_jbf*(qq_jpl(3))  
          q_source_jbf(4) = -sigma_jbf*(qq_jpl(4))  
          q_source_jbf(5) = -sigma_jbf*(qq_jpl(5)-p_target/(gamma-1.0d0))  
          q_source_jbf(6) = -sigma_jbf*(qq_jpl(6)-rho_target)  

  
          ! ARTIFICIAL CONVECTION
          q_source_jbf(1:6) = q_source_jbf(1:6) - Usponge*dQdx(1:6) 



             return
            end subroutine 

!------------------------------------------------------------------------

        subroutine sponge_ymin(gamma,yy,qq_jpl,dQdy,q_source_jbf)


       implicit none

       real(kind=8) :: gamma, qq_jpl(1:6), q_source_jbf(1:6), yy
       real(kind=8) :: dQdy(1:6)

       real(kind=8), parameter :: rho_target =  0.12d0
       real(kind=8), parameter :: p_target = 0.1d5

       real(kind=8), parameter :: ysponge = -0.04d0   
       real(kind=8), parameter :: ymin = -0.0525d0   
       real(kind=8), parameter :: sigma_sponge = 1.0d5 
       real(kind=8), parameter :: beta_sponge = 3.0d0

       real(kind=8) :: term, c_inf, Uconvective 
       real(kind=8) :: sigma_jbf, Usponge



        c_inf = sqrt(gamma*p_target/rho_target)
        Uconvective = -1.3d0*c_inf 
        ! Uconvective has a minus sign for ymin sponge bc


        q_source_jbf = 0.0d0 

 
          term = (yy-ysponge)/(ymin-ysponge)
          if(term.lt.0.0d0) then
           write(*,*) 'bug in sponge bc ', term, yy
           stop 
          endif


          sigma_jbf = sigma_sponge*(term**beta_sponge) 

          Usponge = Uconvective*(term**beta_sponge)


          ! DAMPING TERMS
          q_source_jbf(1) = -sigma_jbf*(qq_jpl(1)-rho_target)   
          q_source_jbf(2) = -sigma_jbf*(qq_jpl(2))  
          q_source_jbf(3) = -sigma_jbf*(qq_jpl(3))  
          q_source_jbf(4) = -sigma_jbf*(qq_jpl(4))  
          q_source_jbf(5) = -sigma_jbf*(qq_jpl(5)-p_target/(gamma-1.0d0))  
          q_source_jbf(6) = -sigma_jbf*(qq_jpl(6)-rho_target)  

  
          ! ARTIFICIAL CONVECTION
          q_source_jbf(1:6) = q_source_jbf(1:6) - Usponge*dQdy(1:6) 



             return
            end subroutine 

!------------------------------------------------------------------------

        subroutine sponge_ymax(gamma,yy,qq_jpl,dQdy,q_source_jbf)


       implicit none

       real(kind=8) :: gamma, qq_jpl(1:6), q_source_jbf(1:6), yy
       real(kind=8) :: dQdy(1:6)

       real(kind=8), parameter :: rho_target =  0.12d0
       real(kind=8), parameter :: p_target = 0.1d5

       real(kind=8), parameter :: ysponge = 0.04d0   
       real(kind=8), parameter :: ymax = 0.0525d0   
       real(kind=8), parameter :: sigma_sponge = 1.0d5 
       real(kind=8), parameter :: beta_sponge = 3.0d0

       real(kind=8) :: term, c_inf, Uconvective 
       real(kind=8) :: sigma_jbf, Usponge



        c_inf = sqrt(gamma*p_target/rho_target)
        Uconvective = 1.3d0*c_inf 


        q_source_jbf = 0.0d0 

 
          term = (yy-ysponge)/(ymax-ysponge)
          if(term.lt.0.0d0) then
           write(*,*) 'bug in sponge bc ', term, yy
           stop 
          endif


          sigma_jbf = sigma_sponge*(term**beta_sponge) 

          Usponge = Uconvective*(term**beta_sponge)


          ! DAMPING TERMS
          q_source_jbf(1) = -sigma_jbf*(qq_jpl(1)-rho_target)   
          q_source_jbf(2) = -sigma_jbf*(qq_jpl(2))  
          q_source_jbf(3) = -sigma_jbf*(qq_jpl(3))  
          q_source_jbf(4) = -sigma_jbf*(qq_jpl(4))  
          q_source_jbf(5) = -sigma_jbf*(qq_jpl(5)-p_target/(gamma-1.0d0))  
          q_source_jbf(6) = -sigma_jbf*(qq_jpl(6)-rho_target)  

  
          ! ARTIFICIAL CONVECTION
          q_source_jbf(1:6) = q_source_jbf(1:6) - Usponge*dQdy(1:6) 



             return
            end subroutine 

!------------------------------------------------------------------------
