#if defined (VARIABLE_GAMMA)
#define _EULER_GAMMA	Gamma_Gas, Gamma 
#else
#define _EULER_GAMMA	Gamma
#endif

  Subroutine Hybrid_Flux (Q, D_Flux, Flag, Flag_WENO)

  USE PseudoPack 
  ! KAUSHIK ADDED THIS USE STATEMENT


  implicit none

  REALTYPE            :: D_F_Epsilon

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, D_Flux
  integer , dimension(N0:N5,M0:M5,K0:K5,3) , OPTIONAL :: Flag_WENO

  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) , OPTIONAL :: Flag

  integer  :: Index

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P, H, C, T
  REALTYPE, dimension(M0:M5,K0:K5,NV)       :: Lambda_x
  REALTYPE, dimension(N0:N5,K0:K5,NV)       :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV)       :: Lambda_z


  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:NV)    :: xflux, yflux, zflux
  REALTYPE  :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr, el, er 
  integer   :: ii, jj, kk
  
  REALTYPE, dimension(1:3) :: sw_x, sw_y, sw_z  
  REALTYPE, dimension(N0:N5) :: switchx, switchx1  
  REALTYPE, dimension(M0:M5) :: switchy, switchy1  
  REALTYPE, dimension(K0:K5) :: switchz, switchz1  

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: entr, int_energy, Yk
  REALTYPE :: thetai, thetaip1, term1, term2
    
  REALTYPE, PARAMETER :: eps_sw =  2.0d-2    
  REALTYPE, dimension(1:6,1:7) :: uu 
  integer, parameter :: n_cells = 5 

  REALTYPE, dimension(1:NV)           :: dQdx, dQdy, dQdz
 

  REALTYPE  :: xx, yy, zz
  REALTYPE, dimension(1:5) :: zzc 
  REALTYPE, dimension(1:5,1:3) :: vel3 
  INTEGER :: ll
  INTEGER, dimension(1:5) :: ibelong
  character ( len = 100 ) filename


#if defined (VARIABLE_GAMMA)
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: Gamma

  call PS_Specific_Heat_Gamma (Gamma_Formulation, Q, _EULER_GAMMA)
#endif

#include "Adjust_D_Flux_Epsilon.i"


 
!-------------------------------------------------------------------------

  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, K0, K5, NV,  &
                                    Gamma, Q, U, V, W, P, H, C)



!--------------------------------------------------------------
     ! EXAMPLE: 0-24 interior cells, 6 ghost cells either side

             ! N0 = -6
             ! N1 = -1
             ! N2 =  0
             ! N3 =  24
             ! N4 =  25
             ! N5 =  30
!--------------------------------------------------------------
!--------------------------------------------------------------

           D_Flux = 0.0d0


           do kk = K0, K5
            do jj = M0, M5
             do ii = N0, N5

              entr(ii,jj,kk) = P(ii,jj,kk)/(Q(ii,jj,kk,1)**gamma) 

              int_energy(ii,jj,kk) = P(ii,jj,kk)/Q(ii,jj,kk,1)/(gamma-1.0d0)

               if(int_energy(ii,jj,kk).le.0.0d0) then
                write(*,*) 'e < 0 in Hybrid_Flux ', int_energy(ii,jj,kk), ii, jj, kk
                stop
               endif


              Yk(ii,jj,kk) = Q(ii,jj,kk,6)/Q(ii,jj,kk,1)

                if(Yk(ii,jj,kk).lt.0.0d0.or.Yk(ii,jj,kk).gt.1.0d0) then
                  write(*,*) 'Yk out of bounds', Yk(ii,jj,kk), ii, jj, kk
                  stop
                endif

             enddo
            enddo
           enddo
!--------------------------------------------------



           ! I-flux
           xflux = 0.0d0


           do kk = K2, K3
            do jj = M2, M3

             switchx = 0.0d0

             do ii = N2-1, N3+1
                
              term1 = Q(ii+1,jj,kk,1) - 2.0d0*Q(ii,jj,kk,1) + Q(ii-1,jj,kk,1)
              term2 = Q(ii+1,jj,kk,1) + 2.0d0*Q(ii,jj,kk,1) + Q(ii-1,jj,kk,1)
              sw_x(1) = abs(term1/term2) 

              term1 = p(ii+1,jj,kk) - 2.0d0*p(ii,jj,kk) + p(ii-1,jj,kk)
              term2 = p(ii+1,jj,kk) + 2.0d0*p(ii,jj,kk) + p(ii-1,jj,kk)
              sw_x(2) = abs(term1/term2) 

              term1 = entr(ii+1,jj,kk) - 2.0d0*entr(ii,jj,kk) + entr(ii-1,jj,kk)
              term2 = entr(ii+1,jj,kk) + 2.0d0*entr(ii,jj,kk) + entr(ii-1,jj,kk)
              sw_x(3) = abs(term1/term2) 


              if(maxval(sw_x(1:3)).ge.eps_sw) then
               switchx(ii+1) = 1.0d0 
              else
               switchx(ii+1) = 0.0d0
              endif

             enddo



             switchx1 = 0.0d0
             do ii = N2-1, N3+1
              if(switchx(ii+1).eq.1.0d0) then
               switchx1(ii+1-n_cells:ii+1+n_cells) = 1.0d0
              endif
             enddo



             do ii = N2-1, N3+1

              uu(1:6,1) = Q(ii-2:ii+3,jj,kk,1)
              uu(1:6,2) = U(ii-2:ii+3,jj,kk)
              uu(1:6,3) = V(ii-2:ii+3,jj,kk)
              uu(1:6,4) = W(ii-2:ii+3,jj,kk)
              uu(1:6,5) = int_energy(ii-2:ii+3,jj,kk)
              uu(1:6,6) = P(ii-2:ii+3,jj,kk)
              uu(1:6,7) = Yk(ii-2:ii+3,jj,kk)

           if(switchx1(ii+1).eq.1.0d0) then
               !call weno5_BS(gamma,1,uu,xflux(ii+1,jj,kk,1:6)) 
               call WENO_cons_flux(gamma,1,uu,xflux(ii+1,jj,kk,1:6))
           else
                call central_flux(gamma,1,uu,xflux(ii+1,jj,kk,1:6))
               !call weno5_BS(gamma,1,uu,xflux(ii+1,jj,kk,1:6)) 
           endif

             enddo                     



            enddo    
           enddo


           do kk = K2, K3
            do jj = M2, M3
             do ii = N2, N3
              D_Flux(ii,jj,kk,:) = (xflux(ii,jj,kk,:) - xflux(ii+1,jj,kk,:))/dx
             enddo
            enddo    
           enddo


!--------------------------------------------------


           ! J-flux
           yflux = 0.0d0


           do kk = K2, K3
            do ii = N2, N3

             switchy = 0.0d0

             do jj = M2-1, M3+1

              term1 = Q(ii,jj+1,kk,1) - 2.0d0*Q(ii,jj,kk,1) + Q(ii,jj-1,kk,1)
              term2 = Q(ii,jj+1,kk,1) + 2.0d0*Q(ii,jj,kk,1) + Q(ii,jj-1,kk,1)
              sw_y(1) = abs(term1/term2) 

              term1 = p(ii,jj+1,kk) - 2.0d0*p(ii,jj,kk) + p(ii,jj-1,kk)
              term2 = p(ii,jj+1,kk) + 2.0d0*p(ii,jj,kk) + p(ii,jj-1,kk)
              sw_y(2) = abs(term1/term2) 

              term1 = entr(ii,jj+1,kk) - 2.0d0*entr(ii,jj,kk) + entr(ii,jj-1,kk)
              term2 = entr(ii,jj+1,kk) + 2.0d0*entr(ii,jj,kk) + entr(ii,jj-1,kk)
              sw_y(3) = abs(term1/term2) 


              if(maxval(sw_y(1:3)).ge.eps_sw) then
               switchy(jj+1) = 1.0d0 
              else
               switchy(jj+1) = 0.0d0
              endif

             enddo



             switchy1 = 0.0d0
             do jj = M2-1, M3+1
              if(switchy(jj+1).eq.1.0d0) then
               switchy1(jj+1-n_cells:jj+1+n_cells) = 1.0d0
              endif
             enddo



             do jj = M2-1, M3+1

              uu(1:6,1) = Q(ii,jj-2:jj+3,kk,1)
              uu(1:6,2) = U(ii,jj-2:jj+3,kk)
              uu(1:6,3) = V(ii,jj-2:jj+3,kk)
              uu(1:6,4) = W(ii,jj-2:jj+3,kk)
              uu(1:6,5) = int_energy(ii,jj-2:jj+3,kk)
              uu(1:6,6) = P(ii,jj-2:jj+3,kk)
              uu(1:6,7) = Yk(ii,jj-2:jj+3,kk)

           if(switchy1(jj+1).eq.1.0d0) then
               !call weno5_BS(gamma,2,uu,yflux(ii,jj+1,kk,1:6)) 
               call WENO_cons_flux(gamma,2,uu,yflux(ii,jj+1,kk,1:6)) 
           else
                call central_flux(gamma,2,uu,yflux(ii,jj+1,kk,1:6))
               !call weno5_BS(gamma,2,uu,yflux(ii,jj+1,kk,1:6)) 
           endif

             enddo                     



            enddo    
           enddo


           do kk = K2, K3
            do jj = M2, M3
             do ii = N2, N3
              D_Flux(ii,jj,kk,:) = D_Flux(ii,jj,kk,:) & 
          + (yflux(ii,jj,kk,:) - yflux(ii,jj+1,kk,:))/dy
             enddo
            enddo    
           enddo

 
!--------------------------------------------------


           ! K-flux
           zflux = 0.0d0


           do ii = N2, N3
            do jj = M2, M3

             switchz = 0.0d0

             do kk = K2-1, K3+1

              term1 = Q(ii,jj,kk+1,1) - 2.0d0*Q(ii,jj,kk,1) + Q(ii,jj,kk-1,1)
              term2 = Q(ii,jj,kk+1,1) + 2.0d0*Q(ii,jj,kk,1) + Q(ii,jj,kk-1,1)
              sw_z(1) = abs(term1/term2) 

              term1 = p(ii,jj,kk+1) - 2.0d0*p(ii,jj,kk) + p(ii,jj,kk-1)
              term2 = p(ii,jj,kk+1) + 2.0d0*p(ii,jj,kk) + p(ii,jj,kk-1)
              sw_z(2) = abs(term1/term2) 

              term1 = entr(ii,jj,kk+1) - 2.0d0*entr(ii,jj,kk) + entr(ii,jj,kk-1)
              term2 = entr(ii,jj,kk+1) + 2.0d0*entr(ii,jj,kk) + entr(ii,jj,kk-1)
              sw_z(3) = abs(term1/term2) 


              if(maxval(sw_z(1:3)).ge.eps_sw) then
               switchz(kk+1) = 1.0d0 
              else
               switchz(kk+1) = 0.0d0
              endif

             enddo



             switchz1 = 0.0d0
             do kk = K2-1, K3+1
              if(switchz(kk+1).eq.1.0d0) then
               switchz1(kk+1-n_cells:kk+1+n_cells) = 1.0d0
              endif
             enddo



             do kk = K2-1, K3+1

              uu(1:6,1) = Q(ii,jj,kk-2:kk+3,1)
              uu(1:6,2) = U(ii,jj,kk-2:kk+3)
              uu(1:6,3) = V(ii,jj,kk-2:kk+3)
              uu(1:6,4) = W(ii,jj,kk-2:kk+3)
              uu(1:6,5) = int_energy(ii,jj,kk-2:kk+3)
              uu(1:6,6) = P(ii,jj,kk-2:kk+3)
              uu(1:6,7) = Yk(ii,jj,kk-2:kk+3)

           if(switchz1(kk+1).eq.1.0d0) then
               !call weno5_BS(gamma,3,uu,zflux(ii,jj,kk+1,1:6)) 
               call WENO_cons_flux(gamma,3,uu,zflux(ii,jj,kk+1,1:6)) 
           else
                call central_flux(gamma,3,uu,zflux(ii,jj,kk+1,1:6))
               !call weno5_BS(gamma,3,uu,zflux(ii,jj,kk+1,1:6)) 
           endif

             enddo                     



            enddo    
           enddo


           do kk = K2, K3
            do jj = M2, M3
             do ii = N2, N3
              D_Flux(ii,jj,kk,:) = D_Flux(ii,jj,kk,:)  & 
          + (zflux(ii,jj,kk,:) - zflux(ii,jj,kk+1,:))/dz
             enddo
            enddo    
           enddo

!--------------------------------------------------
  
!#if defined (NAVIER_STOKES)
!  if (Navier_Stokes) then

    call PS_Temperature (N0, N5, M0, M5, K0, K5, NV, Q, P, T)


!    call PS_CD_NS_Flux             (       N0, N5, N2, N3,                    &
!                                           M0, M5, M2, M3,                    &
!                                           K0, K5, K2, K3,                    &
!                                           NV, NF, NS,                        &
!                                    Q, U, V, W, T, D_F, D_F_Epsilon,          &
!                                    MPI_Comm_Type=Topology%MPI_Communicator)




        ! KAUSHIK CREATED THE FOLLOWING SUBROUTINE

      D_F = 0.0d0

      call viscous_flux(gamma,mu_0,kappa_0,dx,dy,dz, &
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,Q,U,V,W,P,T,D_F)


    D_Flux = D_Flux + D_F




        ! KAUSHIK CREATED THE FOLLOWING SUBROUTINE

      D_F = 0.0d0
 

      call les_smagorinsky(gamma,dx,dy,dz, &
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,Q,U,V,W,P,T,D_F)

    D_Flux = D_Flux + D_F


!  endif
!#endif

!-------------------------------------------------------------------------
!      KAUSHIK CREATED THIS FOR SPONGE BC IN Z

         D_F = 0.0d0



      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(z(ii,jj,kk).ge.0.25d0) then

          ! use backward differences for gradients of Q 
          if(kk.eq.K0) then
           dQdz = 0.0d0
          else 
           dQdz(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii,jj,kk-1,1:NV))/dz
          endif

          call sponge_zmax(gamma,z(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdz,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo



          D_Flux = D_Flux + D_F



!-------------------------------------------------------------------------
!      KAUSHIK CREATED THIS FOR SPONGE BC IN X

         D_F = 0.0d0

      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(x(ii,jj,kk).le.-0.04d0) then

          if(ii.eq.N5) then
           dQdx(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii-1,jj,kk,1:NV))/dx
          else
           dQdx(1:NV) = (Q(ii+1,jj,kk,1:NV)-Q(ii,jj,kk,1:NV))/dx
          endif

          call sponge_xmin(gamma,x(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdx,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo

          D_Flux = D_Flux + D_F

           
          D_F = 0.0d0

      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(x(ii,jj,kk).ge.0.04d0) then

          if(ii.eq.N0) then
           dQdx(1:NV) = (Q(ii+1,jj,kk,1:NV)-Q(ii,jj,kk,1:NV))/dx
          else
           dQdx(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii-1,jj,kk,1:NV))/dx
          endif

          call sponge_xmax(gamma,x(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdx,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo

          D_Flux = D_Flux + D_F 



!------------------------------------------------------------------------- 

!      KAUSHIK CREATED THIS FOR SPONGE BC IN Y

         D_F = 0.0d0

      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(y(ii,jj,kk).le.-0.04d0) then

          if(jj.eq.M5) then
           dQdy(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii,jj-1,kk,1:NV))/dy
          else
           dQdy(1:NV) = (Q(ii,jj+1,kk,1:NV)-Q(ii,jj,kk,1:NV))/dy
          endif

          call sponge_ymin(gamma,y(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdy,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo

          D_Flux = D_Flux + D_F

           
          D_F = 0.0d0

      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(y(ii,jj,kk).ge.0.04d0) then

          if(jj.eq.M0) then
           dQdy(1:NV) = (Q(ii,jj+1,kk,1:NV)-Q(ii,jj,kk,1:NV))/dy
          else
           dQdy(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii,jj-1,kk,1:NV))/dy
          endif

          call sponge_ymax(gamma,y(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdy,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo

          D_Flux = D_Flux + D_F 



!------------------------------------------------------------------------- 

         ! COMPUTE VELOCITY AT PRESCRIBED LOCATION

      zzc = (/0.125d0, 0.15d0, 0.175d0, 0.2d0, 0.225d0/)

      vel3 = 0.0d0 

      ibelong = 0         

 
      do ii = N2, N3
       do jj = M2, M3
        do kk = K2, K3

          xx = x(ii,jj,kk) 
          yy = y(ii,jj,kk) 
          zz = z(ii,jj,kk) 

          if(xx.gt.-dx/2.0d0.and.xx.le.dx/2.0d0.and. & 
             yy.gt.-dy/2.0d0.and.yy.le.dy/2.0d0) then
             do ll = 1, 5
              if(zz-dz/2.0d0.le.zzc(ll).and.zz+dz/2.0d0.gt.zzc(ll)) then
               ibelong(ll) = 1
               vel3(ll,1) = U(ii,jj,kk)
               vel3(ll,2) = V(ii,jj,kk)
               vel3(ll,3) = W(ii,jj,kk)
              endif 
             enddo
          endif
          
        enddo
       enddo
      enddo




        do ll = 1, 5 
         if(ibelong(ll).eq.1) then
          write(filename,'("vel3/output_",I2.2,".dat")'),ll
          OPEN (451, FILE=Filename, STATUS='UNKNOWN', &
                    POSITION='APPEND', IOSTAT=IOS) 
           write(451,452) time, vel3(ll,1:3) 
          close(451) 
         endif
        enddo  

452        format(d16.8,3d13.5)

!------------------------------------------------------------------------- 
  
            return

!--------------------------------------------------------------


  END Subroutine Hybrid_Flux

#if defined (_EULER_GAMMA)
#undef _EULER_GAMMA
#endif


!-----------------------------------------------------------------------------------------------------
