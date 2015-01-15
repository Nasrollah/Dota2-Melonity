!------------------------------------------------------------------
!                  KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!------------------------------------------------------------------


       subroutine viscous_flux(gamma,mu,kappa,dx,dy,dz, & 
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,Q,U,V,W,P,T,dflux)

        implicit none

        integer :: N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, dflux
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P, T
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3,1:3)    :: tau
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)    :: q_cond
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)    :: tau_u
        real(kind=8) :: gamma, mu, kappa, dx, dy, dz
        integer :: ii, jj, kk, l, m  
        real(kind=8) :: taux, tauy, tauz, qx, qy, qz


         dflux = 0.0d0

         tau = 0.0d0
         q_cond = 0.0d0



             call stress_tensor(N0,N5,M0,M5,K0,K5,mu,dx,dy,dz,U,V,W,tau) 



         do ii = N2-2, N3+2 
          do jj = M2-2, M3+2
           do kk = K2-2, K3+2


             tau_u(ii,jj,kk,1) = tau(ii,jj,kk,1,1)*U(ii,jj,kk) + & 
               tau(ii,jj,kk,1,2)*V(ii,jj,kk) + tau(ii,jj,kk,1,3)*W(ii,jj,kk)
             tau_u(ii,jj,kk,2) = tau(ii,jj,kk,2,1)*U(ii,jj,kk) + & 
               tau(ii,jj,kk,2,2)*V(ii,jj,kk) + tau(ii,jj,kk,2,3)*W(ii,jj,kk)
             tau_u(ii,jj,kk,3) = tau(ii,jj,kk,3,1)*U(ii,jj,kk) + & 
               tau(ii,jj,kk,3,2)*V(ii,jj,kk) + tau(ii,jj,kk,3,3)*W(ii,jj,kk)


               q_cond(ii,jj,kk,1) = (T(ii-2,jj,kk) - 8.0d0*T(ii-1,jj,kk) + & 
                        8.0d0*T(ii+1,jj,kk) - T(ii+2,jj,kk))/(12.0d0*dx)
              
               q_cond(ii,jj,kk,2) = (T(ii,jj-2,kk) - 8.0d0*T(ii,jj-1,kk) + & 
                        8.0d0*T(ii,jj+1,kk) - T(ii,jj+2,kk))/(12.0d0*dy)
              
               q_cond(ii,jj,kk,3) = (T(ii,jj,kk-2) - 8.0d0*T(ii,jj,kk-1) + & 
                        8.0d0*T(ii,jj,kk+1) - T(ii,jj,kk+2))/(12.0d0*dz)
              
           enddo 
          enddo    
         enddo
         

!---------------------------------------

                ! q_cond in the above loop is just dT/dx, dT/dy, dT/dz
                ! multiply it by -kappa to get q, thermal conduction

                q_cond = -q_cond*kappa

!-----------------------------------------
!  df/dx(i) = (f(i-2)-8f(i-1)+8f(i+1)-f(i+2))/(12*dx)
!-----------------------------------------


         do ii = N2, N3
          do jj = M2, M3
           do kk = K2, K3

            ! viscous flux, x-direction
            l = 1

            m = 1
            taux = (tau(ii-2,jj,kk,l,m) - 8.0d0*tau(ii-1,jj,kk,l,m) + & 
               8.0d0*tau(ii+1,jj,kk,l,m) - tau(ii+2,jj,kk,l,m))/(12.0d0*dx)   
            m = 2  
            tauy = (tau(ii,jj-2,kk,l,m) - 8.0d0*tau(ii,jj-1,kk,l,m) + & 
               8.0d0*tau(ii,jj+1,kk,l,m) - tau(ii,jj+2,kk,l,m))/(12.0d0*dy)   
            m = 3
            tauz = (tau(ii,jj,kk-2,l,m) - 8.0d0*tau(ii,jj,kk-1,l,m) + & 
               8.0d0*tau(ii,jj,kk+1,l,m) - tau(ii,jj,kk+2,l,m))/(12.0d0*dz)   

  
            dflux(ii,jj,kk,2) = dflux(ii,jj,kk,2) + taux + tauy + tauz   

 

            ! viscous flux, y-direction
            l = 2

            m = 1
            taux = (tau(ii-2,jj,kk,l,m) - 8.0d0*tau(ii-1,jj,kk,l,m) + & 
               8.0d0*tau(ii+1,jj,kk,l,m) - tau(ii+2,jj,kk,l,m))/(12.0d0*dx)   
            m = 2  
            tauy = (tau(ii,jj-2,kk,l,m) - 8.0d0*tau(ii,jj-1,kk,l,m) + & 
               8.0d0*tau(ii,jj+1,kk,l,m) - tau(ii,jj+2,kk,l,m))/(12.0d0*dy)   
            m = 3
            tauz = (tau(ii,jj,kk-2,l,m) - 8.0d0*tau(ii,jj,kk-1,l,m) + & 
               8.0d0*tau(ii,jj,kk+1,l,m) - tau(ii,jj,kk+2,l,m))/(12.0d0*dz)   
 
 
            dflux(ii,jj,kk,3) = dflux(ii,jj,kk,3) + taux + tauy + tauz               


 
            ! viscous flux, z-direction
            l = 3
  
            m = 1
            taux = (tau(ii-2,jj,kk,l,m) - 8.0d0*tau(ii-1,jj,kk,l,m) + & 
               8.0d0*tau(ii+1,jj,kk,l,m) - tau(ii+2,jj,kk,l,m))/(12.0d0*dx)   
            m = 2  
            tauy = (tau(ii,jj-2,kk,l,m) - 8.0d0*tau(ii,jj-1,kk,l,m) + & 
               8.0d0*tau(ii,jj+1,kk,l,m) - tau(ii,jj+2,kk,l,m))/(12.0d0*dy)   
            m = 3
            tauz = (tau(ii,jj,kk-2,l,m) - 8.0d0*tau(ii,jj,kk-1,l,m) + & 
               8.0d0*tau(ii,jj,kk+1,l,m) - tau(ii,jj,kk+2,l,m))/(12.0d0*dz)   
 

            dflux(ii,jj,kk,4) = dflux(ii,jj,kk,4) + taux + tauy + tauz 



            ! viscous dissipation

            m = 1 
            taux = (tau_u(ii-2,jj,kk,m) - 8.0d0*tau_u(ii-1,jj,kk,m) + & 
               8.0d0*tau_u(ii+1,jj,kk,m) - tau_u(ii+2,jj,kk,m))/(12.0d0*dx)   

            m = 2 
            tauy = (tau_u(ii,jj-2,kk,m) - 8.0d0*tau_u(ii,jj-1,kk,m) + & 
               8.0d0*tau_u(ii,jj+1,kk,m) - tau_u(ii,jj+2,kk,m))/(12.0d0*dy)   

            m = 3 
            tauz = (tau_u(ii,jj,kk-2,m) - 8.0d0*tau_u(ii,jj,kk-1,m) + & 
               8.0d0*tau_u(ii,jj,kk+1,m) - tau_u(ii,jj,kk+2,m))/(12.0d0*dz)   

            dflux(ii,jj,kk,5) = dflux(ii,jj,kk,5) + taux + tauy + tauz 



            ! thermal conduction           

            m = 1
            qx = (q_cond(ii-2,jj,kk,m) - 8.0d0*q_cond(ii-1,jj,kk,m) + & 
               8.0d0*q_cond(ii+1,jj,kk,m) - q_cond(ii+2,jj,kk,m))/(12.0d0*dx)   

            m = 2
            qy = (q_cond(ii,jj-2,kk,m) - 8.0d0*q_cond(ii,jj-1,kk,m) + & 
               8.0d0*q_cond(ii,jj+1,kk,m) - q_cond(ii,jj+2,kk,m))/(12.0d0*dy)   

            m = 3
            qz = (q_cond(ii,jj,kk-2,m) - 8.0d0*q_cond(ii,jj,kk-1,m) + & 
               8.0d0*q_cond(ii,jj,kk+1,m) - q_cond(ii,jj,kk+2,m))/(12.0d0*dz)   

  
            ! NOTE: conduction will have a "-" sign  
            dflux(ii,jj,kk,5) = dflux(ii,jj,kk,5) - (qx + qy + qz) 



           enddo
          enddo
         enddo 


                   


         return

        end subroutine 

!-------------------------------------------------------------


             subroutine stress_tensor(N0,N5,M0,M5,K0,K5,mu,dx,dy,dz,U,V,W,tau) 

          implicit none

          integer :: N0,N5,M0,M5,K0,K5, ii, jj, kk
          real(kind=8) :: mu, dx, dy, dz
          real(kind=8), dimension(N0:N5,M0:M5,K0:K5) :: U, V, W 
          real(kind=8), dimension(N0:N5,M0:M5,K0:K5,1:3,1:3) :: tau 
          real(kind=8) :: div_vel      
          real(kind=8) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz    

               tau = 0.0d0



!-------------------------------------------------


       do ii = N0+2, N5-2 
       do jj = M0+2, M5-2 
       do kk = K0+2, K5-2 


            ! compute velocity gradients

         dudx = (U(ii-2,jj,kk) - 8.0d0*U(ii-1,jj,kk) + 8.0d0*U(ii+1,jj,kk) - U(ii+2,jj,kk))/(12.0d0*dx)
         dvdx = (V(ii-2,jj,kk) - 8.0d0*V(ii-1,jj,kk) + 8.0d0*V(ii+1,jj,kk) - V(ii+2,jj,kk))/(12.0d0*dx)
         dwdx = (W(ii-2,jj,kk) - 8.0d0*W(ii-1,jj,kk) + 8.0d0*W(ii+1,jj,kk) - W(ii+2,jj,kk))/(12.0d0*dx)

         dudy = (U(ii,jj-2,kk) - 8.0d0*U(ii,jj-1,kk) + 8.0d0*U(ii,jj+1,kk) - U(ii,jj+2,kk))/(12.0d0*dy)
         dvdy = (V(ii,jj-2,kk) - 8.0d0*V(ii,jj-1,kk) + 8.0d0*V(ii,jj+1,kk) - V(ii,jj+2,kk))/(12.0d0*dy)
         dwdy = (W(ii,jj-2,kk) - 8.0d0*W(ii,jj-1,kk) + 8.0d0*W(ii,jj+1,kk) - W(ii,jj+2,kk))/(12.0d0*dy)
             
         dudz = (U(ii,jj,kk-2) - 8.0d0*U(ii,jj,kk-1) + 8.0d0*U(ii,jj,kk+1) - U(ii,jj,kk+2))/(12.0d0*dz)
         dvdz = (V(ii,jj,kk-2) - 8.0d0*V(ii,jj,kk-1) + 8.0d0*V(ii,jj,kk+1) - V(ii,jj,kk+2))/(12.0d0*dz)
         dwdz = (W(ii,jj,kk-2) - 8.0d0*W(ii,jj,kk-1) + 8.0d0*W(ii,jj,kk+1) - W(ii,jj,kk+2))/(12.0d0*dz)
             

!-------------------------------------------------

            ! compute divergence of velocity

            div_vel = dudx + dvdy + dwdz

!-------------------------------------------------

            ! compute tau matrix

            tau(ii,jj,kk,1,1) = 2.0d0*dudx - 2.0d0/3.0d0*div_vel 
            tau(ii,jj,kk,1,2) = dudy + dvdx
            tau(ii,jj,kk,1,3) = dudz + dwdx 

            tau(ii,jj,kk,2,1) = tau(ii,jj,kk,1,2)
            tau(ii,jj,kk,2,2) = 2.0d0*dvdy - 2.0d0/3.0d0*div_vel 
            tau(ii,jj,kk,2,3) = dvdz + dwdy

            tau(ii,jj,kk,3,1) = tau(ii,jj,kk,1,3)
            tau(ii,jj,kk,3,2) = tau(ii,jj,kk,2,3)
            tau(ii,jj,kk,3,3) = 2.0d0*dwdz - 2.0d0/3.0d0*div_vel 

 
       enddo
       enddo
       enddo


           ! The above tau only has velocity gradients; multiply by mu
            tau = tau*mu 



                    return
                    end subroutine

!-------------------------------------------------
!-------------------------------------------------

          subroutine find_gradient(u,dx,dudx)

       implicit none

            real(kind=8) :: u(1:5), dx, dudx

            ! u(1) = u(i-2)
            ! u(2) = u(i-1)
            ! u(3) = u(i)
            ! u(4) = u(i+1)
            ! u(5) = u(i+2)

            ! gradient at i
            dudx = (u(1) - 8.0d0*u(2) + 8.0d0*u(4) - u(5))/(12.0d0*dx)

          return

          end subroutine
!-------------------------------------------------
