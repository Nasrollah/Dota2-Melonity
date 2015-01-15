!-----------------------------------------------------------------------------------------------------
!       KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!-----------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------
!       input: primitive variables at (i-2,i-1,i,i+1,i+2,i+3)
!       output: flux at i+1/2   
!-----------------------------------------------------------------------------------------------------

               ! WENO BASED ON BALSARA AND SHU, JCP 160, 2000

!-----------------------------------------------------------------------------------------------------

              subroutine weno5_BS(gamma,index,uu,flux)


           implicit none


           real(kind=8) :: gamma
           real(kind=8) :: uu(1:6,1:7), flux(1:6)
           real(kind=8), dimension(1:6) :: rho, u, v, w, e, p, yk
           real(kind=8), dimension(1:6) :: h, Ht, c, ek, vn

           integer :: i, j, s

           real(kind=8) :: omega0, omega1, omega2
           real(kind=8) :: q0, q1, q2, term1, term2

           real(kind=8) :: rho_roe, u_roe, v_roe, w_roe, q_roe
           real(kind=8) :: e_roe, c_roe, p_roe, Ht_roe, yk_roe, h_roe

           real(kind=8), dimension(1:6,1:6) :: Rev, Lev 
           real(kind=8) :: ff(1:6,1:6), cons(1:6,1:6)
           real(kind=8) :: lff(1:6,1:6), fs(1:6,1:6)

           real(kind=8) :: lamtil, ke, twoc2, c2
           real(kind=8), dimension(1:6) :: fluxp, fluxm
           integer :: index
           real(kind=8) :: nx, ny, nz
           real(kind=8), dimension(1:5) :: lam  

           real(kind=8), parameter :: chi = 1.1d0

           ! roe_half = 1 for Roe; = 2 for half
           integer, parameter :: roe_half = 2

!--------------------------------------------

               ! i-2 = 1
               ! i-1 = 2
               ! i   = 3
               ! i+1 = 4
               ! i+2 = 5
               ! i+3 = 6


!--------------------------------------------


            if(index.eq.1) then
              nx = 1.0d0; ny = 0.0d0; nz = 0.0d0
            else if(index.eq.2) then
              nx = 0.0d0; ny = 1.0d0; nz = 0.0d0
            else if(index.eq.3) then
              nx = 0.0d0; ny = 0.0d0; nz = 1.0d0
            else
              write(*,*) 'wrong index ', index
            endif


           rho(:) = uu(:,1)
           u(:) = uu(:,2)
           v(:) = uu(:,3)
           w(:) = uu(:,4)
           e(:) = uu(:,5)
           p(:) = uu(:,6)
           yk(:) = uu(:,7)


           do i = 1, 6
             if(e(i).le.0.0d0) then
               write(*,*) 'e <0 in WENO5_BS ', e, index
               stop 
             endif
           enddo   



           do i = 1, 6
            h(i) = gamma*e(i)
            c(i) = sqrt(gamma*p(i)/rho(i))
            ek(i) = 0.5d0*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
            vn(i) = u(i)*nx + v(i)*ny + w(i)*nz 
            Ht(i) = h(i) + ek(i)
           enddo 





           do i = 1, 6
            ff(i,1) = rho(i)*vn(i)
            ff(i,2) = rho(i)*vn(i)*u(i) + p(i)*nx
            ff(i,3) = rho(i)*vn(i)*v(i) + p(i)*ny
            ff(i,4) = rho(i)*vn(i)*w(i) + p(i)*nz
            ff(i,5) = rho(i)*vn(i)*(e(i) + ek(i)) + p(i)*vn(i)
            ff(i,6) = rho(i)*vn(i)*yk(i)

            cons(i,1) = rho(i)
            cons(i,2) = rho(i)*u(i)
            cons(i,3) = rho(i)*v(i)
            cons(i,4) = rho(i)*w(i)
            cons(i,5) = rho(i)*(e(i) + ek(i))
            cons(i,6) = rho(i)*yk(i)
           enddo


!--------------------------------------------


         if(roe_half.eq.1) then     
 
              ! compute Roe average of i and i+1

           rho_roe = sqrt(rho(3)*rho(4))
           u_roe = (sqrt(rho(3))*u(3) + sqrt(rho(4))*u(4))/(sqrt(rho(3)) + sqrt(rho(4)))
           v_roe = (sqrt(rho(3))*v(3) + sqrt(rho(4))*v(4))/(sqrt(rho(3)) + sqrt(rho(4)))
           w_roe = (sqrt(rho(3))*w(3) + sqrt(rho(4))*w(4))/(sqrt(rho(3)) + sqrt(rho(4)))
           Ht_roe = (sqrt(rho(3))*Ht(3) + sqrt(rho(4))*Ht(4))/(sqrt(rho(3)) + sqrt(rho(4)))
           yk_roe = (sqrt(rho(3))*yk(3) + sqrt(rho(4))*yk(4))/(sqrt(rho(3)) + sqrt(rho(4)))
 
           ke = 0.5d0*(u_roe**2.0d0 + v_roe**2.0d0 + w_roe**2.0d0)
           e_roe = (Ht_roe - ke)/gamma

         else if(roe_half.eq.2) then
 
              ! compute arithmetic average of i and i+1

           rho_roe = 0.5d0*(rho(3) + rho(4))
           u_roe = 0.5d0*(u(3) + u(4))
           v_roe = 0.5d0*(v(3) + v(4)) 
           w_roe = 0.5d0*(w(3) + w(4))
           Ht_roe = 0.5d0*(Ht(3) + Ht(4)) 
           yk_roe = 0.5d0*(yk(3) + yk(4))

           ke = 0.5d0*(u_roe**2.0d0 + v_roe**2.0d0 + w_roe**2.0d0)
           e_roe = (Ht_roe - ke)/gamma
 
         endif 


           p_roe = (gamma-1.0d0)*rho_roe*e_roe
           c_roe = sqrt(gamma*p_roe/rho_roe)

           q_roe = u_roe*nx + v_roe*ny + w_roe*nz 


           if(e_roe.le.0.0d0) then
            write(*,*) 'Roe failed ', e_roe, index
            write(*,*) 'left ', rho(3), u(3), v(3), w(3), e(3)
            write(*,*) 'right ', rho(4), u(4), v(4), w(4), e(4)
            STOP
           endif

           h_roe = gamma*e_roe


!--------------------------------------------

      ! compute right eigenvectors


            Rev(1,1) = 1.0d0
            Rev(2,1) = u_roe - c_roe*nx
            Rev(3,1) = v_roe - c_roe*ny
            Rev(4,1) = w_roe - c_roe*nz
            Rev(5,1) = h_roe + ke - c_roe*q_roe
            Rev(6,1) = yk_roe

            Rev(1,2) = 1.0d0
            Rev(2,2) = u_roe + c_roe*nx
            Rev(3,2) = v_roe + c_roe*ny
            Rev(4,2) = w_roe + c_roe*nz
            Rev(5,2) = h_roe + ke + c_roe*q_roe
            Rev(6,2) = yk_roe   

            Rev(1,3) = 1.0d0
            Rev(2,3) = u_roe 
            Rev(3,3) = v_roe 
            Rev(4,3) = w_roe 
            Rev(5,3) = ke 
            Rev(6,3) = 0.0d0 
            
            Rev(1,4) = 0.0d0
            Rev(2,4) = ny 
            Rev(3,4) = -nx 
            Rev(4,4) = 0.0d0 
            Rev(5,4) = u_roe*ny - v_roe*nx 
            Rev(6,4) = 0.0d0
 
            Rev(1,5) = 0.0d0
            Rev(2,5) = -nz 
            Rev(3,5) = 0.0d0 
            Rev(4,5) = nx 
            Rev(5,5) = w_roe*nx - u_roe*nz 
            Rev(6,5) = 0.0d0

            Rev(1,6) = 0.0d0
            Rev(2,6) = 0.0d0 
            Rev(3,6) = 0.0d0 
            Rev(4,6) = 0.0d0
            Rev(5,6) = 0.0d0
            Rev(6,6) = 1.0d0

!-------------------------------------------------

      ! compute left eigenvectors


            twoc2 = 2.0d0*c_roe*c_roe
            c2 = c_roe*c_roe  


            Lev(1,1) = ((gamma-1.0d0)*ke+c_roe*q_roe)/twoc2
            Lev(1,2) = ((1.0d0-gamma)*u_roe-c_roe*nx)/twoc2
            Lev(1,3) = ((1.0d0-gamma)*v_roe-c_roe*ny)/twoc2
            Lev(1,4) = ((1.0d0-gamma)*w_roe-c_roe*nz)/twoc2
            Lev(1,5) = (gamma-1.0d0)/twoc2
            Lev(1,6) = 0.0d0

            Lev(2,1) = ((gamma-1.0d0)*ke-c_roe*q_roe)/twoc2
            Lev(2,2) = ((1.0d0-gamma)*u_roe+c_roe*nx)/twoc2
            Lev(2,3) = ((1.0d0-gamma)*v_roe+c_roe*ny)/twoc2
            Lev(2,4) = ((1.0d0-gamma)*w_roe+c_roe*nz)/twoc2
            Lev(2,5) = (gamma-1.0d0)/twoc2
            Lev(2,6) = 0.0d0   

            Lev(3,1) = (c_roe*c_roe-(gamma-1.0d0)*ke)/c2
            Lev(3,2) = (gamma-1.0d0)*u_roe/c2
            Lev(3,3) = (gamma-1.0d0)*v_roe/c2
            Lev(3,4) = (gamma-1.0d0)*w_roe/c2
            Lev(3,5) = (1.0d0-gamma)/c2
            Lev(3,6) = 0.0d0 
            
            Lev(4,1) = v_roe*nx - u_roe*ny + u_roe*nz
            Lev(4,2) = ny - nz 
            Lev(4,3) = -nx 
            Lev(4,4) = 0.0d0 
            Lev(4,5) = 0.0d0 
            Lev(4,6) = 0.0d0
 
            Lev(5,1) = -w_roe*nx + w_roe*ny - v_roe*nz
            Lev(5,2) = 0.0d0 
            Lev(5,3) = nz 
            Lev(5,4) = nx - ny 
            Lev(5,5) = 0.0d0
            Lev(5,6) = 0.0d0

            Lev(6,1) = q_roe*q_roe/2.0d0 - h_roe
            Lev(6,2) = -u_roe*(nx+1.0d0)/2.0d0 
            Lev(6,3) = -v_roe*(ny+1.0d0)/2.0d0 
            Lev(6,4) = -w_roe*(nz+1.0d0)/2.0d0
            Lev(6,5) = 1.0d0
            Lev(6,6) = 0.0d0

!-------------------------------------------------

           ! REFER TO Eqn. (2.5) of BALSARA AND SHU, JCP 160, 2000


           flux = 0.0d0

           lam(1) = abs(vn(2)) + c(2) 
           lam(2) = abs(vn(3)) + c(3) 
           lam(3) = abs(vn(4)) + c(4) 
           lam(4) = abs(vn(5)) + c(5) 
           lam(5) = abs(q_roe) + c_roe

           lamtil = chi*maxval(lam(:))

           
           ! s = 1 -> q-c
           ! s = 2 -> q+c
           ! s = 3 -> q
           ! s = 4 -> q
           ! s = 5 -> q
           ! s = 6 -> q

!-----------------------------


           !-------------------------------------------------
           ! FLUX+
           !-------------------------------------------------


          do s = 1, 6

           do i = 1, 5
            do j = 1, 6
              lff(i,j) = Lev(s,j)*0.5d0*(ff(i,j) + lamtil*cons(i,j))
            enddo
           enddo
           
           do j = 1, 6
            call compute_q(lff(1:5,j),q0,q1,q2)
            call compute_omega(lff(1:5,j),omega0,omega1,omega2)
            fs(s,j) = omega0*q0 + omega1*q1 + omega2*q2
           enddo

          enddo

          
           fluxp = 0.0d0
           
           do j = 1, 6
            do s = 1, 6
             fluxp(j) = fluxp(j) + fs(s,j)*Rev(j,s) 
            enddo
           enddo

!-----------------------------


           !-------------------------------------------------
           ! FLUX-
           !-------------------------------------------------



           

         do s = 1, 6

           do i = 1, 5
            do j = 1, 6
              lff(6-i,j) = Lev(s,j)*0.5d0*(ff(i+1,j) - lamtil*cons(i+1,j))
            enddo
           enddo


           do j = 1, 6
            call compute_q(lff(1:5,j),q0,q1,q2)
            call compute_omega(lff(1:5,j),omega0,omega1,omega2)
            fs(s,j) = omega0*q0 + omega1*q1 + omega2*q2
           enddo

         enddo

          
           fluxm = 0.0d0

          do j = 1, 6
           do s = 1, 6
            fluxm(j) = fluxm(j) + fs(s,j)*Rev(j,s) 
           enddo
          enddo

!-----------------------------


           !-------------------
           !  FLUX
           !-------------------


           do j = 1, 6
            flux(j) = fluxp(j) + fluxm(j)
           enddo

!-----------------------------


          return

          end subroutine




!---------------------------------------------------------------------------

        subroutine compute_q(lf,q0,q1,q2)


      implicit none


          real(kind=8) :: lf(1:5), q0, q1, q2

          q0 = 1.0d0/3.0d0*lf(1) - 7.0d0/6.0d0*lf(2) + 11.0d0/6.0d0*lf(3)
          q1 = -1.0d0/6.0d0*lf(2) + 5.0d0/6.0d0*lf(3) + 1.0d0/3.0d0*lf(4)
          q2 = 1.0d0/3.0d0*lf(3) + 5.0d0/6.0d0*lf(4) - 1.0d0/6.0d0*lf(5)


            return



          end subroutine


!----------------------------------------------------------------------------

          subroutine compute_omega(lf,omega0,omega1,omega2)


     implicit none


            real(kind=8) :: lf(1:5), omega0, omega1, omega2
            real(kind=8) :: term1, term2, term3, beta0, beta1, beta2
            real(kind=8) :: alpha0, alpha1, alpha2

            real(kind=8), parameter :: d0 = 1.0d0/10.0d0
            real(kind=8), parameter :: d1 = 6.0d0/10.0d0
            real(kind=8), parameter :: d2 = 3.0d0/10.0d0
            real(kind=8), parameter :: eps = 1.0d-12


            term1 = 13.0d0/12.0d0*(lf(1) - 2.0d0*lf(2) + lf(3))**2.0d0
            term2 = 1.0d0/4.0d0*(lf(1) - 4.0d0*lf(2) + 3.0d0*lf(3))**2.0d0
            beta0 = term1 + term2

            term1 = 13.0d0/12.0d0*(lf(2) - 2.0d0*lf(3) + lf(4))**2.0d0
            term2 = 1.0d0/4.0d0*(lf(2) - lf(4))**2.0d0
            beta1 = term1 + term2

            term1 = 13.0d0/12.0d0*(lf(3) - 2.0d0*lf(4) + lf(5))**2.0d0
            term2 = 1.0d0/4.0d0*(3.0d0*lf(3) - 4.0d0*lf(4) + lf(5))**2.0d0
            beta2 = term1 + term2


            alpha0 = d0/(eps + beta0)**2.0d0
            alpha1 = d1/(eps + beta1)**2.0d0
            alpha2 = d2/(eps + beta2)**2.0d0

            term1 = alpha0 + alpha1 + alpha2
            omega0 = alpha0/term1
            omega1 = alpha1/term1
            omega2 = alpha2/term1


       return


       end subroutine

!---------------------------------------------------------------------         
