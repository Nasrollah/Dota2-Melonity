!-----------------------------------------------------------------------------------------------------
!       KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!-----------------------------------------------------------------------------------------------------

               subroutine WENO_cons_flux(gamma,index,uu,flux) 
 
         implicit none

            integer, parameter :: nprim = 7

            real(kind=8) :: gamma, rhol,rhor,ul,ur,vl,vr,wl,wr,el,er,pl,pr         
            real(kind=8) :: ykl, ykr      
            real(kind=8) :: uu(1:6,1:nprim), uiphm(1:nprim), uiphp(1:nprim)
            real(kind=8) :: uuu(1:2,1:nprim)
            real(kind=8) :: flux(1:6), cl, cr, kel, ker
            integer :: index, fail
            real(kind=8) :: nx, ny, nz, ql, qr 
            real(kind=8), dimension(1:nprim) :: left, right, leftstar, rightstar, godunov  
            real(kind=8) :: pstar, ustar, laml, lamr, lamlstar, lamrstar
            real(kind=8) :: gam, gaml, gamr, ke, hl, hr, clstar, crstar
            real(kind=8) :: rholstar, rhorstar, rhoelstar, rhoerstar 
            real(kind=8) :: WWl, WWr, term1, alpha, Vshockl, Vshockr
            real(kind=8) :: sigma

           real(kind=8) :: dlam1, dlam2, dlam3, nuvisc, dcons(1:6) 
           real(kind=8), parameter :: cvisc = 0.0d0  
           integer, parameter :: iriemann = 3  


             ! iriemann = 1 for Gottlieb-Groth Riemann solver
             ! iriemann = 2 for Rusanov Riemann solver (very dissipative)
             ! iriemann = 3 for HLLC Riemann solver
 
  
 


             
              ! uu is primitive variable
              ! 1 = rho              
              ! 2 = u              
              ! 3 = v              
              ! 4 = w              
              ! 5 = e              
              ! 6 = p              
              ! 7 = Y              


!---------------------------------------------------------------

            call weno5_gas(gamma,nprim,uu,uiphm,uiphp)
            


            if(index.eq.1) then
              nx = 1.0d0; ny = 0.0d0; nz = 0.0d0
            else if(index.eq.2) then
              nx = 0.0d0; ny = 1.0d0; nz = 0.0d0
            else if(index.eq.3) then
              nx = 0.0d0; ny = 0.0d0; nz = 1.0d0
            else
              write(*,*) 'wrong index ', index
            endif  


           rhol = uiphm(1); rhor = uiphp(1)
           ul = uiphm(2); ur = uiphp(2)
           vl = uiphm(3); vr = uiphp(3)
           wl = uiphm(4); wr = uiphp(4)
           !el = uiphm(5); er = uiphp(5)
           pl = uiphm(6); pr = uiphp(6)
           ykl = uiphm(7); ykr = uiphp(7)


           ! interpolation may have violated EOS; ensure energy is conserved 
           el = pl/rhol/(gamma-1.0d0)
           er = pr/rhor/(gamma-1.0d0)
           uiphm(5) = el
           uiphp(5) = er


!-------------------------------------------------

              uuu(1,1:nprim) = uiphm
              uuu(2,1:nprim) = uiphp

!-------------------------------------------------


              if(rhol.le.0.0d0.or.rhor.le.0.0d0) then
               write(*,*) 'rho bug ', rhol, rhor
               stop
              endif   
              if(pl.le.0.0d0.or.pr.le.0.0d0) then
               write(*,*) 'p bug ', pl, pr
              ! use first order + Rusanov here
              uuu(1,1:nprim) = uu(3,:)
              uuu(2,1:nprim) = uu(4,:)
              call rusanov_flux(gamma,index,uuu,flux)
              return 
              print*, 'should not be here '
               stop
              endif   




               ql = ul*nx + vl*ny + wl*nz
               qr = ur*nx + vr*ny + wr*nz


             cl = sqrt(gamma*pl/rhol)
             cr = sqrt(gamma*pr/rhor)



!-------------------------------------------------
             dlam1 = ql-qr 
             dlam2 = (ql-cl)-(qr-cr) 
             dlam3 = (ql+cl)-(qr+cr) 
             nuvisc = cvisc*max(dlam1,dlam2,dlam3,0.0d0)
 
             dcons(1) = rhor - rhol
             dcons(2) = rhor*ur - rhol*ul
             dcons(3) = rhor*vr - rhol*vl
             dcons(4) = rhor*wr - rhol*wl
           kel = 0.5d0*(ul**2 + vl**2 + wl**2)
           ker = 0.5d0*(ur**2 + vr**2 + wr**2)
             dcons(5) = rhor*(er+ker) - rhol*(el+kel)
             dcons(6) = rhor*ykr - rhol*ykl

!-------------------------------------------------


           if(iriemann.eq.1) then

                   ! GOTTLIEB-GROTH RIEMANN SOLVER 

           call Gottlieb_Groth_flux(gamma,index,uuu,flux)

           endif



!-------------------------------------------------

           if(iriemann.eq.2) then              

                   ! RUSANOV FLUX 

           call rusanov_flux(gamma,index,uuu,flux)

           endif 

!-------------------------------------------------

            if(iriemann.eq.3) then

                   ! HLLC FLUX

           call HLLC_flux(gamma,index,uuu,flux)

           endif

!-------------------------------------------------


                         ! artificial viscosity
                    flux(1:6) = flux(1:6) - nuvisc*dcons(1:6)

!-------------------------------------------------


           return

          end subroutine

!---------------------------------------------------------------------------

            subroutine weno5_gas(gamma,nprim,uu,uiphm,uiphp)

           implicit none

           integer :: nprim  
           real(kind=8) :: uu(1:6,1:nprim), uiphm(1:nprim), uiphp(1:nprim)
           real(kind=8) :: d0, d1, d2, dt2, dt1, dt0
           real(kind=8) :: beta0, beta1, beta2, alpha0, alpha1, alpha2, & 
                               omega0, omega1, omega2
           real(kind=8) :: alphat0, alphat1, alphat2
           real(kind=8) :: sum1, p0, p1, p2, eps
           integer :: i, l

           real(kind=8) :: phi, phi1, phi2, ke, gamma


              ! uu is primitive variable
              ! 1 = rho              
              ! 2 = u              
              ! 3 = v              
              ! 4 = w              
              ! 5 = e              
              ! 6 = p              
              ! 7 = Y              


!---------------------------------------------------              
                  ! first order
                    !uiphm = uu(3,:)
                    !uiphp = uu(4,:)
                     !return
                     !write(*,*) 'should not be here in weno'
                     !stop
!---------------------------------------------------              



                  d0 = 3.0d0/10.0d0
                  d1 = 6.0d0/10.0d0
                  d2 = 1.0d0/10.0d0
                  dt2 = 3.0d0/10.0d0
                  dt1 = 6.0d0/10.0d0
                  dt0 = 1.0d0/10.0d0

                  eps = 1.0d-12


                  uiphm = 0.0d0
                  uiphp = 0.0d0

!-----------------------

              ! compute uu(i+1/2)-

      do l = 1, nprim

                  i = 3
!             if(uu(i,l).ne.uu(i-1,l)) then                 

           beta0 = 13.0d0/12.0d0*(uu(i,l) - 2.0d0*uu(i+1,l) + &
            uu(i+2,l))**2.0d0 + 1.0d0/4.0d0*(3.0d0*uu(i,l) - 4.0d0*uu(i+1,l) + uu(i+2,l))**2.0d0

           beta1 = 13.0d0/12.0d0*(uu(i-1,l) - 2.0d0*uu(i,l) + &
                uu(i+1,l))**2.0d0 + 1.0d0/4.0d0*(uu(i-1,l) - uu(i+1,l))**2.0d0

           beta2 = 13.0d0/12.0d0*(uu(i-2,l) - 2.0d0*uu(i-1,l) + uu(i,l))**2.0d0 + &
                 1.0d0/4.0d0*(uu(i-2,l) - 4.0d0*uu(i-1,l) + 3.0d0*uu(i,l))**2.0d0

           alpha0 = d0/(eps + beta0)**2.0d0
           alpha1 = d1/(eps + beta1)**2.0d0
           alpha2 = d2/(eps + beta2)**2.0d0

           sum1 = alpha0 + alpha1 + alpha2

           omega0 = alpha0/sum1
           omega1 = alpha1/sum1
           omega2 = alpha2/sum1


           p0 = 1.0d0/3.0d0*uu(i,l) + 5.0d0/6.0d0*uu(i+1,l) - 1.0d0/6.0d0*uu(i+2,l)
           p1 = -1.0d0/6.0d0*uu(i-1,l) + 5.0d0/6.0d0*uu(i,l) + 1.0d0/3.0d0*uu(i+1,l)
           p2 = 1.0d0/3.0d0*uu(i-2,l) - 7.0d0/6.0d0*uu(i-1,l) + 11.0d0/6.0d0*uu(i,l)


           uiphm(l) = omega0*p0 + omega1*p1 + omega2*p2


                      ! slope limiter

!               phi1 = (uu(i+1,l)-uu(i,l))/(uu(i,l)-uu(i-1,l))
!               phi2 = (uiphm(l)-uu(i,l))/(uu(i,l)-uu(i-1,l))       
!               phi = max(0.0d0,min(2.0d0,2.0d0*phi1,2.0d0*phi2))
                                       
!               uiphm(l) = uu(i,l) + 0.5d0*(uu(i,l) - uu(i-1,l))*phi
!             else
!               uiphm(l) = uu(i,l)
!             endif

      enddo

!-----------------------

                  ! compute uu(i+1/2)+


      do l = 1, nprim

                  i = 4
!             if(uu(i,l).ne.uu(i+1,l)) then               

           beta0 = 13.0d0/12.0d0*(uu(i,l) - 2.0d0*uu(i+1,l) + &
            uu(i+2,l))**2.0d0 + 1.0d0/4.0d0*(3.0d0*uu(i,l) - 4.0d0*uu(i+1,l) + uu(i+2,l))**2.0d0

           beta1 = 13.0d0/12.0d0*(uu(i-1,l) - 2.0d0*uu(i,l) + &
                uu(i+1,l))**2.0d0 + 1.0d0/4.0d0*(uu(i-1,l) - uu(i+1,l))**2.0d0

           beta2 = 13.0d0/12.0d0*(uu(i-2,l) - 2.0d0*uu(i-1,l) + uu(i,l))**2.0d0 + &
                 1.0d0/4.0d0*(uu(i-2,l) - 4.0d0*uu(i-1,l) + 3.0d0*uu(i,l))**2.0d0

           alphat0 = dt0/(eps + beta0)**2.0d0
           alphat1 = dt1/(eps + beta1)**2.0d0
           alphat2 = dt2/(eps + beta2)**2.0d0

           sum1 = alphat0 + alphat1 + alphat2

           omega0 = alphat0/sum1
           omega1 = alphat1/sum1
           omega2 = alphat2/sum1

           p0 = 11.0d0/6.0d0*uu(i,l) - 7.0d0/6.0d0*uu(i+1,l) + 1.0d0/3.0d0*uu(i+2,l)
           p1 = 1.0d0/3.0d0*uu(i-1,l) + 5.0d0/6.0d0*uu(i,l) - 1.0d0/6.0d0*uu(i+1,l)
           p2 = -1.0d0/6.0d0*uu(i-2,l) + 5.0d0/6.0d0*uu(i-1,l) + 1.0d0/3.0d0*uu(i,l)


           uiphp(l) = omega0*p0 + omega1*p1 + omega2*p2




                      ! slope limiter
 
 
!               phi1 = (uu(i-1,l)-uu(i,l))/(uu(i,l)-uu(i+1,l))
!               phi2 = (uiphp(l)-uu(i,l))/(uu(i,l)-uu(i+1,l))        
!               phi = max(0.0d0,min(2.0d0,2.0d0*phi1,2.0d0*phi2))
                        
!               uiphp(l) = uu(i,l) - 0.5d0*(uu(i+1,l) - uu(i,l))*phi
!             else
!               uiphp(l) = uu(i,l)
!             endif 
                      

      enddo

!-----------------------

              ! uu is primitive variable
              ! 1 = rho              
              ! 2 = u              
              ! 3 = v              
              ! 4 = w              
              ! 5 = e              
              ! 6 = p              
              ! 7 = Y              


           uiphm(5) = uiphm(6)/uiphm(1)/(gamma-1.0d0)
           uiphp(5) = uiphp(6)/uiphp(1)/(gamma-1.0d0)

            
           return
           end subroutine

!---------------------------------------------------------------------------
