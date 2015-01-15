!-----------------------------------------------------------------------------------------------------
!       KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!-----------------------------------------------------------------------------------------------------

          subroutine Gottlieb_Groth_flux(gamma,index,uu,flux) 
 
         implicit none

            integer, parameter :: nprim = 7
 
            real(kind=8) :: gamma, rhol,rhor,ul,ur,vl,vr,wl,wr,el,er,pl,pr         
            real(kind=8) :: ykl, ykr      
            real(kind=8), dimension(1:2,1:7) :: uu
            real(kind=8) :: flux(1:6), cl, cr, lamtil, kel, ker
            real(kind=8) :: fl, fr, consl, consr   
            integer :: index, fail  
            real(kind=8) :: nx, ny, nz, ql, qr 

            real(kind=8), dimension(1:nprim) :: left, right, leftstar, rightstar, godunov
            real(kind=8) :: pstar, ustar, laml, lamr, lamlstar, lamrstar
            real(kind=8) :: gam, gaml, gamr, ke, hl, hr, clstar, crstar
            real(kind=8) :: rholstar, rhorstar, rhoelstar, rhoerstar
            real(kind=8) :: WWl, WWr, term1, alpha, Vshockl, Vshockr
            real(kind=8) :: sigma



            
            if(index.eq.1) then
              nx = 1.0d0; ny = 0.0d0; nz = 0.0d0
            else if(index.eq.2) then
              nx = 0.0d0; ny = 1.0d0; nz = 0.0d0
            else if(index.eq.3) then
              nx = 0.0d0; ny = 0.0d0; nz = 1.0d0
            else
              write(*,*) 'wrong index ', index
            endif  


           rhol = uu(1,1); rhor = uu(2,1)
           ul = uu(1,2); ur = uu(2,2)
           vl = uu(1,3); vr = uu(2,3)
           wl = uu(1,4); wr = uu(2,4)
           el = uu(1,5); er = uu(2,5)
           pl = uu(1,6); pr = uu(2,6)
           ykl = uu(1,7); ykr = uu(2,7)


               ql = ul*nx + vl*ny + wl*nz
               qr = ur*nx + vr*ny + wr*nz


             cl = sqrt(gamma*pl/rhol)
             cr = sqrt(gamma*pr/rhor)



              gaml = cl*cl*rhol/pl

              gamr = cr*cr*rhor/pr

 

                 fail = 0

            ! RIEMANN SOLVER OF GOTTLIEB AND GROTH (JCP, 1988)

        call gottlieb_groth(gaml,gamr,rhol,rhor,ql,qr,pl,pr, &
              ustar,pstar,rholstar,rhorstar,Vshockl,Vshockr,fail)


              if(pstar.le.0.0d0.or.fail.eq.1) then
               !write(*,*) 'pstar bug ', pstar, index
               !write(*,*) 'left ', rhol, pl, gaml, ul, vl, wl, cl  
               !write(*,*) 'right ', rhor, pr, gamr, ur, vr, wr, cr  
               write(*,*) 'Gottlieb-Groth Riemann solver failed '
               stop
              endif 


              clstar = sqrt(gaml*pstar/rholstar)  
              crstar = sqrt(gamr*pstar/rhorstar)  

              hl = el + pl/rhol
              rhoelstar = rhol*el + (pstar-pl)*hl/cl/cl
              hr = er + pr/rhor
              rhoerstar = rhor*er + (pstar-pr)*hr/cr/cr

             
              ! uu is primitive variable
              ! 1 = rho              
              ! 2 = u              
              ! 3 = v              
              ! 4 = w              
              ! 5 = e              
              ! 6 = p              


              ! #5 is now rho*e

              left = uu(1,:)
              left(5) = rhol*el
              right = uu(2,:)
              right(5) = rhor*er


              leftstar(1) = rholstar
              rightstar(1) = rhorstar

              if(index.eq.1) then
               leftstar(2) = ustar
               leftstar(3) = vl
               leftstar(4) = wl

               rightstar(2) = ustar
               rightstar(3) = vr
               rightstar(4) = wr
              else if(index.eq.2) then
               leftstar(2) = ul
               leftstar(3) = ustar
               leftstar(4) = wl

               rightstar(2) = ur
               rightstar(3) = ustar
               rightstar(4) = wr
              else if(index.eq.3) then
               leftstar(2) = ul
               leftstar(3) = vl
               leftstar(4) = ustar

               rightstar(2) = ur
               rightstar(3) = vr
               rightstar(4) = ustar
              else
                write(*,*) 'index wrong ', index
                stop
              endif

              leftstar(5) = rhoelstar
              rightstar(5) = rhoerstar

              leftstar(6) = (gaml-1.0d0)*rhoelstar
              rightstar(6) = (gamr-1.0d0)*rhoerstar

              leftstar(7) = ykl
              rightstar(7) = ykr




              godunov = 0.0d0


              if(ustar.gt.0.0d0) then
                laml = ql - cl
                lamlstar = ustar - clstar
                if(pstar.gt.pl) then
                 ! shock to the left
                 sigma = Vshockl !0.5d0*(laml+lamlstar)
                 if(sigma.gt.0.0d0) then
                  godunov = left
                 else
                  godunov = leftstar
                 endif
                else
                  ! rarefaction to the left
                  if(laml.le.0.0d0.and.lamlstar.le.0.0d0) then
                   godunov = leftstar
                  else if(laml.ge.0.0d0.and.lamlstar.ge.0.0d0) then
                   godunov = left
                  else
                   ! transonic rarefaction
                   alpha = laml/(laml-lamlstar)
                   godunov = alpha*leftstar + (1.0d0-alpha)*left

                      if(laml.eq.lamlstar) then
                       print*, 'transonic rar bug-l ', laml, lamlstar
                      endif

                  endif
                endif
              endif


              if(ustar.lt.0.0d0) then
                lamr = qr + cr
                lamrstar = ustar + crstar
                if(pstar.gt.pr) then
                 ! shock to the right
                 sigma = Vshockr !0.5d0*(lamr+lamrstar)
                 if(sigma.lt.0.0d0) then
                  godunov = right
                 else
                  godunov = rightstar
                 endif
                else
                  ! rarefaction to the right
                  if(lamr.ge.0.0d0.and.lamrstar.ge.0.0d0) then
                   godunov = rightstar
                  else if(lamr.le.0.0d0.and.lamrstar.le.0.0d0) then
                   godunov = right
                  else
                   ! transonic rarefaction
                   alpha = lamr/(lamr-lamrstar)
                   godunov = alpha*rightstar + (1.0d0-alpha)*right

                      if(lamr.eq.lamrstar) then
                       print*, 'transonic rar bug-r ', lamr, lamrstar
                      endif

                  endif
                endif
              endif



              if(ustar.eq.0.0d0) then
               godunov = 0.5d0*(leftstar+rightstar)
              endif


               ! assume gamma is a constant
               gam = 0.5d0*(gaml+gamr)
               godunov(6) = (gam-1.0d0)*godunov(5)

               if(godunov(1).le.0.0d0.or.godunov(6).le.0.0d0) then
                write(*,*) 'godunov bug ', godunov
                write(*,*) 'ustar ', ustar, index
                write(*,*) 'pstar ', pstar, pstar/pl, pstar/pr
                stop
               endif


               ustar = godunov(2)*nx + godunov(3)*ny + godunov(4)*nz

               flux(1) = godunov(1)*ustar

               flux(2) = godunov(1)*ustar*godunov(2) + godunov(6)*nx
               flux(3) = godunov(1)*ustar*godunov(3) + godunov(6)*ny
               flux(4) = godunov(1)*ustar*godunov(4) + godunov(6)*nz

               ke = 0.5d0*(godunov(2)**2.0d0 + godunov(3)**2.0d0 + godunov(4)**2.0d0)
               flux(5) = ustar*(godunov(5) + godunov(1)*ke + godunov(6))

               flux(6) = godunov(1)*ustar*godunov(7)




           return

          end subroutine

!---------------------------------------------------------------------------

        subroutine gottlieb_groth(gaml,gamr,rhol,rhor,ul,ur,pl,pr, &
                   ustar,pstar,rholstar,rhorstar,Vshockl,Vshockr,fail)

          implicit none

             real(kind=8) :: gaml,gamr,ul,ur,pl,pr,ustar,pstar,rhol,rhor
             real(kind=8) :: al, ar, Cl, Cr, Wl, Wr, term1, term2
             real(kind=8) :: sigma, u0, p0, ultil, urtil, alstar, arstar
             real(kind=8) :: dpdul, dpdur, rholstar, rhorstar, plstar, prstar
             integer :: fail, iter
             real(kind=8) :: z, Vshockl, Vshockr


              fail = 0

              Vshockl = 0.0d0
              Vshockr = 0.0d0


              if(gaml.le.0.0d0.or.gamr.le.0.0d0.or.pl.le.0.0d0.or.pr.le.0.0d0 &
                        .or.rhol.le.0.0d0.or.rhor.le.0.0d0) then
               write(*,*) 'gottlieb-groth error ', gaml, gamr, pl, pr, rhol, rhor
               stop
              endif



              al = sqrt(gaml*pl/rhol)
              ar = sqrt(gamr*pr/rhor)

              Cl = gaml*pl/al
              Cr = gamr*pr/ar


              ! initial guess

             ultil = ul + 2.0d0/(gaml-1.0d0)*al
             urtil = ur - 2.0d0/(gamr-1.0d0)*ar

             if(pl.ge.pr) then
               sigma = gaml
             else
               sigma = gamr
             endif


             z = (gaml-1.0d0)/(gamr-1.0d0)*ar/al*((pl/pr)**((sigma-1.0d0)/2.0d0/sigma))

             u0 = (ultil*z + urtil)/(1.0d0 + z)

             term1 = (gaml-1.0d0)/2.0d0*(ultil-urtil)/al/(1.0d0+z)
             p0 = pl*term1**(2.0d0*gaml/(gaml-1.0d0))



             ! iterate

             iter = 0

             ustar = u0
             pstar = p0

385          continue

             iter = iter + 1

             if(ustar.le.ul)then
              ! left moving shock
              term1 = 1.0d0 + ((gaml+1.0d0)/4.0d0*(ustar-ul)/al)**2.0d0
              Wl = (gaml+1.0d0)/4.0d0*(ustar-ul)/al - sqrt(term1)

              plstar = pl + Cl*(ustar-ul)*Wl
              dpdul = 2.0d0*Cl*Wl*Wl*Wl/(1.0d0+Wl*Wl)

                    if(plstar.le.0.0d0) then
                     write(*,*) 'plstar gottlieb-groth ', plstar
                     write(*,*) 'lr', pl, pr, rhol, rhor, ul, ur
                     stop
                    endif 


              term1 = (gaml+1.0d0) + (gaml-1.0d0)*plstar/pl
              term2 = (gaml+1.0d0) + (gaml-1.0d0)*pl/plstar
              alstar = al*sqrt(term1/term2)
              Vshockl = ul + al*Wl

             else
              ! left moving rarefaction
              alstar = al - (gaml-1.0d0)/2.0d0*(ustar-ul)
              plstar = pl*((alstar/al)**(2.0d0*gaml/(gaml-1.0d0)))
              dpdul = -gaml*plstar/alstar
             endif

             if(ustar.ge.ur)then
              ! right moving shock
              term1 = 1.0d0 + ((gamr+1.0d0)/4.0d0*(ustar-ur)/ar)**2.0d0
              Wr = (gamr+1.0d0)/4.0d0*(ustar-ur)/ar + sqrt(term1)

              prstar = pr + Cr*(ustar-ur)*Wr
              dpdur = 2.0d0*Cr*Wr*Wr*Wr/(1.0d0+Wr*Wr)

                    if(prstar.le.0.0d0) then
                     write(*,*) 'prstar gottlieb-groth ', prstar
                     write(*,*) 'lr', pl, pr, rhol, rhor, ul, ur
                     stop
                    endif

              term1 = (gamr+1.0d0) + (gamr-1.0d0)*prstar/pr
              term2 = (gamr+1.0d0) + (gamr-1.0d0)*pr/prstar
              arstar = ar*sqrt(term1/term2)
              Vshockr = ur + ar*Wr

             else
              ! right moving rarefaction
              arstar = ar + (gamr-1.0d0)/2.0d0*(ustar-ur)
              prstar = pr*((arstar/ar)**(2.0d0*gamr/(gamr-1.0d0)))
              dpdur = gamr*prstar/arstar
             endif

              if(dpdul.ne.dpdur) then
               ustar = ustar - (plstar-prstar)/(dpdul-dpdur)
              else
               if(abs(1.0d0-plstar/prstar).gt.1.0d-6) then
                fail = 1
               else
                pstar = 0.5d0*(plstar+prstar)
                rholstar = gaml*pstar/alstar/alstar
                rhorstar = gamr*pstar/arstar/arstar
               endif
               return
              endif


               if(iter.ge.100) then
                ! iteration did not converge
                fail = 1
                return
               endif


               if(abs(1.0d0-plstar/prstar).gt.1.0d-6) then
                goto 385
               endif


               pstar = 0.5d0*(plstar+prstar)
               rholstar = gaml*pstar/alstar/alstar
               rhorstar = gamr*pstar/arstar/arstar



               return
               end subroutine

!---------------------------------------------------------------------------



