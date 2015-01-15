!-----------------------------------------------------------------------------------------------------
!       KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!-----------------------------------------------------------------------------------------------------

       subroutine HLLC_flux(gamma,index,uu,flux) 
 
         implicit none

            real(kind=8) :: gamma, rhol,rhor,ul,ur,vl,vr,wl,wr,el,er,pl,pr         
            real(kind=8) :: ykl, ykr      
            real(kind=8), dimension(1:2,1:7) :: uu
            real(kind=8) :: flux(1:6), cl, cr, kel, ker
            real(kind=8), dimension(1:6) :: fl, fr, flstar, frstar   
            integer :: index  
            real(kind=8) :: nx, ny, nz, ql, qr 
            real(kind=8) :: Sl, Sr, ustar, pstar 
            real(kind=8), dimension(1:6) :: lstar, rstar, consl, consr
            real(kind=8) :: rr, util, vtil, wtil, ketil, ctil, Hl, Hr, Htil 
            real(kind=8) :: qtil, term1, term2, energyl, energyr            
            real(kind=8) :: rholstar, rhorstar   

            
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

           kel = 0.5d0*(ul*ul + vl*vl + wl*wl)
           ker = 0.5d0*(ur*ur + vr*vr + wr*wr)

!-----------------------------------------

               ! WAVE SPEEDS


           rr = sqrt(rhor/rhol) 
           util = (ul + ur*rr)/(1.0d0+rr) 
           vtil = (vl + vr*rr)/(1.0d0+rr) 
           wtil = (wl + wr*rr)/(1.0d0+rr) 

           Hl = el + pl/rhol + kel
           Hr = er + pr/rhor + ker

           Htil = (Hl + Hr*rr)/(1.0d0+rr)
           ketil = 0.5d0*(util**2.0d0 + vtil**2.0d0 + wtil**2.0d0)

           ctil = sqrt((gamma-1.0d0)*(Htil - ketil))
           qtil = util*nx + vtil*ny + wtil*nz
            

           Sl = min(ql-cl,qtil-ctil) 
           Sr = max(qr+cr,qtil+ctil) 



                 if(ctil.le.0.0d0) then
                  write(*,*) 'ctil error in HLLC ', ctil
                  write(*,*) 'left ', rhol, ul, vl, wl, pl 
                  write(*,*) 'right ', rhor, ur, vr, wr, pr 
                  stop 
                 endif 


!-----------------------------------------


            term1 = rhor*qr*(Sr-qr) - rhol*ql*(Sl-ql) + pl-pr
            term2 = rhor*(Sr-qr) - rhol*(Sl-ql)

            ustar = term1/term2

 
            pstar = pl + rhol*(ql-Sl)*(ql-ustar)


           rholstar = rhol*(Sl-ql)/(Sl-ustar)  
           rhorstar = rhor*(Sr-qr)/(Sr-ustar)  

!-----------------------------------------

           lstar(1) = rholstar
           rstar(1) = rhorstar


           term1 = (Sl-ql)*rhol*ul + (pstar-pl)*nx  
           lstar(2) = term1/(Sl-ustar)         
 
           term1 = (Sl-ql)*rhol*vl + (pstar-pl)*ny  
           lstar(3) = term1/(Sl-ustar)         
 
           term1 = (Sl-ql)*rhol*wl + (pstar-pl)*nz  
           lstar(4) = term1/(Sl-ustar)         

           energyl = rhol*(el + kel)
           term1 = (Sl-ql)*energyl -pl*ql + pstar*ustar  
           lstar(5) = term1/(Sl-ustar)         

           term1 = (Sl-ql)*rhol*ykl   
           lstar(6) = term1/(Sl-ustar)         


           term1 = (Sr-qr)*rhor*ur + (pstar-pr)*nx
           rstar(2) = term1/(Sr-ustar)

           term1 = (Sr-qr)*rhor*vr + (pstar-pr)*ny
           rstar(3) = term1/(Sr-ustar)

           term1 = (Sr-qr)*rhor*wr + (pstar-pr)*nz
           rstar(4) = term1/(Sr-ustar)

           energyr = rhor*(er + ker)
           term1 = (Sr-qr)*energyr -pr*qr + pstar*ustar
           rstar(5) = term1/(Sr-ustar)

           term1 = (Sr-qr)*rhor*ykr               
           rstar(6) = term1/(Sr-ustar)


!-----------------------------------------

           consl(1) = rhol 
           consl(2) = rhol*ul 
           consl(3) = rhol*vl 
           consl(4) = rhol*wl 
           consl(5) = rhol*(el+kel) 
           consl(6) = rhol*ykl 

           consr(1) = rhor 
           consr(2) = rhor*ur 
           consr(3) = rhor*vr 
           consr(4) = rhor*wr 
           consr(5) = rhor*(er+ker) 
           consr(6) = rhor*ykr 

!-----------------------------------------

           fl(1) = rhol*ql
           fl(2) = rhol*ql*ul + pl*nx
           fl(3) = rhol*ql*vl + pl*ny
           fl(4) = rhol*ql*wl + pl*nz
           fl(5) = rhol*ql*(el+kel) + pl*ql 
           fl(6) = rhol*ql*ykl  

           fr(1) = rhor*qr
           fr(2) = rhor*qr*ur + pr*nx
           fr(3) = rhor*qr*vr + pr*ny
           fr(4) = rhor*qr*wr + pr*nz
           fr(5) = rhor*qr*(er+ker) + pr*qr 
           fr(6) = rhor*qr*ykr  


!-----------------------------------------

           flstar(:) = fl(:) + Sl*(lstar(:) - consl(:))  
           frstar(:) = fr(:) + Sr*(rstar(:) - consr(:))  

!-----------------------------------------

           ! HLLC FLUX

           if(Sl.gt.0.0d0) then
             flux = fl   
           else if(Sl.le.0.0d0.and.ustar.gt.0.0d0) then   
             flux = flstar
           else if(ustar.le.0.0d0.and.Sr.ge.0.0d0) then   
             flux = frstar 
           else if(Sr.lt.0.0d0) then   
             flux = fr
           else
             write(*,*) 'hllc bug ', Sl, Sr, ustar
             stop
           endif



           return

          end subroutine

!---------------------------------------------------------------------------
