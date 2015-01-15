!-----------------------------------------------------------------------------------------------------
!       KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!-----------------------------------------------------------------------------------------------------

       subroutine rusanov_flux(gamma,index,uu,flux) 
 
         implicit none

            real(kind=8) :: gamma, rhol,rhor,ul,ur,vl,vr,wl,wr,el,er,pl,pr         
            real(kind=8) :: ykl, ykr      
            real(kind=8), dimension(1:2,1:7) :: uu
            real(kind=8) :: flux(1:6), cl, cr, lamtil, kel, ker
            real(kind=8) :: fl, fr, consl, consr   
            integer :: index  
            real(kind=8) :: nx, ny, nz, ql, qr 


            
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



           lamtil = max(abs(ql)+cl,abs(qr)+cr)


           fl = rhol*ql
           fr = rhor*qr
           consl = rhol
           consr = rhor
           flux(1) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           fl = rhol*ql*ul + pl*nx
           fr = rhor*qr*ur + pr*nx
           consl = rhol*ul
           consr = rhor*ur
           flux(2) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           fl = rhol*ql*vl + pl*ny
           fr = rhor*qr*vr + pr*ny
           consl = rhol*vl
           consr = rhor*vr
           flux(3) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           fl = rhol*ql*wl + pl*nz
           fr = rhor*qr*wr + pr*nz
           consl = rhol*wl
           consr = rhor*wr
           flux(4) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           kel = 0.5d0*(ul**2 + vl**2 + wl**2)
           ker = 0.5d0*(ur**2 + vr**2 + wr**2)
           fl = rhol*ql*(el + kel) + pl*ql
           fr = rhor*qr*(er + ker) + pr*qr
           consl = rhol*(el + kel)
           consr = rhor*(er + ker)
           flux(5) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           fl = rhol*ql*ykl
           fr = rhor*qr*ykr
           consl = rhol*ykl
           consr = rhor*ykr
           flux(6) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)


           return

          end subroutine

!---------------------------------------------------------------------------
