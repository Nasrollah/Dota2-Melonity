                   subroutine central_flux(gamma,index,uu,flux)

          implicit none

          real(kind=8) :: gamma
          real(kind=8) :: uu(1:6,1:7), flux(1:6), ff(1:6,1:6) 
          integer :: i, index
          real(kind=8), dimension(1:6) :: rho, u, v, w, e, Yk
          real(kind=8), dimension(1:6) :: p, Et, q
          real(kind=8) :: ke
          real(kind=8) :: nx, ny, nz
          real(kind=8) :: ql, qr
          real(kind=8) :: rhol, rhor, ul, ur, vl, vr, wl, wr
          real(kind=8) :: cl, cr, pl, pr, el, er

          integer, parameter :: order = 4




            if(index.eq.1) then
              nx = 1.0d0; ny = 0.0d0; nz = 0.0d0
            else if(index.eq.2) then
              nx = 0.0d0; ny = 1.0d0; nz = 0.0d0
            else if(index.eq.3) then
              nx = 0.0d0; ny = 0.0d0; nz = 1.0d0
            else
              write(*,*) 'wrong index ', index
            endif



           rho(1:6) = uu(1:6,1)
           u(1:6) = uu(1:6,2)
           v(1:6) = uu(1:6,3)             
           w(1:6) = uu(1:6,4)
           e(1:6) = uu(1:6,5)
           P(1:6) = uu(1:6,6)  
           Yk(1:6) = uu(1:6,7)
                        



!-------------------------------------------------


           do i = 1, 6
             q(i) = u(i)*nx + v(i)*ny + w(i)*nz
             ke = 0.5d0*(u(i)**2.0d0 + v(i)**2.0d0 + w(i)**2.0d0)
             Et(i) = rho(i)*(e(i) + ke)
           enddo    


           do i = 1, 6
             ff(i,1) = rho(i)*q(i)
             ff(i,2) = rho(i)*q(i)*u(i) + p(i)*nx
             ff(i,3) = rho(i)*q(i)*v(i) + p(i)*ny
             ff(i,4) = rho(i)*q(i)*w(i) + p(i)*nz
             ff(i,5) = q(i)*(Et(i) + p(i))
             ff(i,6) = rho(i)*q(i)*Yk(i)
           enddo 


           ! i-2 = 1 
           ! i-1 = 2 
           ! i = 3 
           ! i+1 = 4 
           ! i+2 = 5 
           ! i+3 = 6 

            i = 3 

             if(order.eq.6) then

              flux(:) = (ff(i-2,:) - 8.0d0*ff(i-1,:) + 37.0d0*ff(i,:) & 
              + 37.0d0*ff(i+1,:) - 8.0d0*ff(i+2,:) + ff(i+3,:))/60.0d0               

             else if(order.eq.4) then

              flux(:) = (-ff(i-1,:) + 7.0d0*ff(i,:) + & 
                    7.0d0*ff(i+1,:) - ff(i+2,:))/12.0d0  

             else
               write(*,*) 'wrong order ', order
             endif 




                    return
                    end subroutine
