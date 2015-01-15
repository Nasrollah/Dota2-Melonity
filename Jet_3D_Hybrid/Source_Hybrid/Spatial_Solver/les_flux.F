!------------------------------------------------------------------
!                  KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!------------------------------------------------------------------


       subroutine les_smagorinsky(gamma,dx,dy,dz, & 
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,Q,U,V,W,P,T,dflux)

      implicit none

        integer :: N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, dflux
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U, V, W, P, T, H, RHOH
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3,1:3)    :: tau_sgs
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)    :: tau_u, q_sgs
        real(kind=8) :: uu(N0:N5,M0:M5,K0:K5,1:3) 
        real(kind=8) :: Sbar, Skk
        real(kind=8), dimension(N0:N5,M0:M5,K0:K5,1:3,1:3) :: S
        real(kind=8), dimension(1:3,1:3) :: S1, S2

        real(kind=8) :: dx, dy, dz, gamma
        real(kind=8) :: taukk, taux, tauy, tauz
        integer :: ii, jj, kk, l, m  
        real(kind=8) :: hx, hy, hz, hfact, qx, qy, qz 
        
        real(kind=8), parameter :: CS = 0.16d0 
        real(kind=8), parameter :: CI = 0.09d0 
        real(kind=8), parameter :: Pr_sgs = 0.7d0 

        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: CS_smag, CI_smag, Pr_smag
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:NV)    :: Q_test
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: U_test, V_test, W_test, rhoh_test, H_test
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: e_test
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)    :: grad_h_test, grad_h

        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3,1:3)    :: mat_ij, mat1_ij, Lij, Mij
        REALTYPE, dimension(1:3,1:3)    :: mat2_ij

        real(kind=8) :: term, term1, term2
        real(kind=8) :: ke, delbar, del_hat

        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: alpha, beta, alpha_hat
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3,1:3)    :: alpha_ij, beta_ij, alpha_ij_hat

        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)    :: kappa_j, theta_j
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)    :: mat_pr_1, mat_pr_2, mat_pr_3
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)    :: theta1, theta2

        real(kind=8) :: dynsmag1, dynsmag2, dynsmag3, nentries_ds 

        integer, parameter :: les_model = 2
        ! les_model = 1 for constant coefficient Smagorinsky  
        ! les_model = 2 for dynamic Smagorinsky  


!-------------------------------------------------
          
         ! Initial values

         CS_smag = CS
         CI_smag = CI
         Pr_smag = Pr_sgs
         Q_test = Q 

         mat_ij = 0.0d0
         mat1_ij = 0.0d0
         mat2_ij = 0.0d0
         Lij = 0.0d0 
         Mij = 0.0d0

         alpha = 0.0d0
         beta = 0.0d0
         alpha_hat = 0.0d0 

         alpha_ij = 0.0d0
         beta_ij = 0.0d0
         alpha_ij_hat = 0.0d0 

         mat_pr_1 = 0.0d0 
         mat_pr_2 = 0.0d0 

         kappa_j = 0.0d0
         theta_j = 0.0d0 
         theta1 = 0.0d0
         theta2 = 0.0d0

!-------------------------------------------------

         delbar = sqrt(dx*dx + dy*dy + dz*dz)            
         del_hat = 2.0d0*delbar
 

         dflux = 0.0d0
         tau_sgs = 0.0d0 
         q_sgs = 0.0d0
         S = 0.0d0


!-------------------------------------------------

         ! Compute enthalpy
 

         do ii = N0, N5
          do jj = M0, M5
           do kk = K0, K5
             H(ii,jj,kk) = gamma/(gamma-1.0d0)*P(ii,jj,kk)/Q(ii,jj,kk,1) 
             RHOH(ii,jj,kk) = Q(ii,jj,kk,1)*H(ii,jj,kk)
           enddo 
          enddo 
         enddo 

!-------------------------------------------------

          ! DYNAMIC SMAGORINSKY  
 

            if(les_model.eq.2) then

          ! Compute test filter quantities 

         do ii = N0, N5
          do jj = M0, M5
           do kk = K0, K5

              do l = 1, 3
              do m = 1, 3 
               mat_ij(ii,jj,kk,l,m) = Q(ii,jj,kk,1+l)*Q(ii,jj,kk,1+m)/Q(ii,jj,kk,1)
              enddo 
              enddo 

              do l = 1, 3
               mat_pr_1(ii,jj,kk,l) = Q(ii,jj,kk,1+l)*RHOH(ii,jj,kk)/Q(ii,jj,kk,1)
              enddo
    
           enddo 
          enddo 
         enddo 


        
                 
        do l = 1, NV
         call filter(N0,N5,M0,M5,K0,K5,Q(:,:,:,l),Q_test(:,:,:,l))
        enddo

         call filter(N0,N5,M0,M5,K0,K5,RHOH(:,:,:),rhoh_test(:,:,:))
         call filter(N0,N5,M0,M5,K0,K5,H(:,:,:),h_test(:,:,:))

        do l = 1, 3 
         call filter(N0,N5,M0,M5,K0,K5,mat_pr_1(:,:,:,l),mat_pr_2(:,:,:,l))
        enddo

        do l = 1, 3
         do m = 1, 3
          call filter(N0,N5,M0,M5,K0,K5,mat_ij(:,:,:,l,m),mat1_ij(:,:,:,l,m))
         enddo
        enddo


         do ii = N0+1, N5-1
          do jj = M0+1, M5-1
           do kk = K0+1, K5-1

             do l = 1, 3 
              kappa_j(ii,jj,kk,l) = mat_pr_2(ii,jj,kk,l) - Q_test(ii,jj,kk,1+l)*rhoh_test(ii,jj,kk)/Q_test(ii,jj,kk,1)
             enddo

           enddo 
          enddo 
         enddo 


            Q_test(N0,:,:,:) = Q_test(N0+1,:,:,:)
            Q_test(N5,:,:,:) = Q_test(N5-1,:,:,:)
            Q_test(:,M0,:,:) = Q_test(:,M0+1,:,:)
            Q_test(:,M5,:,:) = Q_test(:,M5-1,:,:)
            Q_test(:,:,K0,:) = Q_test(:,:,K0+1,:)
            Q_test(:,:,K5,:) = Q_test(:,:,K5-1,:)

            h_test(N0,:,:) = h_test(N0+1,:,:)
            h_test(N5,:,:) = h_test(N5-1,:,:)
            h_test(:,M0,:) = h_test(:,M0+1,:)
            h_test(:,M5,:) = h_test(:,M5-1,:)
            h_test(:,:,K0) = h_test(:,:,K0+1)
            h_test(:,:,K5) = h_test(:,:,K5-1)




         do ii = N0, N5
          do jj = M0, M5
           do kk = K0, K5
             U_test(ii,jj,kk) = Q_test(ii,jj,kk,2)/Q_test(ii,jj,kk,1)
             V_test(ii,jj,kk) = Q_test(ii,jj,kk,3)/Q_test(ii,jj,kk,1)
             W_test(ii,jj,kk) = Q_test(ii,jj,kk,4)/Q_test(ii,jj,kk,1)
           enddo 
          enddo 
         enddo 
            

         do ii = N0+2, N5-2
          do jj = M0+2, M5-2
           do kk = K0+2, K5-2

            grad_h_test(ii,jj,kk,1) = (h_test(ii-2,jj,kk) - 8.0d0*h_test(ii-1,jj,kk) + & 
                          8.0d0*h_test(ii+1,jj,kk) - h_test(ii+2,jj,kk))/(12.0d0*dx)

            grad_h_test(ii,jj,kk,2) = (h_test(ii,jj-2,kk) - 8.0d0*h_test(ii,jj-1,kk) + & 
                          8.0d0*h_test(ii,jj+1,kk) - h_test(ii,jj+2,kk))/(12.0d0*dy)

            grad_h_test(ii,jj,kk,3) = (h_test(ii,jj,kk-2) - 8.0d0*h_test(ii,jj,kk-1) + & 
                          8.0d0*h_test(ii,jj,kk+1) - h_test(ii,jj,kk+2))/(12.0d0*dz)

            grad_h(ii,jj,kk,1) = (H(ii-2,jj,kk) - 8.0d0*H(ii-1,jj,kk) + & 
                          8.0d0*H(ii+1,jj,kk) - H(ii+2,jj,kk))/(12.0d0*dx) 

            grad_h(ii,jj,kk,2) = (H(ii,jj-2,kk) - 8.0d0*H(ii,jj-1,kk) + & 
                          8.0d0*H(ii,jj+1,kk) - H(ii,jj+2,kk))/(12.0d0*dy) 

            grad_h(ii,jj,kk,3) = (H(ii,jj,kk-2) - 8.0d0*H(ii,jj,kk-1) + & 
                          8.0d0*H(ii,jj,kk+1) - H(ii,jj,kk+2))/(12.0d0*dz) 

           enddo 
          enddo 
         enddo 

!--------------

         ! Compute alpha_ij using grid-filtered quantities

        uu(:,:,:,1) = U(:,:,:)  
        uu(:,:,:,2) = V(:,:,:)  
        uu(:,:,:,3) = W(:,:,:)  

        S = 0.0d0  
        call compute_Sij(N0,N5,M0,M5,K0,K5,N0+2,N5-2,M0+2,M5-2,K0+2,K5-2,dx,dy,dz,uu,S)


!----------------------


         do ii = N0+2, N5-2
          do jj = M0+2, M5-2
           do kk = K0+2, K5-2

             Sbar = 0.0d0

             do l = 1, 3
              do m = 1, 3
               Sbar = Sbar + 2.0d0*S(ii,jj,kk,l,m)*S(ii,jj,kk,l,m)
              enddo
             enddo
             Sbar = sqrt(Sbar)

             Skk = S(ii,jj,kk,1,1) + S(ii,jj,kk,2,2) + S(ii,jj,kk,3,3)

             alpha(ii,jj,kk) = 2.0d0*(delbar**2.0d0)*Q(ii,jj,kk,1)*(Sbar**2.0d0)


             do l = 1, 3
             do m = 1, 3

              term = S(ii,jj,kk,l,m)
              if(l.eq.m) then
               term = term - Skk/3.0d0
              endif

              alpha_ij(ii,jj,kk,l,m) = -2.0d0*(delbar**2.0d0)*Q(ii,jj,kk,1)*Sbar*term
             enddo
             enddo

             do l = 1, 3
              mat_pr_3(ii,jj,kk,l) = (delbar**2.0d0)*Q(ii,jj,kk,1)*Sbar*grad_h(ii,jj,kk,l)
             enddo

           enddo
          enddo
         enddo

!--------------
 
         call filter(N0+2,N5-2,M0+2,M5-2,K0+2,K5-2,alpha(:,:,:),alpha_hat(:,:,:))

         do l = 1, 3
          do m = 1, 3
           call filter(N0+2,N5-2,M0+2,M5-2,K0+2,K5-2,alpha_ij(:,:,:,l,m),alpha_ij_hat(:,:,:,l,m))
          enddo
         enddo

         do l = 1, 3
          call filter(N0+2,N5-2,M0+2,M5-2,K0+2,K5-2,mat_pr_3(:,:,:,l),theta2(:,:,:,l))
         enddo


!--------------

         ! Compute beta_ij using test-filtered quantities

        uu(:,:,:,1) = U_test(:,:,:)  
        uu(:,:,:,2) = V_test(:,:,:)  
        uu(:,:,:,3) = W_test(:,:,:)  

        S = 0.0d0    
        call compute_Sij(N0,N5,M0,M5,K0,K5,N0+2,N5-2,M0+2,M5-2,K0+2,K5-2,dx,dy,dz,uu,S)


         do ii = N0+2, N5-2
          do jj = M0+2, M5-2
           do kk = K0+2, K5-2

             Sbar = 0.0d0
             do l = 1, 3
              do m = 1, 3
               Sbar = Sbar + 2.0d0*S(ii,jj,kk,l,m)*S(ii,jj,kk,l,m)
              enddo
             enddo 
             Sbar = sqrt(Sbar)   
             
             Skk = S(ii,jj,kk,1,1) + S(ii,jj,kk,2,2) + S(ii,jj,kk,3,3)  

             beta(ii,jj,kk) = 2.0d0*(del_hat**2.0d0)*Q_test(ii,jj,kk,1)*(Sbar**2.0d0)


             do l = 1, 3
             do m = 1, 3
        
              mat2_ij(l,m) = Q_test(ii,jj,kk,1+l)*Q_test(ii,jj,kk,1+m)/Q_test(ii,jj,kk,1)    
              Lij(ii,jj,kk,l,m) = mat1_ij(ii,jj,kk,l,m) - mat2_ij(l,m)

              term = S(ii,jj,kk,l,m)
              if(l.eq.m) then
               term = term - Skk/3.0d0 
              endif 

              beta_ij(ii,jj,kk,l,m) = -2.0d0*(del_hat**2.0d0)*Q_test(ii,jj,kk,1)*Sbar*term
              Mij(ii,jj,kk,l,m) = beta_ij(ii,jj,kk,l,m) - alpha_ij_hat(ii,jj,kk,l,m)
             enddo 
             enddo 


             do l = 1, 3
              theta1(ii,jj,kk,l) = -(del_hat**2.0d0)*Q_test(ii,jj,kk,1)*Sbar*grad_h_test(ii,jj,kk,l)
              theta_j(ii,jj,kk,l) = theta1(ii,jj,kk,l) + theta2(ii,jj,kk,l)  
             enddo

             
           enddo 
          enddo 
         enddo 


!------------------

         ! Compute dynamic Smagorinsky coefficients

         do ii = N2-3, N3+3 
          do jj = M2-3, M3+3
           do kk = K2-3, K3+3

             term1 = 0.0d0
             term2 = 0.0d0
             do l = 1, 3
             do m = 1, 3
              term1 = term1 + Lij(ii,jj,kk,l,m)*Mij(ii,jj,kk,l,m)
              term2 = term2 + Mij(ii,jj,kk,l,m)*Mij(ii,jj,kk,l,m)
             enddo 
             enddo 
             
             if(abs(term2).gt.1.0d-6) then
              if(term1/term2.gt.0.0d0) then
               CS_smag(ii,jj,kk) = sqrt(term1/term2)
              else
               CS_smag(ii,jj,kk) = CS
              endif 
             else
              CS_smag(ii,jj,kk) = CS
             endif

             term1 = Lij(ii,jj,kk,1,1) + Lij(ii,jj,kk,2,2) + Lij(ii,jj,kk,3,3) 
             term2 = beta(ii,jj,kk) - alpha_hat(ii,jj,kk)
             if(abs(term2).gt.1.0d-6) then
               CI_smag(ii,jj,kk) = term1/term2
             else
              CI_smag(ii,jj,kk) = CI
             endif 

             term1 = 0.0d0
             term2 = 0.0d0
             do l = 1, 3
              term1 = term1 + theta_j(ii,jj,kk,l)**2.0d0
              term2 = term2 + theta_j(ii,jj,kk,l)*kappa_j(ii,jj,kk,l)
             enddo 
             if(abs(term2).gt.1.0d-6) then
               Pr_smag(ii,jj,kk) = (CS_smag(ii,jj,kk)**2.0d0)*term1/term2
             else
              Pr_smag(ii,jj,kk) = Pr_sgs
             endif  


             ! WE DON'T WANT VALUES WAY OFF FROM "OPTIMUM" VALUES 
             if(CS_smag(ii,jj,kk).gt.2.0d0*CS.or.CS_smag(ii,jj,kk).lt.0.0d0) then
              CS_smag(ii,jj,kk) = CS
             endif
             if(CI_smag(ii,jj,kk).gt.2.0d0*CI.or.CI_smag(ii,jj,kk).lt.0.0d0) then
              CI_smag(ii,jj,kk) = CI
             endif
             if(Pr_smag(ii,jj,kk).gt.2.0d0*Pr_sgs.or.Pr_smag(ii,jj,kk).lt.0.0d0) then
              Pr_smag(ii,jj,kk) = Pr_sgs
             endif


           enddo 
          enddo 
         enddo 


!----------------------------------------------

             ! AVERAGE LOCALLY

        do kk = K2-3, K3+3
 
         dynsmag1 = 0.0d0          
         dynsmag2 = 0.0d0          
         dynsmag3 = 0.0d0          
         nentries_ds = 0.0d0

         do ii = N2-3, N3+3 
          do jj = M2-3, M3+3

            dynsmag1 = dynsmag1 + CS_smag(ii,jj,kk)
            dynsmag2 = dynsmag2 + CI_smag(ii,jj,kk)
            dynsmag3 = dynsmag3 + Pr_smag(ii,jj,kk)
            nentries_ds = nentries_ds + 1.0d0 
           
          enddo 
         enddo

         CS_smag(:,:,kk) = dynsmag1/nentries_ds 
         CI_smag(:,:,kk) = dynsmag2/nentries_ds 
         Pr_smag(:,:,kk) = dynsmag3/nentries_ds 
 
        enddo


            endif 


!-------------------------------------------------
!-------------------------------------------------
!-------------------------------------------------

        uu(:,:,:,1) = U(:,:,:)  
        uu(:,:,:,2) = V(:,:,:)  
        uu(:,:,:,3) = W(:,:,:)  

        S = 0.0d0    
        call compute_Sij(N0,N5,M0,M5,K0,K5,N2-2,N3+2,M2-2,M3+2,K2-2,K3+2,dx,dy,dz,uu,S)


         do ii = N2-2, N3+2 
          do jj = M2-2, M3+2
           do kk = K2-2, K3+2

             Sbar = 0.0d0
             do l = 1, 3
              do m = 1, 3
               Sbar = Sbar + 2.0d0*S(ii,jj,kk,l,m)*S(ii,jj,kk,l,m)
              enddo
             enddo 
             Sbar = sqrt(Sbar)   
             
             Skk = S(ii,jj,kk,1,1) + S(ii,jj,kk,2,2) + S(ii,jj,kk,3,3)  

             S1 = 0.0d0
             S2 = 0.0d0

             S1(:,:) = S(ii,jj,kk,:,:)
             do l = 1, 3
              S1(l,l) = S1(l,l) - Skk/3.0d0
             enddo
             S2(:,:) = -S1(:,:)*Sbar*Q(ii,jj,kk,1)*((delbar*CS_smag(ii,jj,kk))**2.0d0)*2.0d0   

             taukk = CI_smag(ii,jj,kk)*2.0d0*Q(ii,jj,kk,1)*((delbar*Sbar)**2.0d0)

             tau_sgs(ii,jj,kk,:,:) = S2(:,:)
             do l = 1, 3
              tau_sgs(ii,jj,kk,l,l) = tau_sgs(ii,jj,kk,l,l) + taukk/3.0d0
             enddo

             tau_u(ii,jj,kk,1) = tau_sgs(ii,jj,kk,1,1)*U(ii,jj,kk) + &
                                 tau_sgs(ii,jj,kk,1,2)*V(ii,jj,kk) + &
                                 tau_sgs(ii,jj,kk,1,3)*W(ii,jj,kk)
             tau_u(ii,jj,kk,2) = tau_sgs(ii,jj,kk,2,1)*U(ii,jj,kk) + &
                                 tau_sgs(ii,jj,kk,2,2)*V(ii,jj,kk) + &
                                 tau_sgs(ii,jj,kk,2,3)*W(ii,jj,kk)
             tau_u(ii,jj,kk,3) = tau_sgs(ii,jj,kk,3,1)*U(ii,jj,kk) + &
                                 tau_sgs(ii,jj,kk,3,2)*V(ii,jj,kk) + &
                                 tau_sgs(ii,jj,kk,3,3)*W(ii,jj,kk)

            


            hx = (H(ii-2,jj,kk) - 8.0d0*H(ii-1,jj,kk) + & 
                          8.0d0*H(ii+1,jj,kk) - H(ii+2,jj,kk))/(12.0d0*dx) 

            hy = (H(ii,jj-2,kk) - 8.0d0*H(ii,jj-1,kk) + & 
                          8.0d0*H(ii,jj+1,kk) - H(ii,jj+2,kk))/(12.0d0*dy) 

            hz = (H(ii,jj,kk-2) - 8.0d0*H(ii,jj,kk-1) + & 
                          8.0d0*H(ii,jj,kk+1) - H(ii,jj,kk+2))/(12.0d0*dz) 


            hfact = -((delbar*CS_smag(ii,jj,kk))**2.0d0)*Q(ii,jj,kk,1)*Sbar/Pr_smag(ii,jj,kk)
            q_sgs(ii,jj,kk,1) = hx*hfact
            q_sgs(ii,jj,kk,2) = hy*hfact
            q_sgs(ii,jj,kk,3) = hz*hfact

           enddo 
          enddo    
         enddo
         

!---------------------------------------

!-----------------------------------------
!  df/dx(i) = (f(i-2)-8f(i-1)+8f(i+1)-f(i+2))/(12*dx)
!-----------------------------------------


         do ii = N2, N3
          do jj = M2, M3
           do kk = K2, K3

            ! sgs flux, x-direction
            l = 1

            m = 1
            taux = (tau_sgs(ii-2,jj,kk,l,m) - 8.0d0*tau_sgs(ii-1,jj,kk,l,m) + & 
                          8.0d0*tau_sgs(ii+1,jj,kk,l,m) - tau_sgs(ii+2,jj,kk,l,m))/(12.0d0*dx) 

            m = 2
            tauy = (tau_sgs(ii,jj-2,kk,l,m) - 8.0d0*tau_sgs(ii,jj-1,kk,l,m) + & 
                          8.0d0*tau_sgs(ii,jj+1,kk,l,m) - tau_sgs(ii,jj+2,kk,l,m))/(12.0d0*dy) 

            m = 3
            tauz = (tau_sgs(ii,jj,kk-2,l,m) - 8.0d0*tau_sgs(ii,jj,kk-1,l,m) + & 
                          8.0d0*tau_sgs(ii,jj,kk+1,l,m) - tau_sgs(ii,jj,kk+2,l,m))/(12.0d0*dz) 

  
            dflux(ii,jj,kk,2) = dflux(ii,jj,kk,2) - taux - tauy - tauz   

 

            ! sgs flux, y-direction
            l = 2

            m = 1
            taux = (tau_sgs(ii-2,jj,kk,l,m) - 8.0d0*tau_sgs(ii-1,jj,kk,l,m) + & 
                          8.0d0*tau_sgs(ii+1,jj,kk,l,m) - tau_sgs(ii+2,jj,kk,l,m))/(12.0d0*dx) 

            m = 2
            tauy = (tau_sgs(ii,jj-2,kk,l,m) - 8.0d0*tau_sgs(ii,jj-1,kk,l,m) + & 
                          8.0d0*tau_sgs(ii,jj+1,kk,l,m) - tau_sgs(ii,jj+2,kk,l,m))/(12.0d0*dy) 

            m = 3
            tauz = (tau_sgs(ii,jj,kk-2,l,m) - 8.0d0*tau_sgs(ii,jj,kk-1,l,m) + & 
                          8.0d0*tau_sgs(ii,jj,kk+1,l,m) - tau_sgs(ii,jj,kk+2,l,m))/(12.0d0*dz) 


            dflux(ii,jj,kk,3) = dflux(ii,jj,kk,3) - taux - tauy - tauz               


 
            ! sgs flux, z-direction
            l = 3
  
            m = 1
            taux = (tau_sgs(ii-2,jj,kk,l,m) - 8.0d0*tau_sgs(ii-1,jj,kk,l,m) + & 
                          8.0d0*tau_sgs(ii+1,jj,kk,l,m) - tau_sgs(ii+2,jj,kk,l,m))/(12.0d0*dx) 

            m = 2
            tauy = (tau_sgs(ii,jj-2,kk,l,m) - 8.0d0*tau_sgs(ii,jj-1,kk,l,m) + & 
                          8.0d0*tau_sgs(ii,jj+1,kk,l,m) - tau_sgs(ii,jj+2,kk,l,m))/(12.0d0*dy) 

            m = 3
            tauz = (tau_sgs(ii,jj,kk-2,l,m) - 8.0d0*tau_sgs(ii,jj,kk-1,l,m) + & 
                          8.0d0*tau_sgs(ii,jj,kk+1,l,m) - tau_sgs(ii,jj,kk+2,l,m))/(12.0d0*dz) 


            dflux(ii,jj,kk,4) = dflux(ii,jj,kk,4) - taux - tauy - tauz 



            ! sgs viscous dissipation

            m = 1
            taux = (tau_u(ii-2,jj,kk,m) - 8.0d0*tau_u(ii-1,jj,kk,m) + &
                          8.0d0*tau_u(ii+1,jj,kk,m) - tau_u(ii+2,jj,kk,m))/(12.0d0*dx)
 
            m = 2 
            tauy = (tau_u(ii,jj-2,kk,m) - 8.0d0*tau_u(ii,jj-1,kk,m) + & 
                          8.0d0*tau_u(ii,jj+1,kk,m) - tau_u(ii,jj+2,kk,m))/(12.0d0*dy) 

            m = 3 
            tauz = (tau_u(ii,jj,kk-2,m) - 8.0d0*tau_u(ii,jj,kk-1,m) + & 
                          8.0d0*tau_u(ii,jj,kk+1,m) - tau_u(ii,jj,kk+2,m))/(12.0d0*dz) 

            dflux(ii,jj,kk,5) = dflux(ii,jj,kk,5) - taux - tauy - tauz 



            ! sgs thermal conduction

            m = 1
            qx = (q_sgs(ii-2,jj,kk,m) - 8.0d0*q_sgs(ii-1,jj,kk,m) + &
                          8.0d0*q_sgs(ii+1,jj,kk,m) - q_sgs(ii+2,jj,kk,m))/(12.0d0*dx) 

            m = 2
            qy = (q_sgs(ii,jj-2,kk,m) - 8.0d0*q_sgs(ii,jj-1,kk,m) + &
                          8.0d0*q_sgs(ii,jj+1,kk,m) - q_sgs(ii,jj+2,kk,m))/(12.0d0*dy) 

            m = 3
            qz = (q_sgs(ii,jj,kk-2,m) - 8.0d0*q_sgs(ii,jj,kk-1,m) + &
                          8.0d0*q_sgs(ii,jj,kk+1,m) - q_sgs(ii,jj,kk+2,m))/(12.0d0*dz) 

            dflux(ii,jj,kk,5) = dflux(ii,jj,kk,5) - (qx + qy + qz)


           enddo
          enddo
         enddo 


                   


         return

        end subroutine 

!-------------------------------------------------------------


        subroutine compute_Sij(N0,N5,M0,M5,K0,K5,i1,i2,j1,j2,k1,k2,dx,dy,dz,uu,Sij)

          implicit none

          integer :: N0, N5, M0, M5, K0, K5
          integer :: i1, i2, j1, j2, k1, k2, ii, jj, kk
             
          real(kind=8) :: dx, dy, dz
          real(kind=8) :: uu(N0:N5,M0:M5,K0:K5,1:3), Sij(N0:N5,M0:M5,K0:K5,1:3,1:3) 
          real(kind=8), dimension(N0:N5,M0:M5,K0:K5) :: u, v, w
          real(kind=8) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz    

               Sij = 0.0d0



           u = uu(:,:,:,1)
           v = uu(:,:,:,2)             
           w = uu(:,:,:,3)
                        

       do ii = i1, i2 
       do jj = j1, j2 
       do kk = k1, k2 

!-------------------------------------------------


            ! compute velocity gradients


         dudx = (u(ii-2,jj,kk) - 8.0d0*u(ii-1,jj,kk) + 8.0d0*u(ii+1,jj,kk) - u(ii+2,jj,kk))/(12.0d0*dx)
         dvdx = (v(ii-2,jj,kk) - 8.0d0*v(ii-1,jj,kk) + 8.0d0*v(ii+1,jj,kk) - v(ii+2,jj,kk))/(12.0d0*dx)
         dwdx = (w(ii-2,jj,kk) - 8.0d0*w(ii-1,jj,kk) + 8.0d0*w(ii+1,jj,kk) - w(ii+2,jj,kk))/(12.0d0*dx)

         dudy = (u(ii,jj-2,kk) - 8.0d0*u(ii,jj-1,kk) + 8.0d0*u(ii,jj+1,kk) - u(ii,jj+2,kk))/(12.0d0*dy)
         dvdy = (v(ii,jj-2,kk) - 8.0d0*v(ii,jj-1,kk) + 8.0d0*v(ii,jj+1,kk) - v(ii,jj+2,kk))/(12.0d0*dy)
         dwdy = (w(ii,jj-2,kk) - 8.0d0*w(ii,jj-1,kk) + 8.0d0*w(ii,jj+1,kk) - w(ii,jj+2,kk))/(12.0d0*dy)

         dudz = (u(ii,jj,kk-2) - 8.0d0*u(ii,jj,kk-1) + 8.0d0*u(ii,jj,kk+1) - u(ii,jj,kk+2))/(12.0d0*dz)
         dvdz = (v(ii,jj,kk-2) - 8.0d0*v(ii,jj,kk-1) + 8.0d0*v(ii,jj,kk+1) - v(ii,jj,kk+2))/(12.0d0*dz)
         dwdz = (w(ii,jj,kk-2) - 8.0d0*w(ii,jj,kk-1) + 8.0d0*w(ii,jj,kk+1) - w(ii,jj,kk+2))/(12.0d0*dz)


!-------------------------------------------------

            ! compute Sij matrix

            Sij(ii,jj,kk,1,1) = dudx
            Sij(ii,jj,kk,1,2) = 0.5d0*(dudy+dvdx)
            Sij(ii,jj,kk,1,3) = 0.5d0*(dudz+dwdx)
            
            Sij(ii,jj,kk,2,1) = Sij(ii,jj,kk,1,2)
            Sij(ii,jj,kk,2,2) = dvdy
            Sij(ii,jj,kk,2,3) = 0.5d0*(dvdz+dwdy)
            
            Sij(ii,jj,kk,3,1) = Sij(ii,jj,kk,1,3)
            Sij(ii,jj,kk,3,2) = Sij(ii,jj,kk,2,3)
            Sij(ii,jj,kk,3,3) = dwdz



       enddo
       enddo
       enddo


                    return
                    end subroutine

!-------------------------------------------------

        subroutine filter(i1,i2,j1,j2,k1,k2,f33,fbar)

       implicit none

       integer :: i1, i2, j1, j2, k1, k2   
       real(kind=8), dimension(i1:i2,j1:j2,k1:k2) :: f33, fbar
       real(kind=8) :: termx, termy, termz   
       integer :: ii, jj, kk
           

       fbar = 0.0d0
  

       do ii = i1+1, i2-1
       do jj = j1+1, j2-1
       do kk = k1+1, k2-1   

       termx = 0.25d0*f33(ii-1,jj,kk) + 0.5d0*f33(ii,jj,kk) + 0.25d0*f33(ii+1,jj,kk)
       termy = 0.25d0*f33(ii,jj-1,kk) + 0.5d0*f33(ii,jj,kk) + 0.25d0*f33(ii,jj+1,kk)
       termz = 0.25d0*f33(ii,jj,kk-1) + 0.5d0*f33(ii,jj,kk) + 0.25d0*f33(ii,jj,kk+1)
       fbar(ii,jj,kk) = (termx + termy + termz)/3.0d0

       enddo 
       enddo 
       enddo 


            return
        end subroutine

!-------------------------------------------------
