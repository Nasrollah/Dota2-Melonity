        do r = 0,k
          V(r) = SUM(C(r,:)*F(i,r-k:r,l))

          b(r) = F(i,r-4,l)*SUM(S(0,0:k,r)*F(i,r-4:r,l)) + &
                 F(i,r-3,l)*SUM(S(1,1:k,r)*F(i,r-3:r,l)) + &
                 F(i,r-2,l)*SUM(S(2,2:k,r)*F(i,r-2:r,l)) + &
                 F(i,r-1,l)*SUM(S(3,3:k,r)*F(i,r-1:r,l)) + &
                 F(i,r-0,l)*SUM(S(4,4:k,r)*F(i,r-0:r,l))
        enddo
