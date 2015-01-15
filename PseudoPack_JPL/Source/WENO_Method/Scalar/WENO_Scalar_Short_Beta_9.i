        do r = 0,k
          V(r) = SUM(C(r,:)*F(i+r-k:i+r))

          b(r) = F(i+r-4)*SUM(S(0,0:k,r)*F(i+r-4:i+r)) + &
                 F(i+r-3)*SUM(S(1,1:k,r)*F(i+r-3:i+r)) + &
                 F(i+r-2)*SUM(S(2,2:k,r)*F(i+r-2:i+r)) + &
                 F(i+r-1)*SUM(S(3,3:k,r)*F(i+r-1:i+r)) + &
                 F(i+r-0)*SUM(S(4,4:k,r)*F(i+r-0:i+r))
        enddo
