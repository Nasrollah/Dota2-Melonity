      do r = 0,k
          V(r) = SUM(C(r,:)*F(i,r-k:r,l))

          b(r) = F(i,r-0,l)*SUM(S(k-0,k-0:k,r)*F(i,r-0:r,l))
        do j = 1,k
          b(r) = F(i,r-j,l)*SUM(S(k-j,k-j:k,r)*F(i,r-j:r,l)) + b(r)
        enddo
      enddo
