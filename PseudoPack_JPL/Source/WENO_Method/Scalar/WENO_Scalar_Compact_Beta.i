    do r = 0,k
      j = i+r 

        V(r) = SUM(C(r,:)*F(j-k:j))

        b(r) = F(j-0)*SUM(S(k-0,k-0:k,r)*F(j-0:j))
      do m = 1,k
        b(r) = F(j-m)*SUM(S(k-m,k-m:k,r)*F(j-m:j)) + b(r)
      enddo
    enddo
