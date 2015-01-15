  integer                    :: Order, N0, N5, i0, i1

  REALTYPE, OPTIONAL         :: WENO_Epsilon
  REALTYPE                   :: Epsilon_WENO

  REALTYPE, dimension(N0:N5,-(Order-1)/2:(Order-1)/2,2) :: F
  REALTYPE, dimension(N0:N5                           ) :: F_Half

  REALTYPE, dimension(0:(Order-1)/2                            ) :: V, a, b
  REALTYPE, dimension(0:(Order-1)/2                            ) :: d, Omega
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2              ) :: C
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2,0:(Order-1)/2) :: S

  integer                     :: i, l
  REALTYPE                    :: Sum_a
