  integer                    :: Order, N0, N5, i0, i1
  REALTYPE, dimension(N0:N5) :: F, F_Half

  REALTYPE, OPTIONAL         :: WENO_Epsilon
  REALTYPE                   :: Epsilon_WENO

  REALTYPE, dimension(0:(Order-1)/2                            ) :: V, a, b
  REALTYPE, dimension(0:(Order-1)/2                            ) :: d, Omega
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2              ) :: C
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2,0:(Order-1)/2) :: S

  integer  :: i, j, k, m, r
  REALTYPE :: Sum_a, Tau
