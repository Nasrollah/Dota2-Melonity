  Subroutine Geometry

  integer  :: i, j, k

  x0 = Left_x ; x1 = Right_x 
  y0 = Left_y ; y1 = Right_y 
  z0 = Left_z ; z1 = Right_z 

  call PS_WENO_Uniform_Grid (N0_g,N2_g, N3_g,N5_g, x0, x1, dx,  xi, Symmetry_x)
  call PS_WENO_Uniform_Grid (M0_g,M2_g, M3_g,M5_g, y0, y1, dy, eta, Symmetry_y)
  call PS_WENO_Uniform_Grid (K0_g,K2_g, K3_g,K5_g, z0, z1, dz,  fi, Symmetry_z) 

  do k = K0,K5
    do j = M0,M5
      do i = N0,N5
        x(i,j,k) =  xi(i+n_0)
        y(i,j,k) = eta(j+m_0)
        z(i,j,k) =  fi(k+k_0)
      enddo
    enddo
  enddo

  x_1D = xi ; y_1D = eta ; z_1D = fi

  Normal_N0(:,:,:,1) = -ONE ; Normal_N0(:,:,:,2) = ZERO
                              Normal_N0(:,:,:,3) = ZERO
  Normal_N5(:,:,:,1) =  ONE ; Normal_N5(:,:,:,2) = ZERO
                              Normal_N5(:,:,:,3) = ZERO

  Normal_M0(:,:,:,1) = ZERO ; Normal_M0(:,:,:,2) = -ONE
                              Normal_M0(:,:,:,3) = ZERO
  Normal_M5(:,:,:,1) = ZERO ; Normal_M5(:,:,:,2) =  ONE
                              Normal_M5(:,:,:,3) = ZERO

  Normal_K0(:,:,:,1) = ZERO ; Normal_K0(:,:,:,2) = ZERO
                              Normal_K0(:,:,:,3) = -ONE
  Normal_K5(:,:,:,1) = ZERO ; Normal_K5(:,:,:,2) = ZERO
                              Normal_K5(:,:,:,3) =  ONE

  End Subroutine Geometry
