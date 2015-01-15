write (6,*) '------------------------------------------'
write (6,*) 'PS_WENO euler flux general in y'
write (6,*) Index_y, Distributed_y, Order_y, dy
write (6,*) 'n   : ',N0,N2,N3,N5
write (6,*) 'm   : ',M0,M2,M3,M5
write (6,*) 'nv  : ',NV
write (6,*) 'rho : ',MINVAL(Q(:,:,1)),MAXVAL(Q(:,:,1))
write (6,*) 'u   : ',MINVAL(U),MAXVAL(U)
write (6,*) 'v   : ',MINVAL(V),MAXVAL(V)
write (6,*) 'p   : ',MINVAL(P),MAXVAL(P)
write (6,*) 'c   : ',MINVAL(C),MAXVAL(C)
write (6,*) 'F 1 : ',MINVAL(Flux(:,:,1)),MAXVAL(Flux(:,:,1))
write (6,*) 'F 2 : ',MINVAL(Flux(:,:,2)),MAXVAL(Flux(:,:,2))
write (6,*) 'F 3 : ',MINVAL(Flux(:,:,3)),MAXVAL(Flux(:,:,3))
write (6,*) 'F 4 : ',MINVAL(Flux(:,:,4)),MAXVAL(Flux(:,:,4))
write (6,*) 'F 5 : ',MINVAL(Flux(:,:,5)),MAXVAL(Flux(:,:,5))
write (6,*) 'DF1 : ',MINVAL(D_G(:,:,1)),MAXVAL(D_G(:,:,1))
write (6,*) 'DF2 : ',MINVAL(D_G(:,:,2)),MAXVAL(D_G(:,:,2))
write (6,*) 'DF3 : ',MINVAL(D_G(:,:,3)),MAXVAL(D_G(:,:,3))
write (6,*) 'DF4 : ',MINVAL(D_G(:,:,4)),MAXVAL(D_G(:,:,4))
write (6,*) 'DF5 : ',MINVAL(D_G(:,:,5)),MAXVAL(D_G(:,:,5))
write (6,*) 'gam : ',MINVAL(Gamma),MAXVAL(Gamma)
write (6,*) 'gam : ',gamma_gas(1,:)
write (6,*) 'gam : ',gamma_gas(2,:)
if (PRESENT(D_F_Epsilon)) write (6,*) 'eps : ',D_F_Epsilon

write (6,*) '------------------------------------------'
