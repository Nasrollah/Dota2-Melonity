#if defined (VFFT)

#include "pseudopack.h"


MODULE VFFT_FFT

implicit NONE

PRIVATE

#if defined (LOCAL_VFFT)
PUBLIC  :: vcosqi, vcosqf, vcosqb
PUBLIC  ::         vsinqf, vsinqb
PUBLIC  :: vcosti, vcost
PUBLIC  :: vsinti, vsint
PUBLIC  :: vrffti, vrfftf, vrfftb
#endif

PUBLIC  :: VFFT_QWV_Init
PUBLIC  :: VFFT_COS_Init
PUBLIC  :: VFFT_SIN_Init

CONTAINS

  Subroutine VFFT_QWV_Init (N1, QWave)

  integer                :: N1, N
  REALTYPE, dimension(*) :: QWave

  N = N1-1
  call vcosqi (N, QWave)

  END Subroutine VFFT_QWV_Init 
!
! -----------------------------------------------------------------------
!
  Subroutine VFFT_SIN_Init (N0, SWave)

  integer                :: N0, N
  REALTYPE, dimension(*) :: SWave

  N = N0-1
  call vsinti (N, SWave)

  END Subroutine VFFT_SIN_Init 
!
! -----------------------------------------------------------------------
!
  Subroutine VFFT_COS_Init (N0, CWave)

  integer                :: N0, N
  REALTYPE, dimension(*) :: CWave

  N = N0+1
  call vcosti (N, CWave)

  END Subroutine VFFT_COS_Init 

!
! -----------------------------------------------------------------------
!
#if defined (LOCAL_VFFT)
  Subroutine vsinqf (m,n,x,xt,ldx,WSave)

  integer                     :: m, n, ldx
  REALTYPE, dimension(ldx,*)  :: x, xt
  REALTYPE, dimension(*)      :: WSave

  integer                     :: j, k, ns2, kc
  REALTYPE                    :: x_Temp

  if ((m <= 0) .OR. (n <= 1)) RETURN

!     ... preprocessing

  ns2 = n/2

  do k = 1,ns2
    kc = n-k
    do j = 1,m
      x_Temp    = x(j,k)
      x(j,k   ) = x(j,kc+1)
      x(j,kc+1) = x_Temp
    enddo
  enddo
 
!     ... cosine quarter wave transform
 
  call vcosqf (m,n,x,xt,ldx,WSave)
 
!     ... postprocessing
 
  do k = 2,n,2
    do j = 1,m
      x(j,k) = -x(j,k)
    enddo
  enddo

  END Subroutine vsinqf 
!
! -----------------------------------------------------------------------
!
  Subroutine vsinqb (m,n,x,xt,ldx,WSave)

  integer                     :: m, n, ldx
  REALTYPE, dimension(ldx,*)  :: x, xt
  REALTYPE, dimension(*)      :: WSave

  integer                     :: j, k, kc, ns2
  REALTYPE                    :: x_Temp

  if ((m <= 0) .OR. (n <= 1)) RETURN
 
!     ... preprocessing
 
  ns2 = n/2

  do k = 2,n,2
    do j = 1,m
      x(j,k) = -x(j,k)
    enddo
  enddo
 
!     ... cosine quarter wave transform
 
  call vcosqb (m,n,x,xt,ldx,WSave)
 
!     ... postprocessing
 
  do k = 1,ns2
    kc = n-k
    do j = 1,m
      x_Temp    = x(j,k)
      x(j,k   ) = x(j,kc+1)
      x(j,kc+1) = x_Temp
    enddo
  enddo

  END Subroutine vsinqb 
!
! -----------------------------------------------------------------------
!
  Subroutine vcosqf (m,n,x,xt,ldx,WSave)

  integer                     :: m, n, ldx
  REALTYPE, dimension(ldx,*)  :: x, xt
  REALTYPE, dimension(2*n+15) :: WSave

  integer                     :: i, j, k, kc, ns2, np2, modn
  REALTYPE                    :: x_Temp, sqrt2, scale

  if ((m <= 0) .OR. (n <  2)) RETURN
 
!  case  n = 2
 
  if (n == 2) then
    sqrt2 = SQRT(TWO)
    scale = HALF*SQRT_HALF
    do j = 1,m
       x_Temp = sqrt2* x(j,2)
       x(j,2) = scale*(x(j,1)-x_Temp)
       x(j,1) = scale*(x(j,1)+x_Temp)
    enddo
    RETURN
  endif
 
!  case  n > 2
 
!     ... preprocessing
 
  ns2 = (n+1)/2
  np2 =  n+2

  do k = 2,ns2
    kc = np2-k
    do j = 1,m
      xt(j,k ) = x(j,k)+x(j,kc)
      xt(j,kc) = x(j,k)-x(j,kc)
    enddo
  enddo

  modn = MOD(n,2)
  if (modn == 0) then
    do j = 1,m
      xt(j,ns2+1) = x(j,ns2+1)+x(j,ns2+1)
    enddo
  endif

  do k = 2,ns2
    kc = np2-k
    do j = 1,m
      x(j,k ) = WSave(k-1)*xt(j,kc)+WSave(kc-1)*xt(j,k )
      x(j,kc) = WSave(k-1)*xt(j,k )-WSave(kc-1)*xt(j,kc)
    enddo
  enddo

  if (modn == 0) then
    do j = 1,m
      x(j,ns2+1) = WSave(ns2)*xt(j,ns2+1)
    enddo
  endif
 
!     ... real, periodic transform
 
  call vrfftf (m,n,x,xt,ldx,WSave(n+1))
 
!     ... postprocessing
 
  do i = 3,n,2
    do j = 1,m
      x_Temp   = x(j,i-1)-x(j,i)
      x(j,i  ) = x(j,i-1)+x(j,i)
      x(j,i-1) = x_Temp
    enddo
  enddo
 
!     ... normalization
 
  scale = SQRT_HALF
  do i = 1,n
    do j = 1,m
      x(j,i) = scale*x(j,i)
    enddo
  enddo

  END Subroutine vcosqf 
!
! -----------------------------------------------------------------------
!
  Subroutine vcosqb (m,n,x,xt,ldx,WSave)

  integer                     :: m, n, ldx
  REALTYPE, dimension(ldx,*)  :: x, xt
  REALTYPE, dimension(*)      :: WSave

  integer                     :: i, j, k, kc, ns2, np2, modn
  REALTYPE                    :: x_Temp, scale

  if ((m <= 0) .OR. (n <= 1)) RETURN
 
!  case  n = 2
 
  if (n == 2) then
    scale = TWO*SQRT_HALF
    do j = 1,m
      x_Temp = scale*(x(j,1)+x(j,2))
      x(j,2) =        x(j,1)-x(j,2)
      x(j,1) =        x_Temp
    enddo
    RETURN
  endif
 
!  case n > 2
 
!     ... preprocessing
 
  ns2 = (n+1)/2
  np2 =  n+2

  do i = 3,n,2
    do j = 1,m
      x_Temp   = x(j,i-1)+x(j,i  )
      x(j,i  ) = x(j,i  )-x(j,i-1)
      x(j,i-1) = x_Temp
    enddo
  enddo

  do j = 1,m
    x(j,1) = x(j,1)+x(j,1)
  enddo

  modn = MOD(n,2)
  if (modn == 0) then
    do j = 1,m
      x(j,n) = x(j,n)+x(j,n)
    enddo
  endif
 
!     ... real, periodic transform
 
  call vrfftb (m,n,x,xt,ldx,WSave(n+1))
 
!     ... postprocessing
 
  do k = 2,ns2
    kc = np2-k
    do j = 1,m
      xt(j,k ) = WSave(k-1)*x(j,kc)+WSave(kc-1)*x(j,k )
      xt(j,kc) = WSave(k-1)*x(j,k )-WSave(kc-1)*x(j,kc)
    enddo
  enddo

  if (modn == 0) then
    do j = 1,m
      x(j,ns2+1) = WSave(ns2)*(x(j,ns2+1)+x(j,ns2+1))
    enddo
  endif

  do k = 2,ns2
    kc = np2-k
    do j = 1,m
      x(j,k ) = xt(j,k)+xt(j,kc)
      x(j,kc) = xt(j,k)-xt(j,kc)
    enddo
  enddo

  do j = 1,m
    x(j,1) = x(j,1)+x(j,1)
  enddo
 
!     ... normalization
 
  scale = HALF*SQRT_HALF
  do i = 1,n
    do j = 1,m
      x(j,i) = scale*x(j,i)
    enddo
  enddo

  END Subroutine vcosqb 
!
! -----------------------------------------------------------------------
!
  Subroutine vcosqi (n,WSave)

  integer                :: n
  REALTYPE, dimension(*) :: WSave

  integer                :: k
  REALTYPE               :: dt, fk

  if (n <= 2) RETURN

  dt = HALF*PI/n

  fk = ZERO
  do k = 1,n
    fk = fk + ONE
    WSave(k) = COS(fk*dt)
  enddo

  call vrffti (n,WSave(n+1))

  END Subroutine vcosqi 
!
! -----------------------------------------------------------------------
!
  Subroutine vcost  (m,n,x,xt,ldx,WSave)

  integer                     :: m, n, ldx
  REALTYPE, dimension(ldx,n)  :: x, xt
  REALTYPE, dimension(3*n+15) :: WSave
 
  integer                     :: i, j, ic, ns2, modn, nm1, np1
  REALTYPE                    :: x1h, x1p3, tx2, scale
  REALTYPE                    :: t1, t2, factor, xi

  if ((m <= 0) .OR. (n <= 1)) RETURN
 
  if (n == 2) then
    scale = SQRT_HALF
    do j = 1,m
         x1h = scale*(x(j,1)+x(j,2))
      x(j,2) = scale*(x(j,1)-x(j,2))
      x(j,1) = x1h
    enddo
    RETURN
  endif

  if (n == 3) then
    scale = HALF
    do j = 1,m
        x1p3 = x(j,1)+x(j,3)
         tx2 = x(j,2)+x(j,2)
      x(j,2) = scale*(x(j,1)-x(j,3))
      x(j,1) = scale*(x1p3+tx2)
      x(j,3) = scale*(x1p3-tx2)
    enddo
    RETURN
  endif
 
!  case  n >  3
 
!     ... preprocessing
 
  nm1 = n-1
  np1 = n+1
  ns2 = n/2

  do j = 1,m
    xt(j,1) = x(j,1)-x(j,n)
     x(j,1) = x(j,1)+x(j,n)
  enddo

  do i = 2,ns2
    ic = np1-i
    do j = 1,m
      t1 = x(j,i)+x(j,ic)
      t2 = x(j,i)-x(j,ic)
      xt(j,1) = xt(j,1)+WSave(ic)*t2

      t2 = WSave(i)*t2
      x(j,i ) = t1-t2
      x(j,ic) = t1+t2
    enddo
  enddo

  modn = MOD(n,2)
  if (modn /= 0) then
    do j = 1,m
      x(j,ns2+1) = x(j,ns2+1)+x(j,ns2+1)
    enddo
  endif

  do j = 1,m
    x(j,n) = xt(j,1)
  enddo
 
!     ... real periodic transform
 
  call vrfftf (m,nm1,x,xt,ldx,WSave(np1))
 
!     ... postprocessing
 
  factor = SQRT(ONE/nm1) 
  do j = 1,m
    xt(j,1) =        x(j,2)
     x(j,2) = factor*x(j,n)
  enddo

  do i = 4,n,2
    do j = 1,m
      xi        =  x(j,i  )
       x(j,i  ) =  x(j,i-2)-x(j,i-1)
       x(j,i-1) = xt(j,1  )
      xt(j,1  ) = xi
    enddo
  enddo

  if (modn /= 0) then
    do j = 1,m
      x(j,n) = xt(j,1)
    enddo
  endif
 
!     ... normalization
 
  scale = SQRT_HALF
  do i = 1,n
    do j = 1,m
      x(j,i) = scale*x(j,i)
    enddo
  enddo

  END Subroutine vcost  
!
! -----------------------------------------------------------------------
!
  Subroutine vcosti (n,WSave)

  integer                :: n
  REALTYPE, dimension(*) :: WSave

  integer                :: k, kc, nm1, np1, ns2
  REALTYPE               :: dt, fk

  if (n <= 3) RETURN

  nm1 = n-1
  np1 = n+1
  ns2 = n/2
  dt  = PI/nm1

  fk = ZERO
  do k = 2,ns2
    fk = fk + ONE
    WSave(k) = TWO*SIN(fk*dt)
  enddo

  fk = ZERO
  do k = 2,ns2
    kc = np1-k
    fk = fk+ONE
    WSave(kc) = TWO*COS(fk*dt)
  enddo

  call vrffti (nm1,WSave(n+1))

  END Subroutine vcosti 



!
! -----------------------------------------------------------------------
!
  Subroutine vsint (m,n,x,xt,ldx,WSave)
      
  integer                     :: m, n, ldx
  REALTYPE, dimension(ldx,*)  :: x, xt
  REALTYPE, dimension(*)      :: WSave

  integer                     :: i, j, k, kc, np1, ns2, modn, nf
  REALTYPE                    :: xh, sqrth, scale
  REALTYPE                    :: t1, t2

  if ((m <= 0) .OR. (n <= 1)) RETURN

  if (n == 2) then
    sqrth = SQRT_HALF
    do j = 1,m
          xh = sqrth*(x(j,1)+x(j,2))
      x(j,2) = sqrth*(x(j,1)-x(j,2))
      x(j,1) = xh
    enddo
    RETURN
  endif
 
!     ... preprocessing
 
  np1 = n+1
  ns2 = n/2

  do j = 1,m
    xt(j,1) = ZERO
  enddo

  do k = 1,ns2
    kc = np1-k
    do j = 1,m
      t1 =  x(j,k)-x(j,kc)
      t2 = (x(j,k)+x(j,kc))*WSave(k)

      xt(j,k +1) = t1+t2
      xt(j,kc+1) = t2-t1
    enddo
  enddo

  modn = MOD(n,2)
  if (modn /= 0) then
    do j = 1,m
      xt(j,ns2+2) = FOUR*x(j,ns2+1)
    enddo
  endif
 
!     ... real periodic transform
 
  nf = ns2+1
  call vrfftf (m,np1,xt,x,ldx,WSave(nf))
 
!     ... postprocessing
 
  do j = 1,m
    x(j,1) = HALF*xt(j,1)
  enddo

  do i = 3,n,2
    do j = 1,m
      x(j,i-1) = -xt(j,i  )
      x(j,i  ) =  xt(j,i-1)+x(j,i-2)
    enddo
  enddo

  if (modn == 0) then
    do j = 1,m
      x(j,n) = -xt(j,n+1)
    enddo
  endif
 
!     ... normalization
 
  scale = SQRT_HALF
  do i = 1,n
    do j = 1,m
      x(j,i) = scale*x(j,i)
    enddo
  enddo

  END Subroutine vsint 
!
! -----------------------------------------------------------------------
!
  Subroutine vsinti (n,WSave)
      
  integer                :: n
  REALTYPE, dimension(*) :: WSave

  integer                :: k
  REALTYPE               :: dt, fk

  if (n <= 1) RETURN

  dt = PI/(n+1)

  fk = ZERO
  do k = 1,n/2
          fk = fk + ONE
    WSave(k) = TWO*sin(fk*dt)
  enddo

  call vrffti (n+1,WSave(n/2+1))

  END Subroutine vsinti 

!
! -----------------------------------------------------------------------
!
  Subroutine vrfftf (m,n,x,xt,ldx,WSave)

  integer                     :: m, n, ldx
  REALTYPE, dimension(ldx,n)  :: x, xt
  REALTYPE, dimension(n+15)   :: WSave

  if (n == 1) RETURN

  call vrftf1 (m,n,x,xt,ldx,WSave(1),WSave(n+1))

  END Subroutine vrfftf 
!
! -----------------------------------------------------------------------
!
  Subroutine vrftf1 (m,n,c,ch,ldc,wa,fac)

  integer                    :: m, n, ldc
  REALTYPE, dimension(ldc,n) :: c, ch
  REALTYPE, dimension(n)     :: wa
  REALTYPE, dimension(15)    :: fac

  integer                    :: i, j
  integer                    :: nf, na, l2, iw, kh, ip, l1, ido, idl1
  integer                    :: k1, ix2, ix3, ix4
  REALTYPE                   :: scale

  nf = fac(2)
  na = 1
  l2 = n
  iw = n

  do k1 = 1,nf
    kh   = nf-k1
    ip   = fac(kh+3)
    l1   = l2/ip
    ido  = n/l2
    idl1 = ido*l1
    iw   = iw-(ip-1)*ido
    na   = 1-na

    if (ip /= 4) goto 102
    ix2 = iw +ido
    ix3 = ix2+ido
    if (na /= 0) then
      call vradf4 (m,ido,l1,ch,c,ldc,wa(iw),wa(ix2),wa(ix3))
    else
      call vradf4 (m,ido,l1,c,ch,ldc,wa(iw),wa(ix2),wa(ix3))
    endif
    goto 110

 102 continue
    if (ip /= 2) goto 104
    if (na /= 0) then
      call vradf2 (m,ido,l1,ch,c,ldc,wa(iw))
    else
      call vradf2 (m,ido,l1,c,ch,ldc,wa(iw))
    endif
    goto 110

 104 continue
    if (ip /= 3) goto 106
    ix2 = iw +ido
    if (na /= 0) then
      call vradf3 (m,ido,l1,ch,c,ldc,wa(iw),wa(ix2))
    else
      call vradf3 (m,ido,l1,c,ch,ldc,wa(iw),wa(ix2))
    endif
    goto 110

 106 continue
    if (ip /= 5) goto 108
    ix2 = iw +ido
    ix3 = ix2+ido
    ix4 = ix3+ido
    if (na /= 0) then
      call vradf5 (m,ido,l1,ch,c,ldc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
    else
      call vradf5 (m,ido,l1,c,ch,ldc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
    endif
    goto 110

 108 continue
    if (ido == 1) na = 1-na
    if (na /= 0) then
      call vradfg (m,ido,ip,l1,idl1,ch,ch,ch,c ,c ,ldc,wa(iw))
      na = 0
    else
      call vradfg (m,ido,ip,l1,idl1,c ,c ,c ,ch,ch,ldc,wa(iw))
      na = 1
    endif

 110 continue
    l2 = l1
  enddo

  scale=SQRT(ONE/n)

  if (na /= 1) then
    do j = 1,n
      do i = 1,m
        c(i,j) = scale*ch(i,j)
      enddo
    enddo
  else
    do j = 1,n
      do i = 1,m
        c(i,j) = scale*c(i,j)
      enddo
    enddo
  endif

  END Subroutine vrftf1 


!
! -----------------------------------------------------------------------
!
  Subroutine vradf2 (mp,ido,l1,cc,ch,ldc,wa1)

  integer                           :: mp, ido, l1, ldc
  REALTYPE, dimension(ldc,ido,2,l1) :: ch
  REALTYPE, dimension(ldc,ido,l1,2) :: cc
  REALTYPE, dimension(ido)          :: wa1

  integer                           :: i, k, m, idp2, ic

  do k = 1,l1
    do m = 1,mp
      ch(m,1  ,1,k) = cc(m,1,k,1)+cc(m,1,k,2)
      ch(m,ido,2,k) = cc(m,1,k,1)-cc(m,1,k,2)
    enddo
  enddo

  if (ido-2 <  0) RETURN

  if (ido-2 >  0) then
    idp2 = ido+2

    do k = 1,l1
      do i = 3,ido,2
        ic = idp2-i
        do m = 1,mp
          ch(m,i   ,1,k) =  (wa1(i-2)*cc(m,i  ,k,2)-  &
                             wa1(i-1)*cc(m,i-1,k,2))+cc(m,i  ,k,1)
          ch(m,ic  ,2,k) =  (wa1(i-2)*cc(m,i  ,k,2)-  &
                             wa1(i-1)*cc(m,i-1,k,2))-cc(m,i  ,k,1)
          ch(m,i-1 ,1,k) =  (wa1(i-2)*cc(m,i-1,k,2)+  &
                             wa1(i-1)*cc(m,i  ,k,2))+cc(m,i-1,k,1)
          ch(m,ic-1,2,k) = -(wa1(i-2)*cc(m,i-1,k,2)+  &
                             wa1(i-1)*cc(m,i  ,k,2))+cc(m,i-1,k,1)
        enddo
      enddo
    enddo

    if (MOD(ido,2) == 1) RETURN
  endif

  do k = 1,l1
    do m = 1,mp
      ch(m,1  ,2,k) = -cc(m,ido,k,2)
      ch(m,ido,1,k) =  cc(m,ido,k,1)
    enddo
  enddo

  END Subroutine vradf2 
!
! -----------------------------------------------------------------------
!
  Subroutine vradf3 (mp,ido,l1,cc,ch,ldc,wa1,wa2)

  integer                           :: mp, ido, l1, ldc
  REALTYPE, dimension(ldc,ido,3,l1) :: ch
  REALTYPE, dimension(ldc,ido,l1,3) :: cc
  REALTYPE, dimension(ido)          :: wa1, wa2

  integer                           :: i, k, m
  integer                           :: idp2, ic
  REALTYPE                          :: arg, taur, taui

   arg=TWO*PI/THREE
  taur=COS(arg)
  taui=SIN(arg)

  do k = 1,l1
    do m = 1,mp
      ch(m,1  ,1,k) = cc(m,1,k,1)+     (cc(m,1,k,2)+cc(m,1,k,3))
      ch(m,1  ,3,k) =             taui*(cc(m,1,k,3)-cc(m,1,k,2))
      ch(m,ido,2,k) = cc(m,1,k,1)+taur*(cc(m,1,k,2)+cc(m,1,k,3))
    enddo
  enddo

  if (ido == 1) RETURN

  idp2 = ido+2
  do k = 1,l1
    do i = 3,ido,2
      ic = idp2-i
      do m = 1,mp
        ch(m,i-1 ,1,k) =    cc(m,i-1,k,1)+  &
                 ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3)))

        ch(m,i   ,1,k) =    cc(m,i  ,k,1)+  &
                 ((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3)))

        ch(m,i-1,3,k) =     cc(m,i-1,k,1)+  &
            taur*((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3)))+ &
            taui*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))-  &
                  (wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3)))

        ch(m,ic-1,2,k) =    cc(m,i-1,k,1)+  &
            taur*((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3)))- &
            taui*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))-  &
                  (wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3)))

        ch(m,i   ,3,k) =    cc(m,i  ,k,1)+  &
            taur*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3)))+ &
            taui*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))-  &
                  (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2)))

        ch(m,ic  ,2,k) =   -cc(m,i  ,k,1)-  &
            taur*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3)))+ &
            taui*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))-  &
                  (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2)))
      enddo
    enddo
  enddo

  END Subroutine vradf3 
!
! -----------------------------------------------------------------------
!
  Subroutine vradf4 (mp,ido,l1,cc,ch,ldc,wa1,wa2,wa3)

  integer                           :: mp, ido, l1, ldc
  REALTYPE, dimension(ldc,ido,4,l1) :: ch
  REALTYPE, dimension(ldc,ido,l1,4) :: cc
  REALTYPE, dimension(ido)          :: wa1, wa2, wa3

  integer                           :: i, k, m
  integer                           :: idp2, ic
  REALTYPE                          :: hsqt2

  do k = 1,l1
    do m = 1,mp
      ch(m,1  ,1,k) = (cc(m,1,k,2)+cc(m,1,k,4))  &
                     +(cc(m,1,k,1)+cc(m,1,k,3))
      ch(m,ido,4,k) = (cc(m,1,k,1)+cc(m,1,k,3))  &
                     -(cc(m,1,k,2)+cc(m,1,k,4))
      ch(m,ido,2,k) =  cc(m,1,k,1)-cc(m,1,k,3) 
      ch(m,1  ,3,k) =  cc(m,1,k,4)-cc(m,1,k,2)
    enddo
  enddo

  if (ido-2 <  0) RETURN

  if (ido-2 >  0) then
    idp2 = ido+2
    do k = 1,l1
      do i = 3,ido,2
        ic = idp2-i
        do m = 1,mp
          ch(m,i-1 ,1,k) =  cc(m,i-1,k,1)+  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))+  &
                  (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4)) 

          ch(m,ic-1,4,k) =  cc(m,i-1,k,1)+  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))-  &
                 ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4)))

          ch(m,i   ,1,k) =  cc(m,i  ,k,1)+  &
                  (wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))+  &
                 ((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4)))

          ch(m,ic  ,4,k) = -cc(m,i  ,k,1)-  &
                  (wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))+  &
                 ((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4)))

          ch(m,i-1 ,3,k) =  cc(m,i-1,k,1)-  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))+  &
                 ((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))-  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4)))

          ch(m,ic-1,2,k) =  cc(m,i-1,k,1)-  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))-  &
                 ((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))-  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4)))

          ch(m,i   ,3,k) =  cc(m,i  ,k,1)-  &
                  (wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))-  &
                 ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))-  &
                  (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4)))

          ch(m,ic  ,2,k) = -cc(m,i  ,k,1)+  &
                  (wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))-  &
                 ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))-  &
                  (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4)))
        enddo
      enddo
    enddo

    if (MOD(ido,2) == 1) RETURN
  endif

  hsqt2 = SQRT_HALF

  do k = 1,l1
    do m = 1,mp
      ch(m,ido,1,k) =  hsqt2*(cc(m,ido,k,2)-cc(m,ido,k,4))+  &
                              cc(m,ido,k,1)
      ch(m,ido,3,k) = -hsqt2*(cc(m,ido,k,2)-cc(m,ido,k,4))+  &
                              cc(m,ido,k,1)      
      ch(m,  1,2,k) = -hsqt2*(cc(m,ido,k,2)+cc(m,ido,k,4))-  &
                              cc(m,ido,k,3)
      ch(m,  1,4,k) = -hsqt2*(cc(m,ido,k,2)+cc(m,ido,k,4))+  &
                              cc(m,ido,k,3)
    enddo
  enddo

  END Subroutine vradf4 
!
! -----------------------------------------------------------------------
!
  Subroutine vradf5 (mp,ido,l1,cc,ch,ldc,wa1,wa2,wa3,wa4)

  integer                           :: mp, ido, l1, ldc
  REALTYPE, dimension(ldc,ido,5,l1) :: ch
  REALTYPE, dimension(ldc,ido,l1,5) :: cc
  REALTYPE, dimension(ido)          :: wa1, wa2, wa3, wa4

  integer                           :: i, k, m
  integer                           :: idp2, ic
  REALTYPE                          :: arg, tr11, ti11, tr12, ti12

  arg  = TWO*PI/FIVE
  tr11 = COS(arg)
  ti11 = SIN(arg)
  tr12 = COS(TWO*arg)
  ti12 = SIN(TWO*arg)

  do k = 1,l1
    do m = 1,mp
      ch(m,1  ,1,k) =      (cc(m,1,k,5)+cc(m,1,k,2))+  &
                           (cc(m,1,k,4)+cc(m,1,k,3))+cc(m,1,k,1)
      ch(m,ido,2,k) = tr11*(cc(m,1,k,5)+cc(m,1,k,2))+  &
                      tr12*(cc(m,1,k,4)+cc(m,1,k,3))+cc(m,1,k,1)
      ch(m,1  ,3,k) = ti11*(cc(m,1,k,5)-cc(m,1,k,2))+  &
                      ti12*(cc(m,1,k,4)-cc(m,1,k,3))
      ch(m,ido,4,k) = tr12*(cc(m,1,k,5)+cc(m,1,k,2))+  &
                      tr11*(cc(m,1,k,4)+cc(m,1,k,3))+cc(m,1,k,1)
      ch(m,1  ,5,k) = ti12*(cc(m,1,k,5)-cc(m,1,k,2))-  &
                      ti11*(cc(m,1,k,4)-cc(m,1,k,3))
    enddo
  enddo

  if (ido == 1) RETURN

  idp2 = ido+2
  do k = 1,l1
    do i = 3,ido,2
      ic = idp2-i
      do m = 1,mp
        ch(m,i-1 ,1,k) =   (cc(m,i-1,k,1)+  &
                 ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i  ,k,5)))+ &
                 ((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))+  &
                  (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4))))

        ch(m,i   ,1,k) =   (cc(m,i  ,k,1)+  &
                 ((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa4(i-2)*cc(m,i  ,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ &
                 ((wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))+  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4))))

        ch(m,i-1 ,3,k) =   (cc(m,i-1,k,1)+  &
            tr11*((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i  ,k,5)))+ &
            tr12*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))+  &
                  (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4))))+&
            ti11*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))-  &
                  (wa4(i-2)*cc(m,i  ,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ &
            ti12*((wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3)-   &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4))))

        ch(m,ic-1,2,k) =   (cc(m,i-1,k,1)+  &
            tr11*((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i  ,k,5)))+ &
            tr12*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))+  &
                  (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4))))-&
           (ti11*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))-  &
                  (wa4(i-2)*cc(m,i  ,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ &
            ti12*((wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))-  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4))))

        ch(m,i   ,3,k) =   (cc(m,i  ,k,1)+  &
            tr11*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa4(i-2)*cc(m,i  ,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ &
            tr12*((wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))+  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4))))+&
           (ti11*((wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i  ,k,5))-  &
                  (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2)))+ &
            ti12*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4))-  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))))

        ch(m,ic  ,2,k) =  -(cc(m,i  ,k,1)+  &
            tr11*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa4(i-2)*cc(m,i  ,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ &
            tr12*((wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))+  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4))))+&
           (ti11*((wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i  ,k,5))-  &
                  (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2)))+ &
            ti12*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4))-  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))))

        ch(m,i-1 ,5,k) =   (cc(m,i-1,k,1)+  &
            tr12*((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i  ,k,5)))+ &
            tr11*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))+  &
                  (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4))))+&
           (ti12*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))-  &
                  (wa4(i-2)*cc(m,i  ,k,5)-wa4(i-1)*cc(m,i-1,k,5)))- &
            ti11*((wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))-  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4))))

        ch(m,ic-1,4,k) =   (cc(m,i-1,k,1)+  &
            tr12*((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2))+  &
                  (wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i  ,k,5)))+ &
            tr11*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))+  &
                  (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4))))-&
           (ti12*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))-  &
                  (wa4(i-2)*cc(m,i  ,k,5)-wa4(i-1)*cc(m,i-1,k,5)))- &
            ti11*((wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))-  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4))))

        ch(m,i   ,5,k) =   (cc(m,i  ,k,1)+  &
            tr12*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa4(i-2)*cc(m,i  ,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ &
            tr11*((wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))+  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4))))+&
           (ti12*((wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i  ,k,5))-  &
                  (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2)))- &
            ti11*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4))-  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))))

        ch(m,ic  ,4,k) =  -(cc(m,i  ,k,1)+  &
            tr12*((wa1(i-2)*cc(m,i  ,k,2)-wa1(i-1)*cc(m,i-1,k,2))+  &
                  (wa4(i-2)*cc(m,i  ,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ &
            tr11*((wa2(i-2)*cc(m,i  ,k,3)-wa2(i-1)*cc(m,i-1,k,3))+  &
                  (wa3(i-2)*cc(m,i  ,k,4)-wa3(i-1)*cc(m,i-1,k,4))))+&
           (ti12*((wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i  ,k,5))-  &
                  (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i  ,k,2)))- &
            ti11*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i  ,k,4))-  &
                  (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i  ,k,3))))
      enddo
    enddo
  enddo

  END Subroutine vradf5 
!
! -----------------------------------------------------------------------
!
  Subroutine vradfg (mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,ldc,wa)

  integer                            :: mp, ido, ip, l1, idl1, ldc
  REALTYPE, dimension(ldc,ido,ip,l1) :: cc
  REALTYPE, dimension(ldc,ido,l1,ip) :: ch, c1
  REALTYPE, dimension(ldc,idl1  ,ip) :: ch2, c2
  REALTYPE, dimension(ido)           :: wa

  integer                            :: i, j, k, l, m
  integer                            :: ipph, ipp2, idp2, nbd
  integer                            :: ik, is, ic, jc, lc, idij, j2
  REALTYPE                           :: arg, dcp, dsp, ar1, ai1, ar1h
  REALTYPE                           :: dc2, ds2, ar2, ai2, ar2h


  arg = TWO*PI/ip
  dcp = COS(arg)
  dsp = SIN(arg)

  ipph = (ip+1)/2
  ipp2 =  ip+2
  idp2 = ido+2
  nbd  = (ido-1)/2

  if (ido == 1) then
    do ik = 1,idl1
      do m = 1,mp
        c2(m,ik,1) = ch2(m,ik,1)
      enddo
    enddo
  else
    do ik = 1,idl1
      do m = 1,mp
        ch2(m,ik,1) = c2(m,ik,1)
      enddo
    enddo

    do j = 2,ip
      do k = 1,l1
        do m = 1,mp
          ch(m,1,k,j) = c1(m,1,k,j)
        enddo
      enddo
    enddo

    if (nbd >  l1) then
      is = -ido
      do j = 2,ip
        is = is+ido
        do k = 1,l1
          idij = is
          do i = 3,ido,2
            idij = idij+2
            do m = 1,mp
              ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j)+  &
                              wa(idij  )*c1(m,i  ,k,j)
              ch(m,i  ,k,j) = wa(idij-1)*c1(m,i  ,k,j)-  &
                              wa(idij  )*c1(m,i-1,k,j)
            enddo
          enddo
        enddo
      enddo
    else
      is = -ido
      do j = 2,ip
        is = is+ido
        idij = is
        do i = 3,ido,2
          idij = idij+2
          do k = 1,l1
            do m = 1,mp
              ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j)+  &
                              wa(idij  )*c1(m,i  ,k,j)
              ch(m,i  ,k,j) = wa(idij-1)*c1(m,i  ,k,j)-  &
                              wa(idij  )*c1(m,i-1,k,j)
            enddo
          enddo
        enddo
      enddo
    endif

    if (nbd <  l1) then
      do j = 2,ipph
        jc = ipp2-j
        do i = 3,ido,2
          do k = 1,l1
            do m = 1,mp
              c1(m,i-1,k,j ) = ch(m,i-1,k,j )+ch(m,i-1,k,jc)
              c1(m,i-1,k,jc) = ch(m,i  ,k,j )-ch(m,i  ,k,jc)
              c1(m,i  ,k,j ) = ch(m,i  ,k,j )+ch(m,i  ,k,jc)
              c1(m,i  ,k,jc) = ch(m,i-1,k,jc)-ch(m,i-1,k,j )
            enddo
          enddo
        enddo
      enddo
    else
      do j = 2,ipph
        jc = ipp2-j
        do k = 1,l1
          do i = 3,ido,2
            do m = 1,mp
              c1(m,i-1,k,j ) = ch(m,i-1,k,j )+ch(m,i-1,k,jc)
              c1(m,i-1,k,jc) = ch(m,i  ,k,j )-ch(m,i  ,k,jc)
              c1(m,i  ,k,j ) = ch(m,i  ,k,j )+ch(m,i  ,k,jc)
              c1(m,i  ,k,jc) = ch(m,i-1,k,jc)-ch(m,i-1,k,j )
            enddo
          enddo
        enddo
      enddo
    endif
  endif

  do j = 2,ipph
    jc = ipp2-j
    do k = 1,l1
      do m = 1,mp
        c1(m,1,k,j ) = ch(m,1,k,j )+ch(m,1,k,jc)
        c1(m,1,k,jc) = ch(m,1,k,jc)-ch(m,1,k,j )
      enddo
    enddo
  enddo
 
  ar1 = ONE
  ai1 = ZERO
  do l = 2,ipph
      lc = ipp2-l
    ar1h = dcp*ar1-dsp*ai1
    ai1  = dcp*ai1+dsp*ar1
    ar1  = ar1h
    do ik = 1,idl1
      do m = 1,mp
        ch2(m,ik,l ) = c2(m,ik,1)+ar1*c2(m,ik,2 )
        ch2(m,ik,lc) =            ai1*c2(m,ik,ip)
      enddo
    enddo

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1
    do j = 3,ipph
        jc = ipp2-j
      ar2h = dc2*ar2-ds2*ai2
      ai2  = dc2*ai2+ds2*ar2
      ar2  = ar2h
      do ik = 1,idl1
        do m = 1,mp
          ch2(m,ik,l ) = ch2(m,ik,l )+ar2*c2(m,ik,j )
          ch2(m,ik,lc) = ch2(m,ik,lc)+ai2*c2(m,ik,jc)
        enddo
      enddo
    enddo
  enddo

  do j = 2,ipph
    do ik = 1,idl1
      do m = 1,mp
        ch2(m,ik,1) = ch2(m,ik,1)+c2(m,ik,j)
      enddo
    enddo
  enddo
 
  if (ido <  l1) then
    do i = 1,ido
      do k = 1,l1
        do m = 1,mp
          cc(m,i,1,k) = ch(m,i,k,1)
        enddo
      enddo
    enddo
  else
    do k = 1,l1
      do i = 1,ido
        do m = 1,mp
          cc(m,i,1,k) = ch(m,i,k,1)
        enddo
      enddo
    enddo
  endif

  do j = 2,ipph
    jc = ipp2-j
    j2 = j+j
    do k = 1,l1
      do m = 1,mp
        cc(m,ido,j2-2,k) = ch(m,1,k,j )
        cc(m,1  ,j2-1,k) = ch(m,1,k,jc)
      enddo
    enddo
  enddo

  if (ido == 1) RETURN

  if (nbd <  l1) then
    do j = 2,ipph
      jc = ipp2-j
      j2 = j+j
      do i = 3,ido,2
        ic = idp2-i
        do k = 1,l1
          do m = 1,mp
            cc(m,i-1 ,j2-1,k) = ch(m,i-1,k,j )+ch(m,i-1,k,jc)
            cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j )-ch(m,i-1,k,jc)
            cc(m,i   ,j2-1,k) = ch(m,i  ,k,j )+ch(m,i  ,k,jc)
            cc(m,ic  ,j2-2,k) = ch(m,i  ,k,jc)-ch(m,i  ,k,j )
          enddo
        enddo
      enddo
    enddo
  else
    do j = 2,ipph
      jc = ipp2-j
      j2 = j+j
      do k = 1,l1
        do i = 3,ido,2
          ic = idp2-i
          do m = 1,mp
            cc(m,i-1 ,j2-1,k) = ch(m,i-1,k,j )+ch(m,i-1,k,jc)
            cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j )-ch(m,i-1,k,jc)
            cc(m,i   ,j2-1,k) = ch(m,i  ,k,j )+ch(m,i  ,k,jc)
            cc(m,ic  ,j2-2,k) = ch(m,i  ,k,jc)-ch(m,i  ,k,j )
          enddo
        enddo
      enddo
    enddo
  endif

  END Subroutine vradfg 


!
! -----------------------------------------------------------------------
!
  Subroutine vrfftb (m,n,x,xt,ldx,WSave)
 
  integer                     :: m, n, ldx
  REALTYPE, dimension(ldx,n)  :: x, xt
  REALTYPE, dimension(n+15)   :: WSave

  if (n == 1) RETURN

  call vrftb1 (m,n,x,xt,ldx,WSave(1),WSave(n+1))

  END Subroutine vrfftb 
!
! -----------------------------------------------------------------------
!
  Subroutine vrftb1 (m,n,c,ch,ldc,wa,fac)

  integer                    :: m, n, ldc
  REALTYPE, dimension(ldc,n) :: ch, c
  REALTYPE, dimension(n)     :: wa
  REALTYPE, dimension(15)    :: fac

  integer                    :: i, j, k1
  integer                    :: nf, na, l1, iw, ip, l2, ido, idl1
  integer                    :: ix2, ix3, ix4
  REALTYPE                   :: scale

  nf = fac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1,nf
    ip   = fac(k1+2)
    l2   = ip*l1
    ido  = n/l2
    idl1 = ido*l1

    if (ip /= 4) go to 103
    ix2 = iw +ido
    ix3 = ix2+ido
    if (na /= 0) then
      call vradb4 (m,ido,l1,ch,c,ldc,wa(iw),wa(ix2),wa(ix3))
    else
      call vradb4 (m,ido,l1,c,ch,ldc,wa(iw),wa(ix2),wa(ix3))
    endif
    na = 1-na
    go to 115

 103 continue
    if (ip /= 2) go to 106
    if (na /= 0) then
      call vradb2 (m,ido,l1,ch,c,ldc,wa(iw))
    else
      call vradb2 (m,ido,l1,c,ch,ldc,wa(iw))
    endif
    na = 1-na
    go to 115

 106 continue
    if (ip /= 3) go to 109
    ix2 = iw+ido
    if (na /= 0) then
      call vradb3 (m,ido,l1,ch,c,ldc,wa(iw),wa(ix2))
    else
      call vradb3 (m,ido,l1,c,ch,ldc,wa(iw),wa(ix2))
    endif
    na = 1-na
    go to 115

 109 continue
    if (ip /= 5) go to 112
    ix2 = iw +ido
    ix3 = ix2+ido
    ix4 = ix3+ido
    if (na /= 0) then
      call vradb5 (m,ido,l1,ch,c,ldc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
    else
      call vradb5 (m,ido,l1,c,ch,ldc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
    endif
    na = 1-na
    go to 115

 112 continue
    if (na /= 0) then
      call vradbg (m,ido,ip,l1,idl1,ch,ch,ch,c ,c ,ldc,wa(iw))
    else
      call vradbg (m,ido,ip,l1,idl1,c ,c ,c ,ch,ch,ldc,wa(iw))
    endif
    if (ido == 1) na = 1-na

 115 continue
    l1 = l2
    iw = iw+(ip-1)*ido
  enddo

  scale=SQRT(ONE/n)

  if (na == 0) then
    do j = 1,n
      do i = 1,m
        c(i,j)=scale*c(i,j)
      enddo
    enddo
  else
    do j = 1,n
      do i = 1,m
        c(i,j) = scale*ch(i,j)
      enddo
    enddo
  endif

  END Subroutine vrftb1 
!
! -----------------------------------------------------------------------
!
  Subroutine vradb2 (mp,ido,l1,cc,ch,ldc,wa1)

  integer                           :: mp, ido, l1, ldc
  REALTYPE, dimension(ldc,ido,2,l1) :: cc
  REALTYPE, dimension(ldc,ido,l1,2) :: ch
  REALTYPE, dimension(ido)          :: wa1

  integer                           :: i, k, m, idp2, ic

  do k = 1,l1
    do m = 1,mp
      ch(m,1,k,1) = cc(m,1,1,k)+cc(m,ido,2,k)
      ch(m,1,k,2) = cc(m,1,1,k)-cc(m,ido,2,k)
    enddo
  enddo

  if (ido-2 <  0) RETURN

  if (ido-2 >  0) then
    idp2 = ido+2
    do k = 1,l1
      do i = 3,ido,2
        ic = idp2-i
        do m = 1,mp
          ch(m,i-1,k,1) =  cc(m,i-1,1,k)+cc(m,ic-1,2,k)
          ch(m,i  ,k,1) =  cc(m,i  ,1,k)-cc(m,ic  ,2,k)
          ch(m,i-1,k,2) = (cc(m,i-1,1,k)-cc(m,ic-1,2,k))*wa1(i-2) &
                         -(cc(m,i  ,1,k)+cc(m,ic  ,2,k))*wa1(i-1)
          ch(m,i  ,k,2) = (cc(m,i  ,1,k)+cc(m,ic  ,2,k))*wa1(i-2) &
                         +(cc(m,i-1,1,k)-cc(m,ic-1,2,k))*wa1(i-1)
        enddo
      enddo
    enddo

    if (MOD(ido,2) == 1) RETURN
  endif

  do k = 1,l1
    do m = 1,mp
      ch(m,ido,k,1) =   cc(m,ido,1,k)+cc(m,ido,1,k)
      ch(m,ido,k,2) = -(cc(m,1  ,2,k)+cc(m,1  ,2,k))
    enddo
  enddo

  END Subroutine vradb2 
!
! -----------------------------------------------------------------------
!
  Subroutine vradb3 (mp,ido,l1,cc,ch,ldc,wa1,wa2)

  integer                           :: mp, ido, l1, ldc
  REALTYPE, dimension(ldc,ido,3,l1) :: cc
  REALTYPE, dimension(ldc,ido,l1,3) :: ch
  REALTYPE, dimension(ido)          :: wa1, wa2

  integer                           :: i, k, m
  integer                           :: idp2, ic
  REALTYPE                          :: arg, taur, taui

  arg  = TWO*PI/THREE
  taur = COS(arg)
  taui = SIN(arg)

  do k = 1,l1
    do m = 1,mp
       ch(m,1,k,1) = cc(m,1,1,k)+ TWO      *cc(m,ido,2,k)
       ch(m,1,k,2) = cc(m,1,1,k)+(TWO*taur)*cc(m,ido,2,k)  &
                                -(TWO*taui)*cc(m,1  ,3,k)
       ch(m,1,k,3) = cc(m,1,1,k)+(TWO*taur)*cc(m,ido,2,k)  &
                                +(TWO*taui)*cc(m,1  ,3,k)
    enddo
  enddo

  if (ido == 1) RETURN

  idp2 = ido+2
  do k = 1,l1
    do i = 3,ido,2
      ic = idp2-i
      do m = 1,mp
        ch(m,i-1,k,1) = cc(m,i-1,1,k)+(cc(m,i-1,3,k)+cc(m,ic-1,2,k))
        ch(m,i  ,k,1) = cc(m,i  ,1,k)+(cc(m,i  ,3,k)-cc(m,ic  ,2,k))

        ch(m,i-1,k,2) = wa1(i-2)*((cc(m,i-1,1,k)+  &
                            taur* (cc(m,i-1,3,k)+cc(m,ic-1,2,k)))-  &
                           (taui* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))))  &
                       -wa1(i-1)*((cc(m,i  ,1,k)+                   &
                            taur* (cc(m,i  ,3,k)-cc(m,ic  ,2,k)))+  &
                           (taui* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))))

        ch(m,i  ,k,2) = wa1(i-2)*((cc(m,i  ,1,k)+  &
                            taur* (cc(m,i  ,3,k)-cc(m,ic  ,2,k)))+  &
                           (taui* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))))  &
                       +wa1(i-1)*((cc(m,i-1,1,k)+                   &
                            taur* (cc(m,i-1,3,k)+cc(m,ic-1,2,k)))-  &
                           (taui* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))))

        ch(m,i-1,k,3) = wa2(i-2)*((cc(m,i-1,1,k)+  &
                            taur* (cc(m,i-1,3,k)+cc(m,ic-1,2,k)))+  &
                           (taui* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))))  &
                       -wa2(i-1)*((cc(m,i  ,1,k)+                   &
                            taur* (cc(m,i  ,3,k)-cc(m,ic  ,2,k)))-  &
                           (taui* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))))

        ch(m,i  ,k,3) = wa2(i-2)*((cc(m,i  ,1,k)+  &
                            taur* (cc(m,i  ,3,k)-cc(m,ic  ,2,k)))-  &
                           (taui* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))))  &
                       +wa2(i-1)*((cc(m,i-1,1,k)+                   &
                            taur* (cc(m,i-1,3,k)+cc(m,ic-1,2,k)))+  &
                           (taui* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))))
      enddo
    enddo
  enddo

  END Subroutine vradb3 
!
! -----------------------------------------------------------------------
!
  Subroutine vradb4 (mp,ido,l1,cc,ch,ldc,wa1,wa2,wa3)

  integer                           :: mp, ido, l1, ldc
  REALTYPE, dimension(ldc,ido,4,l1) :: cc
  REALTYPE, dimension(ldc,ido,l1,4) :: ch
  REALTYPE, dimension(ido)          :: wa1, wa2, wa3

  integer                           :: i, k, m
  integer                           :: idp2, ic
  REALTYPE                          :: sqrt2

  do k = 1,l1
    do m = 1,mp
      ch(m,1,k,3) = (cc(m,1  ,1,k)+cc(m,ido,4,k))  &
                   -(cc(m,ido,2,k)+cc(m,ido,2,k))
      ch(m,1,k,1) = (cc(m,1  ,1,k)+cc(m,ido,4,k))  &
                   +(cc(m,ido,2,k)+cc(m,ido,2,k))
      ch(m,1,k,4) = (cc(m,1  ,1,k)-cc(m,ido,4,k))  &
                   +(cc(m,1  ,3,k)+cc(m,1  ,3,k))
      ch(m,1,k,2) = (cc(m,1  ,1,k)-cc(m,ido,4,k))  &
                   -(cc(m,1  ,3,k)+cc(m,1  ,3,k))
    enddo
  enddo

  if (ido-2 <  0) RETURN

  if (ido-2 >  0) then
    idp2 = ido+2
    do k = 1,l1
      do i = 3,ido,2
        ic = idp2-i
        do m = 1,mp
          ch(m,i-1,k,1) =           (cc(m,i-1,1,k)+cc(m,ic-1,4,k))  &
                                   +(cc(m,i-1,3,k)+cc(m,ic-1,2,k))

          ch(m,i  ,k,1) =           (cc(m,i  ,1,k)-cc(m,ic  ,4,k))  &
                                   +(cc(m,i  ,3,k)-cc(m,ic  ,2,k))

          ch(m,i-1,k,2) = wa1(i-2)*((cc(m,i-1,1,k)-cc(m,ic-1,4,k))- &
                                    (cc(m,i  ,3,k)+cc(m,ic  ,2,k))) &
                         -wa1(i-1)*((cc(m,i  ,1,k)+cc(m,ic  ,4,k))+ &
                                    (cc(m,i-1,3,k)-cc(m,ic-1,2,k)))

          ch(m  ,i,k,2) = wa1(i-2)*((cc(m,i  ,1,k)+cc(m,ic  ,4,k))+ &
                                    (cc(m,i-1,3,k)-cc(m,ic-1,2,k))) &
                         +wa1(i-1)*((cc(m,i-1,1,k)-cc(m,ic-1,4,k))- &
                                    (cc(m,i  ,3,k)+cc(m,ic  ,2,k)))

          ch(m,i-1,k,3) = wa2(i-2)*((cc(m,i-1,1,k)+cc(m,ic-1,4,k))- &
                                    (cc(m,i-1,3,k)+cc(m,ic-1,2,k))) &
                         -wa2(i-1)*((cc(m,i  ,1,k)-cc(m,ic  ,4,k))- &
                                    (cc(m,i  ,3,k)-cc(m,ic  ,2,k)))

          ch(m,i  ,k,3) = wa2(i-2)*((cc(m,i  ,1,k)-cc(m,ic  ,4,k))- &
                                    (cc(m,i  ,3,k)-cc(m,ic  ,2,k))) &
                         +wa2(i-1)*((cc(m,i-1,1,k)+cc(m,ic-1,4,k))- &
                                    (cc(m,i-1,3,k)+cc(m,ic-1,2,k)))

          ch(m,i-1,k,4) = wa3(i-2)*((cc(m,i-1,1,k)-cc(m,ic-1,4,k))+ &
                                    (cc(m,i  ,3,k)+cc(m,ic  ,2,k))) &
                         -wa3(i-1)*((cc(m,i  ,1,k)+cc(m,ic  ,4,k))- &
                                    (cc(m,i-1,3,k)-cc(m,ic-1,2,k)))

          ch(m,i  ,k,4) = wa3(i-2)*((cc(m,i  ,1,k)+cc(m,ic  ,4,k))- &
                                    (cc(m,i-1,3,k)-cc(m,ic-1,2,k))) &
                         +wa3(i-1)*((cc(m,i-1,1,k)-cc(m,ic-1,4,k))+ &
                                    (cc(m,i  ,3,k)+cc(m,ic  ,2,k)))
        enddo
      enddo
    enddo

    if (MOD(ido,2) == 1) RETURN
  endif

  sqrt2 = SQRT(TWO)

  do k = 1,l1
    do m = 1,mp
       ch(m,ido,k,1) =         (cc(m,ido,1,k)+cc(m,ido,3,k))  &
                              +(cc(m,ido,1,k)+cc(m,ido,3,k))
       ch(m,ido,k,2) =  sqrt2*((cc(m,ido,1,k)-cc(m,ido,3,k))  &
                              -(cc(m,1  ,2,k)+cc(m,1  ,4,k)))
       ch(m,ido,k,3) =         (cc(m,1  ,4,k)-cc(m,1  ,2,k))  &
                              +(cc(m,1  ,4,k)-cc(m,1  ,2,k))
       ch(m,ido,k,4) = -sqrt2*((cc(m,ido,1,k)-cc(m,ido,3,k))  &
                              +(cc(m,1  ,2,k)+cc(m,1  ,4,k)))
    enddo
  enddo

  END Subroutine vradb4 
!
! -----------------------------------------------------------------------
!
  Subroutine vradb5 (mp,ido,l1,cc,ch,ldc,wa1,wa2,wa3,wa4)

  integer                           :: mp, ido, l1, ldc
  REALTYPE, dimension(ldc,ido,5,l1) :: cc
  REALTYPE, dimension(ldc,ido,l1,5) :: ch
  REALTYPE, dimension(ido)          :: wa1, wa2, wa3, wa4

  integer                           :: i, k, m
  integer                           :: idp2, ic
  REALTYPE                          :: arg, tr11, ti11, tr12, ti12

  arg  = TWO*PI/FIVE
  tr11 = COS(arg)
  ti11 = SIN(arg)
  tr12 = COS(TWO*arg)
  ti12 = SIN(TWO*arg)

  do k = 1,l1
    do m = 1,mp
      ch(m,1,k,1) = (cc(m,1,1,k)+  &
                    TWO*(     cc(m,ido,2,k)+     cc(m,ido,4,k)))

      ch(m,1,k,2) = (cc(m,1,1,k)+  &
                    TWO*(tr11*cc(m,ido,2,k)+tr12*cc(m,ido,4,k))) &
                   -TWO*(ti11*cc(m,1  ,3,k)+ti12*cc(m,1  ,5,k))

      ch(m,1,k,3) = (cc(m,1,1,k)+  &
                    TWO*(tr12*cc(m,ido,2,k)+tr11*cc(m,ido,4,k))) &
                   -TWO*(ti12*cc(m,1  ,3,k)-ti11*cc(m,1  ,5,k))

      ch(m,1,k,4) = (cc(m,1,1,k)+  &
                    TWO*(tr12*cc(m,ido,2,k)+tr11*cc(m,ido,4,k))) &
                   +TWO*(ti12*cc(m,1  ,3,k)-ti11*cc(m,1  ,5,k))

      ch(m,1,k,5) = (cc(m,1,1,k)+  &
                    TWO*(tr11*cc(m,ido,2,k)+tr12*cc(m,ido,4,k))) &
                   +TWO*(ti11*cc(m,1  ,3,k)+ti12*cc(m,1  ,5,k))
    enddo
  enddo

  if (ido == 1) RETURN

  idp2 = ido+2
  do k = 1,l1
    do i = 3,ido,2
      ic = idp2-i
      do m = 1,mp
        ch(m,i-1,k,1) =            cc(m,i-1,1,k)+  &
                                  (cc(m,i-1,3,k)+cc(m,ic-1,2,k))  &
                       +          (cc(m,i-1,5,k)+cc(m,ic-1,4,k))

        ch(m,i  ,k,1) =            cc(m,i  ,1,k)+  &
                                  (cc(m,i  ,3,k)-cc(m,ic  ,2,k))  &
                       +          (cc(m,i  ,5,k)-cc(m,ic  ,4,k))

        ch(m,i-1,k,2) = wa1(i-2)*((cc(m,i-1,1,k)+  &
                            tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+  &
                            tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))  &
                          -(ti11* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))+  &
                            ti12* (cc(m,i  ,5,k)+cc(m,ic  ,4,k)))) &
                       -wa1(i-1)*((cc(m,i  ,1,k)+  &
                            tr11* (cc(m,i  ,3,k)-cc(m,ic  ,2,k))+  &
                            tr12* (cc(m,i  ,5,k)-cc(m,ic  ,4,k)))  &
                          +(ti11* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))+  &
                            ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4,k))))

        ch(m,i  ,k,2) = wa1(i-2)*((cc(m,i  ,1,k)+  &
                            tr11* (cc(m,i  ,3,k)-cc(m,ic  ,2,k))+  &
                            tr12* (cc(m,i  ,5,k)-cc(m,ic  ,4,k)))  &
                          +(ti11* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))+  &
                            ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4,k)))) &
                       +wa1(i-1)*((cc(m,i-1,1,k)+  &
                            tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+  &
                            tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))  &
                          -(ti11* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))+  &
                            ti12* (cc(m,i  ,5,k)+cc(m,ic  ,4,k))))

        ch(m,i-1,k,3) = wa2(i-2)*((cc(m,i-1,1,k)+  &
                            tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+  &
                            tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))  &
                          -(ti12* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))-  &
                            ti11* (cc(m,i  ,5,k)+cc(m,ic  ,4,k)))) &
                       -wa2(i-1)*((cc(m,i  ,1,k)+  &
                            tr12* (cc(m,i  ,3,k)-cc(m,ic  ,2,k))+  &
                            tr11* (cc(m,i  ,5,k)-cc(m,ic  ,4,k)))  &
                          +(ti12* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))-  &
                            ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4,k))))

        ch(m,i  ,k,3) = wa2(i-2)*((cc(m,i  ,1,k)+  &
                            tr12* (cc(m,i  ,3,k)-cc(m,ic  ,2,k))+  &
                            tr11* (cc(m,i  ,5,k)-cc(m,ic  ,4,k)))  &
                          +(ti12* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))-  &
                            ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4,k)))) &
                       +wa2(i-1)*((cc(m,i-1,1,k)+  &
                            tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+  &
                            tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))  &
                          -(ti12* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))-  &
                            ti11* (cc(m,i  ,5,k)+cc(m,ic  ,4,k))))

        ch(m,i-1,k,4) = wa3(i-2)*((cc(m,i-1,1,k)+  &
                            tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+  &
                            tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))  &
                          +(ti12* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))-  &
                            ti11* (cc(m,i  ,5,k)+cc(m,ic  ,4,k)))) &
                       -wa3(i-1)*((cc(m,i  ,1,k)+  &
                            tr12* (cc(m,i  ,3,k)-cc(m,ic  ,2,k))+  &
                            tr11* (cc(m,i  ,5,k)-cc(m,ic  ,4,k)))  &
                          -(ti12* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))-  &
                            ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4,k))))

        ch(m,i  ,k,4) = wa3(i-2)*((cc(m,i  ,1,k)+  &
                            tr12* (cc(m,i  ,3,k)-cc(m,ic  ,2,k))+  &
                            tr11* (cc(m,i  ,5,k)-cc(m,ic  ,4,k)))  &
                          -(ti12* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))-  &
                            ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4,k)))) &
                       +wa3(i-1)*((cc(m,i-1,1,k)+  &
                            tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+  &
                            tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))  &
                          +(ti12* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))-  &
                            ti11* (cc(m,i  ,5,k)+cc(m,ic  ,4,k))))

        ch(m,i-1,k,5) = wa4(i-2)*((cc(m,i-1,1,k)+  &
                            tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+  &
                            tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))  &
                          +(ti11* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))+  &
                            ti12* (cc(m,i  ,5,k)+cc(m,ic  ,4,k)))) &
                       -wa4(i-1)*((cc(m,i  ,1,k)+  &
                            tr11* (cc(m,i  ,3,k)-cc(m,ic  ,2,k))+  &
                            tr12* (cc(m,i  ,5,k)-cc(m,ic  ,4,k)))  &
                          -(ti11* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))+  &
                            ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4,k))))

        ch(m,i  ,k,5) = wa4(i-2)*((cc(m,i  ,1,k)+  &
                            tr11* (cc(m,i  ,3,k)-cc(m,ic  ,2,k))+  &
                            tr12* (cc(m,i  ,5,k)-cc(m,ic  ,4,k)))  &
                          -(ti11* (cc(m,i-1,3,k)-cc(m,ic-1,2,k))+  &
                            ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4,k)))) &
                       +wa4(i-1)*((cc(m,i-1,1,k)+  &
                            tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+  &
                            tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))  &
                          +(ti11* (cc(m,i  ,3,k)+cc(m,ic  ,2,k))+  &
                            ti12* (cc(m,i  ,5,k)+cc(m,ic  ,4,k))))
      enddo
    enddo
  enddo

  END Subroutine vradb5 
!
! -----------------------------------------------------------------------
!
  Subroutine vradbg (mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2, ldc,wa)

  integer                            :: mp, ido, ip, l1, idl1, ldc
  REALTYPE, dimension(ldc,ido,ip,l1) :: cc
  REALTYPE, dimension(ldc,ido,l1,ip) :: ch , c1
  REALTYPE, dimension(ldc,idl1  ,ip) :: ch2, c2
  REALTYPE, dimension(ido)           :: wa

  integer                            :: i, j, k, l, m
  integer                            :: ipph, ipp2, idp2, nbd
  integer                            :: ik, is, ic, jc, lc, idij, j2
  REALTYPE                           :: arg, dcp, dsp, ar1, ai1, ar1h
  REALTYPE                           :: dc2, ds2, ar2, ai2, ar2h

  arg = TWO*PI/ip
  dcp = COS(arg)
  dsp = SIN(arg)

  idp2 = ido+2
  nbd = (ido-1)/2
  ipp2 =  ip+2
  ipph = (ip+1)/2

  if (ido <  l1) then
    do i = 1,ido
      do k = 1,l1
        do m = 1,mp
          ch(m,i,k,1) = cc(m,i,1,k)
        enddo
      enddo
    enddo
  else
    do k = 1,l1
      do i = 1,ido
        do m = 1,mp
          ch(m,i,k,1) = cc(m,i,1,k)
        enddo
      enddo
    enddo
  endif

  do j = 2,ipph
    jc = ipp2-j
    j2 = j+j
    do k = 1,l1
      do m = 1,mp
        ch(m,1,k,j ) = cc(m,ido,j2-2,k)+cc(m,ido,j2-2,k)
        ch(m,1,k,jc) = cc(m,1  ,j2-1,k)+cc(m,1  ,j2-1,k)
      enddo
    enddo
  enddo

  if (ido /= 1) then
    if (nbd <  l1) then
      do j = 2,ipph
        jc = ipp2-j
        do i = 3,ido,2
          ic = idp2-i
          do k = 1,l1
            do m = 1,mp
              ch(m,i-1,k,j ) = cc(m,i-1,2*j-1,k)+cc(m,ic-1,2*j-2,k)
              ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k)-cc(m,ic-1,2*j-2,k)
              ch(m,i  ,k,j ) = cc(m,i  ,2*j-1,k)-cc(m,ic  ,2*j-2,k)
              ch(m,i  ,k,jc) = cc(m,i  ,2*j-1,k)+cc(m,ic  ,2*j-2,k)
            enddo
          enddo
        enddo
      enddo
    else
      do j = 2,ipph
        jc = ipp2-j
        do k = 1,l1
          do i = 3,ido,2
            ic = idp2-i
            do m = 1,mp
              ch(m,i-1,k,j ) = cc(m,i-1,2*j-1,k)+cc(m,ic-1,2*j-2,k)
              ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k)-cc(m,ic-1,2*j-2,k)
              ch(m,i  ,k,j ) = cc(m,i  ,2*j-1,k)-cc(m,ic  ,2*j-2,k)
              ch(m,i  ,k,jc) = cc(m,i  ,2*j-1,k)+cc(m,ic  ,2*j-2,k)
            enddo
          enddo
        enddo
      enddo
    endif
  endif

  ar1 = ONE
  ai1 = ZERO
  do l = 2,ipph
      lc = ipp2-l
    ar1h = dcp*ar1-dsp*ai1
    ai1  = dcp*ai1+dsp*ar1
    ar1  = ar1h
    do ik = 1,idl1
      do m = 1,mp
        c2(m,ik,l ) = ch2(m,ik,1)+ar1*ch2(m,ik,2 )
        c2(m,ik,lc) =             ai1*ch2(m,ik,ip)
      enddo
    enddo

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1
    do j = 3,ipph
      jc = ipp2-j
      ar2h = dc2*ar2-ds2*ai2
      ai2  = dc2*ai2+ds2*ar2
      ar2  = ar2h
      do ik = 1,idl1
        do m = 1,mp
          c2(m,ik,l ) = c2(m,ik,l )+ar2*ch2(m,ik,j )
          c2(m,ik,lc) = c2(m,ik,lc)+ai2*ch2(m,ik,jc)
        enddo
      enddo
    enddo
  enddo

  do j = 2,ipph
    do ik = 1,idl1
      do m = 1,mp
        ch2(m,ik,1) = ch2(m,ik,1)+ch2(m,ik,j)
      enddo
    enddo
  enddo

  do j = 2,ipph
    jc = ipp2-j
    do k = 1,l1
      do m = 1,mp
        ch(m,1,k,j ) = c1(m,1,k,j)-c1(m,1,k,jc)
        ch(m,1,k,jc) = c1(m,1,k,j)+c1(m,1,k,jc)
      enddo
    enddo
  enddo

  if (ido /= 1) then
    if (nbd <  l1) then
      do j = 2,ipph
        jc = ipp2-j
        do i = 3,ido,2
          do k = 1,l1
            do m = 1,mp
              ch(m,i-1,k,j ) = c1(m,i-1,k,j)-c1(m,i  ,k,jc)
              ch(m,i-1,k,jc) = c1(m,i-1,k,j)+c1(m,i  ,k,jc)
              ch(m,i  ,k,j ) = c1(m,i  ,k,j)+c1(m,i-1,k,jc)
              ch(m,i  ,k,jc) = c1(m,i  ,k,j)-c1(m,i-1,k,jc)
            enddo
          enddo
        enddo
      enddo
    else
      do j = 2,ipph
        jc = ipp2-j
        do k = 1,l1
          do i = 3,ido,2
            do m = 1,mp
              ch(m,i-1,k,j ) = c1(m,i-1,k,j)-c1(m,i  ,k,jc)
              ch(m,i-1,k,jc) = c1(m,i-1,k,j)+c1(m,i  ,k,jc)
              ch(m,i  ,k,j ) = c1(m,i  ,k,j)+c1(m,i-1,k,jc)
              ch(m,i  ,k,jc) = c1(m,i  ,k,j)-c1(m,i-1,k,jc)
            enddo
          enddo
        enddo
      enddo
    endif
  endif

  if (ido == 1) RETURN

  do ik = 1,idl1
    do m = 1,mp
      c2(m,ik,1) = ch2(m,ik,1)
    enddo
  enddo

  do j = 2,ip
    do k = 1,l1
      do m = 1,mp
        c1(m,1,k,j) = ch(m,1,k,j)
      enddo
    enddo
  enddo

  if (nbd >  l1) then
    is = -ido
    do j = 2,ip
      is = is+ido
      do k = 1,l1
        idij = is
        do i = 3,ido,2
          idij = idij+2
          do m = 1,mp
            c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j)-  &
                            wa(idij  )*ch(m,i  ,k,j)
            c1(m,i  ,k,j) = wa(idij-1)*ch(m,i  ,k,j)+  &
                            wa(idij  )*ch(m,i-1,k,j)
          enddo
        enddo
      enddo
    enddo
  else
    is = -ido
    do j = 2,ip
      is = is+ido
      idij = is
      do i = 3,ido,2
        idij = idij+2
        do k = 1,l1
          do m = 1,mp
            c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j)-  &
                            wa(idij  )*ch(m,i  ,k,j)
            c1(m,i  ,k,j) = wa(idij-1)*ch(m,i  ,k,j)+  &
                            wa(idij  )*ch(m,i-1,k,j)
          enddo
        enddo
      enddo
    enddo
  endif

  END Subroutine vradbg 


!
! -----------------------------------------------------------------------
!
  Subroutine vrffti (n,WSave)

  integer                   :: n
  REALTYPE, dimension(n+15) :: WSave

  if (n == 1) RETURN

  call vrfti1 (n,WSave(1),WSave(n+1))

  END Subroutine vrffti 
!
! -----------------------------------------------------------------------
!
  Subroutine vrfti1 (n,wa,fac)

  integer                 :: n
  REALTYPE, dimension(n)  :: wa
  REALTYPE, dimension(15) :: fac
  REALTYPE, dimension(4)  :: ntryh

  integer                 :: i, j, nl, nf, ntry, nq, nr, ib, is, nfm1
  integer                 :: l1, k1, ip, ld, l2, ido, ipm, ii 
  REALTYPE                :: argld, fi, arg, argh

  ntryh(1) = 4
  ntryh(2) = 2
  ntryh(3) = 3
  ntryh(4) = 5

  nl = n
  nf = 0
  j  = 0

 101 continue
  j = j+1

  if (j-4 <= 0) then
    ntry = ntryh(j)
  else
    ntry = ntry+2
  endif

 104 continue

  nq = nl/ntry
  nr = nl-ntry*nq
  if (nr /= 0) goto 101

  nf = nf+1
  fac(nf+2) = ntry
  nl = nq

  if ((ntry == 2) .AND. (nf /= 1)) then
    do i = 2,nf
      ib = nf-i+2
      fac(ib+2) = fac(ib+1)
    enddo
    fac(3) = 2
  endif

  if (nl /= 1) goto 104

  fac(1) = n
  fac(2) = nf

  argh = TWO*PI/n
  is   = 0
  nfm1 = nf-1
  l1   = 1
  if (nfm1 == 0) RETURN

  do k1 = 1,nfm1
    ip = fac(k1+2)
    ld = 0
    l2 = l1*ip
    ido = n/l2
    ipm = ip-1
    do j = 1,ipm
      ld = ld+l1
      i  = is
      argld = ld*argh
      fi = ZERO
      do ii = 3,ido,2
        i   = i+2
        fi  = fi+ONE
        arg = fi*argld

        wa(i-1) = COS(arg)
        wa(i  ) = SIN(arg)
      enddo
      is = is+ido
    enddo
    l1 = l2
  enddo

  END Subroutine vrfti1 

! =======================================================================
!
!  Subroutine VRFFTF computes the Fourier coefficients (forward
!  transform) of a number of real periodic sequences.  Specifically,
!  for each sequence the subroutine claculates the independent
!  Fourier coefficients described below at output parameter R.
! 
!  The array WSAVE which is used by subroutine VRFFTF must be
!  initialized by calling subroutine VRFFTI(N,WSAVE).
! 
! 
!  Input Parameters
! 
!  M       the number of sequences to be transformed.
! 
!  N       the length of the sequences to be transformed.  The method
!          is most efficient when N is a product of small primes,
!          however n may be any positive integer.
! 
!  R       areal two-dimensional array of size MDIMX x N containing the
!          the sequences to be transformed.  The sequences are stored
!          in the ROWS of R.  Thus, the I-th sequence to be transformed,
!          X(I,J), J=0,1,...,N-1, is stored as
! 
!               R(I,J) = X(I,J-1) , J=1, 2, . . . , N.
! 
!  RT      a real two-dimensional work array of size MDIMX x N.
! 
!  MDIMR   the row (or first) dimension of the arrays R and RT exactly
!          as they appear in the calling program.  This parameter is
!          used to specify the variable dimension of these arrays.
! 
!  WSAVE   a real one-dimensional work array which must be dimensioned
!          at least N+15.  The WSAVE array must be initialized by
!          calling subroutine VRFFTI.  A different WSAVE array must be
!          used for each different value of N.  This initialization does
!          not have to be repeated so long as N remains unchanged.  The
!          same WSAVE array may be used by VRFFTF and VRFFTB.
! 
!  Output Parameters
! 
!  R       contains the Fourier coefficients F(K) for each of the M
!          input sequences.  Specifically, row I of R, R(I,J),
!          J=1,2,..,N, contains the independent Fourier coefficients
!          F(I,K), for the I-th input sequence stored as
! 
!             R(I,1) = REAL( F(I,0) ),
!                    = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ],
! 
!             R(I,2*K) = REAL( F(I,K) )
!                      = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)]
! 
!             R(I,2*K+1) = IMAG( F(I,K) )
!                        =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)]
! 
!                   for K = 1, 2, . . . , M-1,
! 
!              and, when N is even,
! 
!              R(I,N) = REAL( F(I,N/2) ).
!                     = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ].
! 
!  WSAVE   contains results which must not be destroyed between calls
!          to VRFFTF or VRFFTB.
! 
!  -----------------------------------------------------------------
! 
!  NOTE  -  A call of VRFFTF followed immediately by a call of
!           of VRFFTB will return the original sequences R.  Thus,
!           VRFFTB is the correctly normalized inverse of VRFFTF.
!
! =======================================================================
! 
!  Subroutine VRFFTB computes the synthesis (backward transform) of a
!  number of real periodic sequences from their Fourier coefficients.
!  Specifically, for each set of independent Fourier coefficients
!  F(K), the corresponding real periodic sequence is computed.
! 
!  The array WSAVE which is used by subroutine VRFFTB must be
!  initialized by calling subroutine VRFFTI(N,WSAVE).
! 
! 
!  Input Parameters
! 
!  M       the number of sets of coefficients.
! 
!  N       the length of the sequences of coefficients to be
!          transformed.  The method is most efficient when N is a
!          product of small primes, however n may be any positive
!          integer.
! 
!  R       areal two-dimensional array of size MDIMX x N containing the
!          coefficients to be transformed.  Each set of coefficients
!          F(K), K\0,1,..,N-1, is stored as a ROW of R.  Specifically,
!          the I-th set of independent Fourier coefficients is stored
! 
!                R(I,1) = REAL( F(I,0) ),
! 
!                R(I,2*K) = REAL( F(I,K) )
! 
!                R(I,2*K+1) = IMAG( F(I,K) )
! 
!                   for K = 1, 2, . . . , M-1,
! 
!                and, when N is even,
! 
!                R(I,N) = REAL( F(I,N/2) ).
! 
!  RT      a real two-dimensional work array of size MDIMX x N.
! 
!  MDIMR   the row (or first) dimension of the arrays R and RT exactly
!          as they appear in the calling program.  This parameter is
!          used to specify the variable dimension of these arrays.
! 
!  WSAVE   a real one-dimensional work array which must be dimensioned
!          at least N+15.  The WSAVE array must be initialized by
!          calling subroutine VRFFTI.  A different WSAVE array must be
!          used for each different value of N.  This initialization does
!          not have to be repeated so long as N remains unchanged.  The
!          same WSAVE array may be used by VRFFTB and VRFFTB.
! 
!  Output Parameters
! 
!  R       contains M real periodic sequences corresponding to the given
!          coefficients.  Specifically, the I-th row of R contains the
!          real periodic sequence corresponding to the I-th set of
!          independent Fourier coefficients F(I,K) stored as
! 
!               R(I,J) = X(I,J-1) ,   J = 1, 2, . . . , N, where
! 
!               X(I,J) = SQRT(1/N)* F(I,0) + (-1)**J*F(I,N/2)
!                        + 2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
!                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
! 
!                 when N is even, and
! 
!               X(I,J) = SQRT(1/N)* F(I,0) +
!                        2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
!                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
! 
!                 when N is odd.
! 
!  WSAVE   contains results which must not be destroyed between calls
!          to VRFFTF or VRFFTB.
! 
!  -----------------------------------------------------------------
! 
!  NOTE  -  A call of VRFFTF followed immediately by a call of
!           of VRFFTB will return the original sequences R.  Thus,
!           VRFFTB is the correctly normalized inverse of VRFFTF.
! 
!
! =======================================================================
! 
!  Subroutine VCOSQF computes the forward fast Fourier cosine transform
!  of M quarter wave sequences.  That is, cosine series representations
!  with only odd wave numbers.  The transform is defined below at output
!  parameter X.
! 
!  The array WSAVE which is used by subroutine VCOSQF must be
!  initialized by calling subroutine VCOSQI(N,WSAVE).
! 
! 
!  Input Parameters
! 
!  M       the number of sequences to be transformed.
! 
!  N       the length of the sequences to be transformed.  The method
!          is most efficient when N is a product of small primes.
! 
!  X       an array of size at least X(MDIMX,N) which contains the
!          the sequences to be transformed.  The sequences are stored
!          in the ROWS of X.  Thus, the Jth sequence is stored in
!          X(J,I), I=1,..,N.
! 
!  XT      a work array of size at least XT(MDIMX,N).
! 
!  MDIMX   the first dimension of the array X exactly as it appears in
!          the calling program.
! 
!  WSAVE   a work array which must be dimensioned at least 2*N+15
!          in the program that calls VCOSQF.  The WSAVE array must be
!          initialized by calling subroutine VCOSQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.
! 
!  Output Parameters
! 
!  X       For I=1,...,N and J=1,...,M
! 
!               X(I) = ( X(1) + the sum from K=2 to K=N of
! 
!                  2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N)) )/SQRT(4*N)
! 
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of VCOSQF or VCOSQB.
! 
!  -----------------------------------------------------------------
! 
!  NOTE  -  A call of VCOSQF followed immediately by a call of
!           of VCOSQB will return the original sequences X.  Thus,
!           VCOSQB is the correctly normalized inverse VCOSQF.
! 
!  -----------------------------------------------------------------
! 
!  VCOSQF is a straightforward extension of the subprogram COSQF to
!  handle M simultaneous sequences.  COSQF was originally developed
!  by P. N. Swarztrauber of NCAR.
! 
! ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!               pp. 51-83.
!
! =========================================================================
!
!  Subroutine VCOSQB computes the backward fast Fourier cosine transform
!  of M quarter wave sequences.  That is, cosine series representations
!  with only odd wave numbers.  The transform is defined below at output
!  parameter X.
! 
!  The array WSAVE which is used by subroutine VCOSQB must be
!  initialized by calling subroutine VCOSQI(N,WSAVE).
! 
! 
!  Input Parameters
! 
!  M       the number of sequences to be transformed.
! 
!  N       the length of the sequences to be transformed.  The method
!          is most efficient when N is a product of small primes.
! 
!  X       an array of size at least X(MDIMX,N) which contains the
!          the sequences to be transformed.  The sequences are stored
!          in the ROWS of X.  Thus, the Jth sequence is stored in
!          X(J,I), I=1,..,N.
! 
!  XT      a work array of size at least XT(MDIMX,N).
! 
!  MDIMX   the first dimension of the array X exactly as it appears in
!          the calling program.
! 
!  WSAVE   a work array which must be dimensioned at least 2*N+15
!          in the program that calls VCOSQB.  The WSAVE array must be
!          initialized by calling subroutine VCOSQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.
! 
!  Output Parameters
! 
!  X       For I=1,...,N and J=1,...,M
! 
!               X(I)= the sum from K=1 to K=N of
! 
!                 4*X(K)*COS((2*K-1)*(I-1)*PI/(2*N)) /SQRT(4*N)
! 
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of VCOSQF or VCOSQB.
! 
!  -----------------------------------------------------------------
! 
!  NOTE  -  A call of VCOSQF followed immediately by a call of
!           of VCOSQB will return the original sequences X.  Thus,
!           VCOSQB is the correctly normalized inverse VCOSQF.
! 
!  -----------------------------------------------------------------
! 
!  VCOSQB is a straightforward extension of the subprogram COSQB to
!  handle M simultaneous sequences.  COSQB was originally developed
!  by P. N. Swarztrauber of NCAR.
! 
! ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!               pp. 51-83.
!
! =========================================================================
!
!  Subroutine VSINQF computes the forward fast Fourier sine transform
!  of M quarter wave sequences.  That is, sine series representations
!  with only odd wave numbers.  The transform is defined below at output
!  parameter X.
! 
!  The array WSAVE which is used by subroutine VSINQF must be
!  initialized by calling subroutine VSINQI(N,WSAVE).
! 
! 
!  Input Parameters
! 
!  M       the number of sequences to be transformed.
! 
!  N       the length of the sequences to be transformed.  The method
!          is most efficient when N is a product of small primes.
! 
!  X       an array of size at least X(MDIMX,N) which contains the
!          the sequences to be transformed.  The sequences are stored
!          in the ROWS of X.  Thus, the Jth sequence is stored in
!          X(J,I), I=1,..,N.
! 
!  XT      a work array of size at least XT(MDIMX,N).
! 
!  MDIMX   the first dimension of the array X exactly as it appears in
!          the calling program.
! 
!  WSAVE   a work array which must be dimensioned at least 2*N+15
!          in the program that calls VSINQF.  The WSAVE array must be
!          initialized by calling subroutine VSINQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.
! 
!  Output Parameters
! 
!  X       For I=1,...,N and J=1,...,M
! 
!               X(I) = ( (-1)**(I-1)*X(N)
! 
!                        + the sum from K=1 to K=N-1 of
! 
!                        2*X(K)*SIN((2*I-1)*K*PI/(2*N)) )/SQRT(4*N)
! 
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of VSINQF or VSINQB.
! 
!  -----------------------------------------------------------------
! 
!  NOTE  -  A call of VSINQF followed immediately by a call of
!           of VSINQB will return the original sequences X.  Thus,
!           VSINQB is the correctly normalized inverse VSINQF.
! 
!  -----------------------------------------------------------------
! 
!  VSINQF is a straightforward extension of the subprogram SINQF to
!  handle M simultaneous sequences.  SINQF was originally developed
!  by P. N. Swarztrauber of NCAR.
!
! =========================================================================
!
!  Subroutine VSINQB computes the backward fast Fourier sine transform
!  of M quarter wave sequences.  That is, sine series representations
!  with only odd wave numbers.  The transform is defined below at output
!  parameter X.
! 
!  The array WSAVE which is used by subroutine VSINQB must be
!  initialized by calling subroutine VSINQI(N,WSAVE).
! 
! 
!  Input Parameters
! 
!  M       the number of sequences to be transformed.
! 
!  N       the length of the sequences to be transformed.  The method
!          is most efficient when N is a product of small primes.
! 
!  X       an array of size at least X(MDIMX,N) which contains the
!          the sequences to be transformed.  The sequences are stored
!          in the ROWS of X.  Thus, the Jth sequence is stored in
!          X(J,I), I=1,..,N.
! 
!  XT      a work array of size at least XT(MDIMX,N).
! 
!  MDIMX   the first dimension of the array X exactly as it appears in
!          the calling program.
! 
!  WSAVE   a work array which must be dimensioned at least 2*N+15
!          in the program that calls VSINQB.  The WSAVE array must be
!          initialized by calling subroutine VSINQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.
! 
!  Output Parameters
! 
!  X       for I=1,...,N and J=1,...,M
! 
!               X(I)= the sum from K=1 to K=N of
! 
!                 4*X(K)*SIN((2k-1)*I*PI/(2*N)) /SQRT(4*N)
! 
! 
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of VSINQB or VSINQF.
! 
!  -----------------------------------------------------------------
! 
!  NOTE  -  A call of VSINQF followed immediately by a call of
!           of VSINQB will return the original sequences X.  Thus,
!           VSINQB is the correctly normalized inverse VSINQF.
! 
!  -----------------------------------------------------------------
! 
!  VSINQB is a straightforward extension of the subprogram SINQB to
!  handle M simultaneous sequences.  SINQB was originally developed
!  by P. N. Swarztrauber of NCAR.
! 
! ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!               pp. 51-83.
!
! =========================================================================
!
!  Subroutine VCOST computes the discrete Fourier cosine transform
!  of M even sequences X(J,I), J=1,...,M.  The transform is defined
!  below at output parameter X.
! 
!  The array WSAVE which is used by Subroutine VCOST must be
!  initialized by calling Subroutine VCOSTI(N,WSAVE).
! 
!  Input Parameters
! 
!  M       the number of sequences to be transformed.
! 
!  N       the length of the sequence to be transformed.  N must be
!          greater than 1.  The method is most efficient when N-1 is
!          is a product of small primes.
! 
!  X       an array of size at least X(MDIMX,N) which contains the
!          the sequences to be transformed.  The sequences are stored
!          in the ROWS of X.  Thus, the Jth sequence is stored in
!          X(J,I), I=1,..,N.
! 
!  XT      a work array of size at least XT(MDIMX,N-1).
! 
!  MDIMX   the first dimension of the array X exactly as it appears in
!          the calling program.
! 
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls VCOST.  The WSAVE array must be
!          initialized by calling Subroutine VCOSTI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
! 
!  Output Parameters
! 
!  X       For I=1,...,N and J=1,...,M
! 
!             X(J,I) = ( X(J,1)+(-1)**(I-1)*X(J,N)
! 
!               + the sum from K=2 to K=N-1
! 
!                 2*X(J,K)*COS((K-1)*(I-1)*PI/(N-1)) )/SQRT(2*(N-1))
! 
!  WSAVE   contains initialization calculations which must not be
!          destroyed between calls of VCOST.
! 
!  -----------------------------------------------------------------
! 
!  NOTE  -  A call of VCOST followed immediately by another call
!           of VCOST will return the original sequences X.  Thus,
!           VCOST is the correctly normalized inverse of itself.
! 
!  -----------------------------------------------------------------
! 
!  VCOST is a straightforward extension of the subprogram COST to
!  handle M simultaneous sequences.  The scaling of the sequences
!  computed by VCOST is different than that of COST.  COST was
!  originally developed by P. N. Swarztrauber of NCAR.
! 
! ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!               pp. 51-83.
! 
! =====================================================================
! 
! ***begin prologue  vsint
! ***date written   860701   (yymmdd)
! ***revision date  900509   (yymmdd)
! ***category no.  j1a3
! ***keywords  fast fourier transform, sine transform, multiple
!             sequences
! ***author  boisvert, r. f., (nist)
! ***purpose  sine transform of one or more real, odd sequences.
! ***description
! 
!  Subroutine vsint computes the discrete fourier sine transform
!  of m odd sequences x(j,i), j=1,...,m.  the transform is defined
!  below at output parameter x.
! 
!  the array WSave which is used by Subroutine vsint must be
!  initialized by calling Subroutine vsinti(n,WSave).
! 
!  input parameters
! 
!  m       the number of sequences to be transformed.
! 
!  n       the length of the sequence to be transformed.  the method
!          is most efficient when n+1 is the product of small primes.
! 
!  x       an array of size at least x(mdimx,n+1) which contains the
!          the sequences to be transformed.  the sequences are stored
!          in the rows of x.  thus, the jth sequence is stored in
!          x(j,i), i=1,..,n.  the extra column of x is used as work
!          storage.
! 
!  xt      a work array of size at least xt(mdimx,n+1).
! 
!  mdimx   the first dimension of the array x exactly as it appears in
!          the calling program.
! 
!  WSave   a work array with dimension at least int(2.5*n+15)
!          in the program that calls vsint.  the WSave array must be
!          initialized by calling Subroutine vsinti(n,WSave), and a
!          different WSave array must be used for each different
!          value of n.  this initialization does not have to be
!          repeated so long as n remains unchanged.
! 
!  output parameters
! 
!  x       for i=1,...,n and j=1,...,m
! 
!               x(j,i)= the sum from k=1 to k=n
! 
!                    2*x(j,k)*sin(k*i*pi/(n+1))/sqrt(2*(n+1))
! 
!  WSave   contains initialization calculations which must not be
!          destroyed between calls of vsint.
! 
!  -----------------------------------------------------------------
! 
!  note  -  a call of vsint followed immediately by another call
!           of vsint will return the original sequences x.  thus,
!           vsint is the correctly normalized inverse of itself.
! 
!  -----------------------------------------------------------------
! 
!  vsint is a straightforward extension of the subprogram sint to
!  handle m simultaneous sequences.  the scaling of the sequences
!  computed by vsint is different than that of sint.  sint was
!  originally developed by p. n. swarztrauber of ncar.
! 
! ***routines called  vrfftf
! -----------------------------------------------------------------
! ***begin prologue  vsinti
! ***date written   860701   (yymmdd)
! ***revision date  900509   (yymmdd)
! ***category no.  j1a3
! ***keywords  fast fourier transform, sine transform, multiple
!             sequences
! ***author  boisvert, r. f. (nist)
! ***purpose  initialize for vsint.
! ***description
! 
!  Subroutine vsinti initializes the array WSave which is used in
!  Subroutine sint.  the prime factorization of n together with
!  a tabulation of the trigonometric functions are computed and
!  stored in WSave.
! 
!  input parameter
! 
!  n       the length of the sequence to be transformed.  the method
!          is most efficient when n+1 is a product of small primes.
! 
!  output parameter
! 
!  WSave   a work array with at least int(2.5*n+15) locations.
!          different WSave arrays are required for different values
!          of n.  the contents of WSave must not be changed between
!          calls of vsint.
! 
!  -----------------------------------------------------------------
! 
!  vsinti is a straightforward extension of the subprogram sinti to
!  handle m simultaneous sequences.  sinti was originally developed
!  p. n. swarztrauber of ncar.
! 
! ***references  p. n. swarztrauber, vectorizing the ffts, in parallel
!               computations, (g. rodrigue, ed.), academic press, 1982,
!               pp. 51-83.
! ***routines called  vrffti
! ***end prologue  vsinti
#endif

END MODULE VFFT_FFT

#else

MODULE VFFT_FFT

implicit NONE

END MODULE VFFT_FFT

#endif
