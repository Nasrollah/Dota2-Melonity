#if defined (LOCAL_CRAY_FFT)

! JMFFTLIB : A library of portable fourier transform Subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.


#include "pseudopack.h"

MODULE CRAY_FFT_JM

implicit NONE

PRIVATE

PUBLIC  :: ccfft        ! 1D Complex-Complex
PUBLIC  :: ccfft2d      ! 2D Complex-Complex
PUBLIC  :: ccfft3d      ! 3D Complex-Complex
PUBLIC  :: ccfftm       ! 1D Multiple Complex-Complex

PUBLIC  :: csfft        ! 1D Complex-Real
PUBLIC  :: csfft2d      ! 2D Complex-Real
PUBLIC  :: csfft3d      ! 3D Complex-Real
PUBLIC  :: csfftm       ! 1D Multiple Complex-Real

PUBLIC  :: scfft        ! 1D Real-Complex
PUBLIC  :: scfft2d      ! 2D Real-Complex
PUBLIC  :: scfft3d      ! 3D Real-Complex
PUBLIC  :: scfftm       ! 1D Multiple Real-Complex

PUBLIC  :: cfftmlt      ! 1D Multiple Complex-Complex  (Obsolete)
PUBLIC  :: rfftmlt      ! 1D Multiple Complex-Real
PUBLIC  :: cftfax       ! Factorization for cfftmlt
PUBLIC  :: fftfax       ! Factorization for rfftmlt

CONTAINS

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine ccfft(isign,n,scale,x,y,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: n
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:2*n-1) :: x
  REALTYPE, intent(out), dimension(0:2*n-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*n-1) :: table
  REALTYPE, intent(inout), dimension(0:4*n-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i
  integer  :: ioff
  integer  :: ntable, nwork
  integer  :: nfact
  integer , dimension(0:99) :: fact
  character(len=*), parameter :: nomsp = 'CCFFT'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 4*n

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,n)
    return
  else
    nfact = nint(table(0))
    fact(0:nfact-1) = nint(table(0:nfact-1))
  end if

  ! On copie le tableau d'entree dans le tableau de travail
  ! On en profite pour premultiplier et pour tenir compte du signe
  do i = 0,n-1
    work(i)   =       scale* x(2*i)
    work(n+i) = isign*scale* x(2*i+1)
  end do
  ioff = 0

  ! On appelle le sous-programme principal
  call jmccm1d(1,n,fact,100,0,table,ntable,100+0,work,nwork,ioff)

  ! On recopie dans le tableau d'arrivee
  do i = 0,n-1
    y(2*i)   =         work(ioff  +i)
    y(2*i+1) = isign * work(ioff+n+i)
  end do

END Subroutine ccfft

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine ccfft2d(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: n, m, ldx, ldy
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:2*ldx*m-1) :: x
  REALTYPE, intent(out), dimension(0:2*ldy*m-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*(n+m)-1) :: table
  REALTYPE, intent(inout), dimension(0:512*max(n,m)-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i, j
  integer  :: ioff
  integer  :: ntable, nwork
  integer  :: nfact, mfact
  integer , dimension(0:99) :: fact
  integer  :: ideb, ifin, jdeb, jfin, n_temp, m_temp, nwork_temp
  logical  :: debut, fin
  character(len=*), parameter :: nomsp = 'CCFFT2D'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n) call jmerreur2(nomsp,9,ldx,n)
  if (ldy < n) call jmerreur2(nomsp,14,ldy,n)

  ! Gestion de table
  ntable = 100+2*(n+m)

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    call jmfact(m,fact,100,nfact,mfact)
    table(0:mfact-1) = fact(0:mfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0  ,n)
    call jmtable(table,ntable,100+2*n,m)
    return
  else
    nfact = nint(table(0))
    mfact = nint(table(nfact)) + nfact
    fact(0:mfact-1) = nint(table(0:mfact-1))
  end if

  ! Gestion de work
  !nwork = 4*n*m
  !nwork = 512*max(n,m)
  call jmgetnwork(nwork,512*max(n,m),4*max(n,m))

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  debut = .true.
  do

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    call jmdecoup(m,4*n,nwork,debut,.true.,m_temp,jdeb,jfin,nwork_temp,fin)

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! Note : On copie en transposant
    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = jdeb,jfin
        work(j-jdeb+m_temp*i)     =       scale*x(2*i  +j*2*ldx)
        work(j-jdeb+m_temp*(n+i)) = isign*scale*x(2*i+1+j*2*ldx)
      end do
    end do
    ioff = 0

    ! Attention : ioff1 est peut-etre modifie en sortie
    call jmccm1d(m_temp,n,fact,100,0    ,table,ntable,100+0  ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = jdeb,jfin
        y(2*i  +j*2*ldy) = work(ioff+j-jdeb+m_temp*i)
        y(2*i+1+j*2*ldy) = work(ioff+j-jdeb+m_temp*(n+i))
      end do
    end do

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

  ! On fait les T.F. sur l'autre dimension
  debut = .true.
  do

    ! Tronconnage
    call jmdecoup(n,4*m,nwork,debut,.true.,n_temp,ideb,ifin,nwork_temp,fin)

    ! On copie
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
        work(i-ideb+n_temp*j)     = y(2*i  +j*2*ldy)
        work(i-ideb+n_temp*(m+j)) = y(2*i+1+j*2*ldy)
      end do
    end do
    ioff = 0

    call jmccm1d(n_temp,m,fact,100,nfact,table,ntable,100+2*n,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
        y(2*i  +j*2*ldy) =       work(ioff+i-ideb+n_temp*j)
        y(2*i+1+j*2*ldy) = isign*work(ioff+i-ideb+n_temp*(m+j))
      end do
    end do

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

END Subroutine ccfft2d
 
! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine ccfft3d(isign,n,m,l,scale,x,ldx1,ldx2,y,ldy1,ldy2,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: n, m, l, ldx1, ldx2, ldy1, ldy2
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:2*ldx1*ldx2*l-1) :: x
  REALTYPE, intent(out), dimension(0:2*ldy1*ldy2*l-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*(n+m+l)-1) :: table
  REALTYPE, intent(inout), dimension(0:512*max(n,m,l)-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i, j, k
  integer  :: ioff
  integer  :: ntable, nwork
  integer  :: nfact, mfact, lfact
  integer , dimension(0:99) :: fact
  integer  :: ideb, ifin, i1, i2, jdeb, jfin, j1, j2, kdeb, kfin
  integer  :: nwork_temp, nmtemp, mltemp, nltemp, iwork
  logical  :: debut, fini
  character(len=*), parameter :: nomsp = 'CCFFT3D'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (l < 1) call jmerreur1(nomsp,8,l)
  if (ldx1 < n) call jmerreur2(nomsp,11,ldx1,n)
  if (ldx2 < m) call jmerreur2(nomsp,13,ldx2,m)
  if (ldy1 < n) call jmerreur2(nomsp,17,ldy1,n)
  if (ldy2 < m) call jmerreur2(nomsp,20,ldy2,m)

  ! Gestion de table
  ntable = 100+2*(n+m+l)

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    call jmfact(m,fact,100,nfact,mfact)
    call jmfact(l,fact,100,mfact,lfact)
    table(0:lfact-1) = fact(0:lfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,      n)
    call jmtable(table,ntable,100+2*n,    m)
    call jmtable(table,ntable,100+2*(n+m),l)
    return
  else
    nfact = nint(table(0))
    mfact = nint(table(nfact)) + nfact
    lfact = nint(table(mfact)) + mfact
    fact(0:lfact-1) = nint(table(0:lfact-1))
  end if

  ! Gestion de work
  !nwork = 4*n*m*l
  !nwork = 512*max(n,m,l)
  call jmgetnwork(nwork,512*max(n,m,l),4*max(n,m,l))

  ! On fait les T.F. sur la troisieme dimension en tronconnant sur la premiere
  ! et la deuxieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    call jmdecoup3(n,m,4*l,nwork,debut,.true.,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! On prend garde a la gestion des extremites
    do k = 0,l-1
      iwork = 0
      do j = jdeb,jfin
        i1 = 0
        i2 = n-1
        if (j == jdeb) i1 = ideb
        if (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          work(             iwork+k*nmtemp) =       scale*x(2*i  +2*ldx1*j+2*ldx1*ldx2*k)
          work(nwork_temp/4+iwork+k*nmtemp) = isign*scale*x(2*i+1+2*ldx1*j+2*ldx1*ldx2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la troisieme dimension
    ioff = 0
    call jmccm1d(nmtemp,l,fact,100,mfact,table,ntable,100+2*(n+m),work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do k = 0,l-1
      iwork = 0
      do j = jdeb,jfin
        i1 = 0
        i2 = n-1
        if (j == jdeb) i1 = ideb
        if (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          y(2*i  +2*ldy1*j+2*ldy1*ldy2*k) = work(ioff             +iwork+k*nmtemp)
          y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k) = work(ioff+nwork_temp/4+iwork+k*nmtemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

  ! On fait les T.F. sur la deuxieme dimension en tronconnant sur la premiere
  ! et la troisieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    call jmdecoup3(n,l,4*m,nwork,debut,.true.,ideb,ifin,kdeb,kfin,nltemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    do j = 0,m-1
      iwork = 0
      do k = kdeb,kfin
        i1 = 0
        i2 = n-1
        if (k == kdeb) i1 = ideb
        if (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          work(             iwork+j*nltemp) = y(2*i  +2*ldy1*j+2*ldy1*ldy2*k)
          work(nwork_temp/4+iwork+j*nltemp) = y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la deuxieme dimension
    ioff = 0
    call jmccm1d(nltemp,m,fact,100,nfact,table,ntable,100+2*n    ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do j = 0,m-1
      iwork = 0
      do k = kdeb,kfin
        i1 = 0
        i2 = n-1
        if (k == kdeb) i1 = ideb
        if (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          y(2*i  +2*ldy1*j+2*ldy1*ldy2*k) = work(ioff             +iwork+j*nltemp)
          y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k) = work(ioff+nwork_temp/4+iwork+j*nltemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  ! et la troisieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    call jmdecoup3(m,l,4*n,nwork,debut,.true.,jdeb,jfin,kdeb,kfin,mltemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    do i = 0,n-1
      iwork = 0
      do k = kdeb,kfin
        j1 = 0
        j2 = m-1
        if (k == kdeb) j1 = jdeb
        if (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = j1,j2
          work(             iwork+i*mltemp) = y(2*i  +2*ldy1*j+2*ldy1*ldy2*k)
          work(nwork_temp/4+iwork+i*mltemp) = y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la premiere dimension
    ioff = 0
    call jmccm1d(mltemp,n,fact,100,0    ,table,ntable,100+0      ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee en tenant compte du signe
    do i = 0,n-1
      iwork = 0
      do k = kdeb,kfin
        j1 = 0
        j2 = m-1
        if (k == kdeb) j1 = jdeb
        if (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = j1,j2
          y(2*i  +2*ldy1*j+2*ldy1*ldy2*k) =       work(ioff             +iwork+i*mltemp)
          y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k) = isign*work(ioff+nwork_temp/4+iwork+i*mltemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

END Subroutine ccfft3d

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine ccfftm(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: n, m, ldx, ldy
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:2*ldx*m-1) :: x
  REALTYPE, intent(out), dimension(0:2*ldy*m-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*n-1) :: table
  REALTYPE, intent(inout), dimension(0:4*n*m-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i, j
  integer  :: ioff
  integer  :: ntable, nwork
  integer  :: nfact
  integer , dimension(0:99) :: fact
  character(len=*), parameter :: nomsp = 'CCFFTM'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n) call jmerreur2(nomsp,9,ldx,n)
  if (ldy < n) call jmerreur2(nomsp,14,ldy,n)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 4*n*m

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,n)
    return
  else
    nfact = nint(table(0))
    fact(0:nfact-1) = nint(table(0:nfact-1))
  end if

  ! On copie le tableau d'entree dans le tableau de travail
  ! On en profite pour premultiplier et pour tenir compte du signe
  do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
    do j = 0,m-1
      work(j+m*i)     =       scale*x(2*i  +2*ldx*j)
      work(j+m*(n+i)) = isign*scale*x(2*i+1+2*ldx*j)
    end do
  end do

  ! On appelle le sous-programme principal
  ioff = 0
  call jmccm1d(m,n,fact,100,0,table,ntable,100+0,work,nwork,ioff)

  ! On recopie le tableau de travail dans le tableau de sortie
  do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
    do j = 0,m-1
      y(2*i  +2*ldy*j) =       work(ioff+j+m*i)
      y(2*i+1+2*ldy*j) = isign*work(ioff+j+m*(n+i))
    end do
  end do

END Subroutine ccfftm

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

  ! Arguments
  integer , intent(in) :: m, n
  integer , intent(in) :: nfact, ifact
  integer , intent(in), dimension(0:nfact-1) :: fact
  integer , intent(in) :: ntable,itable
  REALTYPE, intent(in), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: nterms
  integer  :: np, pp, lastnp, premier
  integer  :: nprod, nprod1, nprod2
  integer  :: n2, p2, n3, p3, n5, p5
  integer  :: i
  logical , save :: copyright = .false.

  ! Comme tout le monde passe par la, on envoie le copyright (une fois)
  if (.not.copyright) then
    copyright = .true.
    print *,' '
    print *,'************************************************************'
    print *,'* Portable Fourier transforms by JMFFTLIB                  *'
    print *,'* Author : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr) *'
    print *,'************************************************************'
    print *,' '
  end if

  ! On recupere les facteurs
  nterms = fact(ifact)
  n2 = fact(ifact+1)
  p2 = fact(ifact+2)
  n3 = fact(ifact+3)
  p3 = fact(ifact+4)
  n5 = fact(ifact+5)
  p5 = fact(ifact+6)
  nprod = n2*n3*n5
  do i = 7,nterms-1,2
    nprod = nprod*fact(ifact+i+1)**fact(ifact+i)
  end do

  ! On fait n3*n5 T.F. de n2 (qui est en puissances de 2)
  if (n2 /= 1) then
    call jmccm1d2(p2,n2,m*(nprod/n2),table,ntable,itable,n,n/n2,work,nwork,ioff)
  end if

  ! On transpose (on tient compte de ioff) en permutant les deux parties
  ! On en profite pour multiplier par le bon wij
  if (n2 /= 1 .and. nprod /= n2) then
    call jmcctranspcs(m,n,n2,nprod/n2,table,ntable,itable,work,nwork,ioff)
  end if

  ! On fait n5*n2 T.F. de n3 (en puissances de 3)
  if (n3 /= 1) then
    call jmccm1d3(p3,n3,m*(nprod/n3),table,ntable,itable,n,n/n3,work,nwork,ioff)
  end if

  ! On transpose (on tient compte de ioff) en permutant les deux parties
  ! On en profite pour multiplier par le bon wij
  if (n3 /= 1 .and. nprod /= n3) then
    call jmcctranspcs(m*n2,n,n3,nprod/(n2*n3), &
    & table,ntable,itable,work,nwork,ioff)
  end if

  ! On fait n2*n3 T.F. de n5 (en puissances de 5)
  if (n5 /= 1) then
    call jmccm1d5(p5,n5,m*(nprod/n5),table,ntable,itable,n,n/n5,work,nwork,ioff)
  end if

  ! On transpose s'il y a lieu (si on a fait quelque chose et s'il reste des
  ! termes a traiter
  if (n5 /= 1 .and. nprod /= n5 .and. nterms > 7) then
    call jmcctranspcs(m*n2*n3,n,n5,nprod/(n2*n3*n5), &
    & table,ntable,itable,work,nwork,ioff)
  end if
  nprod1 = m*n2*n3
  nprod2 = n2*n3*n5
  lastnp = n5

  ! On passe aux nombres premiers autres que 2, 3 et 5
  do i = 7,nterms-1,2

    pp = fact(ifact+i)
    premier = fact(ifact+i+1)
    np = premier**pp

    call jmccm1dp(premier,pp,m*(nprod/np), &
    & table,ntable,itable,n,n/np,work,nwork,ioff)

    nprod1 = nprod1 * lastnp
    nprod2 = nprod2 * np
    if (np /= 1 .and. nprod /= np .and. nterms > i+1) then
      call jmcctranspcs(nprod1,n,np,nprod/nprod2, &
      & table,ntable,itable,work,nwork,ioff)
    end if
    lastnp = np

  end do

END Subroutine jmccm1d

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmccm1d2(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  ! Arguments
  integer , intent(in) :: p, n, m
  integer , intent(in) :: ntable,itable,ntable2,mtable
  REALTYPE, intent(in), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: k, jl
  integer  :: it1,iu1,it2,iu2
  integer  :: jt1,ju1,jt2,ju2
  REALTYPE :: x1, x2, y1, y2
  REALTYPE :: c, s
  integer  :: ioff1, ioff2

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  if (mod(p,2)==0) then

    ! Si p est pair, on peut travailler entierement en base 4
    call jmccm1d4(p/2,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff1)
    ioff = ioff1

  else

    ! On fait les premieres etapes en base 4
    call jmccm1d4(p/2,n,2*m,table,ntable,itable,ntable2,mtable*2,work,nwork,ioff1)
    ioff2 = nwork/2-ioff1
    ! On fait la derniere etape en base 2
    if (m >= 16 .or. 2**(p-1) < 8) then
      do k = 0,2**(p-1)-1

        ! Les sinus et cosinus
        c = table(itable+        mtable*k)
        s = table(itable+ntable2+mtable*k)

        ! Les indices
        it1 = ioff1        +m*(k*2  )
        iu1 = ioff1+nwork/4+m*(k*2  )
        it2 = ioff1        +m*(k*2+1)
        iu2 = ioff1+nwork/4+m*(k*2+1)
        jt1 = ioff2        +m*( k          )
        ju1 = ioff2+nwork/4+m*( k          )
        jt2 = ioff2        +m*((k+2**(p-1)))
        ju2 = ioff2+nwork/4+m*((k+2**(p-1)))

!dir$ ivdep
!ocl novrec
!cdir nodep
        do jl = 0,m-1
          x1 = work(it1+jl)
          y1 = work(iu1+jl)
          x2 = work(it2+jl)
          y2 = work(iu2+jl)
          work(jt1+jl) = x1 + ( x2*c - y2*s )
          work(ju1+jl) = y1 + ( x2*s + y2*c )
          work(jt2+jl) = x1 - ( x2*c - y2*s )
          work(ju2+jl) = y1 - ( x2*s + y2*c )
        end do
      end do
    else
      do jl = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,2**(p-1)-1
          x1 = work(ioff1+jl        +m*(k*2  ))
          y1 = work(ioff1+jl+nwork/4+m*(k*2  ))
          x2 = work(ioff1+jl        +m*(k*2+1))
          y2 = work(ioff1+jl+nwork/4+m*(k*2+1))
          ! Les sinus et cosinus
          c = table(itable+        mtable*k)
          s = table(itable+ntable2+mtable*k)
          work(ioff2+jl        +m*( k          )) = x1 + ( x2*c - y2*s )
          work(ioff2+jl+nwork/4+m*( k          )) = y1 + ( x2*s + y2*c )
          work(ioff2+jl        +m*((k+2**(p-1)))) = x1 - ( x2*c - y2*s )
          work(ioff2+jl+nwork/4+m*((k+2**(p-1)))) = y1 - ( x2*s + y2*c )
        end do
      end do
    end if

    ioff = ioff2

  end if

END Subroutine jmccm1d2

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmccm1d3(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  ! Arguments
  integer , intent(in) :: p, n, m
  integer , intent(in) :: ntable,itable,ntable2, mtable
  REALTYPE, intent(in), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: i, k, jl
  REALTYPE :: x1, x2, x3, y1, y2, y3, t1, t2, t3, u1, u2, u3
  REALTYPE :: c2, s2, c3, s3
  integer  :: it1,iu1,it2,iu2,it3,iu3
  integer  :: jt1,ju1,jt2,ju2,jt3,ju3
  REALTYPE :: r,s,t,u

  ! Gestion des constantes cosinus
  REALTYPE, save :: ctwopi3, stwopi3
  logical , save :: first = .true.

  integer  :: ioff1, ioff2

  ! On recupere cos et sin de 2*pi/3
  if (first) then
    first = .false.
    ctwopi3 = -HALF
    stwopi3 =  HALF*SQRT(THREE)
  end if

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Boucle sur les etapes
  do i = 0, p-1

    if (m*3**(p-i-1) >= 16 .or. 3**i < 8) then

      do k = 0,3**i-1

        ! Les sinus et cosinus
        c2 = table(itable+        mtable*  3**(p-i-1)*k)
        s2 = table(itable+ntable2+mtable*  3**(p-i-1)*k)
        c3 = table(itable+        mtable*2*3**(p-i-1)*k)
        s3 = table(itable+ntable2+mtable*2*3**(p-i-1)*k)

        ! Les indices
        it1 = ioff1        +m*(k*3**(p-i)             )
        iu1 = ioff1+nwork/4+m*(k*3**(p-i)             )
        it2 = ioff1        +m*(k*3**(p-i)+  3**(p-i-1))
        iu2 = ioff1+nwork/4+m*(k*3**(p-i)+  3**(p-i-1))
        it3 = ioff1        +m*(k*3**(p-i)+2*3**(p-i-1))
        iu3 = ioff1+nwork/4+m*(k*3**(p-i)+2*3**(p-i-1))
        jt1 = ioff2        +m*( k        *3**(p-i-1))
        ju1 = ioff2+nwork/4+m*( k        *3**(p-i-1))
        jt2 = ioff2        +m*((k+  3**i)*3**(p-i-1))
        ju2 = ioff2+nwork/4+m*((k+  3**i)*3**(p-i-1))
        jt3 = ioff2        +m*((k+2*3**i)*3**(p-i-1))
        ju3 = ioff2+nwork/4+m*((k+2*3**i)*3**(p-i-1))

!dir$ ivdep
!ocl novrec
!cdir nodep
        do jl = 0,m*3**(p-i-1)-1

          r = (c2*work(it2+jl))-(s2*work(iu2+jl))
          s = (c2*work(iu2+jl))+(s2*work(it2+jl))
          t = (c3*work(it3+jl))-(s3*work(iu3+jl))
          u = (c3*work(iu3+jl))+(s3*work(it3+jl))
          x1 = work(it1+jl)
          y1 = work(iu1+jl)
          work(jt1+jl) = x1 + r + t
          work(ju1+jl) = y1 + s + u
          work(jt2+jl) = x1 + ctwopi3*(r+t) - stwopi3*(s-u)
          work(ju2+jl) = y1 + ctwopi3*(s+u) + stwopi3*(r-t)
          work(jt3+jl) = x1 + ctwopi3*(r+t) + stwopi3*(s-u)
          work(ju3+jl) = y1 + ctwopi3*(s+u) - stwopi3*(r-t)

        end do

      end do

    else

      do jl = 0,m*3**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,3**i-1

          t1 = work(ioff1+jl        +m*(k*3**(p-i)             ))
          u1 = work(ioff1+jl+nwork/4+m*(k*3**(p-i)             ))
          t2 = work(ioff1+jl        +m*(k*3**(p-i)+  3**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*3**(p-i)+  3**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*3**(p-i)+2*3**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*3**(p-i)+2*3**(p-i-1)))

          ! Les sinus et cosinus
          c2 = table(itable+        mtable*  3**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*  3**(p-i-1)*k)
          c3 = table(itable+        mtable*2*3**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*2*3**(p-i-1)*k)

          ! On premultiplie
          x1 = t1
          y1 = u1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3

          ! Il reste a multiplier par les twopi3
          work(ioff2+jl        +m*( k        *3**(p-i-1))) = &
          & x1 + x2                    + x3
          work(ioff2+jl+nwork/4+m*( k        *3**(p-i-1))) = &
          & y1 + y2                    + y3
          work(ioff2+jl        +m*((k+  3**i)*3**(p-i-1))) = &
          & x1 + ctwopi3*x2-stwopi3*y2 + ctwopi3*x3+stwopi3*y3
          work(ioff2+jl+nwork/4+m*((k+  3**i)*3**(p-i-1))) = &
          & y1 + ctwopi3*y2+stwopi3*x2 + ctwopi3*y3-stwopi3*x3
          work(ioff2+jl        +m*((k+2*3**i)*3**(p-i-1))) = &
          & x1 + ctwopi3*x2+stwopi3*y2 + ctwopi3*x3-stwopi3*y3
          work(ioff2+jl+nwork/4+m*((k+2*3**i)*3**(p-i-1))) = &
          & y1 + ctwopi3*y2-stwopi3*x2 + ctwopi3*y3+stwopi3*x3

        end do

      end do

    end if

    ! On alterne les offsets
    ioff1 = nwork/2-ioff1
    ioff2 = nwork/2-ioff2

  ! Fin boucle sur le nombre d'etapes
  end do

  ioff = ioff1

END Subroutine jmccm1d3

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmccm1d4(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  ! Arguments
  integer , intent(in) :: p, n, m
  integer , intent(in) :: ntable,itable,ntable2, mtable
  REALTYPE, intent(in), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: i, k, jl
  REALTYPE :: x0,x1,x2,x3,y0,y1,y2,y3,t0,t1,t2,t3,u0,u1,u2,u3
  REALTYPE :: x0px2,x0mx2,x1px3,x1mx3
  REALTYPE :: y0py2,y0my2,y1py3,y1my3
  REALTYPE :: c1, s1, c2, s2, c3, s3
  integer  :: ioff1, ioff2
  integer  :: it0,iu0,it1,iu1,it2,iu2,it3,iu3
  integer  :: jt0,ju0,jt1,ju1,jt2,ju2,jt3,ju3

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Boucle sur les etapes
  do i = 0, p-1

    if (m*4**(p-i-1) >= 16 .or. 4**i < 8) then

      do k = 0,4**i-1

        ! Les sinus et cosinus
        c1 = table(itable+        mtable*  4**(p-i-1)*k)
        s1 = table(itable+ntable2+mtable*  4**(p-i-1)*k)
        c2 = table(itable+        mtable*2*4**(p-i-1)*k)
        s2 = table(itable+ntable2+mtable*2*4**(p-i-1)*k)
        c3 = table(itable+        mtable*3*4**(p-i-1)*k)
        s3 = table(itable+ntable2+mtable*3*4**(p-i-1)*k)

        ! Les indices
        it0 = ioff1        +m*(k*4**(p-i)             )
        iu0 = ioff1+nwork/4+m*(k*4**(p-i)             )
        it1 = ioff1        +m*(k*4**(p-i)+  4**(p-i-1))
        iu1 = ioff1+nwork/4+m*(k*4**(p-i)+  4**(p-i-1))
        it2 = ioff1        +m*(k*4**(p-i)+2*4**(p-i-1))
        iu2 = ioff1+nwork/4+m*(k*4**(p-i)+2*4**(p-i-1))
        it3 = ioff1        +m*(k*4**(p-i)+3*4**(p-i-1))
        iu3 = ioff1+nwork/4+m*(k*4**(p-i)+3*4**(p-i-1))
        jt0 = ioff2        +m*( k        *4**(p-i-1))
        ju0 = ioff2+nwork/4+m*( k        *4**(p-i-1))
        jt1 = ioff2        +m*((k+  4**i)*4**(p-i-1))
        ju1 = ioff2+nwork/4+m*((k+  4**i)*4**(p-i-1))
        jt2 = ioff2        +m*((k+2*4**i)*4**(p-i-1))
        ju2 = ioff2+nwork/4+m*((k+2*4**i)*4**(p-i-1))
        jt3 = ioff2        +m*((k+3*4**i)*4**(p-i-1))
        ju3 = ioff2+nwork/4+m*((k+3*4**i)*4**(p-i-1))

!dir$ ivdep
!ocl novrec
!cdir nodep
        do jl = 0,m*4**(p-i-1)-1

          x0px2 = work(it0+jl) + (c2*work(it2+jl)-s2*work(iu2+jl))
          x0mx2 = work(it0+jl) - (c2*work(it2+jl)-s2*work(iu2+jl))
          y0py2 = work(iu0+jl) + (c2*work(iu2+jl)+s2*work(it2+jl))
          y0my2 = work(iu0+jl) - (c2*work(iu2+jl)+s2*work(it2+jl))
          x1px3 = (c1*work(it1+jl)-s1*work(iu1+jl))+(c3*work(it3+jl)-s3*work(iu3+jl))
          x1mx3 = (c1*work(it1+jl)-s1*work(iu1+jl))-(c3*work(it3+jl)-s3*work(iu3+jl))
          y1py3 = (c1*work(iu1+jl)+s1*work(it1+jl))+(c3*work(iu3+jl)+s3*work(it3+jl))
          y1my3 = (c1*work(iu1+jl)+s1*work(it1+jl))-(c3*work(iu3+jl)+s3*work(it3+jl))

          ! Il reste a multiplier par les twopi4
          work(jt0+jl) = (x0px2)+(x1px3)
          work(jt2+jl) = (x0px2)-(x1px3)
          work(ju0+jl) = (y0py2)+(y1py3)
          work(ju2+jl) = (y0py2)-(y1py3)
          work(jt1+jl) = (x0mx2)-(y1my3)
          work(jt3+jl) = (x0mx2)+(y1my3)
          work(ju1+jl) = (y0my2)+(x1mx3)
          work(ju3+jl) = (y0my2)-(x1mx3)

        end do

      end do

    else

      do jl = 0,m*4**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,4**i-1

          t0 = work(ioff1+jl        +m*(k*4**(p-i)             ))
          u0 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)             ))
          t1 = work(ioff1+jl        +m*(k*4**(p-i)+  4**(p-i-1)))
          u1 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)+  4**(p-i-1)))
          t2 = work(ioff1+jl        +m*(k*4**(p-i)+2*4**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)+2*4**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*4**(p-i)+3*4**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)+3*4**(p-i-1)))

          ! Les sinus et cosinus
          c1 = table(itable+        mtable*  4**(p-i-1)*k)
          s1 = table(itable+ntable2+mtable*  4**(p-i-1)*k)
          c2 = table(itable+        mtable*2*4**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*2*4**(p-i-1)*k)
          c3 = table(itable+        mtable*3*4**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*3*4**(p-i-1)*k)

          ! On premultiplie
          x0 = t0
          y0 = u0
          x1 = c1*t1-s1*u1
          y1 = c1*u1+s1*t1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3

          ! Il reste a multiplier par les twopi4
          work(ioff2+jl        +m*( k        *4**(p-i-1))) = x0+x1+x2+x3
          work(ioff2+jl+nwork/4+m*( k        *4**(p-i-1))) = y0+y1+y2+y3
          work(ioff2+jl        +m*((k+  4**i)*4**(p-i-1))) = x0-y1-x2+y3
          work(ioff2+jl+nwork/4+m*((k+  4**i)*4**(p-i-1))) = y0+x1-y2-x3
          work(ioff2+jl        +m*((k+2*4**i)*4**(p-i-1))) = x0-x1+x2-x3
          work(ioff2+jl+nwork/4+m*((k+2*4**i)*4**(p-i-1))) = y0-y1+y2-y3
          work(ioff2+jl        +m*((k+3*4**i)*4**(p-i-1))) = x0+y1-x2-y3
          work(ioff2+jl+nwork/4+m*((k+3*4**i)*4**(p-i-1))) = y0-x1-y2+x3

        end do

      end do

    end if

    ! On alterne les offsets
    ioff1 = nwork/2-ioff1
    ioff2 = nwork/2-ioff2

  ! Fin boucle sur le nombre d'etapes
  end do

  ioff = ioff1

END Subroutine jmccm1d4

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmccm1d5(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  ! Arguments
  integer , intent(in) :: p, n, m
  integer , intent(in) :: ntable,itable,ntable2, mtable
  REALTYPE, intent(in), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: i, k, jl
  REALTYPE :: x0,x1,x2,x3,x4,y0,y1,y2,y3,y4,t0,t1,t2,t3,t4,u0,u1,u2,u3,u4
  REALTYPE :: c1, s1, c2, s2, c3, s3, c4, s4
  integer  :: ioff1, ioff2

  ! Gestion des constantes cosinus
  REALTYPE, save :: twopi5
  REALTYPE, save :: ctwopi51, ctwopi52, ctwopi53, ctwopi54
  REALTYPE, save :: stwopi51, stwopi52, stwopi53, stwopi54
  logical , save :: first = .true.

  ! On recupere cos et sin de 2*pi/5
  if (first) then
    first = .false.
    twopi5   = TWO*PI/FIVE
    ctwopi51 = COS(      twopi5)
    stwopi51 = SIN(      twopi5)
    ctwopi52 = COS(TWO  *twopi5)
    stwopi52 = SIN(TWO  *twopi5)
    ctwopi53 = COS(THREE*twopi5)
    stwopi53 = SIN(THREE*twopi5)
    ctwopi54 = COS(FOUR *twopi5)
    stwopi54 = SIN(FOUR *twopi5)
  end if

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Boucle sur les etapes
  do i = 0, p-1

    if (m*5**(p-i-1) >= 16 .or. 5**i < 8) then

      do k = 0,5**i-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        do jl = 0,m*5**(p-i-1)-1

          t0 = work(ioff1+jl        +m*(k*5**(p-i)             ))
          u0 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)             ))
          t1 = work(ioff1+jl        +m*(k*5**(p-i)+  5**(p-i-1)))
          u1 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+  5**(p-i-1)))
          t2 = work(ioff1+jl        +m*(k*5**(p-i)+2*5**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+2*5**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*5**(p-i)+3*5**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+3*5**(p-i-1)))
          t4 = work(ioff1+jl        +m*(k*5**(p-i)+4*5**(p-i-1)))
          u4 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+4*5**(p-i-1)))

          ! Les sinus et cosinus
          c1 = table(itable+        mtable*  5**(p-i-1)*k)
          s1 = table(itable+ntable2+mtable*  5**(p-i-1)*k)
          c2 = table(itable+        mtable*2*5**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*2*5**(p-i-1)*k)
          c3 = table(itable+        mtable*3*5**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*3*5**(p-i-1)*k)
          c4 = table(itable+        mtable*4*5**(p-i-1)*k)
          s4 = table(itable+ntable2+mtable*4*5**(p-i-1)*k)

          ! On premultiplie
          x0 = t0
          y0 = u0
          x1 = c1*t1-s1*u1
          y1 = c1*u1+s1*t1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3
          x4 = c4*t4-s4*u4
          y4 = c4*u4+s4*t4

          ! Il reste a multiplier par les twopi5
          work(ioff2+jl        +m*( k        *5**(p-i-1))) =   &
          & x0 + x1                    + x2                    &
          &    + x3                    + x4
          work(ioff2+jl+nwork/4+m*( k        *5**(p-i-1))) =   &
          & y0 + y1                    + y2                    &
          &    + y3                    + y4
          work(ioff2+jl        +m*((k+  5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi51*x1 - stwopi51*y1 &
          &    + ctwopi52*x2 - stwopi52*y2 &
          &    + ctwopi53*x3 - stwopi53*y3 &
          &    + ctwopi54*x4 - stwopi54*y4
          work(ioff2+jl+nwork/4+m*((k+  5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi51*y1 + stwopi51*x1 &
          &    + ctwopi52*y2 + stwopi52*x2 &
          &    + ctwopi53*y3 + stwopi53*x3 &
          &    + ctwopi54*y4 + stwopi54*x4
          work(ioff2+jl        +m*((k+2*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi52*x1 - stwopi52*y1 &
          &    + ctwopi54*x2 - stwopi54*y2 &
          &    + ctwopi51*x3 - stwopi51*y3 &
          &    + ctwopi53*x4 - stwopi53*y4
          work(ioff2+jl+nwork/4+m*((k+2*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi52*y1 + stwopi52*x1 &
          &    + ctwopi54*y2 + stwopi54*x2 &
          &    + ctwopi51*y3 + stwopi51*x3 &
          &    + ctwopi53*y4 + stwopi53*x4
          work(ioff2+jl        +m*((k+3*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi53*x1 - stwopi53*y1 &
          &    + ctwopi51*x2 - stwopi51*y2 &
          &    + ctwopi54*x3 - stwopi54*y3 &
          &    + ctwopi52*x4 - stwopi52*y4
          work(ioff2+jl+nwork/4+m*((k+3*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi53*y1 + stwopi53*x1 &
          &    + ctwopi51*y2 + stwopi51*x2 &
          &    + ctwopi54*y3 + stwopi54*x3 &
          &    + ctwopi52*y4 + stwopi52*x4
          work(ioff2+jl        +m*((k+4*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi54*x1 - stwopi54*y1 &
          &    + ctwopi53*x2 - stwopi53*y2 &
          &    + ctwopi52*x3 - stwopi52*y3 &
          &    + ctwopi51*x4 - stwopi51*y4
          work(ioff2+jl+nwork/4+m*((k+4*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi54*y1 + stwopi54*x1 &
          &    + ctwopi53*y2 + stwopi53*x2 &
          &    + ctwopi52*y3 + stwopi52*x3 &
          &    + ctwopi51*y4 + stwopi51*x4

        end do

      end do

    else

      do jl = 0,m*5**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,5**i-1

          t0 = work(ioff1+jl        +m*(k*5**(p-i)             ))
          u0 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)             ))
          t1 = work(ioff1+jl        +m*(k*5**(p-i)+  5**(p-i-1)))
          u1 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+  5**(p-i-1)))
          t2 = work(ioff1+jl        +m*(k*5**(p-i)+2*5**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+2*5**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*5**(p-i)+3*5**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+3*5**(p-i-1)))
          t4 = work(ioff1+jl        +m*(k*5**(p-i)+4*5**(p-i-1)))
          u4 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+4*5**(p-i-1)))

          ! Les sinus et cosinus
          c1 = table(itable+        mtable*  5**(p-i-1)*k)
          s1 = table(itable+ntable2+mtable*  5**(p-i-1)*k)
          c2 = table(itable+        mtable*2*5**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*2*5**(p-i-1)*k)
          c3 = table(itable+        mtable*3*5**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*3*5**(p-i-1)*k)
          c4 = table(itable+        mtable*4*5**(p-i-1)*k)
          s4 = table(itable+ntable2+mtable*4*5**(p-i-1)*k)

          ! On premultiplie
          x0 = t0
          y0 = u0
          x1 = c1*t1-s1*u1
          y1 = c1*u1+s1*t1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3
          x4 = c4*t4-s4*u4
          y4 = c4*u4+s4*t4

          ! Il reste a multiplier par les twopi5
          work(ioff2+jl        +m*( k        *5**(p-i-1))) =   &
          & x0 + x1                    + x2                    &
          &    + x3                    + x4
          work(ioff2+jl+nwork/4+m*( k        *5**(p-i-1))) =   &
          & y0 + y1                    + y2                    &
          &    + y3                    + y4
          work(ioff2+jl        +m*((k+  5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi51*x1 - stwopi51*y1 &
          &    + ctwopi52*x2 - stwopi52*y2 &
          &    + ctwopi53*x3 - stwopi53*y3 &
          &    + ctwopi54*x4 - stwopi54*y4
          work(ioff2+jl+nwork/4+m*((k+  5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi51*y1 + stwopi51*x1 &
          &    + ctwopi52*y2 + stwopi52*x2 &
          &    + ctwopi53*y3 + stwopi53*x3 &
          &    + ctwopi54*y4 + stwopi54*x4
          work(ioff2+jl        +m*((k+2*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi52*x1 - stwopi52*y1 &
          &    + ctwopi54*x2 - stwopi54*y2 &
          &    + ctwopi51*x3 - stwopi51*y3 &
          &    + ctwopi53*x4 - stwopi53*y4
          work(ioff2+jl+nwork/4+m*((k+2*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi52*y1 + stwopi52*x1 &
          &    + ctwopi54*y2 + stwopi54*x2 &
          &    + ctwopi51*y3 + stwopi51*x3 &
          &    + ctwopi53*y4 + stwopi53*x4
          work(ioff2+jl        +m*((k+3*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi53*x1 - stwopi53*y1 &
          &    + ctwopi51*x2 - stwopi51*y2 &
          &    + ctwopi54*x3 - stwopi54*y3 &
          &    + ctwopi52*x4 - stwopi52*y4
          work(ioff2+jl+nwork/4+m*((k+3*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi53*y1 + stwopi53*x1 &
          &    + ctwopi51*y2 + stwopi51*x2 &
          &    + ctwopi54*y3 + stwopi54*x3 &
          &    + ctwopi52*y4 + stwopi52*x4
          work(ioff2+jl        +m*((k+4*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi54*x1 - stwopi54*y1 &
          &    + ctwopi53*x2 - stwopi53*y2 &
          &    + ctwopi52*x3 - stwopi52*y3 &
          &    + ctwopi51*x4 - stwopi51*y4
          work(ioff2+jl+nwork/4+m*((k+4*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi54*y1 + stwopi54*x1 &
          &    + ctwopi53*y2 + stwopi53*x2 &
          &    + ctwopi52*y3 + stwopi52*x3 &
          &    + ctwopi51*y4 + stwopi51*x4

        end do

      end do

    end if

    ! On alterne les offsets
    ioff1 = nwork/2-ioff1
    ioff2 = nwork/2-ioff2

  ! Fin boucle sur le nombre d'etapes
  end do

  ioff = ioff1

END Subroutine jmccm1d5

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmccm1dp(p,q,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  ! On fait m t.f. d'ordre q en base p (p**q)
  ! Note : n n'est pas utilise ici

  ! Arguments
  integer , intent(in) :: p, q, m
  integer , intent(in) :: ntable,itable,ntable2, mtable
  REALTYPE, intent(in), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: i, k, jl, jp, kp
  REALTYPE :: ck, sk, tk, uk, cpjk, spjk
  integer  :: pqq, pi, pqi, pqii
  integer  :: ikpr, ikpi, ijpr, ijpi
  integer  :: itr, iti, jtr, jti
  integer  :: ioff1, ioff2
  REALTYPE :: c11, c12, c21, c22

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Pour le calcul des cos(2*pi/p)
  pqq = p**(q-1)

  ! Boucle sur les etapes
  do i = 0, q-1

    pi   = p**i
    pqi  = p**(q-i)
    pqii = p**(q-i-1)

    do k = 0,pi-1

      do jp = 0,p-1

        do jl = 0,m*pqii-1

          ijpr = ioff2 + jl + m*((k+jp*pi)*pqii)
          ijpi = ijpr + nwork/4

          work(ijpr) = 0
          work(ijpi) = 0

        end do

      end do

      do kp = 0,p-1

        itr = itable+mtable*kp*pqii*k
        iti = itr + ntable2
        ck = table(itr)
        sk = table(iti)

        do jp = 0,p-1

          ! Gymanstique infernale pour recuperer cos(2*pi/p) etc
          jtr = itable+mtable*pqq*mod(jp*kp,p)
          jti = jtr + ntable2
          cpjk = table(jtr)
          spjk = table(jti)
          c11 = (cpjk*ck-spjk*sk)
          c12 = (cpjk*sk+spjk*ck)
          c21 = (cpjk*sk+spjk*ck)
          c22 = (cpjk*ck-spjk*sk)

!dir$ ivdep
!ocl novrec
!cdir nodep
          do jl = 0,m*pqii-1

            ikpr = ioff1+jl+m*(k*pqi+kp*pqii)
            ikpi = ikpr + nwork/4
            tk = work(ikpr)
            uk = work(ikpi)

            ijpr = ioff2+jl+m*((k+jp*pi)*pqii)
            ijpi = ijpr + nwork/4

            work(ijpr) = work(ijpr) + tk*c11-uk*c12
            work(ijpi) = work(ijpi) + tk*c21+uk*c22

          end do

        end do

      end do

    end do

    ! On alterne les offsets
    ioff1 = nwork/2 - ioff1
    ioff2 = nwork/2 - ioff2

  ! Fin boucle sur le nombre d'etapes
  end do

  ioff = ioff1

END Subroutine jmccm1dp

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmcctranspcs(m,n,n2,n3,table,ntable,itable,work,nwork,ioff)

  ! Cette Subroutine permute le contenu du tableau work de la facon suivante
  ! On considere qu'a l'origine ce tableau est en (m,n3,n2)
  ! On doit transposer le terme d'ordre (k,j,i) en (k,i,j)
  ! On en profite pour faire les multiplications par wij
  ! Le role de n est seulement d'attaquer les bonnes valeurs du tableau table
  ! (il y a un ecart de n entre les cos et les sin, et le stride entre
  !  les cos est de n/(n2*n3)
  ! Note : le sens de n2 et n3 ici n'a rien a voir avec celui de jmccm1d

  ! Arguments
  integer , intent(in) :: m, n
  integer , intent(in) :: n2, n3
  integer , intent(in) :: ntable,itable
  REALTYPE, intent(in),  dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout),  dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: i, j, k
  REALTYPE :: t, u, c, s
  integer  :: ioff1, ioff2
  integer  :: is

  ! Gestion des offsets
  ioff1 = ioff
  ioff2 = nwork/2-ioff

  ! Gestion du stride
  is = n/(n2*n3)

  if ( m >= 16 .or. (n2 < 8 .and. n3 < 8) ) then

    do i = 0,n2-1
      do j = 0,n3-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do k = 0,m-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        end do
      end do
    end do

  else if ( n2 >= 16 .or. n3 < 8 ) then

    do j = 0,n3-1
      do k = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n2-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        end do
      end do
    end do

  else

    do i = 0,n2-1
      do k = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,n3-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        end do
      end do
    end do

  end if

  ioff = ioff2

END Subroutine jmcctranspcs

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine cfftmlt(xr,xi,work,trigs,ifax,inc,jump,n,m,isign)

  ! Constantes pour les arguments
  integer , parameter :: nfax = 19

  ! Arguments
  integer , intent(in) :: inc, jump, n, m, isign
  REALTYPE, intent(inout), dimension(0:m*n-1) :: xr, xi
  REALTYPE, intent(out), dimension(0:4*n*m-1) :: work
  REALTYPE, intent(in), dimension(0:2*n-1) :: trigs
  integer , intent(in), dimension(0:nfax-1) :: ifax

  ! Variables locales
  integer  :: ntrigs, nwork
  integer  :: ioff
  integer  :: i, j
  character(len=*), parameter :: nomsp = 'CFFTMLT'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,1,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)

  ! Gestion de table
  ntrigs = 2*n

  ! Gestion de work
  nwork = 4*n*m

  if (isign == 1) then

    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        work(j+m*i)     = xr(jump*j+inc*i)
        work(j+m*(n+i)) = xi(jump*j+inc*i)
      end do
    end do

  else ! isign = -1

    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        work(j+m*i)     =  xr(jump*j+inc*i)
        work(j+m*(n+i)) = -xi(jump*j+inc*i)
      end do
    end do

  end if

  ! On appelle le sous-programme principal
  ioff = 0
  call jmccm1d(m,n,ifax,nfax,0,trigs,ntrigs,0,work,nwork,ioff)

  ! On recopie le tableau de travail dans le tableau de sortie
  if (isign == 1) then

    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        xr(jump*j+inc*i) = work(ioff+j+m*i)
        xi(jump*j+inc*i) = work(ioff+j+m*(n+i))
      end do
    end do

  else

    do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do j = 0,m-1
        xr(jump*j+inc*i) =  work(ioff+j+m*i)
        xi(jump*j+inc*i) = -work(ioff+j+m*(n+i))
      end do
    end do

  end if

END Subroutine cfftmlt

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine cftfax(n,ifax,trigs)

  ! Constantes pour les arguments
  integer , parameter :: nfax = 19

  ! Arguments
  integer , intent(in) :: n
  REALTYPE, dimension(0:2*n-1), intent(out) :: trigs
  integer , dimension(0:nfax-1), intent(out) :: ifax

  ! Variables locales
  integer  :: ntrigs 
  integer  :: ifin
  character(len=*), parameter :: nomsp = 'CFTFAX'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ntrigs = 2*n

  ! Factorisation de n dans ifax
  call jmfact(n,ifax,nfax,0,ifin)

  ! Preparation des tables
  call jmtable(trigs,ntrigs,0,n)

END Subroutine cftfax

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine csfft(isign,n,scale,x,y,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: n
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:2*(n/2)+1) :: x
  REALTYPE, intent(out),  dimension(0:n-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*n-1) :: table
  REALTYPE, intent(inout), dimension(0:2*n-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i
  integer  :: ntable, nwork
  integer  :: nfact
  integer , dimension(0:99) :: fact
  integer  :: dimx, debx, incx, jumpx
  integer  :: dimy, deby, incy, jumpy
  character(len=*), parameter :: nomsp = 'CSFFT'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (mod(n,2) /= 0) call jmerreur1(nomsp,24,n)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 2*n

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,n)
    return
  else
    nfact = nint(table(0))
    fact(0:nfact-1) = nint(table(0:nfact-1))
  end if

  ! On appelle le sous-programme principal
  dimx = 2*(n/2)+2 ; debx = 0 ; incx = 1 ; jumpx = 0
  dimy = n         ; deby = 0 ; incy = 1 ; jumpy = 0
  call jmcsm1dxy(1,n,fact,100,0,table,ntable,100+0,work,nwork, &
  & x,dimx,debx,incx,jumpy,y,dimy,deby,incy,jumpy,isign,scale)

END Subroutine csfft

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine csfft2d(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: m, n, ldx, ldy
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:2*ldx*m-1) :: x
  REALTYPE, intent(out), dimension(0:ldy*m-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*(n+m)-1) :: table
  REALTYPE, intent(inout), dimension(0:512*max(n,m)-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i, j
  integer  :: ntable, nwork, ioff
  integer  :: nfact, mfact
  integer , dimension(0:99) :: fact
  integer  :: ideb, ifin, jdeb, jfin, n_temp, m_temp, nwork_temp
  logical  :: debut, fin
  integer  :: dimy, deby, incy, jumpy
  integer  :: signe
  REALTYPE :: scale_temp
  logical  :: npair, mpair
  character(len=*), parameter :: nomsp = 'CSFFT2D'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Gestion de npair et mpair
  npair = ( mod(n,2) == 0 )
  mpair = ( mod(m,2) == 0 )

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n/2+1) call jmerreur2(nomsp,10,ldx,n/2+1)
  if (ldy < n+2  ) call jmerreur2(nomsp,15,ldy,n+2  )
  if ( .not.npair .and. .not.mpair) &
  & call jmerreur2(nomsp,22,n,m)

  ! Gestion de table
  ntable = 100+2*(n+m)

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    call jmfact(m,fact,100,nfact,mfact)
    table(0:mfact-1) = fact(0:mfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0  ,n)
    call jmtable(table,ntable,100+2*n,m)
    return
  else
    nfact = nint(table(0))
    mfact = nint(table(nfact)) + nfact
    fact(0:mfact-1) = nint(table(0:mfact-1))
  end if

  ! Gestion de work
  !nwork = 2*2*(n/2+1)*m
  !nwork = 512*max(n,m)
  call jmgetnwork(nwork,512*max(n,m),4*max(n,m))

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  debut = .true.
  do

    ! Tronconnage
    call jmdecoup(n/2+1,4*m,nwork,debut,mpair,n_temp,ideb,ifin,nwork_temp,fin)

    ! On copie le tableau d'entree dans le tableau de travail sans permuter
    ! les dimensions (n en premier) pour faire d'abord la tf sur m
    ! On en profite pour premultiplier et pour tenir compte du signe
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
        work(             n_temp*j+i-ideb) =       scale*x(2*i  +2*ldx*j)
        work(nwork_temp/4+n_temp*j+i-ideb) = isign*scale*x(2*i+1+2*ldx*j)
      end do
    end do
    ioff = 0

    ! On fait la FFT complexe -> complexe sur la deuxieme dimension (m)
    call jmccm1d(n_temp,m,fact,100,nfact,table,ntable,100+2*n,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
         y(2*i  +ldy*j) = work(ioff+             n_temp*j+i-ideb)
         y(2*i+1+ldy*j) = work(ioff+nwork_temp/4+n_temp*j+i-ideb)
      end do
    end do

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

  ! On fait les T.F. sur l'autre dimension
  debut = .true.
  do

    ! Tronconnage
    call jmdecoup(m,2*n,nwork,debut,npair,m_temp,jdeb,jfin,nwork_temp,fin)

    ! On fait la FFT complexe -> reel sur le premiere dimension (n)
    dimy = ldy*m   ; deby = jdeb*ldy   ; incy = 1 ; jumpy = ldy
    signe = 1
    scale_temp = ONE
    call jmcsm1dxy(m_temp,n,fact,100,0,table,ntable,100+0, &
    & work,nwork_temp, &
    & y,dimy,deby,incy,jumpy,y,dimy,deby,incy,jumpy,signe,scale_temp)

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

END Subroutine csfft2d

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine csfft3d(isign,n,m,l,scale,x,ldx1,ldx2,y,ldy1,ldy2,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: m, n, l, ldx1, ldx2, ldy1, ldy2
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:2*ldx1*ldx2*l-1) :: x
  REALTYPE, intent(out), dimension(0:ldy1*ldy2*l-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*(n+m+l)-1) :: table
  REALTYPE, intent(inout), dimension(0:512*max(n,m,l)-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i, j, k
  integer  :: ntable, nwork, ioff
  integer  :: nfact, mfact, lfact
  integer , dimension(0:99) :: fact
  integer  :: ideb, ifin, jdeb, jfin, kdeb, kfin, i1, i2, j1, j2
  integer  :: nltemp, nmtemp, mltemp, nwork_temp, iwork
  logical  :: debut, fini
  logical  :: npair, mpair, lpair
  character(len=*), parameter :: nomsp = 'CSFFT3D'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Gestion de npair
  npair = (mod(n,2) == 0)
  mpair = (mod(m,2) == 0)
  lpair = (mod(l,2) == 0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (l < 1) call jmerreur1(nomsp,8,l)
  if (ldx1 < n/2+1) call jmerreur2(nomsp,12,ldx1,n/2+1)
  if (ldy1 < n+2  ) call jmerreur2(nomsp,18,ldy1,n+2  )
  if (ldx2 < m) call jmerreur2(nomsp,13,ldx2,m)
  if (ldy2 < m) call jmerreur2(nomsp,20,ldy2,m)
  if (.not.mpair .and. .not.npair .and. .not.lpair) &
  & call jmerreur3(nomsp,25,n,m,l)

  ! Gestion de table
  ntable = 100+2*(n+m+l)

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    call jmfact(m,fact,100,nfact,mfact)
    call jmfact(l,fact,100,mfact,lfact)
    table(0:lfact-1) = fact(0:lfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0      ,n)
    call jmtable(table,ntable,100+2*n    ,m)
    call jmtable(table,ntable,100+2*(n+m),l)
    return
  else
    nfact = nint(table(0))
    mfact = nint(table(nfact)) + nfact
    lfact = nint(table(mfact)) + mfact
    fact(0:lfact-1) = nint(table(0:lfact-1))
  end if

  ! Gestion de work
  !nwork = 2*2*(n/2+1)*m*l
  !nwork = 512*max(n,m,l)
  call jmgetnwork(nwork,512*max(n,m,l),4*max(n,m,l))

  ! On fait les T.F. sur la troisieme dimension en tronconnant sur la premiere
  ! et la deuxieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    call jmdecoup3(n/2+1,m,4*l,nwork,debut,.true.,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! On prend garde a la gestion des extremites
    do k = 0,l-1
      iwork = 0
      do j = jdeb,jfin
        i1 = 0
        i2 = n/2
        if (j == jdeb) i1 = ideb
        if (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          work(             iwork+k*nmtemp) = &
          &       scale*x(2*i  +2*ldx1*j+2*ldx1*ldx2*k)
          work(nwork_temp/4+iwork+k*nmtemp) = &
          & isign*scale*x(2*i+1+2*ldx1*j+2*ldx1*ldx2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la troisieme dimension
    ioff = 0
    call jmccm1d(nmtemp,l,fact,100,mfact,table,ntable,100+2*(n+m),work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do k = 0,l-1
      iwork = 0
      do j = jdeb,jfin
        i1 = 0
        i2 = n/2
        if (j == jdeb) i1 = ideb
        if (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          y(2*i  +ldy1*j+ldy1*ldy2*k) = work(ioff+             iwork+k*nmtemp)
          y(2*i+1+ldy1*j+ldy1*ldy2*k) = work(ioff+nwork_temp/4+iwork+k*nmtemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

  ! On fait les T.F. sur la deuxieme dimension en tronconnant sur la premiere
  ! et la troisieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    call jmdecoup3(n/2+1,l,4*m,nwork,debut,.true.,ideb,ifin,kdeb,kfin,nltemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    do j = 0,m-1
      iwork = 0
      do k = kdeb,kfin
        i1 = 0
        i2 = n/2
        if (k == kdeb) i1 = ideb
        if (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          work(             iwork+j*nltemp) = &
          & y(2*i  +ldy1*j+ldy1*ldy2*k)
          work(nwork_temp/4+iwork+j*nltemp) = &
          & y(2*i+1+ldy1*j+ldy1*ldy2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la deuxieme dimension
    ioff = 0
    call jmccm1d(nltemp,m,fact,100,nfact,table,ntable,100+2*n    ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do j = 0,m-1
      iwork = 0
      do k = kdeb,kfin
        i1 = 0
        i2 = n/2
        if (k == kdeb) i1 = ideb
        if (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          y(2*i  +ldy1*j+ldy1*ldy2*k) = work(ioff             +iwork+j*nltemp)
          y(2*i+1+ldy1*j+ldy1*ldy2*k) = work(ioff+nwork_temp/4+iwork+j*nltemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  ! et la troisieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    call jmdecoup3(m,l,4*(n/2+1),nwork,debut,npair,jdeb,jfin,kdeb,kfin,mltemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    do i = 0,n/2
      iwork = 0
      do k = kdeb,kfin
        j1 = 0
        j2 = m-1
        if (k == kdeb) j1 = jdeb
        if (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = j1,j2
          work(             iwork+i*mltemp) = y(2*i  +ldy1*j+ldy1*ldy2*k)
          work(nwork_temp/4+iwork+i*mltemp) = y(2*i+1+ldy1*j+ldy1*ldy2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la premiere dimension
    ioff = 0
    call jmcsm1d(mltemp,n,fact,100,0    ,table,ntable,100+0      ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do i = 0,n-1
      iwork = 0
      do k = kdeb,kfin
        j1 = 0
        j2 = m-1
        if (k == kdeb) j1 = jdeb
        if (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = j1,j2
          y(i+ldy1*j+ldy1*ldy2*k) = work(ioff+iwork+i*mltemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

END Subroutine csfft3d

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine csfftm(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: m, n, ldx, ldy
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:2*ldx-1,0:m-1) :: x
  REALTYPE, intent(out), dimension(0:ldy-1,0:m-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*n-1) :: table
  REALTYPE, intent(inout), dimension(0:2*n*m-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i
  integer  :: ntable, nwork
  integer  :: nfact
  integer , dimension(0:99) :: fact
  integer  :: dimx, debx, incx, jumpx
  integer  :: dimy, deby, incy, jumpy
  character(len=*), parameter :: nomsp = 'CSFFTM'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n/2+1) call jmerreur2(nomsp,10,ldx,n/2+1)
  if (ldy < n    ) call jmerreur2(nomsp,14,ldy,n    )
  if (mod(n,2) /= 0 .and. mod(m,2) /= 0) &
  & call jmerreur2(nomsp,22,n,m)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 2*n*m

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,n)
    return
  else
    nfact = nint(table(0))
    fact(0:nfact-1) = nint(table(0:nfact-1))
  end if

  ! On appelle le sous-programme principal
  dimx = 2*ldx*m ; debx = 0 ; incx = 1 ; jumpx = 2*ldx
  dimy = ldy*m   ; deby = 0 ; incy = 1 ; jumpy = ldy
  call jmcsm1dxy(m,n,fact,100,0,table,ntable,100+0,work,nwork, &
  & x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

END Subroutine csfftm

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmcsm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

  ! Arguments
  integer , intent(in) :: m, n
  integer , intent(in) :: nfact, ifact
  integer , intent(inout), dimension(0:nfact-1) :: fact
  integer , intent(in) :: ntable,itable
  REALTYPE, intent(inout), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: ioff1, ioff2
  integer  :: i, j
  REALTYPE :: t, u, v, w, tt, uu, vv, ww
  REALTYPE :: c, s
  integer  :: it

  ! Gestion de work
  ioff1 = ioff
  ioff2 = nwork/2 - ioff1

  ! On doit faire m T.F. complexes -> reelles de longueur n
  ! Si m est pair
  if (mod(m,2) == 0) then

    ! On distribue

    if (m/2 >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          it = n-i
          if (i == 0) it = 0
          t = work(ioff1        +i*m+j    )
          u = work(ioff1+nwork/4+i*m+j    )
          v = work(ioff1        +i*m+j+m/2)
          w = work(ioff1+nwork/4+i*m+j+m/2)
          work(ioff2        + i*m/2+j) = (t-w)
          work(ioff2+nwork/4+ i*m/2+j) = (u+v)
          work(ioff2        +it*m/2+j) = (t+w)
          work(ioff2+nwork/4+it*m/2+j) = (v-u)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2
          it = n-i
          if (i == 0) it = 0
          t = work(ioff1        +i*m+j    )
          u = work(ioff1+nwork/4+i*m+j    )
          v = work(ioff1        +i*m+j+m/2)
          w = work(ioff1+nwork/4+i*m+j+m/2)
          work(ioff2        + i*m/2+j) = (t-w)
          work(ioff2+nwork/4+ i*m/2+j) = (u+v)
          work(ioff2        +it*m/2+j) = (t+w)
          work(ioff2+nwork/4+it*m/2+j) = (v-u)
        end do
      end do

    end if

    ! On fait m/2 t.f. complexes -> complexes de longueur n
    call jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    ioff1 = nwork/2 - ioff2

    ! On reconstitue

    if (m/2 >= 16 .or. n < 8) then

      do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          work(ioff1+i*m+j    ) = work(ioff2        +i*m/2+j)
          work(ioff1+i*m+j+m/2) = work(ioff2+nwork/4+i*m/2+j)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n-1
          work(ioff1+i*m+j    ) = work(ioff2        +i*m/2+j)
          work(ioff1+i*m+j+m/2) = work(ioff2+nwork/4+i*m/2+j)
        end do
      end do

    end if

  ! Si m n'est pas pair mais que n l'est
  else if (mod(n,2) == 0) then

    ! On distribue

    if (m >= 16 .or. n/2 < 8) then

      do i = 0,n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =  work(ioff1        +      i*m+j)
          u = -work(ioff1+nwork/4+      i*m+j)
          v =  work(ioff1        +(n/2-i)*m+j)
          w = -work(ioff1+nwork/4+(n/2-i)*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2
          uu = (u-w)/2
          vv = (c*(t-v)+s*(u+w))/2
          ww = (c*(u+w)-s*(t-v))/2
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(ioff2        +m*i+j) =  2*(tt-ww)
          work(ioff2+nwork/4+m*i+j) = -2*(uu+vv)
        end do
      end do

    else

      do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =  work(ioff1        +      i*m+j)
          u = -work(ioff1+nwork/4+      i*m+j)
          v =  work(ioff1        +(n/2-i)*m+j)
          w = -work(ioff1+nwork/4+(n/2-i)*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2
          uu = (u-w)/2
          vv = (c*(t-v)+s*(u+w))/2
          ww = (c*(u+w)-s*(t-v))/2
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(ioff2        +m*i+j) =  2*(tt-ww)
          work(ioff2+nwork/4+m*i+j) = -2*(uu+vv)
        end do
      end do

    end if

    ! On fait m t.f. complexes de taille n/2
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    call jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1
    ioff1 = nwork/2 - ioff2

    ! On reconstitue

    if (m >= 16 .or. n/2 < 8) then

      do i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0, m-1
          ! Note : le signe - vient de l'inversion
          work(ioff1+m*(2*i  )+j) =  work(ioff2        +m*i+j)
          work(ioff1+m*(2*i+1)+j) = -work(ioff2+nwork/4+m*i+j)
        end do
      end do

    else

      do j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0, n/2-1
          ! Note : le signe - vient de l'inversion
          work(ioff1+m*(2*i  )+j) =  work(ioff2        +m*i+j)
          work(ioff1+m*(2*i+1)+j) = -work(ioff2+nwork/4+m*i+j)
        end do
      end do

    end if

  end if

  ioff = ioff1

END Subroutine jmcsm1d

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

! Variante de jmcsm1d ou on fournit x en entree et en sortie
Subroutine jmcsm1dxy(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

  ! Arguments
  integer , intent(in) :: m, n
  integer , intent(in) :: nfact, ifact
  integer , intent(inout), dimension(0:nfact-1) :: fact
  integer , intent(in) :: ntable,itable
  REALTYPE, intent(in), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(in) :: dimx, debx, incx, jumpx
  integer , intent(in) :: dimy, deby, incy, jumpy
  REALTYPE, intent(in),  dimension(0:dimx-1) :: x
  REALTYPE, intent(out), dimension(0:dimy-1) :: y
  integer , intent(in) :: isign
  REALTYPE, intent(in) :: scale

  ! Variables locales
  integer  :: i, j
  REALTYPE :: t, u, v, w, tt, uu, vv, ww
  REALTYPE :: c, s
  integer  :: it
  integer  :: ioff

  ! On doit faire m T.F. complexes -> reelles de longueur n
  ! Si m est pair
  if (mod(m,2) == 0) then

    ! On distribue

    if (m/2 >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          it = n-i
          if (i == 0) it = 0
          t =       scale*x(debx+incx*(2*i  )+jumpx*(j    ))
          u = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j    ))
          v =       scale*x(debx+incx*(2*i  )+jumpx*(j+m/2))
          w = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j+m/2))
          work(         i*m/2+j) = (t-w)
          work(nwork/4+ i*m/2+j) = (u+v)
          work(        it*m/2+j) = (t+w)
          work(nwork/4+it*m/2+j) = (v-u)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2
          it = n-i
          if (i == 0) it = 0
          t =       scale*x(debx+incx*(2*i  )+jumpx*(j    ))
          u = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j    ))
          v =       scale*x(debx+incx*(2*i  )+jumpx*(j+m/2))
          w = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j+m/2))
          work(         i*m/2+j) = (t-w)
          work(nwork/4+ i*m/2+j) = (u+v)
          work(        it*m/2+j) = (t+w)
          work(nwork/4+it*m/2+j) = (v-u)
        end do
      end do

    end if

    ! On fait m/2 t.f. complexes -> complexes de longueur n
    ioff = 0
    call jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

    ! On reconstitue

    if (m/2 >= 16 .or. n < 8) then

      do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          y(deby+jumpy*(j    )+incy*i) = work(ioff        +i*m/2+j)
          y(deby+jumpy*(j+m/2)+incy*i) = work(ioff+nwork/4+i*m/2+j)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n-1
          y(deby+jumpy*(j    )+incy*i) = work(ioff        +i*m/2+j)
          y(deby+jumpy*(j+m/2)+incy*i) = work(ioff+nwork/4+i*m/2+j)
        end do
      end do

    end if

  ! Si m n'est pas pair mais que n l'est
  else if (mod(n,2) == 0) then

    ! On distribue

    if (m >= 16 .or. n/2 < 8) then

      do i = 0,n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =        scale*x(debx+(2*i)  *incx+j*jumpx)
          u = -isign*scale*x(debx+(2*i+1)*incx+j*jumpx)
          v =        scale*x(debx+(2*(n/2-i)  )*incx+j*jumpx)
          w = -isign*scale*x(debx+(2*(n/2-i)+1)*incx+j*jumpx)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2
          uu = (u-w)/2
          vv = (c*(t-v)+s*(u+w))/2
          ww = (c*(u+w)-s*(t-v))/2
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(        m*i+j) =  2*(tt-ww)
          work(nwork/4+m*i+j) = -2*(uu+vv)
        end do
      end do

    else

      do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =        scale*x(debx+(2*i)  *incx+j*jumpx)
          u = -isign*scale*x(debx+(2*i+1)*incx+j*jumpx)
          v =        scale*x(debx+(2*(n/2-i)  )*incx+j*jumpx)
          w = -isign*scale*x(debx+(2*(n/2-i)+1)*incx+j*jumpx)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2
          uu = (u-w)/2
          vv = (c*(t-v)+s*(u+w))/2
          ww = (c*(u+w)-s*(t-v))/2
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(        m*i+j) =  2*(tt-ww)
          work(nwork/4+m*i+j) = -2*(uu+vv)
        end do
      end do

    end if

    ! On fait m t.f. complexes de taille n/2
    ioff = 0
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    call jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1

    ! On reconstitue

    if (m >= 16 .or. n/2 < 8) then

      do i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0, m-1
          ! Note : le signe - vient de l'inversion
          y(deby+incy*(2*i  )+jumpy*j) =  work(ioff        +m*i+j)
          y(deby+incy*(2*i+1)+jumpy*j) = -work(ioff+nwork/4+m*i+j)
        end do
      end do

    else

!dir$ ivdep
      do j = 0, m-1
!ocl novrec
!cdir nodep
        do i = 0, n/2-1
          ! Note : le signe - vient de l'inversion
          y(deby+incy*(2*i  )+jumpy*j) =  work(ioff        +m*i+j)
          y(deby+incy*(2*i+1)+jumpy*j) = -work(ioff+nwork/4+m*i+j)
        end do
      end do

    end if

  end if

END Subroutine jmcsm1dxy

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

! Pour tronconner de facon a tenir dans le nwork disponible

Subroutine jmdecoup(n,nr,nwork,debut,mpair,n_temp,ideb,ifin,nwork_temp,fin)

  ! Arguments
  integer , intent(in) :: n, nr, nwork
  logical , intent(in) :: debut, mpair
  integer , intent(out) :: n_temp, ideb, nwork_temp
  integer , intent(inout) :: ifin
  logical , intent(out) :: fin

  ! Variables locales
  character(len=*), parameter :: nomsp = 'JMDECOUP'

  ! n*nr est l'espace total qu'il faudrait pour work.
  ! Malheureusement, on n'a que nwork au plus
  ! On va donc decouper n en morceaux pour tenir

  ! Gestion de debut
  if (debut) then
    ideb = 0
  else
    ideb = ifin+1
  end if

  ! Gestion de n_temp et ifin
  n_temp = nwork/nr
  ! Si m impair, on doit eviter que n_temp soit impair (routine cs et sc)
  if (.not.mpair .and. mod(n_temp,2) /= 0) n_temp = n_temp-1
  ifin = min(ideb+n_temp-1,n-1)
  n_temp = ifin-ideb+1
  ! On verifie que n_temp n'est pas nul
  if (n_temp <= 0) then
    call jmerreur3(nomsp,6,n,nr,nwork)
  end if
  nwork_temp = n_temp*nr

  ! Gestion de fin
  if (ifin == n-1) then
    fin = .true.
  else
    fin = .false.
  end if

END Subroutine jmdecoup

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

! Pour tronconner en dimension 3 de facon a tenir dans le nwork disponible

Subroutine jmdecoup3(n,m,nmr,nwork,debut,lpair,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)

  ! Arguments
  integer , intent(in) :: n, m, nmr, nwork
  logical , intent(in) :: debut, lpair
  integer , intent(out) :: nmtemp, nwork_temp
  integer , intent(out)   :: ideb, jdeb
  integer , intent(inout) :: ifin, jfin
  logical , intent(out) :: fini

  ! Variables locales
  integer  :: ijdeb, ijfin
  character(len=*), parameter :: nomsp = 'JMDECOUP3'

  ! n*m*nr est l'espace total qu'il faudrait pour work.
  ! Malheureusement, on n'a que nwork au plus
  ! On va donc decouper n et m en morceaux pour tenir

  ! Gestion de debut
  if (debut) then
    ideb = 0
    jdeb = 0
  else
    if (ifin < n-1) then
      ideb = ifin+1
      jdeb = jfin
    else
      ideb = 0
      jdeb = jfin+1
    end if
  end if

  ! Gestion de nmtemp
  nmtemp = nwork/nmr
  ! Si l impair, on doit eviter que nmtemp soit impair (routine cs et sc)
  if (.not.lpair .and. mod(nmtemp,2) /= 0) nmtemp = nmtemp-1
  ! Pour simplifier, on passe par des indices 2d
  ijdeb = ideb+jdeb*n
  ijfin = min(ijdeb+nmtemp-1,n*m-1)
  nmtemp = ijfin-ijdeb+1
  ! On verifie que nmtemp n'est pas nul
  if (nmtemp <= 0) then
    call jmerreur4(nomsp,6,n,m,nmr,nwork)
  end if
  nwork_temp = nmtemp*nmr

  ! On deduit ifin et jfin de ijfin
  jfin = ijfin/n
  ifin = ijfin-n*jfin

  ! Gestion de fin
  if (ifin == n-1 .and. jfin == m-1) then
    fini = .true.
  else
    fini = .false.
  end if

END Subroutine jmdecoup3

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmerreur1(nomsp,code,var1)

  ! Arguments
  character(len=*), intent(in) :: nomsp
  integer , intent(in) :: code
  integer , intent(in) :: var1

  ! Variables locales
  integer  :: arret
  character(len=80) :: message

  call jmgetstop(arret)
  if (arret == 1) then
    call jmgetmessage(code,message)
    print *,'JMFFT Erreur dans ',trim(nomsp),' : ',trim(message), &
    & ' (',var1,')'
    stop 1
  else
    call jmsetcode(code)
  end if

END Subroutine jmerreur1

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmerreur2(nomsp,code,var1,var2)

  ! Arguments
  character(len=*), intent(in) :: nomsp
  integer , intent(in) :: code
  integer , intent(in) :: var1, var2

  ! Variables locales
  integer  :: arret
  character(len=80) :: message

  call jmgetstop(arret)
  if (arret == 1) then
    call jmgetmessage(code,message)
    print *,'JMFFT Erreur dans ',trim(nomsp),' : ',trim(message), &
    & ' (',var1,var2,')'
    stop 1
  else
    call jmsetcode(code)
  end if

END Subroutine jmerreur2

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmerreur3(nomsp,code,var1,var2,var3)

  ! Arguments
  character(len=*), intent(in) :: nomsp
  integer , intent(in) :: code
  integer , intent(in) :: var1, var2, var3

  ! Variables locales
  integer  :: arret
  character(len=80) :: message

  call jmgetstop(arret)
  if (arret == 1) then
    call jmgetmessage(code,message)
    print *,'JMFFT Erreur dans ',trim(nomsp),' : ',trim(message), &
    & ' (',var1,var2,var3,')'
    stop 1
  else
    call jmsetcode(code)
  end if

END Subroutine jmerreur3

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmerreur4(nomsp,code,var1,var2,var3,var4)

  ! Arguments
  character(len=*), intent(in) :: nomsp
  integer , intent(in) :: code
  integer , intent(in) :: var1, var2, var3,var4

  ! Variables locales
  integer  :: arret
  character(len=80) :: message

  call jmgetstop(arret)
  if (arret == 1) then
    call jmgetmessage(code,message)
    print *,'JMFFT Erreur dans ',trim(nomsp),' : ',trim(message), &
    & ' (',var1,var2,var3,var4,')'
    stop 1
  else
    call jmsetcode(code)
  end if

END Subroutine jmerreur4

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmfact(n,fact,nfact,ideb,ifin)

  ! Arguments
  integer , intent(in) :: n, nfact, ideb
  integer , intent(out) :: ifin
  integer , intent(inout), dimension(0:nfact-1) :: fact

  ! Variables locales
  integer  :: m
  integer  :: n2, p2, n3, p3, n5, p5
  character(len=*), parameter :: nomsp = 'JMFACT'
  ! Nombres premiers
  integer , parameter :: npremiers = 7
  integer , dimension(0:npremiers-1) :: premiers = (/7,11,13,17,19,23,29/)
  integer  :: ip, premier, pp, np

  m = n

  ! Etude des puissances de deux
  p2 = 0
  n2 = 1
  do
    if (mod(m,2) == 0) then
      p2 = p2+1
      n2 = n2*2
      m  = m/2
    else
      exit
    end if
  end do
  ifin = ideb+3
  if (ifin > nfact) &
  & call jmerreur2(nomsp,7,nfact,ifin)
  fact(ideb+1) = n2
  fact(ideb+2) = p2

  ! Etude des puissances de trois
  p3 = 0
  n3 = 1
  do
    if (mod(m,3) == 0) then
      p3 = p3+1
      n3 = n3*3
      m  = m/3
    else
      exit
    end if
  end do
  ifin = ifin+2
  if (ifin > nfact) &
  & call jmerreur2(nomsp,7,nfact,ifin)
  fact(ideb+3) = n3
  fact(ideb+4) = p3

  ! Etude des puissances de cinq
  p5 = 0
  n5 = 1
  do
    if (mod(m,5) == 0) then
      p5 = p5+1
      n5 = n5*5
      m  = m/5
    else
      exit
    end if
  end do
  ifin = ifin+2
  if (ifin > nfact) &
  & call jmerreur2(nomsp,7,nfact,ifin)
  fact(ideb+5) = n5
  fact(ideb+6) = p5

  ! On met a jour le nombre de termes
  fact(ideb) = 7

  ! Si on a fini
  if (n2*n3*n5 == n) return

  ! Il reste maintenant des facteurs premiers bizarres
  ! On va boucler tant qu'on n'a pas fini ou tant qu'on n'a pas epuise la liste

  do ip = 0,npremiers-1

    premier = premiers(ip)

    pp = 0
    np = 1
    do
      if (mod(m,premier) == 0) then
        pp = pp+1
        np = np*premier
        m  = m/premier
      else
        exit
      end if
    end do
    ifin = ifin+2
    if (ifin > nfact) &
    & call jmerreur2(nomsp,7,nfact,ifin)
    fact(ifin-2) = pp
    fact(ifin-1) = premier
    fact(ideb) = fact(ideb) + 2

    ! Si le nombre est completement factorise, inutile de continuer
    if (m == 1) exit

  end do

  ! On regarde si la factorisation est terminee
  if (m == 1) then
    return
  else
    call jmerreur1(nomsp,3,n)
  end if
  
END Subroutine jmfact

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine fftfax(n,ifax,trigs)

  ! Constantes pour les arguments
  integer , parameter :: nfax = 19

  ! Arguments
  integer , intent(in) :: n
  REALTYPE, dimension(0:2*n-1), intent(out) :: trigs
  integer , dimension(0:nfax-1), intent(out) :: ifax

  ! Variables locales
  integer  :: ntrigs 
  integer  :: ifin
  character(len=*), parameter :: nomsp = 'FFTFAX'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ntrigs = 2*n

  ! Factorisation de n dans ifax
  call jmfact(n,ifax,nfax,0,ifin)

  ! Preparation des tables
  call jmtable(trigs,ntrigs,0,n)

END Subroutine fftfax

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmgetcode(code)

  ! Arguments
  integer , intent(out) :: code

  ! Variables locales

  call jmgetsetcode(code,'g')

END Subroutine jmgetcode

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmgeterreur(arret)

  ! Arguments
  logical , intent(out) :: arret

  ! Variables locales

  call jmgetseterreur(arret,'g')

END Subroutine jmgeterreur
 
! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmgetmessage(code,message)

  ! Arguments
  integer , intent(in) :: code
  character(len=*), intent(out) :: message

  ! Variables locales
  integer , parameter :: mm = 26
  character(len=34), dimension(0:mm-1) :: messages = (/ &
  & "Pas d'erreur                     ",                &
  & "Isign doit etre egal a -1 ou 1   ",                &
  & "Isign doit etre egal a 0, -1 ou 1",                &
  & "Nombres premiers trop grands     ",                &
  & "Nwork negatif ou nul             ",                &
  & "Nwork trop petit                 ",                &
  & "Tronconnage impossible           ",                &
  & "Trop de facteurs premiers        ",                &
  & "l doit etre >= 1                 ",                &
  & "ldx doit etre >= n               ",                &
  & "ldx doit etre >= n/2+1           ",                &
  & "ldx1 doit etre >= n              ",                &
  & "ldx1 doit etre >= n/2+1          ",                &
  & "ldx2 doit etre >= m              ",                &
  & "ldy doit etre >= n               ",                &
  & "ldy doit etre >= n+2             ",                &
  & "ldy doit etre >= n/2+1           ",                &
  & "ldy1 doit etre >= n              ",                &
  & "ldy1 doit etre >= n+2            ",                &
  & "ldy1 doit etre >= n/2+1          ",                &
  & "ldy2 doit etre >= m              ",                &
  & "m doit etre >= 1                 ",                &
  & "m ou n doit etre pair            ",                &
  & "n doit etre >= 1                 ",                &
  & "n doit etre pair                 ",                &
  & "n ou m ou l doit etre pair       "                 &
  & /)

  if (code < 0 .or. code >= mm) then
    print *,'JMFFT GETMESSAGE Code invalide : ',code
  end if

  message = messages(code)

END Subroutine jmgetmessage

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmgetnwork(nwork,nwork_def,nwork_min)

  ! On recupere la valeur de nwork si elle a ete augmentee par l'utilisateur
  ! Sinon on prend la valeur par defaut
  ! Il s'agit du nwork des routines 2d et 3d

  ! Arguments
  integer , intent(out) :: nwork
  integer , intent(in)  :: nwork_def, nwork_min

  ! Variables locales
  integer  :: nwork_loc
  character(len=*), parameter :: nomsp = 'JMGETNWORK'

  call jmgetsetnwork(nwork_loc,'g')

  ! Valeur par defaut
  if (nwork_loc == -1) then
    nwork = nwork_def
  ! Valeur invalide (trop petite)
  else if (nwork_loc < nwork_min) then
    call jmerreur2(nomsp,5,nwork_loc,nwork_min)
  ! Valeur correcte
  else
    nwork = nwork_loc
  end if

END Subroutine jmgetnwork

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmgetsetcode(code,type)

  ! Subroutine qui permet de stocker le dernier code de retour obtenu
  ! Ceci evite de recourir a un common ...

  ! Arguments
  integer , intent(inout) :: code
  character(len=1), intent(in) :: type

  ! Variables locales

  ! Variable statique
  integer , save :: code_last = 0

  if (type == 's') then
    code_last = code
  else if (type == 'g') then 
    code = code_last
  end if

END Subroutine jmgetsetcode

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmgetseterreur(arret,type)

  ! Subroutine qui permet de stocker une valeur statique
  ! Ceci evite de recourir a un common ...

  ! Arguments
  logical , intent(inout) :: arret
  character(len=1), intent(in) :: type

  ! Variables locales

  ! Variable statique
  logical , save :: arret_last = .true.

  if (type == 's') then
    arret_last = arret
  else if (type == 'g') then 
    arret = arret_last
  end if

END Subroutine jmgetseterreur

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmgetsetnwork(nwork,type)

  ! Subroutine qui permet de stocker une valeur statique
  ! Ceci evite de recourir a un common ...

  ! Arguments
  integer , intent(inout) :: nwork
  character(len=1), intent(in) :: type

  ! Variables locales

  ! Variable statique
  integer , save :: nwork_last = -1

  if (type == 's') then
    nwork_last = nwork
  else if (type == 'g') then 
    nwork = nwork_last
  end if

END Subroutine jmgetsetnwork

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmgetsetstop(arret,type)

  ! Subroutine qui permet de stocker une valeur statique
  ! Ceci evite de recourir a un common ...

  ! Arguments
  integer , intent(inout) :: arret
  character(len=1), intent(in) :: type

  ! Variables locales

  ! Variable statique
  integer , save :: arret_last = 1

  if (type == 's') then
    arret_last = arret
  else if (type == 'g') then 
    arret = arret_last
  end if

END Subroutine jmgetsetstop

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmgetstop(arret)

  ! Arguments
  integer , intent(out) :: arret

  ! Variables locales

  call jmgetsetstop(arret,'g')

END Subroutine jmgetstop

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine rfftmlt(x,work,trigs,ifax,incx,jumpx,n,m,isign)

  ! Constantes pour les arguments
  integer , parameter :: nfax = 19

  ! Arguments
  integer , intent(in) :: incx, jumpx, n, m, isign
  REALTYPE, intent(inout), dimension(0:m*(n+2)-1) :: x
  REALTYPE, intent(out), dimension(0:2*n*m-1) :: work
  REALTYPE, intent(in), dimension(0:2*n-1) :: trigs
  integer , intent(inout), dimension(0:nfax-1) :: ifax

  ! Variables locales
  integer  :: ntrigs, nwork
  REALTYPE :: scale
  integer  :: dimx, debx
  integer  :: signe
  REALTYPE :: scale_temp
  character(len=*), parameter :: nomsp = 'RFFTMLT'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,1,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (mod(n,2) /= 0 .and. mod(m,2) /= 0) &
  & call jmerreur2(nomsp,22,n,m)

  ! Gestion de trigs
  ntrigs = 2*n

  ! Gestion de work
  nwork = 2*n*m

  if (isign == -1) then

    ! On appelle le sous-programme principal
    scale = ONE/n
    dimx = m*(n+2)
    debx = 0
    signe = -1
    call jmscm1dxy(m,n,ifax,nfax,0,trigs,ntrigs,0,work,nwork, &
    & x,dimx,debx,incx,jumpx,x,dimx,debx,incx,jumpx,signe,scale)

  else

    ! On appelle le sous-programme principal
    dimx = m*(n+2)
    debx = 0
    signe = 1
    scale_temp = ONE
    call jmcsm1dxy(m,n,ifax,nfax,0,trigs,ntrigs,0,work,nwork, &
    & x,dimx,debx,incx,jumpx,x,dimx,debx,incx,jumpx,signe,scale_temp)

  end if

END Subroutine rfftmlt

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine scfft(isign,n,scale,x,y,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: n
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:n-1) :: x
  REALTYPE, intent(out), dimension(0:2*(n/2)+1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*n-1) :: table
  REALTYPE, intent(inout), dimension(0:2*n-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i
  integer  :: ntable, nwork
  integer  :: nfact
  integer , dimension(0:99) :: fact
  integer  :: dimx, debx, incx, jumpx
  integer  :: dimy, deby, incy, jumpy
  character(len=*), parameter :: nomsp = 'SCFFT'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (mod(n,2) /= 0) call jmerreur1(nomsp,24,n)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 2*n

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,n)
    return
  else
    nfact = nint(table(0))
    fact(0:nfact-1) = nint(table(0:nfact-1))
  end if

  ! On appelle le sous-programme principal
  dimx = n         ; debx = 0 ; incx = 1 ; jumpx = 0
  dimy = 2*(n/2)+2 ; deby = 0 ; incy = 1 ; jumpy = 0
  call jmscm1dxy(1,n,fact,100,0,table,ntable,100+0,work,nwork, &
  & x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

END Subroutine scfft

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine scfft2d(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: n, m, ldx, ldy
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:ldx*m-1) :: x
  REALTYPE, intent(out), dimension(0:2*ldy*m-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*(n+m)-1) :: table
  REALTYPE, intent(inout), dimension(0:512*max(n,m)-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i, j
  integer  :: ntable, nwork, ioff
  integer  :: nfact, mfact
  integer , dimension(0:99) :: fact
  integer  :: ideb, ifin, jdeb, jfin, n_temp, m_temp, nwork_temp
  logical  :: debut, fin
  integer  :: dimx, debx, incx, jumpx
  integer  :: dimy, deby, incy, jumpy
  integer  :: signe
  logical  :: npair, mpair
  character(len=*), parameter :: nomsp = 'SCFFT2D'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Gestion de npair et mpair
  npair = ( mod(n,2) == 0 )
  mpair = ( mod(m,2) == 0 )

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n    ) call jmerreur2(nomsp, 9,ldx,n    )
  if (ldy < n/2+1) call jmerreur2(nomsp,16,ldy,n/2+1)
  if ( .not.npair .and. .not.mpair) &
  & call jmerreur2(nomsp,22,n,m)

  ! Gestion de table
  ntable = 100+2*(n+m)

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    call jmfact(m,fact,100,nfact,mfact)
    table(0:mfact-1) = fact(0:mfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0  ,n)
    call jmtable(table,ntable,100+2*n,m)
    return
  else
    nfact = nint(table(0))
    mfact = nint(table(nfact)) + nfact
    fact(0:mfact-1) = nint(table(0:mfact-1))
  end if

  ! Gestion de work
  !nwork = 2*2*(n/2+1)*m
  !nwork = 512*max(n,m)
  call jmgetnwork(nwork,512*max(n,m),4*max(n,m))

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  debut = .true.
  do

    ! Tronconnage
    call jmdecoup(m,2*n,nwork,debut,npair,m_temp,jdeb,jfin,nwork_temp,fin)

    ! On fait les T.F. reelles sur la premiere dimension
    ! Tout se passe comme si on faisait une T.F. 1D multiple (d'ordre m)
    dimx = ldx*m   ; debx = jdeb*ldx   ; incx = 1 ; jumpx = ldx
    dimy = 2*ldy*m ; deby = jdeb*2*ldy ; incy = 1 ; jumpy = 2*ldy
    signe = 1
    call jmscm1dxy(m_temp,n,fact,100,0,table,ntable,100+0, &
    & work,nwork_temp, &
    & x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,signe,scale)

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

  ! On fait les T.F. sur l'autre dimension
  debut = .true.
  do

    ! Tronconnage
    call jmdecoup(n/2+1,4*m,nwork,debut,mpair,n_temp,ideb,ifin,nwork_temp,fin)

    ! On copie
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
        work(             n_temp*j+i-ideb) = y(2*i  +2*ldy*j) 
        work(nwork_temp/4+n_temp*j+i-ideb) = y(2*i+1+2*ldy*j) 
      end do
    end do
    ioff = 0

    ! On fait les T.F. sur l'autre dimension (divisee par deux bien sur)
    ! On va chercher l'autre table des cosinus
    call jmccm1d(n_temp,m,fact,100,nfact,table,ntable,100+2*n,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do i = ideb,ifin
        y(2*i  +2*ldy*j) =       work(ioff             +n_temp*j+i-ideb)
        y(2*i+1+2*ldy*j) = isign*work(ioff+nwork_temp/4+n_temp*j+i-ideb)
      end do
    end do

    ! A-t-on fini ?
    if (fin) then
      exit
    else
      debut = .false.
      cycle
    end if

  end do

END Subroutine scfft2d

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine scfft3d(isign,n,m,l,scale,x,ldx1,ldx2,y,ldy1,ldy2,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: n, m, l, ldx1, ldx2, ldy1, ldy2
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:ldx1*ldx2*l-1) :: x
  REALTYPE, intent(out), dimension(0:2*ldy1*ldy2*l-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*(n+m+l)-1) :: table
  REALTYPE, intent(inout), dimension(0:512*max(n,m,l)-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i, j, k
  integer  :: ntable, nwork, ioff
  integer  :: nfact, mfact, lfact
  integer , dimension(0:99) :: fact
  integer  :: ideb, ifin, jdeb, jfin, kdeb, kfin, i1, i2, j1, j2
  integer  :: nltemp, nmtemp, mltemp, nwork_temp, iwork
  logical  :: debut, fini
  logical  :: npair, mpair, lpair
  character(len=*), parameter :: nomsp = 'SCFFT3D'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Gestion de npair
  npair = (mod(n,2) == 0)
  mpair = (mod(m,2) == 0)
  lpair = (mod(l,2) == 0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (l < 1) call jmerreur1(nomsp,8,l)
  if (ldx1 < n    ) call jmerreur2(nomsp,11,ldx1,n    )
  if (ldy1 < n/2+1) call jmerreur2(nomsp,19,ldy1,n/2+1)
  if (ldx2 < m    ) call jmerreur2(nomsp,13,ldx2,m    )
  if (ldy2 < m    ) call jmerreur2(nomsp,20,ldy2,m    )
  if (.not.mpair .and. .not.npair .and. .not.lpair) &
  & call jmerreur3(nomsp,25,n,m,l)

  ! Gestion de table
  ntable = 100+2*(n+m+l)

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    call jmfact(m,fact,100,nfact,mfact)
    call jmfact(l,fact,100,mfact,lfact)
    table(0:lfact-1) = fact(0:lfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0      ,n)
    call jmtable(table,ntable,100+2*n    ,m)
    call jmtable(table,ntable,100+2*(n+m),l)
    return
  else
    nfact = nint(table(0))
    mfact = nint(table(nfact)) + nfact
    lfact = nint(table(mfact)) + mfact
    fact(0:lfact-1) = nint(table(0:lfact-1))
  end if

  ! Gestion de work
  !nwork = 2*2*(n/2+1)*m*l
  !nwork = 512*max(n,m,l)
  call jmgetnwork(nwork,512*max(n,m,l),4*max(n,m,l))

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  ! et la troisieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    call jmdecoup3(m,l,4*(n/2+1),nwork,debut,npair,jdeb,jfin,kdeb,kfin,mltemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    do i = 0,n-1
      iwork = 0
      do k = kdeb,kfin
        j1 = 0
        j2 = m-1
        if (k == kdeb) j1 = jdeb
        if (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = j1,j2
          work(iwork+i*mltemp) = scale*x(i+ldx1*j+ldx1*ldx2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la premiere dimension
    ioff = 0
    call jmscm1d(mltemp,n,fact,100,0    ,table,ntable,100+0      ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do i = 0,n/2
      iwork = 0
      do k = kdeb,kfin
        j1 = 0
        j2 = m-1
        if (k == kdeb) j1 = jdeb
        if (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = j1,j2
           y(2*i  +2*ldy1*j+2*ldy1*ldy2*k) = work(ioff+             iwork+i*mltemp)
           y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k) = work(ioff+nwork_temp/4+iwork+i*mltemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

  ! On fait les T.F. sur la troisieme dimension en tronconnant sur la premiere
  ! et la deuxieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    call jmdecoup3(n/2+1,m,4*l,nwork,debut,.true.,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! On prend garde a la gestion des extremites
    do k = 0,l-1
      iwork = 0
      do j = jdeb,jfin
        i1 = 0
        i2 = n/2
        if (j == jdeb) i1 = ideb
        if (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          work(             iwork+k*nmtemp) = y(2*i  +2*ldy1*j+2*ldy1*ldy2*k)
          work(nwork_temp/4+iwork+k*nmtemp) = y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la troisieme dimension
    ioff = 0
    call jmccm1d(nmtemp,l,fact,100,mfact,table,ntable,100+2*(n+m),work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do k = 0,l-1
      iwork = 0
      do j = jdeb,jfin
        i1 = 0
        i2 = n/2
        if (j == jdeb) i1 = ideb
        if (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          y(2*i  +2*ldy1*j+2*ldy1*ldy2*k) = work(ioff+             iwork+k*nmtemp)
          y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k) = work(ioff+nwork_temp/4+iwork+k*nmtemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

  ! On fait les T.F. sur la deuxieme dimension en tronconnant sur la premiere
  ! et la troisieme
  debut = .true.
  fini  = .false.
  do while (.not.fini)

    ! Tronconnage
    call jmdecoup3(n/2+1,l,4*m,nwork,debut,.true.,ideb,ifin,kdeb,kfin,nltemp,nwork_temp,fini)
    debut = .false.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    do j = 0,m-1
      iwork = 0
      do k = kdeb,kfin
        i1 = 0
        i2 = n/2
        if (k == kdeb) i1 = ideb
        if (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          work(             iwork+j*nltemp) = y(2*i  +2*ldy1*j+2*ldy1*ldy2*k)
          work(nwork_temp/4+iwork+j*nltemp) = y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k)
          iwork = iwork+1
        end do
      end do
    end do

    ! On fait les T.F. sur la deuxieme dimension
    ioff = 0
    call jmccm1d(nltemp,m,fact,100,nfact,table,ntable,100+2*n    ,work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    do j = 0,m-1
      iwork = 0
      do k = kdeb,kfin
        i1 = 0
        i2 = n/2
        if (k == kdeb) i1 = ideb
        if (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = i1,i2
          y(2*i  +2*ldy1*j+2*ldy1*ldy2*k) = &
          &       work(ioff             +iwork+j*nltemp)
          y(2*i+1+2*ldy1*j+2*ldy1*ldy2*k) = &
          & isign*work(ioff+nwork_temp/4+iwork+j*nltemp)
          iwork = iwork+1
        end do
      end do
    end do

  end do

END Subroutine scfft3d

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine scfftm(isign,n,m,scale,x,ldx,y,ldy,table,work,isys)

  ! Arguments
  integer , intent(in) :: isign
  integer , intent(in) :: m, n, ldx, ldy
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in), dimension(0:ldx*m-1) :: x
  REALTYPE, intent(out), dimension(0:2*ldy*m-1) :: y
  REALTYPE, intent(inout), dimension(0:100+2*n-1) :: table
  REALTYPE, intent(inout), dimension(0:2*n*m-1) :: work
  integer , intent(in) :: isys

  ! Variables locales
  integer  :: i
  integer  :: ntable, nwork
  integer  :: nfact
  integer , dimension(0:99) :: fact
  integer  :: dimx, debx, incx, jumpx
  integer  :: dimy, deby, incy, jumpy
  character(len=*), parameter :: nomsp = 'SCFFTM'

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (m < 1) call jmerreur1(nomsp,21,m)
  if (ldx < n    ) call jmerreur2(nomsp,9,    ldx,n)
  if (ldy < n/2+1) call jmerreur2(nomsp,16,ldy,n/2+1)
  if (mod(n,2) /= 0 .and. mod(m,2) /= 0) &
  & call jmerreur2(nomsp,22,n,m)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 2*n*m

  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,n)
    return
  else
    nfact = nint(table(0))
    fact(0:nfact-1) = nint(table(0:nfact-1))
  end if

  ! On appelle le sous-programme principal
  dimx = ldx*m   ; debx = 0 ; incx = 1 ; jumpx = ldx
  dimy = 2*ldy*m ; deby = 0 ; incy = 1 ; jumpy = 2*ldy
  call jmscm1dxy(m,n,fact,100,0,table,ntable,100+0,work,nwork, &
  & x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

END Subroutine scfftm

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmscm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

  ! Arguments
  integer , intent(in) :: m, n
  integer , intent(in) :: nfact, ifact
  integer , intent(inout), dimension(0:nfact-1) :: fact
  integer , intent(in) :: ntable,itable
  REALTYPE, intent(in), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: ioff1, ioff2
  integer  :: i, j
  REALTYPE :: t, u, v, w
  REALTYPE :: c, s
  integer  :: is, it

  ! Gestion de work
  ioff1 = ioff
  ioff2 = nwork/2 - ioff1

  ! On doit faire m T.F. reelles de longueur n
  ! Si m est pair
  if (mod(m,2) == 0) then

    ! On distribue
    if (m/2 >= 16 .or. n < 8) then

      do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          work(ioff2        +i*m/2+j) = work(ioff1+i*m+j    )
          work(ioff2+nwork/4+i*m/2+j) = work(ioff1+i*m+j+m/2)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n-1
          work(ioff2        +i*m/2+j) = work(ioff1+i*m+j    )
          work(ioff2+nwork/4+i*m/2+j) = work(ioff1+i*m+j+m/2)
        end do
      end do

    end if
        
    ! On fait m/2 t.f. complexes de longueur n
    call jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    ioff1 = nwork/2 - ioff2

    ! On regenere le resultat
    if (m/2 >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          it = n-i
          if (i == 0) it = 0
          t = work(ioff2        + i*m/2+j)
          u = work(ioff2+nwork/4+ i*m/2+j)
          v = work(ioff2        +it*m/2+j)
          w = work(ioff2+nwork/4+it*m/2+j)
          work(ioff1        +i*m+j    ) = (t+v)/2
          work(ioff1+nwork/4+i*m+j    ) = (u-w)/2
          work(ioff1        +i*m+j+m/2) = (u+w)/2
          work(ioff1+nwork/4+i*m+j+m/2) = (v-t)/2
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2
          it = n-i
          if (i == 0) it = 0
          t = work(ioff2        + i*m/2+j)
          u = work(ioff2+nwork/4+ i*m/2+j)
          v = work(ioff2        +it*m/2+j)
          w = work(ioff2+nwork/4+it*m/2+j)
          work(ioff1        +i*m+j    ) = (t+v)/2
          work(ioff1+nwork/4+i*m+j    ) = (u-w)/2
          work(ioff1        +i*m+j+m/2) = (u+w)/2
          work(ioff1+nwork/4+i*m+j+m/2) = (v-t)/2
        end do
      end do

    end if

  ! Si m n'est pas pair mais que n l'est
  else if (mod(n,2) == 0) then

    ! On distribue les indices pairs et impairs selon n

    if (m >= 16 .or. n/2 < 8) then

      do i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0, m-1
          work(ioff2+        m*i+j) = work(ioff1+m*(2*i  )+j)
          work(ioff2+nwork/4+m*i+j) = work(ioff1+m*(2*i+1)+j)
        end do
      end do

    else

      do j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0, n/2-1
          work(ioff2+        m*i+j) = work(ioff1+m*(2*i  )+j)
          work(ioff2+nwork/4+m*i+j) = work(ioff1+m*(2*i+1)+j)
        end do
      end do

    end if

    ! On fait m t.f. complexes de taille n/2
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    call jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1
    ioff1 = nwork/2 - ioff2

    ! Maintenant, il faut reconstituer la t.f. reelle

    if (m >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m-1
          is = i
          it = n/2-i
          if (i == 0 .or. i == n/2) then
            is = 0
            it = 0
          end if
          t = work(ioff2        +is*m+j)
          u = work(ioff2+nwork/4+is*m+j)
          v = work(ioff2        +it*m+j)
          w = work(ioff2+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          work(ioff1        +i*m+j) = (t+v)/2 + c*(u+w)/2 - s*(v-t)/2
          work(ioff1+nwork/4+i*m+j) = (u-w)/2 + c*(v-t)/2 + s*(u+w)/2
        end do
      end do

    else

      do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2
          is = i
          it = n/2-i
          if (i == 0 .or. i == n/2) then
            is = 0
            it = 0
          end if
          t = work(ioff2        +is*m+j)
          u = work(ioff2+nwork/4+is*m+j)
          v = work(ioff2        +it*m+j)
          w = work(ioff2+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          work(ioff1        +i*m+j) = (t+v)/2 + c*(u+w)/2 - s*(v-t)/2
          work(ioff1+nwork/4+i*m+j) = (u-w)/2 + c*(v-t)/2 + s*(u+w)/2
        end do
      end do

    end if

  end if

  ioff = ioff1

END Subroutine jmscm1d

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

! Variante de jmscm1d ou on fournit x en entree et en sortie
Subroutine jmscm1dxy(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

  ! Arguments
  integer , intent(in) :: m, n
  integer , intent(in) :: nfact, ifact
  integer , intent(inout), dimension(0:nfact-1) :: fact
  integer , intent(in) :: ntable,itable
  REALTYPE, intent(in), dimension(0:ntable-1) :: table
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(in) :: dimx, debx, incx, jumpx
  integer , intent(in) :: dimy, deby, incy, jumpy
  REALTYPE, intent(in) :: scale
  REALTYPE, intent(in),  dimension(0:dimx-1) :: x
  REALTYPE, intent(out), dimension(0:dimy-1) :: y
  integer , intent(in) :: isign

  ! Variables locales
  integer  :: i, j
  REALTYPE :: t, u, v, w
  REALTYPE :: c, s
  integer  :: is, it
  integer  :: ioff

  ! On doit faire m T.F. reelles de longueur n
  ! Si m est pair
  if (mod(m,2) == 0) then

    ! On distribue
    if (m/2 >= 16 .or. n < 8) then

      do i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          work(        i*m/2+j) = x(debx+i*incx+(j)    *jumpx)
          work(nwork/4+i*m/2+j) = x(debx+i*incx+(j+m/2)*jumpx)
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n-1
          work(        i*m/2+j) = x(debx+i*incx+(j)    *jumpx)
          work(nwork/4+i*m/2+j) = x(debx+i*incx+(j+m/2)*jumpx)
        end do
      end do

    end if

    ! On fait m/2 t.f. complexes de longueur n
    ioff = 0
    call jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

    ! On regenere le resultat
    if (m/2 >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m/2-1
          it = n-i
          if (i == 0) it = 0
          t = work(ioff        + i*m/2+j)
          u = work(ioff+nwork/4+ i*m/2+j)
          v = work(ioff        +it*m/2+j)
          w = work(ioff+nwork/4+it*m/2+j)
          y(deby+(2*i)  *incy+(j)    *jumpy) =       scale*(t+v)/2
          y(deby+(2*i+1)*incy+(j)    *jumpy) = isign*scale*(u-w)/2
          y(deby+(2*i)  *incy+(j+m/2)*jumpy) =       scale*(u+w)/2
          y(deby+(2*i+1)*incy+(j+m/2)*jumpy) = isign*scale*(v-t)/2
        end do
      end do

    else

      do j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2
          it = n-i
          if (i == 0) it = 0
          t = work(ioff        + i*m/2+j)
          u = work(ioff+nwork/4+ i*m/2+j)
          v = work(ioff        +it*m/2+j)
          w = work(ioff+nwork/4+it*m/2+j)
          y(deby+(2*i)  *incy+(j)    *jumpy) =       scale*(t+v)/2
          y(deby+(2*i+1)*incy+(j)    *jumpy) = isign*scale*(u-w)/2
          y(deby+(2*i)  *incy+(j+m/2)*jumpy) =       scale*(u+w)/2
          y(deby+(2*i+1)*incy+(j+m/2)*jumpy) = isign*scale*(v-t)/2
        end do
      end do

    end if

  ! Si m n'est pas pair mais que n l'est
  else if (mod(n,2) == 0) then

    ! On distribue les indices pairs et impairs selon n

    if (m >= 16 .or. n/2 < 8) then

      do i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0, m-1
          work(        m*i+j) = x(debx+incx*(2*i  )+jumpx*j)
          work(nwork/4+m*i+j) = x(debx+incx*(2*i+1)+jumpx*j)
        end do
      end do

    else

      do j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0, n/2-1
          work(        m*i+j) = x(debx+incx*(2*i  )+jumpx*j)
          work(nwork/4+m*i+j) = x(debx+incx*(2*i+1)+jumpx*j)
        end do
      end do

    end if

    ! On fait m t.f. complexes de taille n/2
    ioff = 0
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    call jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1

    ! Maintenant, il faut reconstituer la t.f. reelle

    if (m >= 16 .or. n/2 < 8) then

      do i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        do j = 0,m-1
          is = i
          it = n/2-i
          if (i == 0 .or. i == n/2) then
            is = 0
            it = 0
          end if
          t = work(ioff        +is*m+j)
          u = work(ioff+nwork/4+is*m+j)
          v = work(ioff        +it*m+j)
          w = work(ioff+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          y(deby+(2*i  )*incy+j*jumpy) = &
          &       scale*((t+v)/2 + c*(u+w)/2 - s*(v-t)/2)
          y(deby+(2*i+1)*incy+j*jumpy) = &
          & isign*scale*((u-w)/2 + c*(v-t)/2 + s*(u+w)/2)
        end do
      end do

    else

      do j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        do i = 0,n/2
          is = i
          it = n/2-i
          if (i == 0 .or. i == n/2) then
            is = 0
            it = 0
          end if
          t = work(ioff        +is*m+j)
          u = work(ioff+nwork/4+is*m+j)
          v = work(ioff        +it*m+j)
          w = work(ioff+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          y(deby+(2*i  )*incy+j*jumpy) = &
          &       scale*((t+v)/2 + c*(u+w)/2 - s*(v-t)/2)
          y(deby+(2*i+1)*incy+j*jumpy) = &
          & isign*scale*((u-w)/2 + c*(v-t)/2 + s*(u+w)/2)
        end do
      end do

    end if

  end if

END Subroutine jmscm1dxy

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmsetcode(code)

  ! Arguments
  integer , intent(in) :: code

  ! Variables locales
  integer  :: errcode

  errcode = code
  call jmgetsetcode(errcode,'s')

END Subroutine jmsetcode

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmseterreur(arret)

  ! Arguments
  logical , intent(in) :: arret

  ! Variables locales
  logical  :: arret2

  arret2 = arret
  call jmgetseterreur(arret2,'s')

END Subroutine jmseterreur

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmsetnwork(nwork)

  ! Subroutine appelee par l'utilisateur pour augmenter le nwork
  ! des routines 2d et 3d

  ! Arguments
  integer , intent(in) :: nwork

  ! Variables locales
  character(len=*), parameter :: nomsp = 'JMSETNWORK'
  integer  :: nwork2

  if (nwork <= 0) then
    call jmerreur1(nomsp,4,nwork)
  end if

  nwork2 = nwork
  call jmgetsetnwork(nwork2,'s')

END Subroutine jmsetnwork

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmsetstop(arret)

  ! Arguments
  integer , intent(in) :: arret

  ! Variables locales
  integer  :: arret2

  arret2 = arret
  call jmgetsetstop(arret2,'s')

END Subroutine jmsetstop

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmtable(table,ntable,itable,n)

  ! Arguments
  integer , intent(in) :: ntable, itable, n
  REALTYPE, intent(out), dimension(0:ntable-1) :: table

  ! Variables locales
  REALTYPE, save :: twopi
  REALTYPE :: temp, temp1, temp2

  integer  :: i

  twopi = TWO*PI

  ! Calcul des sinus et cosinus

  ! Si n est multiple de 4, astuces en serie
  if (mod(n,4) == 0) then
    ! On se debarrasse des cas limite
    table(itable+      0) =  1
    table(itable+n+    0) =  0
    table(itable+    n/4) =  0
    table(itable+n+  n/4) =  1
    table(itable+    n/2) = -1
    table(itable+n+  n/2) =  0
    table(itable+  3*n/4) =  0
    table(itable+n+3*n/4) = -1
    ! Cas general
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = 1,n/4-1
      temp = COS(twopi*i/n)
      table(itable+    i)     =  temp
      table(itable+    n/2-i) = -temp
      table(itable+    n/2+i) = -temp
      table(itable+    n-i)   =  temp
      table(itable+n+  n/4+i) =  temp
      table(itable+n+  n/4-i) =  temp
      table(itable+n+3*n/4+i) = -temp
      table(itable+n+3*n/4-i) = -temp
    end do

  ! Si n est simplement multiple de 2 (sans etre multiple de 4)
  else if (mod(n,2) == 0) then
    ! On se debarrasse des cas limite
    table(itable+    0) =  1
    table(itable+n+  0) =  0
    table(itable+  n/2) = -1
    table(itable+n+n/2) =  0
    ! Cas general
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = 1,n/2-1
      temp1 = COS(twopi*i/n)
      table(itable+      i) =  temp1
      table(itable+  n/2+i) = -temp1
      temp2 = SIN(twopi*i/n)
      table(itable+n+    i) =  temp2
      table(itable+n+n/2+i) = -temp2
    end do

  ! Si n est impair
  else
    ! On se debarrasse des cas limite
    table(itable+  0) =  1
    table(itable+n+0) =  0
!dir$ ivdep
!ocl novrec
!cdir nodep
    do i = 1,n/2
      temp1 = COS(twopi*i/n)
      table(itable+    i) =  temp1
      table(itable+  n-i) =  temp1
      temp2 = SIN(twopi*i/n)
      table(itable+n+  i) =  temp2
      table(itable+n+n-i) = -temp2
    end do

  end if

END Subroutine jmtable

! $Header: /home/moura/PseudoPack_2006/Source/RCS/Module_FFT_CRAY_JM.F,v 1.1 2006/05/18 03:34:28 moura Exp $

Subroutine jmtransp(n,m,l,work,nwork,ioff)

  ! Arguments
  integer , intent(in) :: n, m, l
  integer , intent(in) :: nwork
  REALTYPE, intent(inout), dimension(0:nwork-1) :: work
  integer , intent(inout) :: ioff

  ! Variables locales
  integer  :: ioff1, ioff2
  integer  :: ij, k

  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! On transpose (nm)(l) en (l)(nm) en distinguant les parties reelles et im.
  if (m*n >= 16 .or. l < 8) then

    do k = 0,l-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do ij = 0,m*n-1
        work(ioff2+      ij*l+k) = work(ioff1+      k*n*m+ij)
        work(ioff2+n*m*l+ij*l+k) = work(ioff1+n*m*l+k*n*m+ij)
      end do
    end do

  else

    do ij = 0,m*n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      do k = 0,l-1
        work(ioff2+      ij*l+k) = work(ioff1+      k*n*m+ij)
        work(ioff2+n*m*l+ij*l+k) = work(ioff1+n*m*l+k*n*m+ij)
      end do
    end do

  end if

  ioff = ioff2

END Subroutine jmtransp

END MODULE CRAY_FFT_JM

#else
MODULE CRAY_FFT_JM

implicit NONE

END MODULE CRAY_FFT_JM

#endif
