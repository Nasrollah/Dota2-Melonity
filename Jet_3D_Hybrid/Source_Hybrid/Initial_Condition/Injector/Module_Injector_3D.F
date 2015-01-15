#include "pseudopack.h"


MODULE Injector_3D

  USE PseudoPack

IMPLICIT NONE

INTEGER                :: Jet_Direction = 3

INTEGER                :: Jet_Shape  = 0

REALTYPE, dimension(3) :: Jet_Center = [ 0, 0, 0 ]
REALTYPE               :: Jet_Circle = 0.19450d0
REALTYPE, dimension(2) :: Jet_Square = [ 0.10d0, 0.10d0 ]
REALTYPE               :: Jet_Height = 1.0d-4

REALTYPE               :: Jet_Rho    = 20.90280d0               ! 0.10d0
REALTYPE, dimension(3) :: Jet_U      = [ ZERO, 0.5970d0, ZERO ] ! [ 0, 1, 0 ]
REALTYPE               :: Jet_P      = 5.323320d0               ! 0.250d0
REALTYPE               :: Jet_f1     = ONE
REALTYPE               :: Jet_Gamma  = 1.40d0

REALTYPE               :: Jet_Mach   = 2
REALTYPE               :: Jet_C      = 448.21869660d0          ! m/s
REALTYPE               :: Jet_T      = 500                     ! K

REALTYPE               :: Jet_mu     = 1.510d-5   ! m^2/s  Kinematic
REALTYPE               :: Jet_kappa  = 2.600d-2   ! W/K-m
REALTYPE               :: Jet_zeta   = 1.40d0

REALTYPE               :: Jet_Re     = 2000
REALTYPE               :: Jet_Pr     = 0.72
REALTYPE               :: Jet_R      = 287                     ! J/kgK

INTEGER                :: N9, M9, K9


INTEGER, parameter :: nmaxeddies = 35 
REALTYPE, dimension(1:nmaxeddies,1:3) :: xeddy
REALTYPE, dimension(1:nmaxeddies,1:3) :: epsi
INTEGER :: neddies
REALTYPE :: sigma_sem 
REALTYPE :: xlow, xhigh, ylow, yhigh, zlow, zhigh, vol_eddies 
INTEGER :: generate_eddies = 1
INTEGER :: proc_injector = 0
REALTYPE :: vol_box

INTEGER, parameter :: bogey_forcing = 0 
INTEGER, parameter :: swirl_forcing = 0

 

REALTYPE, dimension(0:2,0:1) :: A_j, St_j, phi_j, psi_j
INTEGER, dimension(0:2,0:1) :: A_i, St_i, phi_i, psi_i


PRIVATE
PUBLIC  :: Jet_Direction
PUBLIC  :: Jet_Shape, Jet_Center, Jet_Circle, Jet_Square, Jet_Height
PUBLIC  :: Jet_Rho, Jet_U, Jet_P, Jet_f1, Jet_Gamma
PUBLIC  :: Jet_mu, Jet_Kappa, Jet_zeta
PUBLIC  :: Jet_Mach, Jet_C, Jet_T, Jet_Re, Jet_Pr, Jet_R

PUBLIC  :: Injector_Input
PUBLIC  :: Injector

CONTAINS
!
! =========================================================================
!
SUBROUTINE Injector_Input (Injector_Input_File)

  CHARACTER(LEN=*)       :: Injector_Input_File
  INTEGER                :: lid, IOS, i

  if (I_Am == First_Processor) &
    write (6,*) 'Reading .... ', TRIM(Injector_Input_File)
  
  lid = 78
  OPEN (UNIT=lid, FILE=TRIM(Injector_Input_File), STATUS='OLD', IOSTAT=IOS)
  
  read (lid,'(//)')
  read (lid,*) Jet_Direction
  read (lid,*)
  read (lid,*) Jet_Shape
  read (lid,*)
  read (lid,*) Jet_Center
  read (lid,*) Jet_Circle
  read (lid,*) Jet_Square
  read (lid,*) Jet_Height
  read (lid,*)
  read (lid,*) Jet_Mach
  read (lid,*) Jet_Re
  read (lid,*)
  read (lid,*) Jet_Rho
  read (lid,*) Jet_U
  read (lid,*) Jet_P
  read (lid,*) Jet_T
  read (lid,*) Jet_f1
  read (lid,*) Jet_Gamma
  read (lid,*) Jet_R
  
  read (lid,*)
  read (lid,*) Jet_mu
  read (lid,*) Jet_kappa
  read (lid,*) Jet_zeta

  CLOSE (lid)
 

   ! KAUSHIK ADDED THE FOLLOWING LINE
      Jet_P = Jet_Rho*Jet_R*Jet_T

 
  Jet_C                = SQRT(Jet_Gamma*Jet_R*Jet_T)
  Jet_U(Jet_Direction) = Jet_Mach*Jet_C

if (Jet_Circle <  ZERO)  &
  Jet_Circle           = Jet_Re*Jet_mu/(Jet_Rho*Jet_U(Jet_Direction))/2

! Jet_Pr               = Jet_mu*Jet_Cp/Jet_Kappa

  if (I_Am == First_Processor) &
    write (6,*) '....... Done ', TRIM(Injector_Input_File)

#include "Echo_Injector.i"


END SUBROUTINE Injector_Input
!
! =========================================================================
!
SUBROUTINE Injector (N0, N5, N2, N3,                                      &
                     M0, M5, M2, M3,                                      &
                     K0, K5, K2, K3,  NV, x, y, z, Q,                     &
                     MPI_Comm_Type,time,dt3)

integer  :: Index
integer  :: N0, N5, N2, N3
integer  :: M0, M5, M2, M3
integer  :: K0, K5, K2, K3,  NV

REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
real(kind=8) :: time, dt3

integer ,                  OPTIONAL :: MPI_Comm_Type
#if defined (PARALLEL_MPI)
integer                             :: MPI_Communicator
integer                             :: MPI_Communicator_1D
integer                             :: N_Dims
logical , dimension(:), ALLOCATABLE :: Sub_Grid
integer :: m_sem, mm_sem
real(kind=8) :: t1, rr, theta, xc, yc, zc, rc


                            MPI_Communicator = MPI_Comm_World
if (PRESENT(MPI_Comm_Type)) MPI_Communicator = MPI_Comm_Type

call MPI_CARTDIM_GET (MPI_Communicator, N_Dims, MPI_Error_Status)

ALLOCATE (Sub_Grid(N_Dims))
#endif


Index = Jet_Direction

#if defined (PARALLEL_MPI)
  call PS_MPI_Processor_Local_Info (I_Am, Last_Processor, &
                                          Index, MPI_Communicator)

  Sub_Grid = .FALSE. ; Sub_Grid(Index) = .TRUE.

  call MPI_CART_SUB (MPI_Communicator, Sub_Grid, MPI_Communicator_1D, &
                                                 MPI_Error_Status)
#endif


             if(generate_eddies.eq.1) then

!---------------------
        A_j = 0.15d0
        St_j = 0.4d0
        phi_j = 0.0d0
        psi_j = 0.0d0
        A_i = 1
        St_i = 1
        phi_i = 1
        psi_i = 1
!---------------------



                 xlow = 1.0d9
                 ylow = 1.0d9
                 xhigh = -1.0d9
                 yhigh = -1.0d9


                 do K9 = K0, K2
                 do M9 = M0, M5
                 do N9 = N0, N5

                   xc = x(N9,M9,K9)-Jet_Center(1)
                   yc = y(N9,M9,K9)-Jet_Center(2)
                   rc = sqrt(xc**2.0d0 + yc**2.0d0)

                   if(rc.le.Jet_Circle) then
                    proc_injector = 1 
                    xlow = min(xlow,x(N9,M9,K9))
                    xhigh = max(xhigh,x(N9,M9,K9))
                    ylow = min(ylow,y(N9,M9,K9))
                    yhigh = max(yhigh,y(N9,M9,K9))
                   endif

                 enddo 
                 enddo 
                 enddo 


                if(proc_injector.eq.1) then
                  if(xlow.eq.1.0d9.or.ylow.eq.1.0d9) then
                   write(*,*) 'x/y low error ', xlow, ylow
                   stop
                  endif
                  if(xhigh.eq.-1.0d9.or.yhigh.eq.-1.0d9) then
                   write(*,*) 'x/y high error ', xhigh, yhigh
                   stop
                  endif
                endif



              zlow = -Jet_Circle
              zhigh = Jet_Circle



              sigma_sem = 0.4d0*Jet_Circle 

              xlow = xlow - sigma_sem
              xhigh = xhigh + sigma_sem
              ylow = ylow - sigma_sem
              yhigh = yhigh + sigma_sem
              zlow = zlow - sigma_sem
              zhigh = zhigh + sigma_sem
        

              vol_eddies = (xhigh-xlow)*(yhigh-ylow)*(zhigh-zlow)
              vol_box = (2.0d0*(Jet_Circle+sigma_sem))**3.0d0


              neddies = min(int(vol_eddies/vol_box*dble(nmaxeddies)),nmaxeddies)


              ! ALL PROCESSORS ABOVE THE INJECTOR PROCESSORS MUST NOT DO SEM
              if(minval(z(:,:,:)).gt.Jet_Circle) then
                proc_injector = 0
                neddies = 0
              endif


              if(proc_injector.eq.1) then 
               write(*,*) 'neddies ', I_Am, neddies, xlow, xhigh, ylow, yhigh, zlow, zhigh
              endif 


              do m_sem = 1, neddies
               call random_number(t1)
               xeddy(m_sem,1) = xlow + (xhigh-xlow)*t1 
               call random_number(t1)
               xeddy(m_sem,2) = ylow + (yhigh-ylow)*t1 
               call random_number(t1) 
               xeddy(m_sem,3) = zlow + (zhigh-zlow)*t1                 

               do mm_sem = 1, 3 
               call randno(t1)
               if(t1.le.0.0d0) then
                 epsi(m_sem,mm_sem) = -1.0d0
               else
                 epsi(m_sem,mm_sem) = 1.0d0 
               endif    
               enddo
              enddo 

              generate_eddies = 0

             endif 



      if (Index == 2) call Injector_y
      if (Index == 3) call Injector_z(time,dt3)

#if defined (PARALLEL_MPI)
  call MPI_COMM_FREE (MPI_Communicator_1D, MPI_Error_Status)

DEALLOCATE (Sub_Grid)

call PS_MPI_Processor (I_Am, Last_Processor)
#endif

CONTAINS

  SUBROUTINE Injector_y

  LOGICAL , dimension(N0:N5,M0:M2,K0:K5) :: Shape

#if defined (PARALLEL_MPI)
  if (I_Am == First_Processor) then
#endif

  SELECT CASE (Jet_Shape)
    CASE DEFAULT ! Circle
      Shape =   ((x(:,M0:M2,:)-Jet_Center(1))**2 +  & 
                 (z(:,M0:M2,:)-Jet_Center(3))**2 ) <= Jet_Circle**2 
     
    CASE (1)    ! Rectangle
      Shape = ABS(x(:,M0:M2,:)-Jet_Center(1)) <= Jet_Square(1) .AND.  &
              ABS(z(:,M0:M2,:)-Jet_Center(3)) <= Jet_Square(3)
 
  END SELECT
 
  WHERE (Shape)
    Q(:,M0:M2,:,1) = Jet_Rho
    Q(:,M0:M2,:,2) = Jet_Rho*Jet_U(1)
    Q(:,M0:M2,:,3) = Jet_Rho*Jet_U(2)
    Q(:,M0:M2,:,4) = Jet_Rho*Jet_U(3)
    Q(:,M0:M2,:,5) = Jet_P/(Jet_Gamma-ONE) + HALF*Jet_Rho*SUM(Jet_U**2)
    Q(:,M0:M2,:,6) = Jet_Rho*Jet_f1
  END WHERE

  do M9 = M2,M5 
    if (y(N2,M9,K2) >= y(N2,M2,K2)+Jet_Height) EXIT 

    Q(:,M9,:,:) = Q(:,M2,:,:)
  enddo
  
#if defined (PARALLEL_MPI)
  endif
#endif

  END SUBROUTINE Injector_y
!
! ======================================================================
!
  SUBROUTINE Injector_z(time,dt3)

  implicit none

  LOGICAL , dimension(N0:N5,M0:M5,K0:K2) :: Shape
  real(kind=8) :: ke
  real(kind=8) :: xc, yc, zc, rc, theta
  
  real(kind=8) :: rj, ujet, t1, t2
  integer :: mm, nn
  real(kind=8), parameter :: pi = acos(-1.0d0) 
  real(kind=8) :: time, dt3, st, rr, rr1, bb
  real(kind=8) :: uprime, dx, dy 
  real(kind=8), dimension(N0:N5,M0:M5,K0:K5,1:3) :: vel
  integer :: idir 
  real(kind=8) :: feddy
  real(kind=8) :: u_axial, u_radial, ua_ring, ur_ring
  real(kind=8) :: bogey_fact, minz



#if defined (PARALLEL_MPI)
  if (I_Am == First_Processor) then
#endif

  SELECT CASE (Jet_Shape)
    CASE DEFAULT ! Circle
      Shape =   ((x(:,:,K0:K2)-Jet_Center(1))**2 +  & 
                 (y(:,:,K0:K2)-Jet_Center(2))**2 ) <= Jet_Circle**2 
     
    CASE (1)    ! Rectangle
      Shape = ABS(x(:,:,K0:K2)-Jet_Center(1)) <= Jet_Square(1) .AND.  &
              ABS(y(:,:,K0:K2)-Jet_Center(2)) <= Jet_Square(2)
 
  END SELECT
 
!  WHERE (Shape)
!    Q(:,:,K0:K2,1) = Jet_Rho
!    Q(:,:,K0:K2,2) = Jet_Rho*Jet_U(1)
!    Q(:,:,K0:K2,3) = Jet_Rho*Jet_U(2)
!    Q(:,:,K0:K2,4) = Jet_Rho*Jet_U(3)
!    Q(:,:,K0:K2,5) = Jet_P/(Jet_Gamma-ONE) + HALF*Jet_Rho*SUM(Jet_U**2)
!    Q(:,:,K0:K2,6) = Jet_Rho*Jet_f1
!  END WHERE


    rj = Jet_Circle
    ujet = Jet_U(3)
    uprime = 0.2d0*ujet

    ! ASSUME UNIFORM GRID 
    dx = x(N0+1,M0,K0) - x(N0,M0,K0)
    dy = y(N0,M0+1,K0) - y(N0,M0,K0)


    vel = 0.0d0

!------------------------------------------------------------------------
!                           SYNTHETIC EDDY METHOD 
!  JARRIN ET AL., INTL J. HEAT AND FLUID FLOW, VOL. 27, P. 585-593, 2006
!------------------------------------------------------------------------

  ! SET THE AXIAL VELOCITY BASED ON FREUND, JFM, VOL. 438, pp. 277-305, 2001
  
!--------------------------


        do mm = 0, 2  
         do nn = 0, 1

          call random_number(rr) 
          if(rr.le.0.05d0) then
            A_i(mm,nn) = A_i(mm,nn)*(-1) 
          endif
          if(A_i(mm,nn).eq.1) then
            A_j(mm,nn) = A_j(mm,nn) + 0.0005d0
          else if(A_i(mm,nn).eq.-1) then
            A_j(mm,nn) = A_j(mm,nn) - 0.0005d0
          endif
          A_j(mm,nn) = max(min(A_j(mm,nn),0.5d0),0.01d0)

          call random_number(rr) 
          if(rr.le.0.05d0) then
            St_i(mm,nn) = St_i(mm,nn)*(-1) 
          endif
          if(St_i(mm,nn).eq.1) then
            St_j(mm,nn) = St_j(mm,nn) + 0.00085d0
          else if(St_i(mm,nn).eq.-1) then
            St_j(mm,nn) = St_j(mm,nn) - 0.00085d0
          endif
          St_j(mm,nn) = max(min(St_j(mm,nn),0.7d0),0.1d0)

          call random_number(rr) 
          if(rr.le.0.05d0) then
            phi_i(mm,nn) = phi_i(mm,nn)*(-1) 
          endif
          if(phi_i(mm,nn).eq.1) then
            phi_j(mm,nn) = phi_j(mm,nn) + 0.00085d0
          else if(phi_i(mm,nn).eq.-1) then
            phi_j(mm,nn) = phi_j(mm,nn) - 0.00085d0
          endif

          call random_number(rr) 
          if(rr.le.0.05d0) then
            psi_i(mm,nn) = psi_i(mm,nn)*(-1) 
          endif
          if(psi_i(mm,nn).eq.1) then
            psi_j(mm,nn) = psi_j(mm,nn) + 0.00085d0
          else if(psi_i(mm,nn).eq.-1) then
            psi_j(mm,nn) = psi_j(mm,nn) - 0.00085d0
          endif

         enddo
        enddo
!--------------------------

  minz = minval(z(:,:,:))


  do K9 = K0, K2 
  do M9 = M0, M5
  do N9 = N0, N5

   xc = x(N9,M9,K9)-Jet_Center(1)
   yc = y(N9,M9,K9)-Jet_Center(2)
   zc = z(N9,M9,K9) 
   rc = sqrt(xc**2.0d0 + yc**2.0d0)


           !if(rc.le.Jet_Circle.and.minz.le.Jet_Circle) then
           if(minz.le.Jet_Circle) then

          theta = atan(yc/(xc+1.0d-6)) 
 
   
        bb = 2.0d0
        do mm = 0, 2  
         do nn = 0, 1

          st = St_j(mm,nn)
          t1 = st*ujet/(2.0d0*rj)*time + phi_j(mm,nn)

          t2 = dble(mm)*theta + psi_j(mm,nn)   

          bb = bb + A_j(mm,nn)*cos(t1)*cos(t2) 
         enddo    
        enddo


    vel(N9,M9,K9,3) = ujet/2.0d0*(1.0d0 - tanh(bb*(rc/rj-rj/rc)))

              
!--------------------------------------
        if(swirl_forcing.eq.1.and.rc.le.Jet_Circle) then

              ! ACTIVE CONTROL OF JET 

              if(zc.le.0.5d0*Jet_Circle) then
               if(xc.le.-0.5d0*rj.and.abs(yc).le.0.5d0*rj) then
                call random_number(rr)
                rr1 = 50.0d0 + 50.0d0*rr  
                vel(N9,M9,K9,1) = rr1 
                vel(N9,M9,K9,2) = -rr1 
               endif
               if(xc.ge.0.5d0*rj.and.abs(yc).le.0.5d0*rj) then
                call random_number(rr)
                rr1 = 50.0d0 + 50.0d0*rr  
                vel(N9,M9,K9,1) = -rr1 
                vel(N9,M9,K9,2) = rr1 
               endif
               if(yc.le.-0.5d0*rj.and.abs(xc).le.0.5d0*rj) then
                call random_number(rr)
                rr1 = 50.0d0 + 50.0d0*rr  
                vel(N9,M9,K9,1) = rr1 
                vel(N9,M9,K9,2) = rr1 
               endif
               if(yc.ge.0.5d0*rj.and.abs(xc).le.0.5d0*rj) then
                call random_number(rr)
                rr1 = 50.0d0 + 50.0d0*rr  
                vel(N9,M9,K9,1) = -rr1 
                vel(N9,M9,K9,2) = -rr1 
               endif
              endif

        endif
!--------------------------------------


           endif

  enddo
  enddo
  enddo

!--------------------------------------

  minz = minval(z(:,:,:))


     if(bogey_forcing.eq.1) then

        ! FORCING JET USING VORTEX RING 
        ! BOGEY & BAILLY, AIAA J., VOL. 43, N 5, PP. 1000-1007 (2005)
        ! BOGEY & BAILLY, PHYS. FLUIDS, VOL. 18, 065101 (2006)


  do K9 = K0, K5 ! NOTE: K5, NOT K2
  do M9 = M0, M5
  do N9 = N0, N5

   xc = x(N9,M9,K9)-Jet_Center(1)
   yc = y(N9,M9,K9)-Jet_Center(2)
   zc = z(N9,M9,K9) 
   rc = sqrt(xc**2.0d0 + yc**2.0d0)


       
              u_axial = 0.0d0
              u_radial = 0.0d0 

              theta = acos(xc/rc)


      if(minz.le.Jet_Circle) then

             if(rc.ne.0.0d0) then
              t2 = ((zc-rj)**2.0d0 + (rc-rj)**2.0d0)/(dx*dx)

              ! DANGER
              !t1 = 2.0d0*rj/(rc*dx)*exp(-log(2.0d0)*t2)
              t1 = 2.0d0*rj/(rc*dx)*exp(-log(1.2d0)*t2)

              ua_ring = t1*(rc-rj) 
              ur_ring = t1*(rj-zc) 

              do mm = 0, 15   
               call random_number(t1)
               rr = 2.0d0*t1 - 1.0d0
               call random_number(t2)
               rr1 = 2.0d0*pi*t2
               
               ! rr is epsilon_i from the paper
               ! rr1 is phi_i from the paper 

               ! DANGER  
               !bogey_fact = 0.01d0*ujet*rr*cos(dble(mm)*theta+rr1)
               bogey_fact = 0.2d0*ujet*rr*cos(dble(mm)*theta+rr1)

               u_axial = u_axial + bogey_fact*ua_ring 
               u_radial = u_radial + bogey_fact*ur_ring
              enddo
             endif


      endif

    vel(N9,M9,K9,1) = vel(N9,M9,K9,1) + u_radial*cos(theta) 
    vel(N9,M9,K9,2) = vel(N9,M9,K9,2) + u_radial*sin(theta) 
    vel(N9,M9,K9,3) = vel(N9,M9,K9,3) + u_axial  


  enddo
  enddo
  enddo


     endif   
!----------------------------------------------


   ! START THE SYNTHETIC EDDY METHOD 



          idir = 1 

73        continue                       


  do K9 = K0, K2
  do M9 = M0, M5
  do N9 = N0, N5

   xc = x(N9,M9,K9)
   yc = y(N9,M9,K9)
   zc = z(N9,M9,K9)

   rc = sqrt(xc**2.0d0 + yc**2.0d0)

     if(rc.le.Jet_Circle) then

       do mm = 1, neddies
          call shape_func(neddies,sigma_sem,xc,yc,zc,xeddy(mm,1:3),feddy)
          vel(N9,M9,K9,idir) = vel(N9,M9,K9,idir) + 1.0d0/sqrt(dble(neddies))*epsi(mm,idir)*feddy 
       enddo

     endif      

  enddo
  enddo
  enddo

 
                if(idir.lt.3) then
                 idir = idir + 1
                 goto 73 
                endif



!----------------------------------------------


             ! CONVECT THE SYNTHETIC EDDIES


! SINCE WE ARE USING RK3, WE NEED TO CONVECT THE SYNTHETIC EDDIES
! BY AN AMOUNT DT/3 EACH TIME THIS SUBROUTINE IS CALLED
! dt3 IS DT/3



         if(proc_injector.eq.1) then

                do mm = 1, neddies
                 xeddy(mm,3) = xeddy(mm,3) + ujet*dt3

                 if(xeddy(mm,3).gt.zhigh) then
                   xeddy(mm,3) = zlow

               call random_number(t1)
               xeddy(mm,1) = xlow + (xhigh-xlow)*t1 
               call random_number(t1)
               xeddy(mm,2) = ylow + (yhigh-ylow)*t1 

                   do nn = 1, 3 
                   call randno(rr)
                   if(rr.le.0.0d0) then
                    epsi(mm,nn) = -1.0d0
                   else
                    epsi(mm,nn) = 1.0d0 
                   endif 
                   enddo

                 endif
                enddo 

         endif


!----------------------------------------------


  do K9 = K0, K2
  do M9 = M0, M5
  do N9 = N0, N5

    ke = 0.5d0*(vel(N9,M9,K9,1)**2.0d0 + vel(N9,M9,K9,2)**2.0d0 + vel(N9,M9,K9,3)**2.0d0)


    ! DANGER
    ! NOTE: THIS IMPLICITLY ASSUMES DENSITY RATIO = 1

    Q(N9,M9,K9,1) = Jet_Rho
    Q(N9,M9,K9,2) = Jet_Rho*vel(N9,M9,K9,1)
    Q(N9,M9,K9,3) = Jet_Rho*vel(N9,M9,K9,2)
    Q(N9,M9,K9,4) = Jet_Rho*vel(N9,M9,K9,3)
    Q(N9,M9,K9,5) = Jet_P/(Jet_Gamma-ONE) + Jet_Rho*ke



   xc = x(N9,M9,K9)
   yc = y(N9,M9,K9)
   zc = z(N9,M9,K9)

   rc = sqrt(xc**2.0d0 + yc**2.0d0)

     ! DANGER 
     if(rc.le.Jet_Circle) then
       Q(N9,M9,K9,6) = 0.0d0 !Jet_Rho*Jet_f1  
     else
       Q(N9,M9,K9,6) = Jet_Rho*1.0d0 
     endif 


  enddo
  enddo
  enddo


!----------------------------------------------

!  do K9 = K2,K5 
!    if (z(N2,M2,K9) >= z(N2,M2,K2)+Jet_Height) EXIT 

!    Q(:,:,K9,:) = Q(:,:,K2,:)
!  enddo
  
#if defined (PARALLEL_MPI)
  endif
#endif

  END SUBROUTINE Injector_z

END SUBROUTINE Injector

END MODULE Injector_3D


!-----------------------------------------------

! KAUSHIK CREATED THIS SUBROUTINE

                  subroutine randno(x)

          implicit none

             real(kind=8) :: x

             call random_number(x)

             
             x = 2.0d0*x - 1.0d0


                 return
                 end subroutine

!--------------------------------------------------

         subroutine shape_func(neddy,sigma,xc,yc,zc,xe,f)

           implicit none

                  integer :: neddy
                  real(kind=8) :: sigma,xc, yc, zc, xe(1:3), f, fx, fy, fz
                  real(kind=8) :: xx, yy, zz 


                  xx = abs((xc - xe(1))/sigma)
                  yy = abs((yc - xe(2))/sigma)
                  zz = abs((zc - xe(3))/sigma)

                  if(xx.lt.1.0d0) then
                    fx = sqrt(3.0d0/2.0d0)*(1.0d0-xx)
                  else
                    fx = 0.0d0
                  endif
                  
                  if(yy.lt.1.0d0) then
                    fy = sqrt(3.0d0/2.0d0)*(1.0d0-yy)
                  else
                    fy = 0.0d0
                  endif

                  if(zz.lt.1.0d0) then
                    fz = sqrt(3.0d0/2.0d0)*(1.0d0-zz)
                  else
                    fz = 0.0d0
                  endif  


                   f = fx*fy*fz

                 
                   f = f*100.0d0


                   return
                   end subroutine     

!--------------------------------------------------

