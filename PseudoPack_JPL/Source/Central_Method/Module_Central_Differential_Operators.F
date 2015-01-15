#include "pseudopack.h"


MODULE Central_Differential_Operators

  USE Central_Derivative

IMPLICIT NONE

INTERFACE PS_Central_D
  MODULE PROCEDURE PS_Central_Derivative_1D
  MODULE PROCEDURE PS_Central_Derivative_2D
  MODULE PROCEDURE PS_Central_Derivative_3D
END INTERFACE

INTERFACE PS_Central_Gradient
  MODULE PROCEDURE PS_Central_Gradient_1D
  MODULE PROCEDURE PS_Central_Gradient_2D
  MODULE PROCEDURE PS_Central_Gradient_3D
END INTERFACE

INTERFACE PS_Central_Divergence
  MODULE PROCEDURE PS_Central_Divergence_1D
  MODULE PROCEDURE PS_Central_Divergence_2D
  MODULE PROCEDURE PS_Central_Divergence_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_Central_D
PUBLIC  :: PS_Central_Gradient
PUBLIC  :: PS_Central_Divergence

CONTAINS
!
! =======================================================================
!
  Subroutine PS_Central_Derivative_1D (                                    &
                        N0, N5, N2, N3, Order_x, dx, D_Order_x, Degree_x,  &
                                       Index, U, Field                     )
  
  integer  :: NU, Index
  integer  :: N0, N5, N2, N3, Order_x, D_Order_x, Degree_x

  REALTYPE :: dx

  REALTYPE, dimension(N0:N5) :: U
  REALTYPE, dimension(N0:N5) :: Field

  integer  :: D_Order, N_D, Degree
  REALTYPE :: delta

  integer  :: n

      n = Index

      call Adjust_Data (Order_x,Degree_x,D_Order_x,dx, D_Order,Degree,N_D,delta)

    call PS_Central_Derivative (n, D_Order, Degree, N_D, delta,                &
                                N0,N5,N2,N3,                                   &
                                U, Field                                       )

  END Subroutine PS_Central_Derivative_1D 
!
! =======================================================================
!
  Subroutine PS_Central_Gradient_1D (                                      &
                        N0, N5, N2, N3, Order_x, dx, D_Order_x, Degree_x,  &
                                        NV, U, Index, NU, Field            )

  integer  :: NV, NU, Index
  integer  :: N0, N5, N2, N3, Order_x, D_Order_x, Degree_x

  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,NV) :: U
  REALTYPE, dimension(N0:N5,NU) :: Field

  integer  :: D_Order, N_D, Degree
  REALTYPE :: delta

  integer  :: n

  do n = 1,NU
      call Adjust_Data (Order_x,Degree_x,D_Order_x,dx, D_Order,Degree,N_D,delta)

    D_Order = 1
    call PS_Central_Derivative (n, D_Order, Degree, N_D, delta,                &
                                N0,N5,N2,N3,                                   &
                                U(:,Index), Field(:,n)                         )
   enddo

  END Subroutine PS_Central_Gradient_1D 
!
! =======================================================================
!
  Subroutine PS_Central_Divergence_1D (                                    &
                        N0, N5, N2, N3, Order_x, dx, D_Order_x, Degree_x,  &
                                    NU, U, Field                           )

  integer  :: NU
  integer  :: N0, N5, N2, N3, Order_x, D_Order_x, Degree_x

  REALTYPE :: dx

  REALTYPE, dimension(N0:N5,NU) :: U
  REALTYPE, dimension(N0:N5)    :: Field, D_U

  integer  :: D_Order, N_D, Degree
  REALTYPE :: delta

  integer  :: n

  n = 1

      call Adjust_Data (Order_x,Degree_x,D_Order_x,dx, D_Order,Degree,N_D,delta)

    D_Order = 1
    call PS_Central_Derivative (n, D_Order, Degree, N_D, delta,                &
                                N0,N5,N2,N3,                                   &
                                U(:,n), Field                                  )

  END Subroutine PS_Central_Divergence_1D  
!
! =======================================================================
!
  Subroutine PS_Central_Derivative_2D (                                    &
                        N0, N5, N2, N3, Order_x, dx, D_Order_x, Degree_x,  &
                        M0, M5, M2, M3, Order_y, dy, D_Order_y, Degree_y,  &
                                       Index, U, Field                     )
  
  integer  :: NU, Index
  integer  :: N0, N5, N2, N3, Order_x, D_Order_x, Degree_x
  integer  :: M0, M5, M2, M3, Order_y, D_Order_y, Degree_y

  REALTYPE :: dx, dy

  REALTYPE, dimension(N0:N5,M0:M5) :: U
  REALTYPE, dimension(N0:N5,M0:M5) :: Field

  integer  :: D_Order, N_D, Degree
  REALTYPE :: delta

  integer  :: n

      n = Index

    if (n == 1) then
      call Adjust_Data (Order_x,Degree_x,D_Order_x,dx, D_Order,Degree,N_D,delta)
    endif
    if (n == 2) then
      call Adjust_Data (Order_y,Degree_y,D_Order_y,dy, D_Order,Degree,N_D,delta)
    endif

    call PS_Central_Derivative (n, D_Order, Degree, N_D, delta,                &
                                N0,N5,N2,N3, M0,M5,M2,M3,                      &
                                U, Field                                       )

  END Subroutine PS_Central_Derivative_2D 
!
! =======================================================================
!
  Subroutine PS_Central_Gradient_2D (                                      &
                        N0, N5, N2, N3, Order_x, dx, D_Order_x, Degree_x,  &
                        M0, M5, M2, M3, Order_y, dy, D_Order_y, Degree_y,  &
                                        NV, U, Index, NU, Field            )

  integer  :: NV, NU, Index
  integer  :: N0, N5, N2, N3, Order_x, D_Order_x, Degree_x
  integer  :: M0, M5, M2, M3, Order_y, D_Order_y, Degree_y

  REALTYPE :: dx, dy

  REALTYPE, dimension(N0:N5,M0:M5,NV) :: U
  REALTYPE, dimension(N0:N5,M0:M5,NU) :: Field

  integer  :: D_Order, N_D, Degree
  REALTYPE :: delta

  integer  :: n

  do n = 1,NU
    if (n == 1) then
      call Adjust_Data (Order_x,Degree_x,D_Order_x,dx, D_Order,Degree,N_D,delta)
    endif
    if (n == 2) then
      call Adjust_Data (Order_y,Degree_y,D_Order_y,dy, D_Order,Degree,N_D,delta)
    endif

    D_Order = 1
    call PS_Central_Derivative (n, D_Order, Degree, N_D, delta,                &
                                N0,N5,N2,N3, M0,M5,M2,M3,                      &
                                U(:,:,Index), Field(:,:,n)                     )
   enddo

  END Subroutine PS_Central_Gradient_2D 
!
! =======================================================================
!
  Subroutine PS_Central_Divergence_2D (                                    &
                        N0, N5, N2, N3, Order_x, dx, D_Order_x, Degree_x,  &
                        M0, M5, M2, M3, Order_y, dy, D_Order_y, Degree_y,  &
                                    NU, U, Field                           )

  integer  :: NU 
  integer  :: N0, N5, N2, N3, Order_x, D_Order_x, Degree_x
  integer  :: M0, M5, M2, M3, Order_y, D_Order_y, Degree_y

  REALTYPE :: dx, dy

  REALTYPE, dimension(N0:N5,M0:M5,NU) :: U
  REALTYPE, dimension(N0:N5,M0:M5)    :: Field, D_U

  integer  :: D_Order, N_D, Degree
  REALTYPE :: delta

  integer  :: n

  Field = ZERO 

  do n = 1,NU
    if (n == 1) then
      call Adjust_Data (Order_x,Degree_x,D_Order_x,dx, D_Order,Degree,N_D,delta)
    endif
    if (n == 2) then
      call Adjust_Data (Order_y,Degree_y,D_Order_y,dy, D_Order,Degree,N_D,delta)
    endif

    D_U = ZERO 

    D_Order = 1
    call PS_Central_Derivative (n, D_Order, Degree, N_D, delta,                &
                                N0,N5,N2,N3, M0,M5,M2,M3,                      &
                                U(:,:,n), D_U                                  )
    Field = Field + D_U
  enddo

  END Subroutine PS_Central_Divergence_2D  
!
! =======================================================================
!
  Subroutine PS_Central_Derivative_3D (                                    &
                        N0, N5, N2, N3, Order_x, dx, D_Order_x, Degree_x,  &
                        M0, M5, M2, M3, Order_y, dy, D_Order_y, Degree_y,  &
                        K0, K5, K2, K3, Order_z, dz, D_Order_z, Degree_z,  &
                                       Index, U, Field                     )
  
  integer  :: NU, Index
  integer  :: N0, N5, N2, N3, Order_x, D_Order_x, Degree_x
  integer  :: M0, M5, M2, M3, Order_y, D_Order_y, Degree_y
  integer  :: K0, K5, K2, K3, Order_z, D_Order_z, Degree_z

  REALTYPE :: dx, dy, dz

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: U
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Field

  integer  :: D_Order, N_D, Degree
  REALTYPE :: delta

  integer  :: n

      n = Index

    if (n == 1) then
      call Adjust_Data (Order_x,Degree_x,D_Order_x,dx, D_Order,Degree,N_D,delta)
    endif
    if (n == 2) then
      call Adjust_Data (Order_y,Degree_y,D_Order_y,dy, D_Order,Degree,N_D,delta)
    endif
    if (n == 3) then
      call Adjust_Data (Order_z,Degree_z,D_Order_z,dz, D_Order,Degree,N_D,delta)
    endif

    write (6,*) n, D_Order, Degree, N_D, delta
    call PS_Central_Derivative (n, D_Order, Degree, N_D, delta,                &
                                N0,N5,N2,N3, M0,M5,M2,M3, K0,K5,K2,K3,         &
                                U, Field                                       )

  END Subroutine PS_Central_Derivative_3D 
!
! =======================================================================
!
  Subroutine PS_Central_Gradient_3D (                                      &
                        N0, N5, N2, N3, Order_x, dx, D_Order_x, Degree_x,  &
                        M0, M5, M2, M3, Order_y, dy, D_Order_y, Degree_y,  &
                        K0, K5, K2, K3, Order_z, dz, D_Order_z, Degree_z,  &
                                        NV, U, Index, NU, Field            )

  integer  :: NV, NU, Index
  integer  :: N0, N5, N2, N3, Order_x, D_Order_x, Degree_x
  integer  :: M0, M5, M2, M3, Order_y, D_Order_y, Degree_y
  integer  :: K0, K5, K2, K3, Order_z, D_Order_z, Degree_z

  REALTYPE :: dx, dy, dz

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: U
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NU) :: Field

  integer  :: D_Order, N_D, Degree
  REALTYPE :: delta

  integer  :: n

  do n = 1,NU
    if (n == 1) then
      call Adjust_Data (Order_x,Degree_x,D_Order_x,dx, D_Order,Degree,N_D,delta)
    endif
    if (n == 2) then
      call Adjust_Data (Order_y,Degree_y,D_Order_y,dy, D_Order,Degree,N_D,delta)
    endif
    if (n == 3) then
      call Adjust_Data (Order_z,Degree_z,D_Order_z,dz, D_Order,Degree,N_D,delta)
    endif

    D_Order = 1
    call PS_Central_Derivative (n, D_Order, Degree, N_D, delta,               &
                                N0,N5,N2,N3, M0,M5,M2,M3, K0,K5,K2,K3,        &
                                U(:,:,:,Index), Field(:,:,:,n)                )
   enddo

  END Subroutine PS_Central_Gradient_3D 
!
! =======================================================================
!
  Subroutine PS_Central_Divergence_3D (                                    &
                        N0, N5, N2, N3, Order_x, dx, D_Order_x, Degree_x,  &
                        M0, M5, M2, M3, Order_y, dy, D_Order_y, Degree_y,  &
                        K0, K5, K2, K3, Order_z, dz, D_Order_z, Degree_z,  &
                                    NU, U, Field                           )

  integer  :: NU
  integer  :: N0, N5, N2, N3, Order_x, D_Order_x, Degree_x
  integer  :: M0, M5, M2, M3, Order_y, D_Order_y, Degree_y
  integer  :: K0, K5, K2, K3, Order_z, D_Order_z, Degree_z

  REALTYPE :: dx, dy, dz

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NU) :: U
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: Field, D_U

  integer  :: D_Order, N_D, Degree
  REALTYPE :: delta

  integer  :: n

  Field = ZERO 

  do n = 1,NU
    if (n == 1) then
      call Adjust_Data (Order_x,Degree_x,D_Order_x,dx, D_Order,Degree,N_D,delta)
    endif
    if (n == 2) then
      call Adjust_Data (Order_y,Degree_y,D_Order_y,dy, D_Order,Degree,N_D,delta)
    endif
    if (n == 3) then
      call Adjust_Data (Order_z,Degree_z,D_Order_z,dz, D_Order,Degree,N_D,delta)
    endif

    D_U = ZERO 

    D_Order = 1
    call PS_Central_Derivative (n, D_Order, Degree, N_D, delta,                &
                                N0,N5,N2,N3, M0,M5,M2,M3, K0,K5,K2,K3,         &
                                U(:,:,:,n), D_U                                )
    Field = Field + D_U
   enddo

  END Subroutine PS_Central_Divergence_3D 
!
!============================================================================
!
  Subroutine Adjust_Data (Order_x, Degree_x, D_Order_x, dx,  &
                          D_Order, Degree, N_D, delta)
 
  integer  :: Order_x, Degree_x, D_Order_x, D_Order, Degree, N_D
  REALTYPE :: dx, delta

  integer  :: M_s_x, O_s_x

  M_s_x = (Order_x+1)/2 ; O_s_x = MIN(M_s_x,(Degree_x+1)/2)

  D_Order = D_Order_x ; Degree = Degree_x ; N_D = 2*O_s_x+1 ; delta = dx

  END Subroutine Adjust_Data 

END MODULE Central_Differential_Operators
