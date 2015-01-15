#include "pseudopack.h"

MODULE MPI_Extend_Domain_Data

  USE Processor
  USE MPI_Data_ReAssemble
  USE MPI_Data_DisAssemble

implicit NONE

INTERFACE PS_MPI_Extend_Domain_Data
  MODULE PROCEDURE MPI_Extended_Domain_Data_1D
  MODULE PROCEDURE MPI_Extended_Domain_Data_2D
  MODULE PROCEDURE MPI_Extended_Domain_Data_3D

  MODULE PROCEDURE     Extended_Domain_Data_1D
  MODULE PROCEDURE     Extended_Domain_Data_2D
  MODULE PROCEDURE     Extended_Domain_Data_3D
END INTERFACE

PRIVATE
PUBLIC  :: PS_FD_Resize_Domain
PUBLIC  :: PS_MPI_Extend_Domain_Data

CONTAINS

  Subroutine PS_FD_Resize_Domain (Adaptive, Adaptive_x0, Adaptive_x1,    &
                                              Extend_x0,   Extend_x1,    &
                                  Symmetry, dx,                          &
                                  x0_Extension, x1_Extension,            & 
                                  x0_Limit, x1_Limit,                    &
                                  N_x, N0_g, N5_g, N0, N5, n_0, x0, x1,  &
                                  L_x, L0_g, L5_g, L0, L5, l_0, r0, r1,  &
                                  Index_Shift)

  logical  :: Adaptive, Adaptive_x0, Adaptive_x1
  logical  ::             Extend_x0,   Extend_x1
  integer  :: Symmetry

  integer  :: N_x, N0_g, N5_g, N0, N5, n_0
  integer  :: L_x, L0_g, L5_g, L0, L5, l_0

  REALTYPE :: dx, x0_Extension, x1_Extension, x0_Limit, x1_Limit
  REALTYPE :: x0, x1, r0, r1
  integer  :: Index_Shift
  integer  :: N_x0, N_x1

  L_x = N_x ; L0_g = N0_g ; L5_g = N5_g ; L0 = N0 ; L5 = N5 ; l_0 = n_0
  
  r0 = x0 ; r1 = x1 ; Index_Shift = 0

  if ((.NOT. Adaptive) .OR. (Symmetry == 1)) RETURN

  N_x0 = 0 
  if (Extend_x0 .AND. Adaptive_x0) then
    N_x0 = ABS(x0_Extension/dx)

    if (x0 - N_x0*dx <  x0_Limit) then
      N_x0 = ABS(x0_Limit-x0)/dx ; Adaptive_x0 = .FALSE.
    endif
  endif
 

  N_x1 = 0 
  if (Extend_x1 .AND. Adaptive_x1) then
    N_x1 = ABS(x1_Extension/dx)

    if (x1 + N_x1*dx >  x1_Limit) then
      N_x1 = ABS(x1_Limit-x1)/dx ; Adaptive_x1 = .FALSE.
    endif
  endif
  
  x0 = x0 - N_x0*dx 
  x1 = x1 + N_x1*dx

  N_x = N_x + N_x0 + N_x1 ; Index_Shift = N_x0
  
  if (.NOT. (Adaptive_x0 .OR. Adaptive_x1)) Adaptive = .FALSE.

  END Subroutine PS_FD_Resize_Domain 
!
!=======================================================================
!
  Subroutine MPI_Extended_Domain_Data_1D (Q_Old, L0_g, L5_g, L0, L5, n_00, &
                                          Q_New, N0_g, N5_g, N0, N5, n_0 , &
                                          Distributed,                     &
                                          Spatial_Order_x, Index_Shift_x   )

  integer  :: L0_g, L5_g, L0, L5, n_00
  integer  :: N0_g, N5_g, N0, N5, n_0

  logical , OPTIONAL :: Distributed
  integer , OPTIONAL :: Spatial_Order_x, Index_Shift_x

  integer  :: i_Shift, M_s_x, is, ie
  logical  :: Distribute

  REALTYPE, dimension(L0:L5) :: Q_Old
  REALTYPE, dimension(N0:N5) :: Q_New

  REALTYPE, dimension(:), ALLOCATABLE :: P_Old, P_New

                                Distribute = .FALSE. 
  if (PRESENT(Distributed)    ) Distribute = Distributed

                                M_s_x = 0 
  if (PRESENT(Spatial_Order_x)) M_s_x = (Spatial_Order_x+1)/2

                                i_Shift = 0 
  if (PRESENT(Index_Shift_x)  ) i_Shift = ABS(Index_Shift_x)

  if (Distribute) then
    if (I_Am == First_Processor)    &
      ALLOCATE (P_Old(L0_g:L5_g), P_New(N0_g:N5_g))

    call PS_MPI_Data_ReAssemble  (Q_Old, L0  , L5,           &
                                  P_Old, L0_g, L5_g, n_00-M_s_x) 

    call PS_MPI_Data_ReAssemble  (Q_New, N0  , N5  ,         &
                                  P_New, N0_g, N5_g, n_0 -M_s_x)

    if (I_Am == First_Processor) then
      is = N0_g+i_Shift ; ie = is+SIZE(P_Old,DIM=1)-1

        P_New(is:ie) = P_Old
    endif

    call PS_MPI_Data_DisAssemble (P_New, N0_g, N5_g,         &
                                  Q_New, N0  , N5  , n_0 -M_s_x)

    if (I_Am == First_Processor) DEALLOCATE (P_Old, P_New)

  else
      is = N0_g+i_Shift ; ie = is+SIZE(Q_Old,DIM=1)-1

        Q_New(is:ie) = Q_Old
  endif
  
  END Subroutine MPI_Extended_Domain_Data_1D
!
!=======================================================================
!
  Subroutine MPI_Extended_Domain_Data_2D (Q_Old, L0_g, L5_g, L0, L5, n_00, &
                                                 I0_g, I5_g, I0, I5, m_00, &
                                          Q_New, N0_g, N5_g, N0, N5, n_0 , &
                                                 M0_g, M5_g, M0, M5, m_0 , &
                                          Distributed,                     &
                                          Spatial_Order_x, Index_Shift_x,  &
                                          Spatial_Order_y, Index_Shift_y   )

  integer  :: L0_g, L5_g, L0, L5, n_00
  integer  :: I0_g, I5_g, I0, I5, m_00
  integer  :: N0_g, N5_g, N0, N5, n_0 
  integer  :: M0_g, M5_g, M0, M5, m_0 

  logical , OPTIONAL :: Distributed
  integer , OPTIONAL :: Spatial_Order_x, Index_Shift_x
  integer , OPTIONAL :: Spatial_Order_y, Index_Shift_y

  integer  :: i_Shift, M_s_x, is, ie
  integer  :: j_Shift, M_s_y, js, je
  logical  :: Distribute

  REALTYPE, dimension(L0:L5,I0:I5) :: Q_Old
  REALTYPE, dimension(N0:N5,M0:M5) :: Q_New

  REALTYPE, dimension(:,:), ALLOCATABLE :: P_Old, P_New

                                Distribute = .FALSE. 
  if (PRESENT(Distributed)    ) Distribute = Distributed

                                M_s_x = 0 
  if (PRESENT(Spatial_Order_x)) M_s_x = (Spatial_Order_x+1)/2

                                M_s_y = 0 
  if (PRESENT(Spatial_Order_y)) M_s_y = (Spatial_Order_y+1)/2

                                i_Shift = 0 
  if (PRESENT(Index_Shift_x)  ) i_Shift = ABS(Index_Shift_x)

                                j_Shift = 0 
  if (PRESENT(Index_Shift_y)  ) j_Shift = ABS(Index_Shift_y)

  if (Distribute) then
    if (I_Am == First_Processor)              &
      ALLOCATE (P_Old(L0_g:L5_g,I0_g:I5_g),   &
                P_New(N0_g:N5_g,M0_g:M5_g))

    call PS_MPI_Data_ReAssemble  (Q_Old, L0  , L5  , I0  , I5  ,  &
                                  P_Old, L0_g, L5_g, I0_g, I5_g,  &
                                  n_00-M_s_x, m_00-M_s_y) 

    call PS_MPI_Data_ReAssemble  (Q_New, N0  , N5  , M0  , M5  ,  &
                                  P_New, N0_g, N5_g, M0_g, M5_g,  &
                                  n_0 -M_s_x, m_0 -M_s_y) 

    if (I_Am == First_Processor) then
      is = N0_g+i_Shift ; ie = is+SIZE(P_Old,DIM=1)-1
      js = M0_g+j_Shift ; je = js+SIZE(P_Old,DIM=2)-1

        P_New(is:ie,js:je) = P_Old
    endif

    call PS_MPI_Data_DisAssemble (P_New, N0_g, N5_g, M0_g, M5_g,  &
                                  Q_New, N0  , N5,   M0  , M5  ,  & 
                                  n_0 -M_s_x, m_0 -M_s_y) 

    if (I_Am == First_Processor) DEALLOCATE (P_Old, P_New)
      
  else
      is = N0_g+i_Shift ; ie = is+SIZE(Q_Old,DIM=1)-1
      js = M0_g+j_Shift ; je = js+SIZE(Q_Old,DIM=2)-1

        Q_New(is:ie,js:je) = Q_Old
  endif
  
  END Subroutine MPI_Extended_Domain_Data_2D
!
!=======================================================================
!
  Subroutine MPI_Extended_Domain_Data_3D (Q_Old, L0_g, L5_g, L0, L5, n_00, &
                                                 I0_g, I5_g, I0, I5, m_00, &
                                                 J0_g, J5_g, J0, J5, k_00, &
                                          Q_New, N0_g, N5_g, N0, N5, n_0 , &
                                                 M0_g, M5_g, M0, M5, m_0 , &
                                                 K0_g, K5_g, K0, K5, k_0 , &
                                          Distributed,                     &
                                          Spatial_Order_x, Index_Shift_x,  &
                                          Spatial_Order_y, Index_Shift_y,  &
                                          Spatial_Order_z, Index_Shift_z   )

  integer  :: L0_g, L5_g, L0, L5, n_00
  integer  :: I0_g, I5_g, I0, I5, m_00
  integer  :: J0_g, J5_g, J0, J5, k_00
  integer  :: N0_g, N5_g, N0, N5, n_0 
  integer  :: M0_g, M5_g, M0, M5, m_0 
  integer  :: K0_g, K5_g, K0, K5, k_0 

  logical , OPTIONAL :: Distributed
  integer , OPTIONAL :: Spatial_Order_x, Index_Shift_x
  integer , OPTIONAL :: Spatial_Order_y, Index_Shift_y
  integer , OPTIONAL :: Spatial_Order_z, Index_Shift_z

  integer  :: i_Shift, M_s_x, is, ie
  integer  :: j_Shift, M_s_y, js, je
  integer  :: k_Shift, M_s_z, ks, ke
  logical  :: Distribute

  REALTYPE, dimension(L0:L5,I0:I5,J0:J5) :: Q_Old
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5) :: Q_New

  REALTYPE, dimension(:,:,:), ALLOCATABLE :: P_Old, P_New

                                Distribute = .FALSE. 
  if (PRESENT(Distributed)    ) Distribute = Distributed

                                M_s_x = 0 
  if (PRESENT(Spatial_Order_x)) M_s_x = (Spatial_Order_x+1)/2

                                M_s_y = 0 
  if (PRESENT(Spatial_Order_y)) M_s_y = (Spatial_Order_y+1)/2

                                M_s_z = 0 
  if (PRESENT(Spatial_Order_z)) M_s_z = (Spatial_Order_z+1)/2

                                i_Shift = 0 
  if (PRESENT(Index_Shift_x)  ) i_Shift = ABS(Index_Shift_x)

                                j_Shift = 0 
  if (PRESENT(Index_Shift_y)  ) j_Shift = ABS(Index_Shift_y)

                                k_Shift = 0 
  if (PRESENT(Index_Shift_z)  ) k_Shift = ABS(Index_Shift_z)

  if (Distribute) then
    if (I_Am == First_Processor)                        &
      ALLOCATE (P_Old(L0_g:L5_g,I0_g:I5_g,J0_g:J5_g),   &
                P_New(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g))

    call PS_MPI_Data_ReAssemble  (Q_Old, L0  , L5  , I0  , I5  , J0  , J5  ,  &
                                  P_Old, L0_g, L5_g, I0_g, I5_g, J0_g, J5_g,  &
                                         n_00-M_s_x, m_00-M_s_y, k_00-M_s_z) 

    call PS_MPI_Data_ReAssemble  (Q_New, N0  , N5  , M0  , M5  , K0  , K5  ,  &
                                  P_New, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g,  &
                                         n_0 -M_s_x, m_0 -M_s_y, k_0 -M_s_z) 

    if (I_Am == First_Processor) then
      is = N0_g+i_Shift ; ie = is+SIZE(P_Old,DIM=1)-1
      js = M0_g+j_Shift ; je = js+SIZE(P_Old,DIM=2)-1
      ks = K0_g+k_Shift ; ke = ks+SIZE(P_Old,DIM=3)-1

        P_New(is:ie,js:je,ks:ke) = P_Old
    endif

    call PS_MPI_Data_DisAssemble (P_New, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g,  &
                                  Q_New, N0  , N5  , M0  , M5  , K0  , K5  ,  &
                                         n_0 -M_s_x, m_0 -M_s_y, k_0 -M_s_z) 

    if (I_Am == First_Processor) DEALLOCATE (P_Old, P_New)

  else
      is = N0_g+i_Shift ; ie = is+SIZE(Q_Old,DIM=1)-1
      js = M0_g+j_Shift ; je = js+SIZE(Q_Old,DIM=2)-1
      ks = K0_g+k_Shift ; ke = ks+SIZE(Q_Old,DIM=3)-1

        Q_New(is:ie,js:je,ks:ke) = Q_Old
  endif
  
  END Subroutine MPI_Extended_Domain_Data_3D
!
!=======================================================================
!
  Subroutine Extended_Domain_Data_1D (Q_Old, L0_g, L5_g, &
                                      Q_New, N0_g, N5_g, &
                                      Index_Shift_x)

  integer  :: L0_g, L5_g
  integer  :: N0_g, N5_g

  integer , OPTIONAL :: Index_Shift_x

  integer  :: i_Shift, is, ie

  REALTYPE, dimension(L0_g:L5_g) :: Q_Old
  REALTYPE, dimension(N0_g:N5_g) :: Q_New

                                i_Shift = 0 
  if (PRESENT(Index_Shift_x)  ) i_Shift = ABS(Index_Shift_x)

      is = N0_g+i_Shift ; ie = is+SIZE(Q_Old,DIM=1)-1

        Q_New(is:ie) = Q_Old

  END Subroutine Extended_Domain_Data_1D
!
!=======================================================================
!
  Subroutine Extended_Domain_Data_2D (Q_Old, L0_g, L5_g, &
                                             I0_g, I5_g, &
                                      Q_New, N0_g, N5_g, &
                                             M0_g, M5_g, &
                                      Index_Shift_x, Index_Shift_y   )

  integer  :: L0_g, L5_g
  integer  :: I0_g, I5_g
  integer  :: N0_g, N5_g
  integer  :: M0_g, M5_g

  integer , OPTIONAL :: Index_Shift_x
  integer , OPTIONAL :: Index_Shift_y

  integer  :: i_Shift, is, ie
  integer  :: j_Shift, js, je

  REALTYPE, dimension(L0_g:L5_g,I0_g:I5_g) :: Q_Old
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g) :: Q_New

                                i_Shift = 0 
  if (PRESENT(Index_Shift_x)  ) i_Shift = ABS(Index_Shift_x)

                                j_Shift = 0 
  if (PRESENT(Index_Shift_y)  ) j_Shift = ABS(Index_Shift_y)

      is = N0_g+i_Shift ; ie = is+SIZE(Q_Old,DIM=1)-1
      js = M0_g+j_Shift ; je = js+SIZE(Q_Old,DIM=2)-1

        Q_New(is:ie,js:je) = Q_Old

  END Subroutine Extended_Domain_Data_2D
!
!=======================================================================
!
  Subroutine Extended_Domain_Data_3D (Q_Old, L0_g, L5_g, &
                                             I0_g, I5_g, &
                                             J0_g, J5_g, &
                                      Q_New, N0_g, N5_g, &
                                             M0_g, M5_g, &
                                             K0_g, K5_g, &
                                      Index_Shift_x,     &
                                      Index_Shift_y,     &
                                      Index_Shift_z   )

  integer  :: L0_g, L5_g
  integer  :: I0_g, I5_g
  integer  :: J0_g, J5_g
  integer  :: N0_g, N5_g
  integer  :: M0_g, M5_g
  integer  :: K0_g, K5_g

  integer , OPTIONAL :: Index_Shift_x
  integer , OPTIONAL :: Index_Shift_y
  integer , OPTIONAL :: Index_Shift_z

  integer  :: i_Shift, is, ie
  integer  :: j_Shift, js, je
  integer  :: k_Shift, ks, ke

  REALTYPE, dimension(L0_g:L5_g,I0_g:I5_g,J0_g:J5_g) :: Q_Old
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g) :: Q_New

                                i_Shift = 0 
  if (PRESENT(Index_Shift_x)  ) i_Shift = ABS(Index_Shift_x)

                                j_Shift = 0 
  if (PRESENT(Index_Shift_y)  ) j_Shift = ABS(Index_Shift_y)

                                k_Shift = 0 
  if (PRESENT(Index_Shift_z)  ) k_Shift = ABS(Index_Shift_z)

      is = N0_g+i_Shift ; ie = is+SIZE(Q_Old,DIM=1)-1
      js = M0_g+j_Shift ; je = js+SIZE(Q_Old,DIM=2)-1
      ks = K0_g+k_Shift ; ke = ks+SIZE(Q_Old,DIM=3)-1

        Q_New(is:ie,js:je,ks:ke) = Q_Old
  
  END Subroutine Extended_Domain_Data_3D

END MODULE MPI_Extend_Domain_Data
