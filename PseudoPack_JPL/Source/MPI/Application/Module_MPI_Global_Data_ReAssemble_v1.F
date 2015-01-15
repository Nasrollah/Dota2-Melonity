#include "pseudopack.h"


MODULE MPI_Global_Data_ReAssemble_v1

  USE MPI_Data_ReAssemble_v1

implicit NONE

INTERFACE PS_MPI_Global_Data_ReAssemble_v1
  MODULE PROCEDURE  Global_Data_ReAssemble_0D
  MODULE PROCEDURE  Global_Data_ReAssemble_1D
  MODULE PROCEDURE  Global_Data_ReAssemble_2D
  MODULE PROCEDURE  Global_Data_ReAssemble_3D
END INTERFACE

PRIVATE

PUBLIC  :: PS_MPI_Global_Data_ReAssemble_v1

CONTAINS
!
!===================================================================
!
  Subroutine Global_Data_ReAssemble_0D                              &
                            (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                             NV, Q, Q_g, MPI_Comm_Type)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: NV
  integer  :: nn_0, k

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0  :N5  ), TARGET  :: Q
  REALTYPE, dimension(N0_g:N5_g), TARGET  :: Q_g

  REALTYPE, dimension(:)           , POINTER :: f, f_g

  nn_0 = n_0-M_s_x

    f => Q ; f_g => Q_g

    call PS_MPI_Data_ReAssemble_v1 (f  , N0  , N5  ,                         &
                                    f_g, N0_g, N5_g,                         &
                                    nn_0, MPI_Comm_Type)

  END Subroutine Global_Data_ReAssemble_0D
!
!===================================================================
!
  Subroutine Global_Data_ReAssemble_1D                              &
                            (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                             NV, Q, Q_g, MPI_Comm_Type)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: NV
  integer  :: nn_0, k

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0  :N5  ,NV), TARGET  :: Q
  REALTYPE, dimension(N0_g:N5_g,NV), TARGET  :: Q_g

  REALTYPE, dimension(:)           , POINTER :: f, f_g

  nn_0 = n_0-M_s_x 

  do k = 1,NV
    f => Q(:,k) ; f_g => Q_g(:,k)

    call PS_MPI_Data_ReAssemble_v1 (f  , N0  , N5  ,                         &
                                    f_g, N0_g, N5_g,                         &
                                    nn_0, MPI_Comm_Type)
  enddo

  END Subroutine Global_Data_ReAssemble_1D
!
!===================================================================
!
  Subroutine Global_Data_ReAssemble_2D                              &
                            (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                             M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                             NV, Q, Q_g, MPI_Comm_Type)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: NV
  integer  :: nn_0, mm_0, k

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,NV), TARGET  :: Q
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,NV), TARGET  :: Q_g

  REALTYPE, dimension(:,:)                   , POINTER :: f, f_g

  nn_0 = n_0-M_s_x ; mm_0 = m_0-M_s_y 

  do k = 1,NV
    f => Q(:,:,k) ; f_g => Q_g(:,:,k)

    call PS_MPI_Data_ReAssemble_v1 (f  , N0  , N5  , M0  , M5  ,             &
                                    f_g, N0_g, N5_g, M0_g, M5_g,             &
                                    nn_0, mm_0, MPI_Comm_Type)
  enddo

  END Subroutine Global_Data_ReAssemble_2D
!
!===================================================================
!
  Subroutine Global_Data_ReAssemble_3D                              &
                            (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                             M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                             K0_g, K5_g, K0, K5, k_0, M_s_z,        &
                             NV, Q, Q_g, MPI_Comm_Type)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: K0_g, K5_g, K0, K5, k_0, M_s_z
  integer  :: NV

  integer  :: nn_0, mm_0, kk_0, k

  integer , OPTIONAL :: MPI_Comm_Type

  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  ,NV), TARGET  :: Q
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g,NV), TARGET  :: Q_g

  REALTYPE, dimension(:,:,:)                           , POINTER :: f, f_g

  nn_0 = n_0-M_s_x ; mm_0 = m_0-M_s_y ; kk_0 = k_0-M_s_z

  do k = 1,NV
    f => Q(:,:,:,k) ; f_g => Q_g(:,:,:,k)

    call PS_MPI_Data_ReAssemble_v1 (f  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                    f_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                    nn_0, mm_0, kk_0, MPI_Comm_Type)
  enddo

  END Subroutine Global_Data_ReAssemble_3D

END MODULE MPI_Global_Data_ReAssemble_v1
