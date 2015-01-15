#include "pseudopack.h"

MODULE MPI_3D_Data_Slice_2D

  USE PS_IO_Unit, ONLY: lid6
  USE Processor
  USE File_OPEN
  USE Filename_Maker
  USE MPI_Data_ReAssemble
  USE MPI_Synchronize

implicit NONE

INTERFACE PS_3D_Data_Slice_2D
  MODULE PROCEDURE Data_Slice_2D_0
  MODULE PROCEDURE Data_Slice_2D_1
END INTERFACE

PRIVATE 
PUBLIC  :: PS_3D_Data_Slice_2D

CONTAINS

  Subroutine Data_Slice_2D_0  (Direction, Step, Time, dt,                      &
                               N0_g, N5_g, N0, N5, n_0, M_s_x,                 &
                               M0_g, M5_g, M0, M5, m_0, M_s_y,                 &
                               K0_g, K5_g, K0, K5, k_0, M_s_z,                 &
                               NV, x, y, z, Q, Filename)

  integer  :: Direction, Step
  REALTYPE :: Time, dt

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: K0_g, K5_g, K0, K5, k_0, M_s_z
  integer  :: NV

  Character(LEN=256) :: Filename

  integer            :: i, k, iq1, iq2
  integer            :: nn_0, mm_0, kk_0
  integer            :: i0, i1, j0, j1, l0, l1
  integer            :: Nx, Ny, Nz

  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )    :: x, y, z
  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  ,NV) :: Q
  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )    :: F
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g)    :: Q_g

  integer            :: lid
  Character(LEN=60)  :: F_Slice, File_Ext, Var_Name, Var_Name_f

  Nx = ABS(N5_g-N0_g+1) ; nn_0 = n_0-M_s_x
  Ny = ABS(M5_g-M0_g+1) ; mm_0 = m_0-M_s_y 
  Nz = ABS(K5_g-K0_g+1) ; kk_0 = k_0-M_s_z

  SELECT CASE (Direction)
    CASE (1)
      i0 = N0_g ; i1 = N5_g  ; iq1 = Ny ; iq2 = Nz ; File_Ext = '_x.m'

    CASE (2)
      i0 = M0_g ; i1 = M5_g  ; iq1 = Nx ; iq2 = Ny ; File_Ext = '_y.m'

    CASE (3)
      i0 = K0_g ; i1 = K5_g  ; iq1 = Nx ; iq2 = Ny ; File_Ext = '_z.m'

  END SELECT

  if (I_Am == First_Processor) then
      write (lid6,*) '          (Step, Time, dt  ) : ',Step, Time, dt

    do i = i0,i1
      call PS_Make_String   (Filename, File_Ext, F_Slice, Step, i)

      lid = 200+i

      call PS_OPEN_File (lid, F_Slice, Position='REWIND', Status='UNKNOWN')

      write (lid,*) 'Step = ',Step
      write (lid,*) 'Time = ',Time
      write (lid,*) '  dt = ',dt

      write (lid6,*) 'OPEN File (Slice, lid, Name) : ',i, lid, TRIM(F_Slice)
    enddo
  endif

  call PS_MPI_Synchronize

  do k = -1,NV
    SELECT CASE (Direction)
      CASE (1)
        if (k == -1) then ; F = y                 ; Var_Name = '  y' ; endif
        if (k ==  0) then ; F = z                 ; Var_Name = '  z' ; endif
 
      CASE (2)
        if (k == -1) then ; F = x                 ; Var_Name = '  x' ; endif
        if (k ==  0) then ; F = z                 ; Var_Name = '  z' ; endif

      CASE (3)
        if (k == -1) then ; F = x                 ; Var_Name = '  x' ; endif
        if (k ==  0) then ; F = y                 ; Var_Name = '  y' ; endif

    END SELECT 

    if (k ==  1) then ; F = Q(:,:,:,k)            ; Var_Name = 'Rho' ; endif
    if (k ==  2) then ; F = Q(:,:,:,k)            ; Var_Name = ' RU' ; endif
    if (k ==  3) then ; F = Q(:,:,:,k)            ; Var_Name = ' RV' ; endif
    if (k ==  4) then ; F = Q(:,:,:,k)            ; Var_Name = ' RW' ; endif
    if (k ==  5) then ; F = Q(:,:,:,k)            ; Var_Name = '  E' ; endif
    if (k >=  6) then 
                        F = Q(:,:,:,k)            ; Var_Name_f = ' Rf' 

      call PS_Make_Filename (Var_Name_f, '', Var_Name, k-5)
    endif

    call PS_MPI_Data_ReAssemble (F  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                 Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                 nn_0, mm_0, kk_0)

    if (I_Am == First_Processor)  then
      do i = i0,i1
        lid = 200+i

        write (lid,100) TRIM(Var_Name)

        SELECT CASE (Direction)
          CASE (1)
            write (lid,101) Q_g(i,:,:)

          CASE (2)
            write (lid,101) Q_g(:,i,:)

          CASE (3)
            write (lid,101) Q_g(:,:,i)

        END SELECT

        write (lid,102)
        write (lid,103) TRIM(Var_Name), TRIM(Var_Name), iq1, iq2
      enddo

    endif
 enddo

  if (I_Am == First_Processor) then
    do i = i0,i1
      lid = 200+i ; CLOSE (lid)
    enddo
  endif

  100 format (1x,a16,'=[')
  101 format ((1x,g20.14))
  102 format (1x,' ];')
  103 format (1x,a16,'=reshape(',a16,',',i10,',',i10,');')

  END Subroutine Data_Slice_2D_0  
!
! ------------------------------------------------------------------------
!
  Subroutine Data_Slice_2D_1  (Direction, Step, Time, dt,                      &
                               N0_g, N5_g, N0, N5, n_0, M_s_x,                 &
                               M0_g, M5_g, M0, M5, m_0, M_s_y,                 &
                               K0_g, K5_g, K0, K5, k_0, M_s_z,                 &
                               NV, x, y, z, Q, P, Filename)

  integer  :: Direction, Step
  REALTYPE :: Time, dt

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: K0_g, K5_g, K0, K5, k_0, M_s_z
  integer  :: NV

  Character(LEN=256) :: Filename

  integer            :: i, k, iq1, iq2
  integer            :: nn_0, mm_0, kk_0
  integer            :: i0, i1, j0, j1, l0, l1
  integer            :: Nx, Ny, Nz

  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )    :: x, y, z
  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  ,NV) :: Q
  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )    :: F, P
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g)    :: Q_g

  integer            :: lid
  Character(LEN=60)  :: F_Slice, File_Ext, Var_Name, Var_Name_f

  Nx = ABS(N5_g-N0_g+1) ; nn_0 = n_0-M_s_x
  Ny = ABS(M5_g-M0_g+1) ; mm_0 = m_0-M_s_y 
  Nz = ABS(K5_g-K0_g+1) ; kk_0 = k_0-M_s_z

  SELECT CASE (Direction)
    CASE (1)
      i0 = N0_g ; i1 = N5_g  ; iq1 = Ny ; iq2 = Nz ; File_Ext = '_x.m'

    CASE (2)
      i0 = M0_g ; i1 = M5_g  ; iq1 = Nx ; iq2 = Ny ; File_Ext = '_y.m'

    CASE (3)
      i0 = K0_g ; i1 = K5_g  ; iq1 = Nx ; iq2 = Ny ; File_Ext = '_z.m'

  END SELECT

  if (I_Am == First_Processor) then
      write (lid6,*) '          (Step, Time, dt  ) : ',Step, Time, dt

    do i = i0,i1
      call PS_Make_String   (Filename, File_Ext, F_Slice, Step, i)

      lid = 200+i

      call PS_OPEN_File (lid, F_Slice, Position='REWIND', Status='UNKNOWN')

      write (lid,*) 'Step = ',Step
      write (lid,*) 'Time = ',Time
      write (lid,*) '  dt = ',dt

      write (lid6,*) 'OPEN File (Slice, lid, Name) : ',i, lid, TRIM(F_Slice)
    enddo
  endif

  call PS_MPI_Synchronize

  do k = -1,NV
    SELECT CASE (Direction)
      CASE (1)
        if (k == -1) then ; F = y                 ; Var_Name = '  y' ; endif
        if (k ==  0) then ; F = z                 ; Var_Name = '  z' ; endif
 
      CASE (2)
        if (k == -1) then ; F = x                 ; Var_Name = '  x' ; endif
        if (k ==  0) then ; F = z                 ; Var_Name = '  z' ; endif

      CASE (3)
        if (k == -1) then ; F = x                 ; Var_Name = '  x' ; endif
        if (k ==  0) then ; F = y                 ; Var_Name = '  y' ; endif

    END SELECT 

    if (k ==  1) then ; F = Q(:,:,:,k)            ; Var_Name = 'Rho' ; endif
    if (k ==  2) then ; F = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  U' ; endif
    if (k ==  3) then ; F = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  V' ; endif
    if (k ==  4) then ; F = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  W' ; endif
    if (k ==  5) then ; F = P(:,:,:)              ; Var_Name = '  P' ; endif
    if (k >=  6) then 
                        F = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name_f = '  f' 

      call PS_Make_Filename (Var_Name_f, '', Var_Name, k-5)
    endif

    call PS_MPI_Data_ReAssemble (F  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                 Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                 nn_0, mm_0, kk_0)

    if (I_Am == First_Processor)  then
      do i = i0,i1
        lid = 200+i

        write (lid,100) TRIM(Var_Name)

        SELECT CASE (Direction)
          CASE (1)
            write (lid,101) Q_g(i,:,:)

          CASE (2)
            write (lid,101) Q_g(:,i,:)

          CASE (3)
            write (lid,101) Q_g(:,:,i)

        END SELECT

        write (lid,102)
        write (lid,103) TRIM(Var_Name), TRIM(Var_Name), iq1, iq2
      enddo

    endif
 enddo

  if (I_Am == First_Processor) then
    do i = i0,i1
      lid = 200+i ; CLOSE (lid)
    enddo
  endif

  100 format (1x,a16,'=[')
  101 format ((1x,g20.14))
  102 format (1x,' ];')
  103 format (1x,a16,'=reshape(',a16,',',i10,',',i10,');')

  END Subroutine Data_Slice_2D_1  

END MODULE MPI_3D_Data_Slice_2D





!
! ------------------------------------------------------------------------
!
#if 0
MODULE MPI_3D_Data_Slice_2D

  USE PseudoPack

!  USE Analysis_Function

implicit NONE

PRIVATE 
PUBLIC  : PS_Data_Slice_2D

CONTAINS

  Subroutine PS_Data_Slice_2D (Direction, Step,                                &
                               N0_g, N2_g, N3_g, N5_g, N0, N2, N3, N5,         &
                                 n_0, Order_x, dx, x,                          &
                               M0_g, M2_g, M3_g, M5_g, M0, M2, M3, M5,         &
                                 m_0, Order_y, dy, y,                          &
                               K0_g, K2_g, K3_g, K5_g, K0, K2, K3, K5,         &
                                 k_0, Order_z, dz, z,                          &
                               NV, Gamma, Q, Filename)

  integer  :: Direction, Step
  integer  :: N0_g, N2_g, N3_g, N5_g, N0, N2, N3, N5, n_0, Order_x
  integer  :: M0_g, M2_g, M3_g, M5_g, M0, M2, M3, M5, m_0, Order_y
  integer  :: K0_g, K2_g, K3_g, K5_g, K0, K2, K3, K5, k_0, Order_z
  integer  :: NV

  REALTYPE           :: dx, dy, dz
  Character(LEN=256) :: Filename

  integer            :: i, k, iq1, iq2
  integer            :: nn_0, mm_0, kk_0
  integer            :: i0, i1, j0, j1, l0, l1
  integer            :: M_s_x, Nx, M_s_y, Ny, M_s_z, Nz

  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )    :: x, y, z
  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  ,NV) :: Q
  REALTYPE, dimension(N0  :N5  ,M0  :M5  ,K0  :K5  )    :: F, Gamma, P
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g)    :: Q_g

  integer            :: lid
  Character(LEN=60)  :: F_Slice, Var_Name

  M_s_x = (Order_x+1)/2 ; Nx = ABS(N5_g-N0_g+1) ; nn_0 = n_0-M_s_x
  M_s_y = (Order_y+1)/2 ; Ny = ABS(M5_g-M0_g+1) ; mm_0 = m_0-M_s_y 
  M_s_z = (Order_z+1)/2 ; Nz = ABS(K5_g-K0_g+1) ; kk_0 = k_0-M_s_z

  call Pressure_P                (N0, N5, N2, N3, Order_x, dx,   &
                                  M0, M5, M2, M3, Order_y, dy,   &
                                  K0, K5, K2, K3, Order_z, dz,   &
                                  NV, Q, Gamma, P)

  SELECT CASE (Direction)
    CASE (1)
      i0 = N0_g ; i1 = N5_g  ; iq1 = Ny ; iq2 = Nz

    CASE (2)
      i0 = M0_g ; i1 = M5_g  ; iq1 = Nx ; iq2 = Ny

    CASE (3)
      i0 = K0_g ; i1 = K5_g  ; iq1 = Nx ; iq2 = Ny

  END SELECT

  if (I_Am == First_Processor) then
    do i = i0,i1
      call PS_Make_Filename (Filename, '.m', F_Slice, Step, i)

      lid = 200+i

      call PS_OPEN_File (lid, F_Slice, Position='REWIND', Status='UNKNOWN')

      write (lid6,*) 'OPEN File (Slice, lid, Name) : ',i, lid, TRIM(F_Slice)
    enddo
  endif

  do k = -1,NV
    SELECT CASE (Direction)
      CASE (1)
        if (k == -1) then ; F = y                 ; Var_Name = '  y' ; endif
        if (k ==  0) then ; F = z                 ; Var_Name = '  z' ; endif
 
      CASE (2)
        if (k == -1) then ; F = x                 ; Var_Name = '  x' ; endif
        if (k ==  0) then ; F = z                 ; Var_Name = '  z' ; endif

      CASE (3)
        if (k == -1) then ; F = x                 ; Var_Name = '  x' ; endif
        if (k ==  0) then ; F = y                 ; Var_Name = '  y' ; endif

    END SELECT 

    if (k ==  1) then ; F = Q(:,:,:,k)            ; Var_Name = 'Rho' ; endif
    if (k ==  2) then ; F = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  U' ; endif
    if (k ==  3) then ; F = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  V' ; endif
    if (k ==  4) then ; F = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  W' ; endif
    if (k ==  5) then ; F = P(:,:,:)              ; Var_Name = '  P' ; endif
    if (k ==  6) then ; F = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  f' ; endif

    call PS_MPI_Data_ReAssemble (F  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                 Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                 nn_0, mm_0, kk_0)

    if (I_Am == First_Processor)  then
      do i = i0,i1
        lid = 200+i

        write (lid,100) TRIM(Var_Name)

        SELECT CASE (Direction)
          CASE (1)
            write (lid,101) Q_g(i,:,:)

          CASE (2)
            write (lid,101) Q_g(:,i,:)

          CASE (3)
            write (lid,101) Q_g(:,:,i)

        END SELECT

        write (lid,102)
        write (lid,103) TRIM(Var_Name), TRIM(Var_Name), iq1, iq2
      enddo

    endif
 enddo

  if (I_Am == First_Processor) then
    do i = i0,i1
      lid = 200+i ; CLOSE (lid)
    enddo
  endif

  100 format (1x,a16,'=[')
  101 format ((1x,g20.14))
  102 format (1x,' ];')
  103 format (1x,a16,'=reshape(',a16,',',i10,',',i10,');')

  END Subroutine PS_Data_Slice_2D 

END MODULE MPI_3D_Data_Slice_2D
#endif
