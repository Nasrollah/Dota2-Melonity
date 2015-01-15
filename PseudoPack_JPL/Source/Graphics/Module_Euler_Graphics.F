#include "pseudopack.h"


MODULE Euler_Graphics

  USE Processor
  USE File_OPEN
  USE Filename_Maker
  USE MPI_Data_ReAssemble

implicit NONE

INTERFACE PS_Euler_Graphics
  MODULE PROCEDURE  PS_Euler_Graphics_1D
  MODULE PROCEDURE  PS_Euler_Graphics_2D
  MODULE PROCEDURE  PS_Euler_Graphics_3D
END INTERFACE

INTERFACE PS_Euler_MATLAB_Graphics
  MODULE PROCEDURE  MatLab_1D
  MODULE PROCEDURE  MatLab_2D
  MODULE PROCEDURE  MatLab_3D

  MODULE PROCEDURE  MatLab_1D_P
  MODULE PROCEDURE  MatLab_2D_P
  MODULE PROCEDURE  MatLab_3D_P
END INTERFACE

INTERFACE PS_Euler_TECPLOT_Graphics
  MODULE PROCEDURE TecPlot_1D
  MODULE PROCEDURE TecPlot_2D
  MODULE PROCEDURE TecPlot_3D

  MODULE PROCEDURE TecPlot_1D_P
  MODULE PROCEDURE TecPlot_2D_P
  MODULE PROCEDURE TecPlot_3D_P
END INTERFACE

INTERFACE PS_Euler_PLOT3D_Graphics
  MODULE PROCEDURE  Plot3D_1D
  MODULE PROCEDURE  Plot3D_2D
  MODULE PROCEDURE  Plot3D_3D
 
  MODULE PROCEDURE  Plot3D_1D_P
  MODULE PROCEDURE  Plot3D_2D_P
  MODULE PROCEDURE  Plot3D_3D_P
END INTERFACE

PRIVATE

PUBLIC  :: PS_Euler_Graphics

PUBLIC  :: PS_Euler_MATLAB_Graphics
PUBLIC  :: PS_Euler_TECPLOT_Graphics
PUBLIC  :: PS_Euler_PLOT3D_Graphics

CONTAINS
!
!===================================================================
!
  Subroutine PS_Euler_Graphics_1D (Graphic_Format, Data_Format,           &
                                   N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                                   Nv, x, Q, Step, Time, dt,              &
                                   Filename, TecPlot_Header, TecPlot_Title)

  integer  :: Graphic_Format, Data_Format
  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: Nv, Step
  REALTYPE :: Time, dt

  logical , OPTIONAL  :: TecPlot_Header
  logical             ::       T_Header

  Character(LEN=*), OPTIONAL :: TecPlot_Title
  Character(LEN=256)         ::       T_Title

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5)    :: x
  REALTYPE, dimension(N0:N5,Nv) :: Q

  SELECT CASE (Graphic_Format)
    CASE DEFAULT   ! TECPLOT
                                   T_Header = .TRUE. 
      if (PRESENT(TecPlot_Header)) T_Header = TecPlot_Header

                                   T_Title  = ' '    
      if (PRESENT(TecPlot_Title )) T_Title  = '1D '//TecPlot_Title
  
      SELECT CASE (Data_Format)
        CASE DEFAULT  ! Distributed
          call PS_Euler_TECPLOT_Graphics (N0, N5,                             &
                                          Nv, x, Q, Step, Time, dt,           &
                                          Filename, T_Header, T_Title)

        CASE (1)      ! Single
          call PS_Euler_TECPLOT_Graphics (N0_g, N5_g, N0, N5, n_0, M_s_x,     &
                                          Nv, x, Q, Step, Time, dt,           &
                                          Filename, T_Header, T_Title)

      END SELECT

      if (PRESENT(TecPlot_Header)) TecPlot_Header = T_Header

    CASE (1)       ! MATLAB
      SELECT CASE (Data_Format)
        CASE DEFAULT  ! Distributed
          call PS_Euler_MATLAB_Graphics  (N0, N5,                             &
                                          Nv, x, Q, Step, Time, dt,           &
                                          Filename)

        CASE (1)      ! Single
          call PS_Euler_MATLAB_Graphics  (N0_g, N5_g, N0, N5, n_0, M_s_x,     &
                                          Nv, x, Q, Step, Time, dt,           &
                                          Filename)

      END SELECT

    CASE (2)       ! Plot3D
                                   T_Header = .TRUE.
      if (PRESENT(TecPlot_Header)) T_Header = TecPlot_Header
 
      SELECT CASE (Data_Format)
        CASE DEFAULT  ! Distributed
          call PS_Euler_PLOT3D_Graphics  (N0, N5,                             &
                                          Nv, x, Q, Step, Time, dt,           &
                                          Filename, T_Header) 
 
        CASE (1)      ! Single
          call PS_Euler_PLOT3D_Graphics  (N0_g, N5_g, N0, N5, n_0, M_s_x,     &
                                          Nv, x, Q, Step, Time, dt,           &
                                          Filename, T_Header) 
 
      END SELECT
 
      if (PRESENT(TecPlot_Header)) TecPlot_Header = T_Header

  END SELECT

  END Subroutine PS_Euler_Graphics_1D 
!
!===================================================================
!
  Subroutine PS_Euler_Graphics_2D (Graphic_Format, Data_Format,           &
                                   N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                                   M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                                   Nv, x, y, Q, Step, Time, dt,           &
                                   Filename, TecPlot_Header, TecPlot_Title)

  integer  :: Graphic_Format, Data_Format
  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: Nv, Step
  REALTYPE :: Time, dt

  logical , OPTIONAL  :: TecPlot_Header
  logical             ::       T_Header

  Character(LEN=*), OPTIONAL :: TecPlot_Title
  Character(LEN=256)         ::       T_Title

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5)    :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,Nv) :: Q

  SELECT CASE (Graphic_Format)
    CASE DEFAULT   ! TECPLOT
                                   T_Header = .TRUE. 
      if (PRESENT(TecPlot_Header)) T_Header = TecPlot_Header

                                   T_Title  = ' '    
      if (PRESENT(TecPlot_Title )) T_Title  = '2D '//TecPlot_Title
  
      SELECT CASE (Data_Format)
        CASE DEFAULT  ! Distributed
          call PS_Euler_TECPLOT_Graphics (N0, N5, M0, M5,                     &
                                          Nv, x, y, Q, Step, Time, dt,        &
                                          Filename, T_Header, T_Title)

        CASE (1)      ! Single
          call PS_Euler_TECPLOT_Graphics (N0_g, N5_g, N0, N5, n_0, M_s_x,     &
                                          M0_g, M5_g, M0, M5, m_0, M_s_y,     &
                                          Nv, x, y, Q, Step, Time, dt,        &
                                          Filename, T_Header, T_Title)

      END SELECT

      if (PRESENT(TecPlot_Header)) TecPlot_Header = T_Header

    CASE (1)       ! MATLAB
      SELECT CASE (Data_Format)
        CASE DEFAULT  ! Distributed
          call PS_Euler_MATLAB_Graphics  (N0, N5, M0, M5,                     &
                                          Nv, x, y, Q, Step, Time, dt,        &
                                          Filename)

        CASE (1)      ! Single
          call PS_Euler_MATLAB_Graphics  (N0_g, N5_g, N0, N5, n_0, M_s_x,     &
                                          M0_g, M5_g, M0, M5, m_0, M_s_y,     &
                                          Nv, x, y, Q, Step, Time, dt,        &
                                          Filename)
      END SELECT

    CASE (2)       ! Plot3D
                                   T_Header = .TRUE.
      if (PRESENT(TecPlot_Header)) T_Header = TecPlot_Header
 
      SELECT CASE (Data_Format)
        CASE DEFAULT  ! Distributed
          call PS_Euler_PLOT3D_Graphics  (N0, N5, M0, M5,                     &
                                          Nv, x, y, Q, Step, Time, dt,        &
                                          Filename, T_Header) 
 
        CASE (1)      ! Single
          call PS_Euler_PLOT3D_Graphics  (N0_g, N5_g, N0, N5, n_0, M_s_x,     &
                                          M0_g, M5_g, M0, M5, m_0, M_s_y,     &
                                          Nv, x, y, Q, Step, Time, dt,        &
                                          Filename, T_Header) 
 
      END SELECT
 
      if (PRESENT(TecPlot_Header)) TecPlot_Header = T_Header

  END SELECT

  END Subroutine PS_Euler_Graphics_2D
!
!===================================================================
!
  Subroutine PS_Euler_Graphics_3D (Graphic_Format, Data_Format,           &
                                   N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                                   M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                                   K0_g, K5_g, K0, K5, k_0, M_s_z,        &
                                   Nv, x, y, z, Q, Step, Time, dt,        &
                                   Filename, TecPlot_Header, TecPlot_Title)

  integer  :: Graphic_Format, Data_Format
  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: K0_g, K5_g, K0, K5, k_0, M_s_z
  integer  :: Nv, Step
  REALTYPE :: Time, dt

  logical , OPTIONAL  :: TecPlot_Header
  logical             ::       T_Header

  Character(LEN=*), OPTIONAL :: TecPlot_Title
  Character(LEN=256)         ::       T_Title

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Nv) :: Q

  SELECT CASE (Graphic_Format)
    CASE DEFAULT   ! TECPLOT
                                   T_Header = .TRUE. 
      if (PRESENT(TecPlot_Header)) T_Header = TecPlot_Header

                                   T_Title  = ' '    
      if (PRESENT(TecPlot_Title )) T_Title  = '3D '//TecPlot_Title
  
      SELECT CASE (Data_Format)
        CASE DEFAULT  ! Distributed
          call PS_Euler_TECPLOT_Graphics (N0, N5, M0, M5, K0, K5,             &
                                          Nv, x, y, z, Q, Step, Time, dt,     &
                                          Filename, T_Header, T_Title)

        CASE (1)      ! Single
          call PS_Euler_TECPLOT_Graphics (N0_g, N5_g, N0, N5, n_0, M_s_x,     &
                                          M0_g, M5_g, M0, M5, m_0, M_s_y,     &
                                          K0_g, K5_g, K0, K5, k_0, M_s_z,     &
                                          Nv, x, y, z, Q, Step, Time, dt,     &
                                          Filename, T_Header, T_Title)

      END SELECT

      if (PRESENT(TecPlot_Header)) TecPlot_Header = T_Header

    CASE (1)       ! MATLAB
      SELECT CASE (Data_Format)
        CASE DEFAULT  ! Distributed
          call PS_Euler_MATLAB_Graphics  (N0, N5, M0, M5, K0, K5,             &
                                          Nv, x, y, z, Q, Step, Time, dt,     &
                                          Filename)

        CASE (1)      ! Single
          call PS_Euler_MATLAB_Graphics  (N0_g, N5_g, N0, N5, n_0, M_s_x,     &
                                          M0_g, M5_g, M0, M5, m_0, M_s_y,     &
                                          K0_g, K5_g, K0, K5, k_0, M_s_z,     &
                                          Nv, x, y, z, Q, Step, Time, dt,     &
                                          Filename)

      END SELECT

    CASE (2)       ! Plot3D
                                   T_Header = .TRUE.
      if (PRESENT(TecPlot_Header)) T_Header = TecPlot_Header
 
      SELECT CASE (Data_Format)
        CASE DEFAULT  ! Distributed
          call PS_Euler_PLOT3D_Graphics  (N0, N5, M0, M5, K0, K5,             &
                                          Nv, x, y, z, Q, Step, Time, dt,     &
                                          Filename, T_Header)
 
        CASE (1)      ! Single
          call PS_Euler_PLOT3D_Graphics  (N0_g, N5_g, N0, N5, n_0, M_s_x,     &
                                          M0_g, M5_g, M0, M5, m_0, M_s_y,     &
                                          K0_g, K5_g, K0, K5, k_0, M_s_z,     &
                                          Nv, x, y, z, Q, Step, Time, dt,     &
                                          Filename, T_Header)
 
      END SELECT

      if (PRESENT(TecPlot_Header)) TecPlot_Header = T_Header

  END SELECT

  END Subroutine PS_Euler_Graphics_3D



!
!===================================================================
!
  Subroutine  MatLab_1D (N0_g, N5_g, N0, N5, n_0, M_s_x,  &
                         Nv, x, Q, Step, Time, dt, Filename)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: Nv, Step
  REALTYPE :: Time, dt

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5)    :: x
  REALTYPE, dimension(N0:N5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  Character(LEN=20 ) :: Var_Name
  integer            :: lid, i0, i1, Nx, k
  integer            :: nn_0

  REALTYPE, dimension(N0:N5)           :: P
  REALTYPE, dimension(:), ALLOCATABLE  :: Q_g

  nn_0 = n_0-M_s_x 

  i0 = N0_g ; i1 = N5_g 

  if (I_Am == First_Processor) then
    lid = 101+I_Am

    call PS_Make_Filename (Filename, '.m', F_Tecplot, Step)
  
    call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

    ALLOCATE (Q_g(N0_g:N5_g))

    Nx = ABS(i1-i0+1) 
  endif

  do k = 1,1
    if (k == 1) then ; P = x ; Var_Name = 'x' ; endif

    call PS_MPI_Data_ReAssemble (P  , N0  , N5  ,  &
                                 Q_g, N0_g, N5_g,  &
                                 nn_0)

    if (I_Am == First_Processor)  then
      write (lid,100) TRIM(Var_Name)
      write (lid,101) Q_g(i0:i1)
      write (lid,102) 

      write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx
    endif
  enddo

  do k = 1,Nv
    if (k == 1 ) then ; P = Q(:,k)        ; Var_Name = 'Rho' ; endif
    if (k == 2 ) then ; P = Q(:,k)/Q(:,1) ; Var_Name = '  U' ; endif
    if (k == Nv) then ; P = Q(:,k)        ; Var_Name = '  E' ; endif
      
    call PS_MPI_Data_ReAssemble (P  , N0  , N5  ,  &
                                 Q_g, N0_g, N5_g,  &
                                 nn_0)

    if (I_Am == First_Processor)  then
      write (lid,100) TRIM(Var_Name)
      write (lid,101) Q_g(i0:i1)
      write (lid,102) 

      write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx
    endif
  enddo

  if (I_Am == First_Processor) then
    DeAllocate (Q_g) ; CLOSE (lid)
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (1x,a16,'=[')
  101 format ((1x,g20.14))
  102 format (1x,' ];')
  103 format (1x,a16,'=reshape(',a16,',',i10,');')

  END Subroutine MatLab_1D 
!
!===================================================================
!
  Subroutine  MatLab_2D (N0_g, N5_g, N0, N5, n_0, M_s_x,           &
                         M0_g, M5_g, M0, M5, m_0, M_s_y,           &
                         Nv, x, y, Q, Step, Time, dt, Filename)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: Nv, Step
  REALTYPE :: Time, dt

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5)    :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  Character(LEN=20 ) :: Var_Name
  integer            :: lid, i0, i1, j0, j1, Nx, Ny, k
  integer            :: nn_0, mm_0

  REALTYPE, dimension(N0:N5,M0:M5)       :: P
  REALTYPE, dimension(:,:), ALLOCATABLE  :: Q_g

  nn_0 = n_0-M_s_x ; mm_0 = m_0-M_s_y  

  i0 = N0_g ; i1 = N5_g ; j0 = M0_g ; j1 = M5_g  

  if (I_Am == First_Processor) then
    lid = 101+I_Am

    call PS_Make_Filename (Filename, '.m', F_Tecplot, Step)
  
    call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

    ALLOCATE (Q_g(N0_g:N5_g,M0_g:M5_g))

    Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) 
  endif

  do k = 1,2
    if (k == 1) then ; P = x ; Var_Name = 'x' ; endif
    if (k == 2) then ; P = y ; Var_Name = 'y' ; endif

    call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  ,  &
                                 Q_g, N0_g, N5_g, M0_g, M5_g,  &
                                 nn_0, mm_0)

    if (I_Am == First_Processor)  then
      write (lid,100) TRIM(Var_Name)
      write (lid,101) Q_g(i0:i1,j0:j1)
      write (lid,102) 

      write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx, Ny
    endif
  enddo

  do k = 1,Nv
    if (k == 1 ) then ; P = Q(:,:,k)          ; Var_Name = 'Rho' ; endif
    if (k == 2 ) then ; P = Q(:,:,k)/Q(:,:,1) ; Var_Name = '  U' ; endif
    if (k == 3 ) then ; P = Q(:,:,k)/Q(:,:,1) ; Var_Name = '  V' ; endif
    if (k == Nv) then ; P = Q(:,:,k)          ; Var_Name = '  E' ; endif
      
    call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  ,  &
                                 Q_g, N0_g, N5_g, M0_g, M5_g,  &
                                 nn_0, mm_0)

    if (I_Am == First_Processor)  then
      write (lid,100) TRIM(Var_Name)
      write (lid,101) Q_g(i0:i1,j0:j1)
      write (lid,102) 

      write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx, Ny
    endif
  enddo

  if (I_Am == First_Processor) then
    DeAllocate (Q_g) ; CLOSE (lid)
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (1x,a16,'=[')
  101 format ((1x,g20.14))
  102 format (1x,' ];')
  103 format (1x,a16,'=reshape(',a16,',',i10,',',i10,');')

  END Subroutine MatLab_2D 
!
!===================================================================
!
  Subroutine  MatLab_3D (N0_g, N5_g, N0, N5, n_0, M_s_x,           &
                         M0_g, M5_g, M0, M5, m_0, M_s_y,           &
                         K0_g, K5_g, K0, K5, k_0, M_s_z,           &
                         Nv, x, y, z, Q, Step, Time, dt, Filename)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: K0_g, K5_g, K0, K5, k_0, M_s_z
  integer  :: Nv, Step
  REALTYPE :: Time, dt

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  Character(LEN=20 ) :: Var_Name
  integer            :: lid, i0, i1, j0, j1, l0, l1, Nx, Ny, Nz, k
  integer            :: nn_0, mm_0, kk_0

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   :: P
  REALTYPE, dimension(:,:,:), ALLOCATABLE  :: Q_g

  nn_0 = n_0-M_s_x ; mm_0 = m_0-M_s_y ; kk_0 = k_0-M_s_z 

  i0 = N0_g ; i1 = N5_g ; j0 = M0_g ; j1 = M5_g ; l0 = K0_g ; l1 = K5_g 

  if (I_Am == First_Processor) then
    lid = 101+I_Am

    call PS_Make_Filename (Filename, '.m', F_Tecplot, Step)
  
    call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

    ALLOCATE (Q_g(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g))

    Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) ; Nz = ABS(l1-l0+1)
  endif

  do k = 1,3
    if (k == 1) then ; P = x ; Var_Name = 'x' ; endif
    if (k == 2) then ; P = y ; Var_Name = 'y' ; endif
    if (k == 3) then ; P = z ; Var_Name = 'z' ; endif

    call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                 Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                 nn_0, mm_0, kk_0)

    if (I_Am == First_Processor)  then
      write (lid,100) TRIM(Var_Name)
      write (lid,101) Q_g(i0:i1,j0:j1,l0:l1)
      write (lid,102) 

      write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx, Ny, Nz
    endif
  enddo

  do k = 1,Nv
    if (k == 1 ) then ; P = Q(:,:,:,k)            ; Var_Name = 'Rho' ; endif
    if (k == 2 ) then ; P = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  U' ; endif
    if (k == 3 ) then ; P = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  V' ; endif
    if (k == 4 ) then ; P = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  W' ; endif
    if (k == Nv) then ; P = Q(:,:,:,k)            ; Var_Name = '  E' ; endif
      
    call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                 Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                 nn_0, mm_0, kk_0)

    if (I_Am == First_Processor)  then
      write (lid,100) TRIM(Var_Name)
      write (lid,101) Q_g(i0:i1,j0:j1,l0:l1)
      write (lid,102) 

      write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx, Ny, Nz
    endif
  enddo

  if (I_Am == First_Processor) then
    DeAllocate (Q_g) ; CLOSE (lid)
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (1x,a16,'=[')
  101 format ((1x,g20.14))
  102 format (1x,' ];')
  103 format (1x,a16,'=reshape(',a16,',',i10,',',i10,',',i10,');')

  END Subroutine MatLab_3D 
!
!===================================================================
!
  Subroutine  MatLab_1D_P (N0, N5,                           &
                           Nv, x, Q, Step, Time, dt, Filename)

  integer  :: N0, N5
  integer  :: Nv, Step
  REALTYPE :: Time, dt

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5)    :: x
  REALTYPE, dimension(N0:N5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  Character(LEN=20 ) :: Var_Name
  integer            :: lid, i0, i1, Nx, k

  REALTYPE, dimension(N0:N5)   :: P

  i0 = N0 ; i1 = N5 

  lid = 101+I_Am

  call PS_Make_Filename (Filename, '.m', F_Tecplot, Step, I_Am)

  call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

  Nx = ABS(i1-i0+1) 

  do k = 1,1
    if (k == 1) then ; P = x ; Var_Name = 'x' ; endif

    write (lid,100) TRIM(Var_Name)
    write (lid,101) P(i0:i1)  
    write (lid,102) 
    write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx
  enddo

  do k = 1,Nv
    if (k == 1 ) then ; P = Q(:,k)        ; Var_Name = 'Rho' ; endif
    if (k == 2 ) then ; P = Q(:,k)/Q(:,1) ; Var_Name = '  U' ; endif
    if (k == Nv) then ; P = Q(:,k)        ; Var_Name = '  E' ; endif

    write (lid,100) TRIM(Var_Name)
    write (lid,101) P(i0:i1)
    write (lid,102) 

    write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx
  enddo

  CLOSE (lid)

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (1x,a16,'=[')
  101 format ((1x,g20.14))
  102 format (1x,' ];')
  103 format (1x,a16,'=reshape(',a16,',',i10,');')

  END Subroutine MatLab_1D_P
!
!===================================================================
!
  Subroutine  MatLab_2D_P (N0, N5, M0, M5,                      &
                           Nv, x, y, Q, Step, Time, dt, Filename)

  integer  :: N0, N5, M0, M5
  integer  :: Nv, Step
  REALTYPE :: Time, dt

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5)    :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  Character(LEN=20 ) :: Var_Name
  integer            :: lid, i0, i1, j0, j1, Nx, Ny, k

  REALTYPE, dimension(N0:N5,M0:M5)   :: P

  i0 = N0 ; i1 = N5 ; j0 = M0 ; j1 = M5 

  lid = 101+I_Am

  call PS_Make_Filename (Filename, '.m', F_Tecplot, Step, I_Am)

  call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

  Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) 

  do k = 1,2
    if (k == 1) then ; P = x ; Var_Name = 'x' ; endif
    if (k == 2) then ; P = y ; Var_Name = 'y' ; endif

    write (lid,100) TRIM(Var_Name)
    write (lid,101) P(i0:i1,j0:j1)  
    write (lid,102) 

    write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx, Ny
  enddo

  do k = 1,Nv
    if (k == 1 ) then ; P = Q(:,:,k)          ; Var_Name = 'Rho' ; endif
    if (k == 2 ) then ; P = Q(:,:,k)/Q(:,:,1) ; Var_Name = '  U' ; endif
    if (k == 3 ) then ; P = Q(:,:,k)/Q(:,:,1) ; Var_Name = '  V' ; endif
    if (k == Nv) then ; P = Q(:,:,k)          ; Var_Name = '  E' ; endif

    write (lid,100) TRIM(Var_Name)
    write (lid,101) P(i0:i1,j0:j1)
    write (lid,102) 

    write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx, Ny
  enddo

  CLOSE (lid)

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (1x,a16,'=[')
  101 format ((1x,g20.14))
  102 format (1x,' ];')
  103 format (1x,a16,'=reshape(',a16,',',i10,',',i10,');')

  END Subroutine MatLab_2D_P
!
!===================================================================
!
  Subroutine  MatLab_3D_P (N0, N5, M0, M5, K0, K5,                 &
                           Nv, x, y, z, Q, Step, Time, dt, Filename)

  integer  :: N0, N5, M0, M5, K0, K5
  integer  :: Nv, Step
  REALTYPE :: Time, dt

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  Character(LEN=20 ) :: Var_Name
  integer            :: lid, i0, i1, j0, j1, l0, l1, Nx, Ny, Nz, k

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   :: P

  i0 = N0 ; i1 = N5 ; j0 = M0 ; j1 = M5 ; l0 = K0 ; l1 = K5

  lid = 101+I_Am

  call PS_Make_Filename (Filename, '.m', F_Tecplot, Step, I_Am)

  call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

  Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) ; Nz = ABS(l1-l0+1)

  do k = 1,3
    if (k == 1) then ; P = x ; Var_Name = 'x' ; endif
    if (k == 2) then ; P = y ; Var_Name = 'y' ; endif
    if (k == 3) then ; P = z ; Var_Name = 'z' ; endif

    write (lid,100) TRIM(Var_Name)
    write (lid,101) P(i0:i1,j0:j1,l0:l1)  
    write (lid,102) 

    write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx, Ny, Nz
  enddo

  do k = 1,Nv
    if (k == 1 ) then ; P = Q(:,:,:,k)            ; Var_Name = 'Rho' ; endif
    if (k == 2 ) then ; P = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  U' ; endif
    if (k == 3 ) then ; P = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  V' ; endif
    if (k == 4 ) then ; P = Q(:,:,:,k)/Q(:,:,:,1) ; Var_Name = '  W' ; endif
    if (k == Nv) then ; P = Q(:,:,:,k)            ; Var_Name = '  E' ; endif

    write (lid,100) TRIM(Var_Name)
    write (lid,101) P(i0:i1,j0:j1,l0:l1)
    write (lid,102) 

    write (lid,103) TRIM(Var_Name), TRIM(Var_Name), Nx, Ny, Nz
  enddo

  CLOSE (lid)

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (1x,a16,'=[')
  101 format ((1x,g20.14))
  102 format (1x,' ];')
  103 format (1x,a16,'=reshape(',a16,',',i10,',',i10,',',i10,');')

  END Subroutine MatLab_3D_P
!
!===================================================================
!
  Subroutine TecPlot_1D (N0_g, N5_g, N0, N5, n_0, M_s_x,  &
                         Nv, x, Q, Step, Time, dt,        &
                         Filename, TecPlot_Header, TecPlot_Title)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: TecPlot_Header

  Character(LEN=*) :: TecPlot_Title
  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5)    :: x
  REALTYPE, dimension(N0:N5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  integer            :: lid, i0, i1, Nx, k
  integer            :: nn_0

  REALTYPE, dimension(N0:N5)           :: P
  REALTYPE, dimension(:), ALLOCATABLE  :: Q_g

  nn_0 = n_0-M_s_x 

  i0 = N0_g ; i1 = N5_g 

  if (I_Am == First_Processor) then
    lid = 101+I_Am

    call PS_Make_Filename (Filename, '.dat', F_Tecplot, Step)

    ALLOCATE (Q_g(N0_g:N5_g))

    Nx = ABS(i1-i0+1) 
  endif

  if (TecPlot_Header) then
    TecPlot_Header = .FALSE.

    if (I_Am == First_Processor) then
      call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

      write (lid,996) TRIM(TecPlot_Title), Step, Time, dt
      write (lid,997)
      write (lid,998) I_Am, Step, Time, dt, Nx
      write (lid,200) I_Am, Step, Time, dt, Nx, Nv, &
                      N0_g, N5_g
    endif

    do k = 1,1
      if (k == 1) P = x 

      call PS_MPI_Data_ReAssemble (P  , N0  , N5  ,  &
                                   Q_g, N0_g, N5_g,  &
                                   nn_0)

      if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1)
    enddo
  else
    if (I_Am == First_Processor) then
      call PS_OPEN_File (lid, F_Tecplot, Position='APPEND', Status='UNKNOWN')

      write (lid,999) I_Am, Step, Time, dt, Nx
      write (lid,200) I_Am, Step, Time, dt, Nx, Nv, &
                      N0_g, N5_g
    endif
  endif

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      P = Q(:,k)
    else
      P = Q(:,k)/Q(:,1)
    endif

    call PS_MPI_Data_ReAssemble (P  , N0  , N5  ,  &
                                 Q_g, N0_g, N5_g,  &
                                 nn_0)

    if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1)
  enddo

  if (I_Am == First_Processor) then
    DeAllocate (Q_g) ; CLOSE (lid)
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  996 format (1x,'title="',(a),': Step=',i10, ' Time=',g20.10,' dt=',g20.10,'"')
  997 format (1x,'variables=x, rho,u,E')
  998 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' F=Block')
  999 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,          &
                 ' F=Block, D=(1)')
  200 format (1x,'AUXDATA Processor_ID = "',i5     ,'"'         /   &
              1x,'AUXDATA Step         = "',i5     ,'"'         /   &
              1x,'AUXDATA Time         = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA dt           = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA Nx           = "',i5     ,'"'         /   &
              1x,'AUXDATA Nv           = "',i5     ,'"'         /   &
              1x,'AUXDATA N0           = "',i5     ,'"'         /   &
              1x,'AUXDATA N5           = "',i5     ,'"')

  END Subroutine TecPlot_1D
!
!===================================================================
!
  Subroutine TecPlot_2D (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                         M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                         Nv, x, y, Q, Step, Time, dt,           &
                         Filename, TecPlot_Header, TecPlot_Title)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: TecPlot_Header

  Character(LEN=*) :: TecPlot_Title
  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5)    :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  integer            :: lid, i0, i1, j0, j1, Nx, Ny, k
  integer            :: nn_0, mm_0

  REALTYPE, dimension(N0:N5,M0:M5)       :: P
  REALTYPE, dimension(:,:), ALLOCATABLE  :: Q_g

  nn_0 = n_0-M_s_x ; mm_0 = m_0-M_s_y  

  i0 = N0_g ; i1 = N5_g ; j0 = M0_g ; j1 = M5_g 

  if (I_Am == First_Processor) then
    lid = 101+I_Am

    call PS_Make_Filename (Filename, '.dat', F_Tecplot, Step)

    ALLOCATE (Q_g(N0_g:N5_g,M0_g:M5_g))

    Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) 
  endif

  if (TecPlot_Header) then
    TecPlot_Header = .FALSE.

    if (I_Am == First_Processor) then
      call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

      write (lid,996) TRIM(TecPlot_Title), Step, Time, dt
      write (lid,997)
      write (lid,998) I_Am, Step, Time, dt, Nx, Ny
      write (lid,200) I_Am, Step, Time, dt, Nx, Ny, Nv, &
                      N0_g, N5_g, M0_g, M5_g
    endif

    do k = 1,2
      if (k == 1) P = x ; if (k == 2) P = y 

      call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  ,  &
                                   Q_g, N0_g, N5_g, M0_g, M5_g,  &
                                   nn_0, mm_0)

      if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1,j0:j1)
    enddo
  else
    if (I_Am == First_Processor) then
      call PS_OPEN_File (lid, F_Tecplot, Position='APPEND', Status='UNKNOWN')

      write (lid,999) I_Am, Step, Time, dt, Nx, Ny
      write (lid,200) I_Am, Step, Time, dt, Nx, Ny, Nv, &
                      N0_g, N5_g, M0_g, M5_g
    endif
  endif

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      P = Q(:,:,k)
    else
      P = Q(:,:,k)/Q(:,:,1)
    endif

    call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  ,  &
                                 Q_g, N0_g, N5_g, M0_g, M5_g,  &
                                 nn_0, mm_0)

    if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1,j0:j1)
  enddo

  if (I_Am == First_Processor) then
    DeAllocate (Q_g) ; CLOSE (lid)
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  996 format (1x,'title="',(a),': Step=',i10, ' Time=',g20.10,' dt=',g20.10,'"')
  997 format (1x,'variables=x,y, rho,u,v,E')
  998 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' J=',i5,' F=Block')
  999 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' J=',i5,          &
                 ' F=Block, D=(1,2)')
  200 format (1x,'AUXDATA Processor_ID = "',i5     ,'"'         /   &
              1x,'AUXDATA Step         = "',i5     ,'"'         /   &
              1x,'AUXDATA Time         = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA dt           = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA Nx           = "',i5     ,'"'         /   &
              1x,'AUXDATA Ny           = "',i5     ,'"'         /   &
              1x,'AUXDATA Nv           = "',i5     ,'"'         /   &
              1x,'AUXDATA N0           = "',i5     ,'"'         /   &
              1x,'AUXDATA N5           = "',i5     ,'"'         /   &
              1x,'AUXDATA M0           = "',i5     ,'"'         /   &
              1x,'AUXDATA M5           = "',i5     ,'"')

  END Subroutine TecPlot_2D
!
!===================================================================
!
  Subroutine TecPlot_3D (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                         M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                         K0_g, K5_g, K0, K5, k_0, M_s_z,        &
                         Nv, x, y, z, Q, Step, Time, dt,        &
                         Filename, TecPlot_Header, TecPlot_Title)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: K0_g, K5_g, K0, K5, k_0, M_s_z
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: TecPlot_Header

  Character(LEN=*) :: TecPlot_Title
  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  integer            :: lid, i0, i1, j0, j1, l0, l1, Nx, Ny, Nz, k
  integer            :: nn_0, mm_0, kk_0

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   :: P
  REALTYPE, dimension(:,:,:), ALLOCATABLE  :: Q_g

  nn_0 = n_0-M_s_x ; mm_0 = m_0-M_s_y ; kk_0 = k_0-M_s_z 

  i0 = N0_g ; i1 = N5_g ; j0 = M0_g ; j1 = M5_g ; l0 = K0_g ; l1 = K5_g

  if (I_Am == First_Processor) then
    lid = 101+I_Am

    call PS_Make_Filename (Filename, '.dat', F_Tecplot, Step)

    ALLOCATE (Q_g(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g))

    Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) ; Nz = ABS(l1-l0+1)
  endif

  if (TecPlot_Header) then
    TecPlot_Header = .FALSE.

    if (I_Am == First_Processor) then
      call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

      write (lid,996) TRIM(TecPlot_Title), Step, Time, dt
      write (lid,997)
      write (lid,998) I_Am, Step, Time, dt, Nx, Ny, Nz
      write (lid,200) I_Am, Step, Time, dt, Nx, Ny, Nz, Nv, &
                      N0_g, N5_g, M0_g, M5_g, K0_g, K5_g
    endif

    do k = 1,3
      if (k == 1) P = x ; if (k == 2) P = y ; if (k == 3) P = z

      call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                   Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                   nn_0, mm_0, kk_0)

      if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1,j0:j1,l0:l1)
    enddo
  else
    if (I_Am == First_Processor) then
      call PS_OPEN_File (lid, F_Tecplot, Position='APPEND', Status='UNKNOWN')

      write (lid,999) I_Am, Step, Time, dt, Nx, Ny, Nz
      write (lid,200) I_Am, Step, Time, dt, Nx, Ny, Nz, Nv, &
                      N0_g, N5_g, M0_g, M5_g, K0_g, K5_g
    endif
  endif

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      P = Q(:,:,:,k)
    else
      P = Q(:,:,:,k)/Q(:,:,:,1)
    endif

    call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                 Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                 nn_0, mm_0, kk_0)

    if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1,j0:j1,l0:l1)
  enddo

  if (I_Am == First_Processor) then
    DeAllocate (Q_g) ; CLOSE (lid)
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  996 format (1x,'title="',(a),': Step=',i10, ' Time=',g20.10,' dt=',g20.10,'"')
  997 format (1x,'variables=x,y,z, rho,u,v,w,E')
  998 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' J=',i5,' K=',i5,' F=Block')
  999 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' J=',i5,' K=',i5,          &
                 ' F=Block, D=(1,2,3)')
  200 format (1x,'AUXDATA Processor_ID = "',i5     ,'"'         /   &
              1x,'AUXDATA Step         = "',i5     ,'"'         /   &
              1x,'AUXDATA Time         = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA dt           = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA Nx           = "',i5     ,'"'         /   &
              1x,'AUXDATA Ny           = "',i5     ,'"'         /   &
              1x,'AUXDATA Nz           = "',i5     ,'"'         /   &
              1x,'AUXDATA Nv           = "',i5     ,'"'         /   &
              1x,'AUXDATA N0           = "',i5     ,'"'         /   &
              1x,'AUXDATA N5           = "',i5     ,'"'         /   &
              1x,'AUXDATA M0           = "',i5     ,'"'         /   &
              1x,'AUXDATA M5           = "',i5     ,'"'         /   &
              1x,'AUXDATA K0           = "',i5     ,'"'         /   &
              1x,'AUXDATA K5           = "',i5     ,'"')

  END Subroutine TecPlot_3D
!
!===================================================================
!
  Subroutine TecPlot_1D_P (N0, N5,                                &
                           Nv, x, Q, Step, Time, dt,              &
                           Filename, TecPlot_Header, TecPlot_Title)

  integer  :: N0, N5
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: TecPlot_Header

  Character(LEN=*) :: TecPlot_Title
  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5)    :: x
  REALTYPE, dimension(N0:N5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  integer            :: lid, i0, i1, Nx, k

  i0 = N0 ; i1 = N5 

  lid = 101+I_Am

  call PS_Make_Filename (Filename, '.dat', F_Tecplot, Step, I_Am)

  Nx = ABS(i1-i0+1) 

  if (TecPlot_Header) then
    TecPlot_Header = .FALSE.

    call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

    write (lid,996) TRIM(TecPlot_Title), Step, Time, dt
    write (lid,997)
    write (lid,998) I_Am, Step, Time, dt, Nx
    write (lid,200) I_Am, Step, Time, dt, Nx, &
                    N0, N5

    write (lid,100) x(i0:i1)
  else
    call PS_OPEN_File (lid, F_Tecplot, Position='APPEND', Status='UNKNOWN')

    write (lid,999) I_Am, Step, Time, dt, Nx
    write (lid,200) I_Am, Step, Time, dt, Nx, &
                    N0, N5
  endif

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      write (lid,100) Q(i0:i1,k)
    else
      write (lid,100) Q(i0:i1,k)/Q(i0:i1,1)
    endif
  enddo

  CLOSE (lid)

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  996 format (1x,'title="',(a),': Step=',i10, ' Time=',g20.10,' dt=',g20.10,'"')
  997 format (1x,'variables=x, rho,u,E')
  998 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' F=Block')
  999 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,          &
                 ' F=Block, D=(1)')
  200 format (1x,'AUXDATA Processor_ID = "',i5     ,'"'         /   &
              1x,'AUXDATA Step         = "',i5     ,'"'         /   &
              1x,'AUXDATA Time         = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA dt           = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA Nx           = "',i5     ,'"'         /   &
              1x,'AUXDATA Nv           = "',i5     ,'"'         /   &
              1x,'AUXDATA N0           = "',i5     ,'"'         /   &
              1x,'AUXDATA N5           = "',i5     ,'"')

  END Subroutine TecPlot_1D_P 
!
!===================================================================
!
  Subroutine TecPlot_2D_P (N0, N5, M0, M5,                        &
                           Nv, x, y, Q, Step, Time, dt,           &
                           Filename, TecPlot_Header, TecPlot_Title)

  integer  :: N0, N5, M0, M5
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: TecPlot_Header

  Character(LEN=*) :: TecPlot_Title
  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5)    :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  integer            :: lid, i0, i1, j0, j1, Nx, Ny, k

  i0 = N0 ; i1 = N5 ; j0 = M0 ; j1 = M5 

  lid = 101+I_Am

  call PS_Make_Filename (Filename, '.dat', F_Tecplot, Step, I_Am)

  Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) 

  if (TecPlot_Header) then
    TecPlot_Header = .FALSE.

    call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

    write (lid,996) TRIM(TecPlot_Title), Step, Time, dt
    write (lid,997)
    write (lid,998) I_Am, Step, Time, dt, Nx, Ny
    write (lid,200) I_Am, Step, Time, dt, Nx, Ny, Nv, &
                    N0, N5, M0, M5

    write (lid,100) x(i0:i1,j0:j1)
    write (lid,100) y(i0:i1,j0:j1)
  else
    call PS_OPEN_File (lid, F_Tecplot, Position='APPEND', Status='UNKNOWN')

    write (lid,999) I_Am, Step, Time, dt, Nx, Ny
    write (lid,200) I_Am, Step, Time, dt, Nx, Ny, Nv, &
                    N0, N5, M0, M5
  endif

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      write (lid,100) Q(i0:i1,j0:j1,k)
    else
      write (lid,100) Q(i0:i1,j0:j1,k)/Q(i0:i1,j0:j1,1)
    endif
  enddo

  CLOSE (lid)

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  996 format (1x,'title="',(a),': Step=',i10, ' Time=',g20.10,' dt=',g20.10,'"')
  997 format (1x,'variables=x,y, rho,u,v,E')
  998 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' J=',i5,' F=Block')
  999 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' J=',i5,          &
                 ' F=Block, D=(1,2)')
  200 format (1x,'AUXDATA Processor_ID = "',i5     ,'"'         /   &
              1x,'AUXDATA Step         = "',i5     ,'"'         /   &
              1x,'AUXDATA Time         = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA dt           = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA Nx           = "',i5     ,'"'         /   &
              1x,'AUXDATA Ny           = "',i5     ,'"'         /   &
              1x,'AUXDATA Nv           = "',i5     ,'"'         /   &
              1x,'AUXDATA N0           = "',i5     ,'"'         /   &
              1x,'AUXDATA N5           = "',i5     ,'"'         /   &
              1x,'AUXDATA M0           = "',i5     ,'"'         /   &
              1x,'AUXDATA M5           = "',i5     ,'"')

  END Subroutine TecPlot_2D_P 
!
!===================================================================
!
  Subroutine TecPlot_3D_P (N0, N5, M0, M5, K0, K5,                &
                           Nv, x, y, z, Q, Step, Time, dt,        &
                           Filename, TecPlot_Header, TecPlot_Title)

  integer  :: N0, N5, M0, M5, K0, K5
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: TecPlot_Header

  Character(LEN=*) :: TecPlot_Title
  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Nv) :: Q

  Character(LEN=256) :: F_Tecplot
  integer            :: lid, i0, i1, j0, j1, l0, l1, Nx, Ny, Nz, k

  i0 = N0 ; i1 = N5 ; j0 = M0 ; j1 = M5 ; l0 = K0 ; l1 = K5

  lid = 101+I_Am

  call PS_Make_Filename (Filename, '.dat', F_Tecplot, Step, I_Am)

  Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) ; Nz = ABS(l1-l0+1)

  if (TecPlot_Header) then
    TecPlot_Header = .FALSE.

    call PS_OPEN_File (lid, F_Tecplot, Position='REWIND', Status='UNKNOWN')

    write (lid,996) TRIM(TecPlot_Title), Step, Time, dt
    write (lid,997)
    write (lid,998) I_Am, Step, Time, dt, Nx, Ny, Nz
    write (lid,200) I_Am, Step, Time, dt, Nx, Ny, Nz, Nv, &
                    N0, N5, M0, M5, K0, K5

    write (lid,100) x(i0:i1,j0:j1,l0:l1)
    write (lid,100) y(i0:i1,j0:j1,l0:l1)
    write (lid,100) z(i0:i1,j0:j1,l0:l1)
  else
    call PS_OPEN_File (lid, F_Tecplot, Position='APPEND', Status='UNKNOWN')

    write (lid,999) I_Am, Step, Time, dt, Nx, Ny, Nz
    write (lid,200) I_Am, Step, Time, dt, Nx, Ny, Nz, Nv, &
                    N0, N5, M0, M5, K0, K5
  endif

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      write (lid,100) Q(i0:i1,j0:j1,l0:l1,k)
    else
      write (lid,100) Q(i0:i1,j0:j1,l0:l1,k)/Q(i0:i1,j0:j1,l0:l1,1)
    endif
  enddo

  CLOSE (lid)

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  996 format (1x,'title="',(a),': Step=',i10, ' Time=',g20.10,' dt=',g20.10,'"')
  997 format (1x,'variables=x,y,z, rho,u,v,w,E')
  998 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' J=',i5,' K=',i5,' F=Block')
  999 format (1x,'zone T="',2i5,2g20.10,'" I=',i5,' J=',i5,' K=',i5,          &
                 ' F=Block, D=(1,2,3)')
  200 format (1x,'AUXDATA Processor_ID = "',i5     ,'"'         /   &
              1x,'AUXDATA Step         = "',i5     ,'"'         /   &
              1x,'AUXDATA Time         = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA dt           = "',EN14.4 ,'"'         /   &
              1x,'AUXDATA Nx           = "',i5     ,'"'         /   &
              1x,'AUXDATA Ny           = "',i5     ,'"'         /   &
              1x,'AUXDATA Nz           = "',i5     ,'"'         /   &
              1x,'AUXDATA Nv           = "',i5     ,'"'         /   &
              1x,'AUXDATA N0           = "',i5     ,'"'         /   &
              1x,'AUXDATA N5           = "',i5     ,'"'         /   &
              1x,'AUXDATA M0           = "',i5     ,'"'         /   &
              1x,'AUXDATA M5           = "',i5     ,'"'         /   &
              1x,'AUXDATA K0           = "',i5     ,'"'         /   &
              1x,'AUXDATA K5           = "',i5     ,'"')

  END Subroutine TecPlot_3D_P 



!
!===================================================================
!
  Subroutine  Plot3D_1D (N0_g, N5_g, N0, N5, n_0, M_s_x,  &
                         Nv, x, Q, Step, Time, dt,        &
                         Filename, Plot3D_Grid)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: Plot3D_Grid

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5)    :: x
  REALTYPE, dimension(N0:N5,Nv) :: Q

  Character(LEN=256) :: F_Plot3D_Grid, F_Plot3D_Data
  integer            :: lid, i0, i1, Nx, k
  integer            :: nn_0

  REALTYPE, dimension(N0:N5)           :: P
  REALTYPE, dimension(:), ALLOCATABLE  :: Q_g

  nn_0 = n_0-M_s_x 

  i0 = N0_g ; i1 = N5_g 

  if (I_Am == First_Processor) then
    lid = 101+I_Am

    call PS_Make_Filename (Filename, '.g', F_Plot3D_Grid, Step)
    call PS_Make_Filename (Filename, '.q', F_Plot3D_Data, Step)

    ALLOCATE (Q_g(N0_g:N5_g))

    Nx = ABS(i1-i0+1) 
  endif

  if (Plot3D_Grid) then
    Plot3D_Grid = .FALSE.

    if (I_Am == First_Processor) then
      call PS_OPEN_File (lid, F_Plot3D_Grid, Position='REWIND', &
                                             Status='UNKNOWN')
      write (lid,998) Nx
    endif

    do k = 1,1
      if (k == 1) P = x 

      call PS_MPI_Data_ReAssemble (P  , N0  , N5  ,  &
                                   Q_g, N0_g, N5_g,  &
                                   nn_0)

      if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1)
    enddo

    if (I_Am == First_Processor) CLOSE (lid)
  endif

  if (I_Am == First_Processor) then
    call PS_OPEN_File (lid, F_Plot3D_Data, Position='REWIND', &
                                           Status='UNKNOWN')

    write (lid,998) Nx
    write (lid,999) ZERO, ZERO, ZERO, Time
  endif

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      P = Q(:,k)
    else
      P = Q(:,k)/Q(:,1)
    endif

    call PS_MPI_Data_ReAssemble (P  , N0  , N5  ,  &
                                 Q_g, N0_g, N5_g,  &
                                 nn_0)

    if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1)
  enddo

  if (I_Am == First_Processor) then
    DeAllocate (Q_g) ; CLOSE (lid)
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  998 format (1x,3(i10,:))
  999 format (1x,4g14.7)

  END Subroutine  Plot3D_1D
!
!===================================================================
!
  Subroutine  Plot3D_2D (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                         M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                         Nv, x, y, Q, Step, Time, dt,           &
                         Filename, Plot3D_Grid)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: Plot3D_Grid

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5)    :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,Nv) :: Q

  Character(LEN=256) :: F_Plot3D_Grid, F_Plot3D_Data
  integer            :: lid, i0, i1, j0, j1, Nx, Ny, k
  integer            :: nn_0, mm_0

  REALTYPE, dimension(N0:N5,M0:M5)       :: P
  REALTYPE, dimension(:,:), ALLOCATABLE  :: Q_g

  nn_0 = n_0-M_s_x ; mm_0 = m_0-M_s_y  

  i0 = N0_g ; i1 = N5_g ; j0 = M0_g ; j1 = M5_g 

  if (I_Am == First_Processor) then
    lid = 101+I_Am

    call PS_Make_Filename (Filename, '.g', F_Plot3D_Grid, Step)
    call PS_Make_Filename (Filename, '.q', F_Plot3D_Data, Step)

    ALLOCATE (Q_g(N0_g:N5_g,M0_g:M5_g))

    Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) 
  endif

  if (Plot3D_Grid) then
    Plot3D_Grid = .FALSE.

    if (I_Am == First_Processor) then
      call PS_OPEN_File (lid, F_Plot3D_Grid, Position='REWIND', &
                                             Status='UNKNOWN')

      write (lid,998) Nx, Ny
    endif

    do k = 1,2
      if (k == 1) P = x ; if (k == 2) P = y 

      call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  ,  &
                                   Q_g, N0_g, N5_g, M0_g, M5_g,  &
                                   nn_0, mm_0)

      if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1,j0:j1)
    enddo

    if (I_Am == First_Processor) CLOSE (lid)
  endif

  if (I_Am == First_Processor) then
    call PS_OPEN_File (lid, F_Plot3D_Data, Position='REWIND', &
                                           Status='UNKNOWN')

    write (lid,998) Nx, Ny
    write (lid,999) ZERO, ZERO, ZERO, Time
  endif

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      P = Q(:,:,k)
    else
      P = Q(:,:,k)/Q(:,:,1)
    endif

    call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  ,  &
                                 Q_g, N0_g, N5_g, M0_g, M5_g,  &
                                 nn_0, mm_0)

    if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1,j0:j1)
  enddo

  if (I_Am == First_Processor) then
    DeAllocate (Q_g) ; CLOSE (lid)
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  998 format (1x,3(i10,:))
  999 format (1x,4g14.7)

  END Subroutine  Plot3D_2D
!
!===================================================================
!
  Subroutine  Plot3D_3D (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                         M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                         K0_g, K5_g, K0, K5, k_0, M_s_z,        &
                         Nv, x, y, z, Q, Step, Time, dt,        &
                         Filename, Plot3D_Grid)

  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: K0_g, K5_g, K0, K5, k_0, M_s_z
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: Plot3D_Grid

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Nv) :: Q

  Character(LEN=256) :: F_Plot3D_Grid, F_Plot3D_Data
  integer            :: lid, i0, i1, j0, j1, l0, l1, Nx, Ny, Nz, k
  integer            :: nn_0, mm_0, kk_0

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)   :: P
  REALTYPE, dimension(:,:,:), ALLOCATABLE  :: Q_g

  nn_0 = n_0-M_s_x ; mm_0 = m_0-M_s_y ; kk_0 = k_0-M_s_z 

  i0 = N0_g ; i1 = N5_g ; j0 = M0_g ; j1 = M5_g ; l0 = K0_g ; l1 = K5_g

  if (I_Am == First_Processor) then
    lid = 101+I_Am

    call PS_Make_Filename (Filename, '.g', F_Plot3D_Grid, Step)
    call PS_Make_Filename (Filename, '.q', F_Plot3D_Data, Step)

    ALLOCATE (Q_g(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g))

    Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) ; Nz = ABS(l1-l0+1)
  endif

  if (Plot3D_Grid) then
    Plot3D_Grid = .FALSE.

    if (I_Am == First_Processor) then
      call PS_OPEN_File (lid, F_Plot3D_Grid, Position='REWIND', &
                                             Status='UNKNOWN')

      write (lid,998) Nx, Ny, Nz
    endif

    do k = 1,3
      if (k == 1) P = x ; if (k == 2) P = y ; if (k == 3) P = z

      call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                   Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                   nn_0, mm_0, kk_0)

      if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1,j0:j1,l0:l1)
    enddo

    if (I_Am == First_Processor) CLOSE (lid)
  endif

  if (I_Am == First_Processor) then
    call PS_OPEN_File (lid, F_Plot3D_Data, Position='REWIND', &
                                           Status='UNKNOWN')

    write (lid,998) Nx, Ny, Nz
    write (lid,999) ZERO, ZERO, ZERO, Time
  endif

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      P = Q(:,:,:,k)
    else
      P = Q(:,:,:,k)/Q(:,:,:,1)
    endif

    call PS_MPI_Data_ReAssemble (P  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                 Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                 nn_0, mm_0, kk_0)

    if (I_Am == First_Processor) write (lid,100) Q_g(i0:i1,j0:j1,l0:l1)
  enddo

  if (I_Am == First_Processor) then
    DeAllocate (Q_g) ; CLOSE (lid)
  endif

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  998 format (1x,3(i10,:))
  999 format (1x,4g14.7)

  END Subroutine  Plot3D_3D
!
!===================================================================
!
  Subroutine  Plot3D_1D_P (N0, N5,                                &
                           Nv, x, Q, Step, Time, dt,              &
                           Filename, Plot3D_Grid)

  integer  :: N0, N5
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: Plot3D_Grid

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5)    :: x
  REALTYPE, dimension(N0:N5,Nv) :: Q

  Character(LEN=256) :: F_Plot3D_Grid, F_Plot3D_Data
  integer            :: lid, i0, i1, Nx, k

  i0 = N0 ; i1 = N5 

  lid = 101+I_Am

  call PS_Make_Filename (Filename, '.g', F_Plot3D_Grid, Step, I_Am)
  call PS_Make_Filename (Filename, '.q', F_Plot3D_Data, Step, I_Am)

  Nx = ABS(i1-i0+1) 

  if (Plot3D_Grid) then
    Plot3D_Grid = .FALSE.

      call PS_OPEN_File (lid, F_Plot3D_Grid, Position='REWIND', &
                                             Status='UNKNOWN')

    write (lid,998) Nx

    write (lid,100) x(i0:i1)

    CLOSE (lid)
  endif

    call PS_OPEN_File (lid, F_Plot3D_Data, Position='REWIND', &
                                           Status='UNKNOWN')

    write (lid,998) Nx
    write (lid,999) ZERO, ZERO, ZERO, Time

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      write (lid,100) Q(i0:i1,k)
    else
      write (lid,100) Q(i0:i1,k)/Q(i0:i1,1)
    endif
  enddo

  CLOSE (lid)

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  998 format (1x,3(i10,:))
  999 format (1x,4g14.7)

  END Subroutine  Plot3D_1D_P 
!
!===================================================================
!
  Subroutine  Plot3D_2D_P (N0, N5, M0, M5,                        &
                           Nv, x, y, Q, Step, Time, dt,           &
                           Filename, Plot3D_Grid)

  integer  :: N0, N5, M0, M5
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: Plot3D_Grid

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5)    :: x, y
  REALTYPE, dimension(N0:N5,M0:M5,Nv) :: Q

  Character(LEN=256) :: F_Plot3D_Grid, F_Plot3D_Data
  integer            :: lid, i0, i1, j0, j1, Nx, Ny, k

  i0 = N0 ; i1 = N5 ; j0 = M0 ; j1 = M5 

  lid = 101+I_Am

  call PS_Make_Filename (Filename, '.g', F_Plot3D_Grid, Step, I_Am)
  call PS_Make_Filename (Filename, '.q', F_Plot3D_Data, Step, I_Am)

  Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) 

  if (Plot3D_Grid) then
    Plot3D_Grid = .FALSE.

      call PS_OPEN_File (lid, F_Plot3D_Grid, Position='REWIND', &
                                             Status='UNKNOWN')

    write (lid,998) Nx, Ny

    write (lid,100) x(i0:i1,j0:j1)
    write (lid,100) y(i0:i1,j0:j1)

    CLOSE (lid)
  endif

    call PS_OPEN_File (lid, F_Plot3D_Data, Position='REWIND', &
                                           Status='UNKNOWN')

    write (lid,998) Nx, Ny
    write (lid,999) ZERO, ZERO, ZERO, Time

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      write (lid,100) Q(i0:i1,j0:j1,k)
    else
      write (lid,100) Q(i0:i1,j0:j1,k)/Q(i0:i1,j0:j1,1)
    endif
  enddo

  CLOSE (lid)

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  998 format (1x,3(i10,:))
  999 format (1x,4g14.7)

  END Subroutine  Plot3D_2D_P 
!
!===================================================================
!
  Subroutine  Plot3D_3D_P (N0, N5, M0, M5, K0, K5,                &
                           Nv, x, y, z, Q, Step, Time, dt,        &
                           Filename, Plot3D_Grid)

  integer  :: N0, N5, M0, M5, K0, K5
  integer  :: Nv, Step
  REALTYPE :: Time, dt
  logical  :: Plot3D_Grid

  Character(LEN=*) :: Filename

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: x, y, z
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,Nv) :: Q

  Character(LEN=256) :: F_Plot3D_Grid, F_Plot3D_Data
  integer            :: lid, i0, i1, j0, j1, l0, l1, Nx, Ny, Nz, k

  i0 = N0 ; i1 = N5 ; j0 = M0 ; j1 = M5 ; l0 = K0 ; l1 = K5

  lid = 101+I_Am

  call PS_Make_Filename (Filename, '.g', F_Plot3D_Grid, Step, I_Am)
  call PS_Make_Filename (Filename, '.q', F_Plot3D_Data, Step, I_Am)

  Nx = ABS(i1-i0+1) ; Ny = ABS(j1-j0+1) ; Nz = ABS(l1-l0+1)

  if (Plot3D_Grid) then
    Plot3D_Grid = .FALSE.

      call PS_OPEN_File (lid, F_Plot3D_Grid, Position='REWIND', &
                                             Status='UNKNOWN')

    write (lid,998) Nx, Ny, Nz

    write (lid,100) x(i0:i1,j0:j1,l0:l1)
    write (lid,100) y(i0:i1,j0:j1,l0:l1)
    write (lid,100) z(i0:i1,j0:j1,l0:l1)

    CLOSE (lid)
  endif

    call PS_OPEN_File (lid, F_Plot3D_Data, Position='REWIND', &
                                           Status='UNKNOWN')

    write (lid,998) Nx, Ny, Nz
    write (lid,999) ZERO, ZERO, ZERO, Time

  do k = 1,Nv
    if ((k == 1) .OR. (k == Nv)) then
      write (lid,100) Q(i0:i1,j0:j1,l0:l1,k)
    else
      write (lid,100) Q(i0:i1,j0:j1,l0:l1,k)/Q(i0:i1,j0:j1,l0:l1,1)
    endif
  enddo

  CLOSE (lid)

#if defined (PARALLEL_MPI)
  call MPI_Barrier (MPI_Comm_World, MPI_Error_Status)
#endif

  100 format (8(1x,g14.7))
  998 format (1x,3(i10,:))
  999 format (1x,4g14.7)

  END Subroutine  Plot3D_3D_P 

END MODULE Euler_Graphics
