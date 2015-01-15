#include "pseudopack.h"


MODULE Read_Write_Restart_File

  USE Processor
  USE PS_IO_Unit, ONLY: lid6, lid101
  USE MPI_Abort_Routine
  USE MPI_Data_ReAssemble
  USE MPI_Data_DisAssemble
  USE File_OPEN
  USE Filename_Maker

implicit NONE

INTERFACE PS_Read_Restart_File
  MODULE PROCEDURE  Read_Restart_File_1D_A
  MODULE PROCEDURE  Read_Restart_File_2D_A
  MODULE PROCEDURE  Read_Restart_File_3D_A

  MODULE PROCEDURE  Read_Restart_File_1D_B
  MODULE PROCEDURE  Read_Restart_File_2D_B
  MODULE PROCEDURE  Read_Restart_File_3D_B
END INTERFACE

INTERFACE PS_Write_Restart_File
  MODULE PROCEDURE Write_Restart_File_1D_A
  MODULE PROCEDURE Write_Restart_File_2D_A
  MODULE PROCEDURE Write_Restart_File_3D_A

  MODULE PROCEDURE Write_Restart_File_1D_B
  MODULE PROCEDURE Write_Restart_File_2D_B
  MODULE PROCEDURE Write_Restart_File_3D_B
END INTERFACE

PRIVATE

PUBLIC  :: PS_Read_Restart_File
PUBLIC  :: PS_Write_Restart_File

CONTAINS
!
!===================================================================
!
  Subroutine  Read_Restart_File_1D_A (Restart, Step, Time, dt, Q, Filename)

  integer  :: Restart, Step
  REALTYPE :: Time, dt

  REALTYPE, dimension(:,:) :: Q

  Character(LEN=*) :: Filename

  integer :: Nx, NV, Mx, MV, M_Processor, lid, n

  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'OLD'

  lid = lid101+I_Am

  call PS_Make_Filename (Filename, '.dat', F_Restart, Restart, I_Am)

  call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                          Position=F_POSITION, Status=F_STATUS)

  Nx = SIZE(Q,DIM=1) 
  NV = SIZE(Q,DIM=2)

  read  (lid) Step, Time, dt
  read  (lid) Mx, MV, M_Processor

  call IO_Check_A (Restart, Step, Nx,  1,  1, NV, N_Processor, &
                                  Mx,  1,  1, MV, M_Processor)

  do n = 1,NV
    read  (lid) Q(:,n)
  enddo

  CLOSE (lid)

  if (I_Am == First_Processor) write (lid6,100) Step, Time, dt

 100 format (1x/1x,72('*')/                                                &
             1x,'Restarting from Step=',i10,' Time= ',g10.3,' dt= ',g10.3/ &
             1x/1x,72('*')/)

  End Subroutine Read_Restart_File_1D_A
!
!===================================================================
!
  Subroutine  Read_Restart_File_2D_A (Restart, Step, Time, dt, Q, Filename)

  integer  :: Restart, Step
  REALTYPE :: Time, dt

  REALTYPE, dimension(:,:,:) :: Q

  Character(LEN=*) :: Filename

  integer :: Nx, Ny, NV, Mx, My, MV, M_Processor, lid, n

  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'OLD'

  lid = lid101+I_Am

  call PS_Make_Filename (Filename, '.dat', F_Restart, Restart, I_Am)

  call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                          Position=F_POSITION, Status=F_STATUS)

  Nx = SIZE(Q,DIM=1) ; Ny = SIZE(Q,DIM=2) 
  NV = SIZE(Q,DIM=3)

  read  (lid) Step, Time, dt
  read  (lid) Mx, My, MV, M_Processor

  call IO_Check_A (Restart, Step, Nx, Ny,  1, NV, N_Processor, &
                                  Mx, My,  1, MV, M_Processor)

  do n = 1,NV
    read  (lid) Q(:,:,n)
  enddo

  CLOSE (lid)

  if (I_Am == First_Processor) write (lid6,100) Step, Time, dt

 100 format (1x/1x,72('*')/                                                &
             1x,'Restarting from Step=',i10,' Time= ',g10.3,' dt= ',g10.3/ &
             1x/1x,72('*')/)

  End Subroutine Read_Restart_File_2D_A
!
!===================================================================
!
  Subroutine  Read_Restart_File_3D_A (Restart, Step, Time, dt, Q, Filename)

  integer  :: Restart, Step
  REALTYPE :: Time, dt

  REALTYPE, dimension(:,:,:,:) :: Q

  Character(LEN=*) :: Filename

  integer :: Nx, Ny, Nz, NV, Mx, My, Mz, MV, M_Processor, lid, n

  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'OLD'

  lid = lid101+I_Am

  call PS_Make_Filename (Filename, '.dat', F_Restart, Restart, I_Am)

  call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                          Position=F_POSITION, Status=F_STATUS)

  Nx = SIZE(Q,DIM=1) ; Ny = SIZE(Q,DIM=2) ; Nz = SIZE(Q,DIM=3)  
  NV = SIZE(Q,DIM=4)

  read  (lid) Step, Time, dt
  read  (lid) Mx, My, Mz, MV, M_Processor

  call IO_Check_A (Restart, Step, Nx, Ny, Nz, NV, N_Processor, &
                                  Mx, My, Mz, MV, M_Processor)

  do n = 1,NV
    read  (lid) Q(:,:,:,n)
  enddo

  CLOSE (lid)

  if (I_Am == First_Processor) write (lid6,100) Step, Time, dt

 100 format (1x/1x,72('*')/                                                &
             1x,'Restarting from Step=',i10,' Time= ',g10.3,' dt= ',g10.3/ &
             1x/1x,72('*')/)

  End Subroutine Read_Restart_File_3D_A
!
!===================================================================
!
  Subroutine Write_Restart_File_1D_A (Restart, Step, Time, dt, Q, Filename)

  integer           :: Restart, Step
  REALTYPE          :: Time, dt

  REALTYPE, dimension(:,:) :: Q

  Character(LEN=*) :: Filename

  integer :: Nx, NV, lid, n

  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'REPLACE'

  lid = lid101+I_Am

  call PS_Make_Filename (Filename, '.dat', F_Restart, Step, I_Am)

  call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                          Position=F_POSITION, Status=F_STATUS)

  Nx = SIZE(Q,DIM=1) 
  NV = SIZE(Q,DIM=2)

  write (lid) Step, Time, dt
  write (lid) Nx, NV, N_Processor

  do n = 1,NV
    write (lid) Q(:,n)
  enddo

  CLOSE (lid)

  END Subroutine Write_Restart_File_1D_A
!
!===================================================================
!
  Subroutine Write_Restart_File_2D_A (Restart, Step, Time, dt, Q, Filename)

  integer           :: Restart, Step
  REALTYPE          :: Time, dt

  REALTYPE, dimension(:,:,:) :: Q

  Character(LEN=*) :: Filename

  integer :: Nx, Ny, NV, lid, n

  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'REPLACE'

  lid = lid101+I_Am

  call PS_Make_Filename (Filename, '.dat', F_Restart, Step, I_Am)

  call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                          Position=F_POSITION, Status=F_STATUS)

  Nx = SIZE(Q,DIM=1) ; Ny = SIZE(Q,DIM=2) 
  NV = SIZE(Q,DIM=3)

  write (lid) Step, Time, dt
  write (lid) Nx, Ny, NV, N_Processor

  do n = 1,NV
    write (lid) Q(:,:,n)
  enddo

  CLOSE (lid)

  END Subroutine Write_Restart_File_2D_A
!
!===================================================================
!
  Subroutine Write_Restart_File_3D_A (Restart, Step, Time, dt, Q, Filename)

  integer           :: Restart, Step
  REALTYPE          :: Time, dt

  REALTYPE, dimension(:,:,:,:) :: Q

  Character(LEN=*) :: Filename

  integer :: Nx, Ny, Nz, NV, lid, n

  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'REPLACE'

  lid = lid101+I_Am

  call PS_Make_Filename (Filename, '.dat', F_Restart, Step, I_Am)

  call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                          Position=F_POSITION, Status=F_STATUS)

  Nx = SIZE(Q,DIM=1) ; Ny = SIZE(Q,DIM=2) ; Nz = SIZE(Q,DIM=3)  
  NV = SIZE(Q,DIM=4)

  write (lid) Step, Time, dt
  write (lid) Nx, Ny, Nz, NV, N_Processor

  do n = 1,NV
    write (lid) Q(:,:,:,n)
  enddo

  CLOSE (lid)

  END Subroutine Write_Restart_File_3D_A
!
!===================================================================
!
  Subroutine Read_Restart_File_1D_B (N0_g, N5_g, N0, N5, n_0, M_s_x, NV,    &
                                     Restart, Step, Time, dt, Q, Filename)
      
  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: NV, Restart, Step             
  REALTYPE :: Time, dt
    
  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: P
    
  REALTYPE, dimension(N0_g:N5_g) :: Q_g
    
  integer :: Nx, Mx, MV, M_Processor, lid, n
  logical :: Error

  Character(LEN=*)  :: Filename
  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'OLD'
                                 
  if (I_Am == First_Processor) then 
    lid = 101

    call PS_Make_Filename (Filename, '.dat', F_Restart, Restart)

    call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                            Position=F_POSITION, Status=F_STATUS)

    read  (lid) Step, Time, dt
    read  (lid) Mx, MV, M_Processor

    Nx = SIZE(Q_g,DIM=1) 

    call IO_Check_B (Restart, Step, Nx,  1,  1, NV, N_Processor, &
                                    Mx,  1,  1, MV, N_Processor, Error)
  endif          

#if defined (PARALLEL_MPI)
  call BCAST (Step, Time, dt, Error)
#endif

  do n = 1,NV
    if (I_Am == First_Processor) read  (lid) Q_g
 
    call PS_MPI_Data_DisAssemble (Q_g, N0_g, N5_g,  &
                                  P  , N0  , N5  ,  &
                                  n_0-M_s_x)

    Q(:,n) = P
  enddo

  if (I_Am == First_Processor) then
     CLOSE (lid) ; write (lid6,100) Step, Time, dt
  endif

 100 format (1x/1x,72('*')/                                                &
             1x,'Restarting from Step=',i10,' Time= ',g10.3,' dt= ',g10.3/ &
             1x/1x,72('*')/)

  End Subroutine Read_Restart_File_1D_B
!
!===================================================================
!
  Subroutine Write_Restart_File_1D_B (N0_g, N5_g, N0, N5, n_0, M_s_x, NV,    &
                                      Restart, Step, Time, dt, Q, Filename)
      
  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: NV, Restart, Step             
  REALTYPE :: Time, dt
    
  REALTYPE, dimension(N0:N5,NV) :: Q
  REALTYPE, dimension(N0:N5)    :: P
    
  REALTYPE, dimension(N0_g:N5_g) :: Q_g
    
  integer :: Nx, lid, n

  Character(LEN=*)  :: Filename
  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'REPLACE'
                                 
  if (I_Am == First_Processor) then 
    lid = 101

    call PS_Make_Filename (Filename, '.dat', F_Restart, Step)

    call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                            Position=F_POSITION, Status=F_STATUS)

    Nx = SIZE(Q_g,DIM=1) 

    write (lid) Step, Time, dt
    write (lid) Nx, NV, N_Processor
  endif          

  do n = 1,NV
    P = Q(:,n)

    call PS_MPI_Data_ReAssemble  (P  , N0  , N5  ,  &
                                  Q_g, N0_g, N5_g,  &
                                  n_0-M_s_x)

    if (I_Am == First_Processor) write (lid) Q_g
  enddo

  if (I_Am == First_Processor) CLOSE (lid)

  END Subroutine Write_Restart_File_1D_B
!
!===================================================================
!
  Subroutine Read_Restart_File_2D_B (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                                     M0_g, M5_g, M0, M5, m_0, M_s_y, NV,    &
                                     Restart, Step, Time, dt, Q, Filename)
      
  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: NV, Restart, Step             
  REALTYPE :: Time, dt
    
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: P
    
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g) :: Q_g
    
  integer :: Nx, Ny, Mx, My, MV, M_Processor, lid, n
  logical :: Error

  Character(LEN=*)  :: Filename
  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'OLD'
                                 
  if (I_Am == First_Processor) then 
    lid = 101

    call PS_Make_Filename (Filename, '.dat', F_Restart, Restart)

    call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                            Position=F_POSITION, Status=F_STATUS)

    read  (lid) Step, Time, dt
    read  (lid) Mx, My, MV, M_Processor

    Nx = SIZE(Q_g,DIM=1) ; Ny = SIZE(Q_g,DIM=2) 

    call IO_Check_B (Restart, Step, Nx, Ny,  1, NV, N_Processor, &
                                    Mx, My,  1, MV, N_Processor, Error)
  endif          

#if defined (PARALLEL_MPI)
  call BCAST (Step, Time, dt, Error)
#endif

  do n = 1,NV
    if (I_Am == First_Processor) read  (lid) Q_g
 
    call PS_MPI_Data_DisAssemble (Q_g, N0_g, N5_g, M0_g, M5_g,  &
                                  P  , N0  , N5  , M0  , M5  ,  &
                                  n_0-M_s_x, m_0-M_s_y)

    Q(:,:,n) = P
  enddo

  if (I_Am == First_Processor) then
     CLOSE (lid) ; write (lid6,100) Step, Time, dt
  endif

 100 format (1x/1x,72('*')/                                                &
             1x,'Restarting from Step=',i10,' Time= ',g10.3,' dt= ',g10.3/ &
             1x/1x,72('*')/)

  End Subroutine Read_Restart_File_2D_B
!
!===================================================================
!
  Subroutine Write_Restart_File_2D_B (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                                      M0_g, M5_g, M0, M5, m_0, M_s_y, NV,    &
                                      Restart, Step, Time, dt, Q, Filename)
      
  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: NV, Restart, Step             
  REALTYPE :: Time, dt
    
  REALTYPE, dimension(N0:N5,M0:M5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5)    :: P
    
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g) :: Q_g
    
  integer :: Nx, Ny, lid, n

  Character(LEN=*)  :: Filename
  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'REPLACE'
                                 
  if (I_Am == First_Processor) then 
    lid = 101

    call PS_Make_Filename (Filename, '.dat', F_Restart, Step)

    call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                            Position=F_POSITION, Status=F_STATUS)

    Nx = SIZE(Q_g,DIM=1) ; Ny = SIZE(Q_g,DIM=2) 

    write (lid) Step, Time, dt
    write (lid) Nx, Ny, NV, N_Processor
  endif          

  do n = 1,NV
    P = Q(:,:,n)

    call PS_MPI_Data_ReAssemble  (P  , N0  , N5  , M0  , M5  ,  &
                                  Q_g, N0_g, N5_g, M0_g, M5_g,  &
                                  n_0-M_s_x, m_0-M_s_y)

    if (I_Am == First_Processor) write (lid) Q_g
  enddo

  if (I_Am == First_Processor) CLOSE (lid)

  END Subroutine Write_Restart_File_2D_B
!
!===================================================================
!
  Subroutine Read_Restart_File_3D_B (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                                     M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                                     K0_g, K5_g, K0, K5, k_0, M_s_z, NV,    &
                                     Restart, Step, Time, dt, Q, Filename)
      
  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: K0_g, K5_g, K0, K5, k_0, M_s_z
  integer  :: NV, Restart, Step             
  REALTYPE :: Time, dt
    
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: P
    
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g) :: Q_g
    
  integer :: Nx, Ny, Nz, Mx, My, Mz, MV, M_Processor, lid, n
  logical :: Error

  Character(LEN=*)  :: Filename
  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'OLD'
                                 
  if (I_Am == First_Processor) then 
    lid = 101

    call PS_Make_Filename (Filename, '.dat', F_Restart, Restart)

    call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                            Position=F_POSITION, Status=F_STATUS)

    read  (lid) Step, Time, dt
    read  (lid) Mx, My, Mz, MV, M_Processor

    Nx = SIZE(Q_g,DIM=1) ; Ny = SIZE(Q_g,DIM=2) ; Nz = SIZE(Q_g,DIM=3)  

    call IO_Check_B (Restart, Step, Nx, Ny, Nz, NV, N_Processor, &
                                    Mx, My, Mz, MV, N_Processor, Error)
  endif          

#if defined (PARALLEL_MPI)
  call BCAST (Step, Time, dt, Error)
#endif

  do n = 1,NV
    if (I_Am == First_Processor) read  (lid) Q_g
 
    call PS_MPI_Data_DisAssemble (Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                  P  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                  n_0-M_s_x, m_0-M_s_y, k_0-M_s_z)

    Q(:,:,:,n) = P
  enddo

  if (I_Am == First_Processor) then
     CLOSE (lid) ; write (lid6,100) Step, Time, dt
  endif

 100 format (1x/1x,72('*')/                                                &
             1x,'Restarting from Step=',i10,' Time= ',g10.3,' dt= ',g10.3/ &
             1x/1x,72('*')/)

  End Subroutine Read_Restart_File_3D_B
!
!===================================================================
!
  Subroutine Write_Restart_File_3D_B (N0_g, N5_g, N0, N5, n_0, M_s_x,        &
                                      M0_g, M5_g, M0, M5, m_0, M_s_y,        &
                                      K0_g, K5_g, K0, K5, k_0, M_s_z, NV,    &
                                      Restart, Step, Time, dt, Q, Filename)
      
  integer  :: N0_g, N5_g, N0, N5, n_0, M_s_x
  integer  :: M0_g, M5_g, M0, M5, m_0, M_s_y
  integer  :: K0_g, K5_g, K0, K5, k_0, M_s_z
  integer  :: NV, Restart, Step             
  REALTYPE :: Time, dt
    
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: P
    
  REALTYPE, dimension(N0_g:N5_g,M0_g:M5_g,K0_g:K5_g) :: Q_g
    
  integer :: Nx, Ny, Nz, lid, n

  Character(LEN=*)  :: Filename
  character(LEN=80) :: F_Restart
  character(LEN=20) :: F_POSITION = 'REWIND' , F_STATUS = 'REPLACE'
                                 
  if (I_Am == First_Processor) then 
    lid = 101

    call PS_Make_Filename (Filename, '.dat', F_Restart, Step)

    call PS_OPEN_File (lid, F_Restart, Form='UNFORMATTED',  &
                            Position=F_POSITION, Status=F_STATUS)

     Nx = SIZE(Q_g,DIM=1) ; Ny = SIZE(Q_g,DIM=2) ; Nz = SIZE(Q_g,DIM=3)  

    write (lid) Step, Time, dt
    write (lid) Nx, Ny, Nz, NV, N_Processor
  endif

  do n = 1,NV
    P = Q(:,:,:,n)

    call PS_MPI_Data_ReAssemble  (P  , N0  , N5  , M0  , M5  , K0  , K5  , &
                                  Q_g, N0_g, N5_g, M0_g, M5_g, K0_g, K5_g, &
                                  n_0-M_s_x, m_0-M_s_y, k_0-M_s_z)

    if (I_Am == First_Processor) write (lid) Q_g
  enddo

  if (I_Am == First_Processor) CLOSE (lid)

  END Subroutine Write_Restart_File_3D_B
!
!===================================================================
!
  Subroutine IO_Check_A (Restart, Step, Nx, Ny, Nz, NV, NP, Mx, My, Mz, MV, MP)

  integer :: Restart, Step
  logical :: Error_n, Error
  integer :: Nx, Ny, Nz, NV, NP, Mx, My, Mz, MV, MP

  Error_n = (MV /= NV) .OR. (Mx /= Nx) .OR. (My /= Ny) .OR. (Mz /= Nz) 
  Error_n = Error_n    .OR. (MP /= NP)
  Error_n = Error_n    .OR. (Restart /= Step)

#if defined (PARALLEL_MPI)
  call MPI_ALLREDUCE (Error_n, Error, 1,                                     &
                      MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, MPI_Error_Status)
#else
  Error = Error_n
#endif

  if (Error) then
    if (I_Am == First_Processor) then
      write (lid6,*) 'STOP!  Mismatch parameter(s) in Read_Restart'
      write (lid6,*) '       Requested Step, Nx, Ny, Nz, NV, N_Processor = ', &
                                    Restart, Nx, Ny, Nz, NV, NP
      write (lid6,*) '       Stored    Step, Nx, Ny, Nz, NV, N_Processor = ', &
                                       Step, Mx, My, Mz, MV, MP
    endif

    call PS_STOP ('Read_Restart_File', ' ')
  endif

  END Subroutine IO_Check_A
!
!===================================================================
!
  Subroutine IO_Check_B (Restart, Step, Nx, Ny, Nz, NV, NP,        &
                                        Mx, My, Mz, MV, MP, Error)

  integer :: Restart, Step
  logical :: Error_n, Error
  integer :: Nx, Ny, Nz, NV, NP, Mx, My, Mz, MV, MP

  Error_n = (MV /= NV) .OR. (Mx /= Nx) .OR. (My /= Ny) .OR. (Mz /= Nz) 
  Error_n = Error_n    .OR. (MP /= NP)
  Error_n = Error_n    .OR. (Restart /= Step)

  Error = Error_n

  if (Error) then
      write (lid6,*) 'STOP!  Mismatch parameter(s) in Read_Restart'
      write (lid6,*) '       Requested Step, Nx, Ny, Nz, NV, N_Processor = ', &
                                    Restart, Nx, Ny, Nz, NV, NP
      write (lid6,*) '       Stored    Step, Nx, Ny, Nz, NV, N_Processor = ', &
                                       Step, Mx, My, Mz, MV, MP
  endif

  END Subroutine IO_Check_B
#if defined (PARALLEL_MPI)
!
!===================================================================
!
  Subroutine BCAST (Step, Time, dt, Error)

  logical  :: Error
  integer  :: Step
  REALTYPE :: Time, dt
  
  call MPI_BCAST (Error, 1, MPI_LOGICAL , First_Processor,     &
                            MPI_Comm_World, MPI_Error_Status)

  call MPI_BCAST (Step , 1, MPI_INTEGER , First_Processor,     &
                            MPI_Comm_World, MPI_Error_Status)

  call MPI_BCAST (Time , 1, MPI_REALTYPE, First_Processor,     &
                            MPI_Comm_World, MPI_Error_Status)

  call MPI_BCAST (dt   , 1, MPI_REALTYPE, First_Processor,     &
                            MPI_Comm_World, MPI_Error_Status)

  if (Error) call PS_STOP ('Read_Restart_File', ' ')

  END Subroutine BCAST 
#endif

END MODULE Read_Write_Restart_File
