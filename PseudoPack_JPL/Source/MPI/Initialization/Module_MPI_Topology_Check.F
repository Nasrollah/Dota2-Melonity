#include "pseudopack.h"


MODULE MPI_Topology_Check

  USE PseudoPack

IMPLICIT NONE

INTERFACE Topology_Check_Ghost
  MODULE PROCEDURE Topology_Check_1D_0
  MODULE PROCEDURE Topology_Check_2D_0
  MODULE PROCEDURE Topology_Check_3D_0
END INTERFACE

INTERFACE Topology_Check_Processor
  MODULE PROCEDURE Topology_Check_1D_1
  MODULE PROCEDURE Topology_Check_2D_1
  MODULE PROCEDURE Topology_Check_3D_1
END INTERFACE

INTERFACE Topology_Check_N
  MODULE PROCEDURE Topology_Check_1D_2
  MODULE PROCEDURE Topology_Check_2D_2
  MODULE PROCEDURE Topology_Check_3D_2
END INTERFACE

PRIVATE
PUBLIC  :: Topology_Check_Ghost
PUBLIC  :: Topology_Check_Processor
PUBLIC  :: Topology_Check_N

CONTAINS
! 
! ================================================================
!
  Subroutine Topology_Check_1D_0 (N2, N3, M_s_x)

  integer  :: N2, N3, M_s_x
 
  if (N3-N2+1 < M_s_x) then
    write (6,*) 'I_Am, N2, N3, M_s_x :', I_Am, N2, N3, M_s_x

    call PS_STOP ('Topology', ' ','1st Dim : N3-N2+1 < M_s_x')
  endif

  END Subroutine Topology_Check_1D_0 
! 
! ================================================================
!
  Subroutine Topology_Check_2D_0 (N2, N3, M_s_x,  &
                                  M2, M3, M_s_y   )

  integer  :: N2, N3, M_s_x
  integer  :: M2, M3, M_s_y
 
  if (N3-N2+1 < M_s_x) then
    write (6,*) 'I_Am, N2, N3, M_s_x :', I_Am, N2, N3, M_s_x

    call PS_STOP ('Topology', ' ','1st Dim : N3-N2+1 < M_s_x')
  endif

  if (M3-M2+1 < M_s_y) then
    write (6,*) 'I_Am, M2, M3, M_s_y :', I_Am, M2, M3, M_s_y

    call PS_STOP ('Topology', ' ','2st Dim : M3-M2+1 < M_s_y')
  endif

  END Subroutine Topology_Check_2D_0
! 
! ================================================================
!
  Subroutine Topology_Check_3D_0 (N2, N3, M_s_x,  &
                                  M2, M3, M_s_y,  &
                                  K2, K3, M_s_z   )

  integer  :: N2, N3, M_s_x
  integer  :: M2, M3, M_s_y
  integer  :: K2, K3, M_s_z
 
  if (N3-N2+1 <  M_s_x) then
    write (6,*) 'I_Am, N2, N3, M_s_x :', I_Am, N2, N3, M_s_x

    call PS_STOP ('Topology', ' ','1st Dim : N3-N2+1 < M_s_x')
  endif

  if (M3-M2+1 <  M_s_y) then
    write (6,*) 'I_Am, M2, M3, M_s_y :', I_Am, M2, M3, M_s_y

    call PS_STOP ('Topology', ' ','2st Dim : M3-M2+1 < M_s_y')
  endif

  if (K3-K2+1 <  M_s_z) then
    write (6,*) 'I_Am, N2, N3, M_s_z :', I_Am, N2, N3, M_s_z

    call PS_STOP ('Topology', ' ','3st Dim : K3-K2+1 < M_s_z')
  endif

  END Subroutine Topology_Check_3D_0
! 
! ================================================================
!
  Subroutine Topology_Check_1D_1 (N, N_Processor, M_s_x)

  integer  :: N, N_Processor, M_s_x

  if (N/N_Processor < M_s_x) then
    write (6,*) 'I_Am, N, N_Processor, M_s_x :', I_Am, N, N_Processor, M_s_x

    call PS_STOP ('Topology', ' ','1st Dim : N/N_Processor < M_s_x')
  endif

  END Subroutine Topology_Check_1D_1
!
! ================================================================
!
  Subroutine Topology_Check_2D_1 (N, N_Processor, M_s_x, &
                                M, M_Processor, M_s_y)

  integer  :: N, N_Processor, M_s_x
  integer  :: M, M_Processor, M_s_y

  if (N/N_Processor < M_s_x) then
    write (6,*) 'I_Am, N, N_Processor, M_s_x :', I_Am, N, N_Processor, M_s_x

    call PS_STOP ('Topology', ' ','1st Dim : N/N_Processor < M_s_x')
  endif

  if (M/M_Processor < M_s_y) then
    write (6,*) 'I_Am, M, M_Processor, M_s_y :', I_Am, M, M_Processor, M_s_y

    call PS_STOP ('Topology', ' ','2nd Dim : M/M_Processor < M_s_y')
  endif

  END Subroutine Topology_Check_2D_1
!
! ================================================================
!
  Subroutine Topology_Check_3D_1 (N, N_Processor, M_s_x, &
                                  M, M_Processor, M_s_y, &
                                  K, K_Processor, M_s_z)

  integer  :: N, N_Processor, M_s_x
  integer  :: M, M_Processor, M_s_y
  integer  :: K, K_Processor, M_s_z

  if (N/N_Processor < M_s_x) then
    write (6,*) 'I_Am, N, N_Processor, M_s_x :', I_Am, N, N_Processor, M_s_x

    call PS_STOP ('Topology', ' ','1st Dim : N/N_Processor < M_s_x')
  endif

  if (M/M_Processor < M_s_y) then
    write (6,*) 'I_Am, M, M_Processor, M_s_y :', I_Am, M, M_Processor, M_s_y

    call PS_STOP ('Topology', ' ','2nd Dim : M/M_Processor < M_s_y')
  endif

  if (K/K_Processor < M_s_z) then
    write (6,*) 'I_Am, K, K_Processor, M_s_z :', I_Am, K, K_Processor, M_s_z

    call PS_STOP ('Topology', ' ','3rd Dim : K/K_Processor < M_s_z')
  endif

  END Subroutine Topology_Check_3D_1

! Added 2-18-2014 for JPL_Jet_3D program (might not needed!)
#if 1 

! 
! ================================================================
!
  Subroutine Topology_Check_1D_2 (N, N_Processor)

  integer  :: N, N_Processor
 
  if (N < N_Processor) then
    write (6,*) 'I_Am, N_Processor :', I_Am, N_Processor

    call PS_STOP ('Topology', ' ','1st Dim : N < N_Processor')
  endif

  END Subroutine Topology_Check_1D_2 
! 
! ================================================================
!
  Subroutine Topology_Check_2D_2 (N, N_Processor, &
                                  M, M_Processor  )

  integer  :: N, N_Processor
  integer  :: M, M_Processor
 
  if (N < N_Processor) then
    write (6,*) 'I_Am, N_Processor :', I_Am, N_Processor

    call PS_STOP ('Topology', ' ','1st Dim : N < N_Processor')
  endif

  if (M < M_Processor) then
    write (6,*) 'I_Am, M_Processor :', I_Am, M_Processor

    call PS_STOP ('Topology', ' ','2nd Dim : M < M_Processor')
  endif

  END Subroutine Topology_Check_2D_2
! 
! ================================================================
!
  Subroutine Topology_Check_3D_2 (N, N_Processor, &
                                  M, M_Processor, &
                                  K, K_Processor  )

  integer  :: N, N_Processor
  integer  :: M, M_Processor
  integer  :: K, K_Processor
 
  if (N < N_Processor) then
    write (6,*) 'I_Am, N_Processor :', I_Am, N_Processor

    call PS_STOP ('Topology', ' ','1st Dim : N < N_Processor')
  endif

  if (M < M_Processor) then
    write (6,*) 'I_Am, M_Processor :', I_Am, M_Processor

    call PS_STOP ('Topology', ' ','2nd Dim : M < M_Processor')
  endif

  if (K < K_Processor) then
    write (6,*) 'I_Am, K_Processor :', I_Am, K_Processor

    call PS_STOP ('Topology', ' ','3rd Dim : K < K_Processor')
  endif

  END Subroutine Topology_Check_3D_2

#endif


END MODULE MPI_Topology_Check
