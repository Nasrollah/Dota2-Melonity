#include "pseudopack.h"


! ------------------------------------------------------------------------
! AUTHOR      : Wai-Sun Don and Bruno Costa  (Rio)
! Description : Uitility subroutine for input data for the Subroutines
!               Setup_Aux_Diff and Setup_Aux_Smth calls.
! 
!  lid     Input file unit (Input)
!
!  For the definition of all other parameters, see Subroutines
!    Setup_Aux_Diff and Setup_Aux_Smth.
!
!  Parameters ended with D are corresponded with those in Setup_Aux_Diff
!  Parameters ended with S are corresponded with those in Setup_Aux_Smth
!
!  The parameters (Index, Method, Algorithm, Symmetry, N, BlockSize)
!  are assumed to be the same for both Setup_Aux_Diff and Setup_Aux_Smth calls.
!
!  Otherwise, used Subroutine PS_Input_Diff and/or PS_Input_Smooth instead.
!
! ------------------------------------------------------------------------

MODULE PseudoPack_Input

  USE PseudoPack_Type

implicit NONE

PRIVATE

INTERFACE PS_Input
! MODULE PROCEDURE PS_Input_N
  MODULE PROCEDURE PS_Input_0
  MODULE PROCEDURE PS_Input_1
  MODULE PROCEDURE PS_Input_Diff
  MODULE PROCEDURE PS_Input_Smooth
END INTERFACE

PUBLIC  :: PS_Property  
PUBLIC  :: PS_Grid_Index
PUBLIC  :: PS_Domain    
PUBLIC  :: PS_Mapping   
PUBLIC  :: PS_Filtering 

PUBLIC  :: PS_Input
PUBLIC  :: PS_Input_D
PUBLIC  :: PS_Input_S

PUBLIC  :: PS_Input_NameList
PUBLIC  :: PS_Output_NameList

CONTAINS

  Subroutine PS_Input_N (lid,                                              &
                         Index, Method, Point_Type, Max_Order_Of_D,        &
                           Algorithm, Symmetry,                            &
                         N, BlockSize,                                     &
                         x0, x1, Map_F, Map, Manual, alpha, beta,          &
                         Filter_Choice_D, Mode_CutOff_D, Omega_D, Order_D, &
                           Smooth_1, Smooth_2,                             &
                         Filter_Choice_S, Mode_CutOff_S, Omega_S, Order_S, &
                         Filter_Choice_R, Mode_CutOff_R, Omega_R, Order_R  )

  integer            :: lid

  integer , OPTIONAL :: Index, Method, Point_Type, Max_Order_Of_D,          &
                          Algorithm, Symmetry,                              &
                        N, BlockSize, Map_F, Map, Manual,                   &
                        Filter_Choice_D, Mode_CutOff_D, Smooth_1, Smooth_2, &
                        Filter_Choice_S, Mode_CutOff_S,                     &
                        Filter_Choice_R, Mode_CutOff_R

  REALTYPE, OPTIONAL :: x0, x1, alpha, beta
  REALTYPE, OPTIONAL :: Omega_D, Order_D
  REALTYPE, OPTIONAL :: Omega_S, Order_S
  REALTYPE, OPTIONAL :: Omega_R, Order_R
 
                                        read (lid,'(//)')
  if (PRESENT(Index           )) then ; read (lid,*) Index 
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Method          )) then ; read (lid,*) Method
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Point_Type      )) then ; read (lid,*) Point_Type
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Max_Order_Of_D  )) then ; read (lid,*) Max_Order_Of_D
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Algorithm       )) then ; read (lid,*) Algorithm
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Symmetry        )) then ; read (lid,*) Symmetry
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(N               )) then ; read (lid,*) N
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(BlockSize       )) then ; read (lid,*) BlockSize
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(x0              )) then ; read (lid,*) x0
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(x1              )) then ; read (lid,*) x1
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Map_F           )) then ; read (lid,*) Map_F
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Map             )) then ; read (lid,*) Map
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Manual          )) then ; read (lid,*) Manual
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(alpha           )) then ; read (lid,*) alpha
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(beta            )) then ; read (lid,*) beta
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(Filter_Choice_D )) then ; read (lid,*) Filter_Choice_D
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Mode_CutOff_D   )) then ; read (lid,*) Mode_CutOff_D
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Omega_D         )) then ; read (lid,*) Omega_D
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Order_D         )) then ; read (lid,*) Order_D
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Smooth_1        )) then ; read (lid,*) Smooth_1
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Smooth_2        )) then ; read (lid,*) Smooth_2
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(Filter_Choice_S )) then ; read (lid,*) Filter_Choice_S
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Mode_CutOff_S   )) then ; read (lid,*) Mode_CutOff_S
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Omega_S         )) then ; read (lid,*) Omega_S
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Order_S         )) then ; read (lid,*) Order_S
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(Filter_Choice_R )) then ; read (lid,*) Filter_Choice_R
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Mode_CutOff_R   )) then ; read (lid,*) Mode_CutOff_R
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Omega_R         )) then ; read (lid,*) Omega_R
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Order_R         )) then ; read (lid,*) Order_R
                                 else ; read (lid,*)                   ; endif

  END Subroutine PS_Input_N
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Input_D (lid,                                              &
                         Index, Method, Point_Type, Max_Order_Of_D,        &
                           Algorithm, Symmetry,                            &
                         N, BlockSize,                                     &
                         x0, x1, Map_F, Map, Manual, alpha, beta,          &
                         Filter_Choice_D, Mode_CutOff_D, Omega_D, Order_D, &
                           Smooth_1, Smooth_2                              )

  integer            :: lid

  integer , OPTIONAL :: Index, Method, Point_Type, Max_Order_Of_D,          &
                          Algorithm, Symmetry,                              &
                        N, BlockSize, Map_F, Map, Manual,                   &
                        Filter_Choice_D, Mode_CutOff_D, Smooth_1, Smooth_2

  REALTYPE, OPTIONAL :: x0, x1, alpha, beta
  REALTYPE, OPTIONAL :: Omega_D, Order_D
 
                                        read (lid,'(//)')
  if (PRESENT(Index           )) then ; read (lid,*) Index 
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Method          )) then ; read (lid,*) Method
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Point_Type      )) then ; read (lid,*) Point_Type
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Max_Order_Of_D  )) then ; read (lid,*) Max_Order_Of_D  
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Algorithm       )) then ; read (lid,*) Algorithm
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Symmetry        )) then ; read (lid,*) Symmetry
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(N               )) then ; read (lid,*) N
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(BlockSize       )) then ; read (lid,*) BlockSize
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(x0              )) then ; read (lid,*) x0
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(x1              )) then ; read (lid,*) x1
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Map_F           )) then ; read (lid,*) Map_F
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Map             )) then ; read (lid,*) Map
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Manual          )) then ; read (lid,*) Manual
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(alpha           )) then ; read (lid,*) alpha
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(beta            )) then ; read (lid,*) beta
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(Filter_Choice_D )) then ; read (lid,*) Filter_Choice_D
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Mode_CutOff_D   )) then ; read (lid,*) Mode_CutOff_D
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Omega_D         )) then ; read (lid,*) Omega_D
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Order_D         )) then ; read (lid,*) Order_D
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Smooth_1        )) then ; read (lid,*) Smooth_1
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Smooth_2        )) then ; read (lid,*) Smooth_2
                                 else ; read (lid,*)                   ; endif

  END Subroutine PS_Input_D
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Input_S (lid,                                              &
                         Index, Method, Point_Type,                        &
                           Algorithm, Symmetry,                            &
                         N, BlockSize,                                     &
                         Filter_Choice_S, Mode_CutOff_S, Omega_S, Order_S  )

  integer            :: lid

  integer , OPTIONAL :: Index, Method, Point_Type,                          &
                          Algorithm, Symmetry,                              &
                        N, BlockSize,                                       &
                        Filter_Choice_S, Mode_CutOff_S

  REALTYPE, OPTIONAL :: Omega_S, Order_S
 
                                        read (lid,'(//)')
  if (PRESENT(Index           )) then ; read (lid,*) Index 
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Method          )) then ; read (lid,*) Method
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Point_Type      )) then ; read (lid,*) Point_Type
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Algorithm       )) then ; read (lid,*) Algorithm
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Symmetry        )) then ; read (lid,*) Symmetry
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(N               )) then ; read (lid,*) N
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(BlockSize       )) then ; read (lid,*) BlockSize
                                 else ; read (lid,*)                   ; endif

                                        read (lid,'(//)')
  if (PRESENT(Filter_Choice_S )) then ; read (lid,*) Filter_Choice_S
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Mode_CutOff_S   )) then ; read (lid,*) Mode_CutOff_S
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Omega_S         )) then ; read (lid,*) Omega_S
                                 else ; read (lid,*)                   ; endif
  if (PRESENT(Order_S         )) then ; read (lid,*) Order_S
                                 else ; read (lid,*)                   ; endif

  END Subroutine PS_Input_S
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Input_0 (Index, Method, Point_Type, Max_Order_Of_D,        &
                           Algorithm, Symmetry,                            &
                         N, BlockSize,                                     &
                         x0, x1,  Map_F, Map, Manual, alpha, beta,         &
                         Filter_Choice_D, Mode_CutOff_D, Omega_D, Order_D, &
                           Smooth_1_D, Smooth_2_D,                         &
                         Filter_Choice_S, Mode_CutOff_S, Omega_S, Order_S, &
                         Filter_Choice_R, Mode_CutOff_R, Omega_R, Order_R, &
                         lid)

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry, &
              N, BlockSize, Map, Manual, Map_F,                               &
              Filter_Choice_D, Mode_CutOff_D, Smooth_1_D, Smooth_2_D,         &
              Filter_Choice_S, Mode_CutOff_S,                                 &
              Filter_Choice_R, Mode_CutOff_R

  REALTYPE :: x0, x1, alpha, beta
  REALTYPE :: Omega_D, Order_D
  REALTYPE :: Omega_S, Order_S
  REALTYPE :: Omega_R, Order_R
 
  integer  :: lid

  read (lid,555)
  read (lid,*) Index
  read (lid,*) Method
  read (lid,*) Point_Type
  read (lid,*) Max_Order_Of_D
  read (lid,*) Algorithm
  read (lid,*) Symmetry

  read (lid,555)
  read (lid,*) N
  read (lid,*) BlockSize

  read (lid,555)
  read (lid,*) x0
  read (lid,*) x1
  read (lid,*) Map_F
  read (lid,*) Map
  read (lid,*) Manual
  read (lid,*) alpha
  read (lid,*) beta

  read (lid,555)
  read (lid,*) Filter_Choice_D
  read (lid,*) Mode_CutOff_D
  read (lid,*) Omega_D
  read (lid,*) Order_D
  read (lid,*) Smooth_1_D
  read (lid,*) Smooth_2_D

  read (lid,555)
  read (lid,*) Filter_Choice_S
  read (lid,*) Mode_CutOff_S
  read (lid,*) Omega_S
  read (lid,*) Order_S

  read (lid,555)
  read (lid,*) Filter_Choice_R
  read (lid,*) Mode_CutOff_R
  read (lid,*) Omega_R
  read (lid,*) Order_R

 555  format (//)

  END Subroutine PS_Input_0
!
! -----------------------------------------------------------------------
!
  Subroutine PS_Input_1 (Index, Method, Point_Type, Max_Order_Of_D,        &
                           Algorithm, Symmetry,                            &
                         N, BlockSize,                                     &
                         x0, x1,  Map_F, Map, Manual, alpha, beta,         &
                         Filter_Choice_D, Mode_CutOff_D, Omega_D, Order_D, &
                           Smooth_1_D, Smooth_2_D,                         &
                         Filter_Choice_S, Mode_CutOff_S, Omega_S, Order_S, &
                         lid)

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry, &
              N, BlockSize, Map, Manual, Map_F,                               &
              Filter_Choice_D, Mode_CutOff_D, Smooth_1_D, Smooth_2_D,         &
              Filter_Choice_S, Mode_CutOff_S

  REALTYPE :: x0, x1, alpha, beta
  REALTYPE :: Omega_D, Order_D
  REALTYPE :: Omega_S, Order_S
 
  integer  :: lid

  read (lid,555)
  read (lid,*) Index
  read (lid,*) Method
  read (lid,*) Point_Type
  read (lid,*) Max_Order_Of_D
  read (lid,*) Algorithm
  read (lid,*) Symmetry

  read (lid,555)
  read (lid,*) N
  read (lid,*) BlockSize

  read (lid,555)
  read (lid,*) x0
  read (lid,*) x1
  read (lid,*) Map_F
  read (lid,*) Map
  read (lid,*) Manual
  read (lid,*) alpha
  read (lid,*) beta

  read (lid,555)
  read (lid,*) Filter_Choice_D
  read (lid,*) Mode_CutOff_D
  read (lid,*) Omega_D
  read (lid,*) Order_D
  read (lid,*) Smooth_1_D
  read (lid,*) Smooth_2_D

  read (lid,555)
  read (lid,*) Filter_Choice_S
  read (lid,*) Mode_CutOff_S
  read (lid,*) Omega_S
  read (lid,*) Order_S

 555  format (//)

  END Subroutine PS_Input_1

! ------------------------------------------------------------------------
! Author  : Wai Sun Don
! AUTHOR  : Wai-Sun Don and Bruno Costa  (Rio)
! Description : Uitility subroutine for input data for the Setup_Aux_Diff ONLY.
! 
!  lid     Input file unit (Input)
!
!  Parameters ended with D are corresponded with those in Setup_Aux_Diff.
!
!  For the definition of all other parameters, see Subroutine Setup_Aux_Diff.
!
! ------------------------------------------------------------------------

  Subroutine PS_Input_Diff (Index, Method, Point_Type, Max_Order_Of_D,        &
                              Algorithm, Symmetry,                            &
                            N, BlockSize,                                     &
                            x0, x1,  Map_F, Map, Manual, alpha, beta,         &
                            Filter_Choice_D, Mode_CutOff_D, Omega_D, Order_D, &
                              Smooth_1_D, Smooth_2_D,                         &
                            lid)

  integer  :: Index, Method, Point_Type, Max_Order_Of_D, Algorithm, Symmetry, &
              N, BlockSize, Map, Manual, Map_F,                               &
              Filter_Choice_D, Mode_CutOff_D, Smooth_1_D, Smooth_2_D

  REALTYPE :: x0, x1, alpha, beta
  REALTYPE :: Omega_D, Order_D

  integer  :: lid

  read (lid,555)
  read (lid,*) Index
  read (lid,*) Method
  read (lid,*) Point_Type
  read (lid,*) Algorithm
  read (lid,*) Symmetry

  read (lid,555)
  read (lid,*) N
  read (lid,*) BlockSize

  read (lid,555)
  read (lid,*) x0
  read (lid,*) x1
  read (lid,*) Map_F
  read (lid,*) Map
  read (lid,*) Manual
  read (lid,*) alpha
  read (lid,*) beta

  read (lid,555)
  read (lid,*) Filter_Choice_D
  read (lid,*) Mode_CutOff_D
  read (lid,*) Omega_D
  read (lid,*) Order_D
  read (lid,*) Smooth_1_D
  read (lid,*) Smooth_2_D

 555  format (//)

  END Subroutine PS_Input_Diff 
! ------------------------------------------------------------------------
! Author  : Wai Sun Don
! AUTHOR  : Wai-Sun Don and Bruno Costa  (Rio)
! Description : Uitility subroutine for input data for the Setup_Aux_Smth ONLY.
! 
!  lid     Input file unit (Input)
!
!  Parameters ended with S are corresponded with those in Setup_Aux_Smth
!
!  For the definition of all other parameters, see Subroutine Setup_Aux_Smth.
!
! ------------------------------------------------------------------------

  Subroutine PS_Input_Smooth (Index, Method, Point_Type,                       &
                                Algorithm, Symmetry,                           &
                              N, BlockSize,                                    &
                              Filter_Choice_S, Mode_CutOff_S, Omega_S, Order_S,&
                              lid)

  integer  :: Index, Method, Point_Type, Algorithm, Symmetry, &
              N, BlockSize,                                   &
              Filter_Choice_S, Mode_CutOff_S

  REALTYPE :: Omega_S, Order_S

  integer  :: lid

  read (lid,555)
  read (lid,*) Index
  read (lid,*) Method
  read (lid,*) Point_Type
  read (lid,*) Algorithm
  read (lid,*) Symmetry

  read (lid,555)
  read (lid,*) N
  read (lid,*) BlockSize

  read (lid,555)
  read (lid,*) Filter_Choice_S
  read (lid,*) Mode_CutOff_S
  read (lid,*) Omega_S
  read (lid,*) Order_S

 555  format (//)

  END Subroutine PS_Input_Smooth 
!
! ------------------------------------------------------------------------
!
  Subroutine PS_Input_NameList (lid, Property, Grid_Index, Domain, Mapping, &
                                     Filtering_D, Filtering_S, Filtering_R  )
                                     

  integer  :: lid

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Grid_Index), OPTIONAL :: Grid_Index
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping
  TYPE (PS_Filtering) , OPTIONAL :: Filtering_D
  TYPE (PS_Filtering) , OPTIONAL :: Filtering_S
  TYPE (PS_Filtering) , OPTIONAL :: Filtering_R

  NAMELIST /List_Property    / Property
  NAMELIST /List_Grid_Index  / Grid_Index
  NAMELIST /List_Domain      / Domain
  NAMELIST /List_Mapping     / Mapping
  NAMELIST /List_Filtering_D / Filtering_D
  NAMELIST /List_Filtering_S / Filtering_S
  NAMELIST /List_Filtering_R / Filtering_R

  if (PRESENT(Property)   ) Property    = Default_Property 
  if (PRESENT(Grid_Index) ) Grid_Index  = Default_Grid_Index 
  if (PRESENT(Domain)     ) Domain      = Default_Domain     
  if (PRESENT(Mapping)    ) Mapping     = Default_Mapping
  if (PRESENT(Filtering_D)) Filtering_D = Default_Filtering  
  if (PRESENT(Filtering_S)) Filtering_S = Default_Filtering  
  if (PRESENT(Filtering_R)) Filtering_R = Default_Filtering  

  if (PRESENT(Property    )) read (lid, NML=List_Property    )
  if (PRESENT(Grid_Index  )) read (lid, NML=List_Grid_Index  )
  if (PRESENT(Domain      )) read (lid, NML=List_Domain      )
  if (PRESENT(Mapping     )) read (lid, NML=List_Mapping     )
  if (PRESENT(Filtering_D )) read (lid, NML=List_Filtering_D )
  if (PRESENT(Filtering_S )) read (lid, NML=List_Filtering_S )
  if (PRESENT(Filtering_R )) read (lid, NML=List_Filtering_R )

  END Subroutine PS_Input_NameList
!
! ------------------------------------------------------------------------
!
  Subroutine PS_Output_NameList (lid, Property, Grid_Index, Domain, Mapping, &
                                      Filtering_D, Filtering_S, Filtering_R  )

  integer  :: lid

  TYPE (PS_Property)  , OPTIONAL :: Property
  TYPE (PS_Grid_Index), OPTIONAL :: Grid_Index
  TYPE (PS_Domain)    , OPTIONAL :: Domain
  TYPE (PS_Mapping)   , OPTIONAL :: Mapping
  TYPE (PS_Filtering) , OPTIONAL :: Filtering_D
  TYPE (PS_Filtering) , OPTIONAL :: Filtering_S
  TYPE (PS_Filtering) , OPTIONAL :: Filtering_R

  NAMELIST /List_Property    / Property
  NAMELIST /List_Grid_Index  / Grid_Index
  NAMELIST /List_Domain      / Domain
  NAMELIST /List_Mapping     / Mapping
  NAMELIST /List_Filtering_D / Filtering_D
  NAMELIST /List_Filtering_S / Filtering_S
  NAMELIST /List_Filtering_R / Filtering_R

  if (PRESENT(Property    )) write (lid, NML=List_Property    )
  if (PRESENT(Grid_Index  )) write (lid, NML=List_Grid_Index  )
  if (PRESENT(Domain      )) write (lid, NML=List_Domain      )
  if (PRESENT(Mapping     )) write (lid, NML=List_Mapping     )
  if (PRESENT(Filtering_D )) write (lid, NML=List_Filtering_D )
  if (PRESENT(Filtering_S )) write (lid, NML=List_Filtering_S )
  if (PRESENT(Filtering_R )) write (lid, NML=List_Filtering_R )

  END Subroutine PS_Output_NameList

END MODULE PseudoPack_Input
