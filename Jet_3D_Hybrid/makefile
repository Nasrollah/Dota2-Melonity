.SILENT :
.IGNORE :
.SUFFIX : .f .F .f90 .F90 .a .c .o .mod .h

#=========================================================================
#
 Silent_Make = -s

MAKE        = gmake -f makefile $(Silent_Make)

#=========================================================================
#  Machine Name :


 machine_id = 3001

 machine_id = 5001
 machine_id = 5002
 machine_id = 6001

 MACHINE_ID = $(machine_id)

#=========================================================================
#               ---- MakeFile Directory ----

MAKE_DIR     = MakeFiles

include $(MAKE_DIR)/makeMachine

#=========================================================================
  HOME_PATH   = ~

ifeq ($(strip $(MACHINE_ID)),6001)
  HOME_PATH    = /home5/kbalakr2/jet_code3D
endif

#=========================================================================
#
  Tar_Dir        =
  Tar_Filename   = Program_2009

#=========================================================================
#  Code Options

 TARGET             = jet_$(MACHINE_NAME)
 DIMENSION          = 3
xDIMENSION          = 3

 PROGRAM_PATH       = ../Richtmyer_Meshkov
 CODE_PATH          = $(HOME_PATH)/HOPE/Test
 HOPE_PATH          = $(HOME_PATH)/HOPE/Test

PARALLEL_OR_SERIAL = Parallel
# PARALLEL_OR_SERIAL = Serial
#PARALLEL_OR_SERIAL = OpenMP

#USE_ESSL_LIBRARY   = on
 USE_ESSL_LIBRARY   = off

# OPTIMIZATION_LEVEL = 3
OPTIMIZATION_LEVEL = 2

 COMPILER_OPTION    = on
 DEBUG              = off
#DEBUG              = on
 PROFILE            = off
 WARNING            = off
 VISUS              = off

#=========================================================================
#  Executable Path and Name

 SRCDIR     = Source_WENO/
 SRCDIR     = Source_Hybrid/

 MODDIR     = Module/
 OBJDIR     = Object/
 EXEDIR     = Executable/

 EXECUTABLE = $(EXEDIR)$(TARGET)_$(DIMENSION)d_$(PARALLEL_OR_SERIAL)

#=========================================================================
#               ---- PseudoPack library ----

 PS_MAIN_PATH = $(HOME_PATH)/PseudoPack_JPL
 PS_VERSION   = dd

#=========================================================================
#               ---- CPP Pre-Processor ----

#--------------
ifeq ($(strip $(PARALLEL_OR_SERIAL)),Parallel)
 CPP_MPI       = PARALLEL_MPI
endif

ifeq ($(strip $(PARALLEL_OR_SERIAL)),OpenMP)
 CPP_OPENMP    = PARALLEL_OPENMP
endif

 CPP_DIMENSION =   ONE_DIMENSION
ifeq ($(strip $(DIMENSION)),2)
 CPP_DIMENSION =   TWO_DIMENSION
endif
ifeq ($(strip $(DIMENSION)),3)
 CPP_DIMENSION = THREE_DIMENSION
endif

 CPP_PSEUDOPACK           =  PSEUDOPACK
 CPP_NAVIER_STOKES        =  NAVIER_STOKES
xCPP_LATINI               =  LATINI
 CPP_SPECIES              =  WITH_SPECIES_EQUATION
 CPP_ECHO_FLOW_CONDITION  =  ECHO_FLOW_CONDITION
xCPP_VARIABLE_GAMMA       =  VARIABLE_GAMMA
 CPP_HYBRID               =  HYBRID
 CPP_MULTI_RESOLUTION     =  MULTI_RESOLUTION
 CPP_MR_OUTPUT            =  MR_OUTPUT
 CPP_WENO_OPTION_A        =  SHU_VERSION
xCPP_WENO_OPTION_B        =  LF_COMPONENT_WISE
xCPP_ADAPTIVE_DOMAIN      =  ADAPTIVE_DOMAIN
xCPP_ADAPTIVE_VOLUME      =  ADAPTIVE_VOLUME
xCPP_ANALYSIS             =  ANALYSIS
xCPP_PARTICLE             =  PARTICLE
xCPP_INTERPOLATION        =  ENO_INTERPOLATION

xCPP_TIMING               =  TIMING

xCPP_DEBUG_A              =  DEBUG_HOPE_CASE_STUDY
xCPP_DEBUG_B              =  DEBUG_HOPE_MOLECULAR_DIFFUSION
xCPP_DEBUG_C              =  DEBUG_HOPE_NS_NU_KAPPA_ZETA
xCPP_DEBUG_D              =  DEBUG_HOPE_WENO_FLUX
xCPP_DEBUG_E              =  DEBUG_HOPE_RK
xCPP_DEBUG_F              =  DEBUG_HOPE_TIME_STEP
xCPP_DEBUG_G              =  DEBUG_HOPE_EXTEND_DOMAIN

#--------------

CPP_DEFINED = $(CPP_DIMENSION)        $(CPP_MACHINE)          $(CPP_PSEUDOPACK)\
	      $(CPP_MPI)              $(CPP_OPENMP)                            \
	      $(CPP_HYBRID)                                                    \
	      $(CPP_NAVIER_STOKES)                                             \
	      $(CPP_LATINI)           $(CPP_SPECIES)                           \
	      $(CPP_VARIABLE_GAMMA)                                            \
	      $(CPP_PARTICLE)                                                  \
	      $(CPP_INTERPOLATION)                                             \
	      $(CPP_ANALYSIS)                                                  \
	      $(CPP_MULTI_RESOLUTION) $(CPP_MR_OUTPUT)                         \
	      $(CPP_ADAPTIVE_DOMAIN)  $(CPP_ADAPTIVE_VOLUME)                   \
	      $(CPP_WENO_OPTION_A)    $(CPP_WENO_OPTION_B)                     \
	      $(CPP_PGI)                                                       \
	      $(CPP_ECHO_FLOW_CONDITION)                                       \
	      $(CPP_DEBUG  )                                                   \
	      $(CPP_DEBUG_A) $(CPP_DEBUG_B) $(CPP_DEBUG_C)                     \
	      $(CPP_DEBUG_D) $(CPP_DEBUG_E) $(CPP_DEBUG_F)                     \
	      $(CPP_DEBUG_G)                                                    

CPP_OPTION = $(CPP_DEFINED:%=$(CPP)%)

#=========================================================================
#               ---- MakeFile Directory ----

MACHINE_PATH = $(MACHINE_TYPE)/$(MACHINE_OS)/$(MACHINE_COMPILER)
MAKE_PATH    = $(MAKE_DIR)/Compilers/$(MACHINE_PATH)

include $(MAKE_DIR)/makeBase
include $(MAKE_PATH)/make$(MACHINE_COMPILER)
include $(MAKE_PATH)/Compiler_Option

#=========================================================================
#               ---- PseudoPack Library ----

 PS_LIB_NAME  = pseudopack
 PS_MOD_NAME  = pseudopack

ifneq ($(strip $(PS_VERSION)),)
 PS_LIBRARY_NAME = $(strip $(PS_LIB_NAME))_$(strip $(PS_VERSION))
 PS_MODULE_NAME  = $(strip $(PS_MOD_NAME))_$(strip $(PS_VERSION))
else
 PS_LIBRARY_NAME = $(strip $(PS_LIB_NAME))
 PS_MODULE_NAME  = $(strip $(PS_MOD_NAME))
endif

PS_SRCDIR  = $(PS_MAIN_PATH)/Source/
PS_LIBDIR  = $(PS_MAIN_PATH)/Library/$(MACHINE_PATH)/$(PARALLEL_OR_SERIAL)
PS_MODDIR  = $(subst Library,Module,$(PS_LIBDIR))
PS_INCDIR  = $(PS_MAIN_PATH)/Include/

PS_LIBRARY  = -L $(PS_LIBDIR) -l$(PS_LIBRARY_NAME)

#=========================================================================
# Source Directory

SRCDIR_1 = Initial_Condition
SRCDIR_2 = WENO
SRCDIR_3 = Navier_Stokes
SRCDIR_4 = Hybrid
SRCDIR_5 = Input_Output

SUBDIR  = Header Input_Output MPI Domain Gas Physics Functions 
SUBDIR += Boundary_Condition Multi_Resolution
SUBDIR += Temporal_Solver Spatial_Solver Spectral Central_Difference CPU_Timing
SUBDIR += $(SRCDIR_3)                                                        \
          $(SRCDIR_3)/Navier_Stokes_Boundary_Condition
SUBDIR += $(SRCDIR_1)                                                        \
          $(SRCDIR_1)/Uniform_Flow                                           \
	  $(SRCDIR_1)/Injector
SUBDIR += $(SRCDIR_2)                                                        
SUBDIR += $(SRCDIR_4)                                                        
SUBDIR += $(SRCDIR_5)                                                        \
	  $(SRCDIR_5)/Input_Operation                                        \
	  $(SRCDIR_5)/Output_Operation                                       \
	  $(SRCDIR_5)/Restart_Operation

#SRC_PATH = -I$(SRCDIR) $(SUBDIR:%=-I$(SRCDIR)%) -I$(MODDIR)

SRC_PATH = -I $(SRCDIR) $(SUBDIR:%=-I $(SRCDIR)%) -I $(PS_SRCDIR) -I $(PS_INCDIR)
ifeq ($(strip $(MACHINE_ID)),2000)
  MOD_PATH = $(MOD_OPTION)=$(MODDIR) -I $(PS_MODDIR) -I $(MODDIR)
else
  MOD_PATH = $(MOD_OPTION) $(MODDIR) -I $(PS_MODDIR)
endif

  VPATH = $(SRCDIR) $(PS_SRCDIR)

#=========================================================================

include $(MAKE_DIR)/makeCompiler_Option

USER_DEFINED_OPTION = 

OPTION_FLAG = $(ADD_OPTION) $(CPP_OPTION) $(SRC_PATH) $(MOD_PATH)
OPTION_FLAG += $(LIBRARY) $(USER_DEFINED_OPTION)

#=========================================================================
# Update Flags and compiler

F_COMPILER_FLAG += $(OPTION_FLAG)
F_COMPILE        = $(FF) $(F_COMPILER_FLAG)

LOADER_FLAG     += $(LIBPATH) $(LIBLOC)
LOADER           = $(FF) $(LOADER_FLAG) 

#=========================================================================
# Object files

Object = Jet_3D_DNS

OBJS = $(Object:%=%.o)
SRCS = $(Object:%=%.F)

OBJECTS = $(OBJS:%=$(OBJDIR)%)
SOURCES = $(SRCS:%=$(SRCDIR)%)

#=========================================================================
# Build targets

HOPEE:
	@$(MAKE) System_Info
	@$(MAKE) Directory_check
	@$(MAKE) $(TARGET)

HOPE:
	@$(MAKE) System_Info
	@$(MAKE) Directory_check
	@$(KILL) $(OBJECTS)
	@$(MAKE) $(TARGET)
	@echo
	@echo 
	@$(MAKE) clean
	@$(MAKE) System_Done

$(TARGET): $(OBJECTS) 
	$(LOADER) $(OBJECTS) $(LIBS) -o $(EXECUTABLE) $(PS_LIBRARY) 

$(OBJECTS): $(SRCDIR)$(Object).F
	$(F_COMPILE) -c $< -o $@

#=========================================================================
include $(MAKE_DIR)/makeSetup
include $(MAKE_DIR)/makeCode

#=========================================================================
#  Utilities

.PHONY: clean

clean:
	$(KILL) core* *.o *.mod *.F log ifc* *.h *.lst

fresh:
	$(KILL) core* *.o *.mod log $(EXECUTABLE) $(MODDIR)* $(OBJDIR)*
	$(KILL) ifc* *.h *.lst

#	$(LOADER) $(OBJECTS) $(LIBS) $(VISUSLIBS) $(LIBRARY) -o $(EXECUTABLE)
#$(OBJECTS): $(SRCDIR)WENO.F 
