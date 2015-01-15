.SILENT :
.IGNORE :
.SUFFIX : .f .F .f90 .F90 .a .c .o .mod .h

#=========================================================================
 Silent_Make = -s

MAKE         = gmake -f makefile $(Silent_Make)

#=========================================================================
#  Machine ID : (See MakeFiles/makeMachine)

 machine_id = 5001
 machine_id = 5002

 machine_id = 6001

 MACHINE_ID = $(machine_id)

#=========================================================================
#
MAKE_DIR     = MakeFiles

include $(MAKE_DIR)/makeMachine

#=========================================================================
  HOME_PATH   = ~

ifeq ($(strip $(MACHINE_ID)),5001)
  HOME_PATH    = /export/home/jacobs/HOPE_2009b
endif

ifeq ($(strip $(MACHINE_ID)),5002)
  HOME_PATH    = /home/wsdon/JPL
endif

ifeq ($(strip $(MACHINE_ID)),6001)
  HOME_PATH    = /home5/kbalakr2/jet_code3D
endif

#=========================================================================
#  Code Options

#PARALLEL_OR_SERIAL = Serial
 PARALLEL_OR_SERIAL = Parallel
#PARALLEL_OR_SERIAL = OpenMP

#USE_ESSL_LIBRARY   = on
 USE_ESSL_LIBRARY   = off

 OPTIMIZATION_LEVEL = 2

 COMPILER_OPTION    = on
 DEBUG              = off
 PROFILE            = off
 WARNING            = on

 LOCAL_LIBRARY      = on
  LOCAL_GEMM        = on
  LOCAL_GAMMA       = on
  LOCAL_GAUSS       = on
  LOCAL_VFFT        = on

 ADD_OPTION        +=

#=========================================================================
#
  PS_Tar_Dir        = 
  PS_Tar_Filename   = PseudoPack

#=========================================================================
#               ---- PseudoPack library ----

 PS_ROOT_DIR  = 

 PS_MAIN_PATH = $(HOME_PATH)/PseudoPack_JPL
 PS_VERSION   = dd

 PS_LIB_NAME  = pseudopack
 PS_MOD_NAME  = pseudopack

ifneq ($(strip $(PS_VERSION)),)
 PS_LIBRARY_NAME = $(strip $(PS_LIB_NAME))_$(strip $(PS_VERSION))
 PS_MODULE_NAME  = $(strip $(PS_MOD_NAME))_$(strip $(PS_VERSION))
else
 PS_LIBRARY_NAME = $(strip $(PS_LIB_NAME))
 PS_MODULE_NAME  = $(strip $(PS_MOD_NAME))
endif

PS_SRCDIR  = $(PS_MAIN_PATH)/Source
PS_HDRDIR  = $(PS_MAIN_PATH)/Source/Header
PS_LIBDIR  = $(PS_MAIN_PATH)/Library/$(MACHINE_PATH)/$(PARALLEL_OR_SERIAL)
PS_MODDIR  = $(subst Library,Module,$(PS_LIBDIR))
PS_LSTDIR  = $(subst Library,Listing,$(PS_LIBDIR))
PS_OBJDIR  = $(subst Library,Object,$(PS_LIBDIR))

#=========================================================================
#               ---- Search Path -----

 MOD_PATH   = $(MOD_OPTION) $(PS_MODDIR)
 SRC_PATH   = -I$(PS_SRCDIR) 
    VPATH   = $(PS_SRCDIR) $(PS_SRC_SUBDIR) $(PS_MODDIR) $(PS_OBJDIR)

include $(MAKE_DIR)/Compilers/make$(MACHINE_TYPE)

 SRC_PATH  += -I$(PS_HDRDIR) 
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Euler
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Euler/Euler_Long_Beta
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Euler/Euler_Short_Beta
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Euler/Euler_Long_HAP
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Euler/Euler_Long_JS
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Euler/Euler_Long_Zico
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Euler/Zico_Tau
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Scalar
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/General
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Boundary_Condition
 SRC_PATH  += -I$(PS_SRCDIR)/WENO_Method/Burgers
 SRC_PATH  += -I$(PS_SRCDIR)/Hybrid_Method/Euler
 SRC_PATH  += -I$(PS_SRCDIR)/Hybrid_Method/Scalar

#=========================================================================
#               ---- CPP Pre-Processor ----

 CPP_MACHINE      = $(MACHINE_TYPE)
#CPP_PRECISION    = SINGLE_PRECISION
 CPP_FOURIER      = DEFAULT_FOURIER

ifeq ($(strip $(DEBUG)),on)
 CPP_DEBUG        = DEBUG
endif

 CPP_FFT          = VFFT

ifeq ($(strip $(MACHINE_TYPE)),CRAY)
 CPP_FFT          = CRAY
endif

ifeq ($(strip $(MACHINE_TYPE)),IBM_SP)
 ifeq ($(strip $(USE_ESSL_LIBRARY)),on)
  CPP_FFT          = ESSL
 ifeq ($(strip $(PARALLEL_OR_SERIAL)),Parallel)
  CPP_PESSL        = PESSL
 endif
 endif
endif

ifeq ($(strip $(PARALLEL_OR_SERIAL)),Parallel)
  CPP_PESSL        = SCALAPACK
endif

ifeq ($(strip $(LOCAL_LIBRARY)),on)

 ifeq ($(strip $(LOCAL_GEMM)),on)
  CPP_GEMM         = LOCAL_GEMM
 endif

 ifeq ($(strip $(LOCAL_GAMMA)),on)
  CPP_GAMMA        = LOCAL_GAMMA
 endif

 ifeq ($(strip $(LOCAL_GAUSS)),on)
  CPP_GAUSS        = LOCAL_GAUSS_QUADRATURE
 endif

  CPP_VFFT         = LOCAL_VFFT
 ifeq ($(strip $(MACHINE_TYPE)),CRAY)
  CPP_VFFT         = LOCAL_CRAY_FFT
 endif

endif

ifeq ($(strip $(PARALLEL_OR_SERIAL)),OpenMP)
 CPP_OPENMP       = PARALLEL_OPENMP
endif

ifeq ($(strip $(PARALLEL_OR_SERIAL)),Parallel)
 CPP_MPI          = PARALLEL_MPI
endif

CPP_WENO = 

CPP_DEFINED = $(CPP_MACHINE) $(CPP_PRECISION)                     \
	      $(CPP_FOURIER) $(CPP_DEBUG)      $(CPP_FFT)         \
	      $(CPP_VFFT)    $(CPP_GEMM)       $(CPP_GAMMA)       \
	      $(CPP_GAUSS)   $(CPP_PARALLEL)   $(CPP_MPI)         \
	      $(CPP_OPENMP)  $(CPP_PESSL)                         \
	      $(CPP_WENO)    $(CPP_PGI)        $(CPP_COMPILER)

CPP_OPTION = $(CPP_DEFINED:%=$(CPP)%)

#=========================================================================
#               ---- Fortran Compiler ----

MACHINE_PATH = $(MACHINE_TYPE)/$(MACHINE_OS)/$(MACHINE_COMPILER)
MAKE_PATH    = $(MAKE_DIR)/Compilers/$(MACHINE_PATH)

include $(MAKE_DIR)/makeBase
include $(MAKE_DIR)/PseudoPack_Source

include $(MAKE_PATH)/make$(MACHINE_COMPILER)
include $(MAKE_PATH)/Compiler_Option

#=========================================================================

USER_DEFINED_OPTION =

include $(MAKE_DIR)/makeCompiler_Option

OPTION_FLAG  = $(ADD_OPTION) $(CPP_OPTION) $(SRC_PATH) $(MOD_PATH) 
OPTION_FLAG += $(LIB) $(LIBRARY) $(USER_DEFINED_OPTION)

#=========================================================================
# Update Flags and compiler

  COMPILER_FLAG := $(F_COMPILER_FLAG)
    LOADER_FLAG := $(F_COMPILER_FLAG)

F_COMPILER_FLAG += $(OPTION_FLAG)
C_COMPILER_FLAG += $(OPTION_FLAG)

F_COMPILE = $(FF) $(F_COMPILER_FLAG) 
C_COMPILE = $(CC) $(C_COMPILER_FLAG)

#=========================================================================
# Build targets

TARGET     = $(PS_LIBDIR)/lib$(PS_LIBRARY_NAME).a

LIB_SOURCE = $(LIBSC:%=%.o)

PseudoPack:
	@$(MAKE) System_Info
	@$(MAKE) Directory_Check
	@touch $(PS_SRCDIR)/WENO_Method/Module_PseudoPack*.F
	@$(MAKE) $(TARGET) 
	@-if (test -f *.mod) ; then \
	    $(MOVE) *.mod $(PS_MODDIR)    ;      \
	  fi
	@$(MAKE) clean
	@$(MAKE) System_Done

$(TARGET) : $(TARGET)($(LIB_SOURCE)) $(PS_LIB_NAME).h
	ranlib  $@
	@-if (test -f $(PS_MODDIR)/$(PS_MOD_NAME).mod) ; then \
	    $(COPY) $(PS_MODDIR)/$(PS_MOD_NAME).mod           \
	            $(PS_MODDIR)/$(PS_MODULE_NAME).mod ;      \
	  fi

.F.a:
	$(F_COMPILE) -c $< -o $(PS_OBJDIR)/$%
	$(ARCHIVE) $@ $(PS_OBJDIR)/$%

#=========================================================================
include $(MAKE_DIR)/PseudoPack_Setup
include $(MAKE_DIR)/makeUpdate

#=========================================================================
#  Utilities

.PHONY : clean fresh

clean:
	$(KILL) core* *.o *.mod *.F log ifc* *.h *.lst

fresh:
	$(KILL) core* *.o *.mod *.F log ifc* *.h *.lst
	$(KILL) $(PS_LIBDIR)/* $(PS_MODDIR)/* $(PS_LSTDIR)/* $(PS_OBJDIR)/*

all:
	echo $(MACHINE_ID) $(MACHINE_COMPILER) $(COMPILER_VERSION) $(CPP_PGI)

