# name of the library -- make sure it is in the object directory
LIBNAME = libthornado.a
LIBDIR = ../../../object
INCDIR = ../../../include

# what the makefile should do -- this is the default thing to do
.PHONY: default
default: $(LIBNAME)

# Import site specific information.
# At setup time, this file will exist in the local library directory.
include ../../../Makefile.h

# these should get defined by invoking the Makefile from FLASH
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB
MICROPHYSICS = WEAKLIB

# Import thornado Makefiles
THORNADO_DIR = ../../
include $(THORNADO_DIR)/Build/Makefile_GPU
include $(THORNADO_DIR)/Build/Makefile_Path
include $(THORNADO_DIR)/Build/Makefile_Defines
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

# BUILDFLAG can be passed to the Makefile to choose specific
# compilation flags (one of OPT,TEST,DEBUG).
ifdef BUILDFLAG
  FFLAGS = $(FFLAGS_$(BUILDFLAG))
else
  FFLAGS = $(FFLAGS_OPT)
endif
FFLAGS += $(FFLAGS_HDF5)
FFLAGS += $(FFLAGS_LAPACK)
ifeq ($(USE_GPU),TRUE)
  FFLAGS += $(FFLAGS_CUDA)
endif
ifeq ($(USE_OACC),TRUE)
  FFLAGS += $(FFLAGS_OACC)
endif
ifeq ($(USE_OMP_OL),TRUE)
  FFLAGS += $(FFLAGS_OMP_OL)
endif

# command for generating a library
MKLIB = ar -r

OBJS = \
       KindModule.o \
       ProgramHeaderModule.o \
       ArrayUtilitiesModule.o \
       SubcellReconstructionModule.o \
       TimersModule.o \
       QuadratureModule.o \
       LagrangePolynomialsModule.o \
       PolynomialBasisModuleX_Lagrange.o \
       LinearAlgebraModule.o \
       DeviceModule.o \
       ReferenceElementModuleE.o \
       ReferenceElementModuleX.o \
       ReferenceElementModule.o \
       PolynomialBasisModule_Lagrange.o \
       ReferenceElementModuleX_Lagrange.o \
       ReferenceElementModuleE_Lagrange.o \
       ReferenceElementModule_Lagrange.o \
       MeshModule.o \
       GeometryFieldsModuleE.o \
       GeometryFieldsModule.o \
       GeometryComputationModule.o \
       GeometryComputationModuleE.o \
       FluidFieldsModule.o \
       RadiationFieldsModule.o \
       EquationOfStateModule_TABLE.o \
       EquationOfStateModule_IDEAL.o \
       TwoMoment_ClosureModule.o \
       TwoMoment_PositivityLimiterModule.o \
       TwoMoment_DiscretizationModule_Streaming.o \
       TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
       TwoMoment_NeutrinoMatterSolverModule.o \
       TwoMoment_BoundaryConditionsModule.o \
       TwoMoment_MeshRefinementModule.o \
       TwoMoment_UtilitiesModule.o \
       TwoMoment_MeshRefinementModule.o \
       NeutrinoOpacitiesModule.o \
       OpacityModule_TABLE.o \
       NeutrinoOpacitiesComputationModule.o \
       PhysicalConstantsModule.o \
       UnitsModule.o \
       ThornadoInitializationModule.o \
       TimeSteppingModule_Flash.o

$(LIBNAME): $(OBJS)
	$(MKLIB) $(LIBNAME) $(OBJS)
	cp $(LIBNAME) $(LIBDIR)
	cp *.mod $(INCDIR)

.SUFFIXES: .f90 .F90
%.o %.mod : %.f90
	$(FCOMP) $(FFLAGS) $(f90FLAGS) $(DEFINES) $< -o $@
%.o %.mod : %.F90
	$(FCOMP) $(FFLAGS) $(F90FLAGS) $(DEFINES) $< -o $@

.PHONY: clean
clean:
	rm -f $(LIBNAME) $(LIBDIR)/$(LIBNAME) *.mod $(INCDIR)/*.mod *.o

TimeSteppingModule_Flash.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TimersModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_DiscretizationModule_Streaming.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
  TwoMoment_PositivityLimiterModule.o \
  PhysicalConstantsModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o  \
  TimeSteppingModule_Flash.F90

ThornadoInitializationModule.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  DeviceModule.o \
  TimersModule.o \
  QuadratureModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleE.o \
  ReferenceElementModule.o \
  PolynomialBasisModuleX_Lagrange.o \
  PolynomialBasisModule_Lagrange.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModuleE_Lagrange.o \
  ReferenceElementModule_Lagrange.o \
  SubcellReconstructionModule.o \
  EquationOfStateModule_TABLE.o \
  EquationOfStateModule_IDEAL.o \
  OpacityModule_TABLE.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  GeometryFieldsModuleE.o \
  GeometryComputationModuleE.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_PositivityLimiterModule.o \
  TwoMoment_MeshRefinementModule.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
  ThornadoInitializationModule.F90
