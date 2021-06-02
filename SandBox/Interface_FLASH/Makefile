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
MOMENT_CLOSURE = MINERBO
MICROPHYSICS = WEAKLIB
TWOMOMENT_ORDER = ORDER_V

# Import thornado Makefiles
THORNADO_DIR = ../../
include $(THORNADO_DIR)/Build/Makefile_GPU
include $(THORNADO_DIR)/Build/Makefile_Path
include $(THORNADO_DIR)/Build/Makefile_Defines
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

# BUILDFLAG can be passed to the Makefile to choose specific
# compilation flags (one of OPT,TEST,DEBUG).

ifeq ($(TWOMOMENT_ORDER),ORDER_V)
  VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV
endif

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
ifeq ($(MICROPHYSICS),WEAKLIB)
  FFLAGS += -I../../../../weaklib/include
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
       TwoMoment_UtilitiesModule.o \
       TwoMoment_MeshRefinementModule.o \
       NeutrinoOpacitiesModule.o \
       OpacityModule_TABLE.o \
       NeutrinoOpacitiesComputationModule.o \
       PhysicalConstantsModule.o \
       UnitsModule.o \
       ThornadoInitializationModule.o \
       TimeSteppingModule_Flash.o

ifeq ($(TWOMOMENT_ORDER),ORDER_V)
OBJS += \
       EquationOfStateModule.o \
       Euler_BoundaryConditionsModule.o \
       Euler_dgDiscretizationModule.o \
       Euler_UtilitiesModule_NonRelativistic.o \
       Euler_SlopeLimiterModule_NonRelativistic_TABLE.o \
       Euler_PositivityLimiterModule_NonRelativistic_TABLE.o \
       Euler_ErrorModule.o \
       InputOutputModuleHDF.o \
       ReferenceElementModuleZ.o \
       PolynomialBasisMappingModule.o \
       UtilitiesModule.o \
       TimersModule_Euler.o \
       TwoMoment_UtilitiesModule_OrderV.o \
       TwoMoment_TroubledCellIndicatorModule.o \
       TwoMoment_SlopeLimiterModule_OrderV.o \
       TwoMoment_TimersModule_OrderV.o \
       TwoMoment_PositivityLimiterModule_OrderV.o \
       TwoMoment_OpacityModule_OrderV.o \
       TwoMoment_DiscretizationModule_Collisions_OrderV.o \
       TwoMoment_NeutrinoMatterSolverModule_OrderV.o \
       TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
       TwoMoment_DiscretizationModule_Streaming_OrderV.o
endif

ifeq ($(USE_GPU),TRUE)
OBJS += \
       CudaModule.o
endif
ifeq ($(USE_OACC),TRUE)
OBJS += \
       OpenACCModule.o
else ifeq ($(USE_OMP_OL),TRUE)
OBJS += \
       OpenMPModule.o
endif
ifeq ($(USE_MAGMA),TRUE)
OBJS += \
       MagmaModule.o CublasModule.o CusolverModule.o CusparseModule.o
else ifeq ($(USE_CUBLAS),TRUE)
OBJS += \
       CublasModule.o CusolverModule.o CusparseModule.o
endif

$(LIBNAME): $(OBJS)
	$(MKLIB) $(LIBNAME) $(OBJS)
	cp $(LIBNAME) $(LIBDIR)
	cp *.mod $(INCDIR)

.SUFFIXES: .f90 .F90
%.o %.mod : %.f90
	$(FCOMP) $(FFLAGS) $(filter-out -qrealsize=8,$(f90FLAGS)) $(DEFINES) $< -o $@
%.o %.mod : %.F90
	$(FCOMP) $(FFLAGS) $(filter-out -qrealsize=8,$(F90FLAGS)) $(DEFINES) $< -o $@

.PHONY: clean
clean:
	rm -f $(LIBNAME) $(LIBDIR)/$(LIBNAME) *.mod $(INCDIR)/*.mod *.o


ifeq ($(TWOMOMENT_ORDER),ORDER_1)
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
  ThornadoInitializationModule.F90
else
TimeSteppingModule_Flash.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TimersModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_TroubledCellIndicatorModule.o \
  PolynomialBasisMappingModule.o \
  TwoMoment_SlopeLimiterModule_OrderV.o \
  TwoMoment_PositivityLimiterModule_OrderV.o \
  Euler_BoundaryConditionsModule.o \
  TwoMoment_DiscretizationModule_Streaming_OrderV.o \
  ReferenceElementModuleZ.o \
  TwoMoment_OpacityModule_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_OrderV.o \
  Euler_SlopeLimiterModule_NonRelativistic_TABLE.o \
  Euler_PositivityLimiterModule_NonRelativistic_TABLE.o \
  InputOutputModuleHDF.o \
  Euler_dgDiscretizationModule.o \
  TwoMoment_NeutrinoMatterSolverModule_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
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
  TwoMoment_TimersModule_OrderV.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_TroubledCellIndicatorModule.o \
  PolynomialBasisMappingModule.o \
  TwoMoment_SlopeLimiterModule_OrderV.o \
  TwoMoment_PositivityLimiterModule_OrderV.o \
  Euler_BoundaryConditionsModule.o \
  TwoMoment_DiscretizationModule_Streaming_OrderV.o \
  ReferenceElementModuleZ.o \
  TwoMoment_OpacityModule_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_OrderV.o \
  Euler_SlopeLimiterModule_NonRelativistic_TABLE.o \
  Euler_PositivityLimiterModule_NonRelativistic_TABLE.o \
  InputOutputModuleHDF.o \
  Euler_dgDiscretizationModule.o \
  TwoMoment_MeshRefinementModule.o \
  ThornadoInitializationModule.F90
endif
