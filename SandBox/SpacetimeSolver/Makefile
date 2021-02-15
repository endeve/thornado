include ../../../Build/Makefile_Build

all: SpacetimeTest

SpacetimeTest: \
	$(thornado) \
	ScalarFieldsModule.o \
	ScalarWave_ProgramInitializationModule.o \
	InitializationModule.o \
	ScalarWave_InputOutputModuleHDF.o \
	ScalarWave_dgDiscretizationModule.o \
	TimeSteppingModule_SSPRK.o \
	ScalarWaveTest.o 
	$(FLINKER) $(FLAGS) -o SpacetimeTest_$(MACHINE) \
	$(thornado) \
	ScalarFieldsModule.o \
	ScalarWave_ProgramInitializationModule.o \
	InitializationModule.o \
	ScalarWave_InputOutputModuleHDF.o \
	ScalarWave_dgDiscretizationModule.o
	TimeSteppingModule_SSPRK.o \
	ScalarWaveTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f SpacetimeTest_$(MACHINE)

ScalarFieldsModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ScalarFieldsModule.F90

ScalarWave_ProgramInitializationModule.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  UtilitiesModule.o \
  QuadratureModule.o \
  ReferenceElementModule.o \
  PolynomialBasisModule_Lagrange.o \
  PolynomialBasisModule_Legendre.o \
  PolynomialBasisMappingModule.o \
  MeshModule.o \
  ScalarFieldsModule.o \
  InitializationModule.F90

InitializationModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  InitializationModule.F90

ScalarWave_InputOutputModuleHDF.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModule.o \
  MeshModule.o \
  InputOutputUtilitiesModule.o \
  ScalarFieldsModule.o \
  ScalarWave_InputOutputModuleHDF.F90

ScalarWave_dgDiscretizationModule.o: \
  KindModule.o \
  ProgramheaderModule.o \
  MeshModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModulex_Lagrange.o \
  ScalarFieldsModule.F90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ScalarFieldsModule.o \
  TimeSteppingModule_SSPRK.f90

ScalarWaveTest.o: \
  KindModule.o \
  ScalarWave_ProgramInitializationModule.o \
  InitializationModule.o \
  ScalarWave_InputOutputModuleHDF.o \
  ScalarWaveTest.f90

