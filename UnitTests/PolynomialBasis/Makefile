include ../../../Build/Makefile_Build

all: PolynomialBasisMappingTest

PolynomialBasisMappingTest: \
	$(thornado) PolynomialBasisMappingTest.o
	$(FLINKER) $(FLAGS) -o PolynomialBasisMappingTest_$(MACHINE) \
	$(thornado) PolynomialBasisMappingTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f PolynomialBasisMappingTest_$(MACHINE)
