include ../../../Build/Makefile_Build

all: SubcellReconstructionTest

SubcellReconstructionTest: \
	$(thornado) SubcellReconstructionTest.o
	$(FLINKER) $(FLAGS) -o SubcellReconstructionTest_$(MACHINE) \
	$(thornado) SubcellReconstructionTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f SubcellReconstructionTest_$(MACHINE)
