include ../../../Build/Makefile_Build

all: dgMatrixTest

dgMatrixTest: \
	$(thornado) dgMatrixTest.o
	$(FLINKER) $(FLAGS) -o dgMatrixTest_$(MACHINE) \
	$(thornado) dgMatrixTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f dgMatrixTest_$(MACHINE)
