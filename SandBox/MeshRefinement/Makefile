include ../../../Build/Makefile_Build

all: MeshRefinementTest_TwoMoment

MeshRefinementTest_TwoMoment: \
	$(thornado) MeshRefinementTest_TwoMoment.o
	$(FLINKER) $(FLAGS) -o MeshRefinementTest_TwoMoment_$(MACHINE) \
	$(thornado) MeshRefinementTest_TwoMoment.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f MeshRefinementTest_TwoMoment_$(MACHINE)
