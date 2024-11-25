include ../../../Build/Makefile_Build

all: ProgenitorReader

ProgenitorReader: \
	$(thornado) ProgenitorReader.o
	$(FLINKER) $(FLAGS) -o ProgenitorReader_$(MACHINE) \
	$(thornado) ProgenitorReader.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ProgenitorReader_$(MACHINE)

ProgenitorReader.o: \
  ProgenitorReader.f90
