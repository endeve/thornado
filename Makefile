MACHINE = $(THORNADO_MACHINE)

OPT_LEVEL = DEBUG
FLAGS     = $(FLAGS_$(OPT_LEVEL))

include ./Build/Makefile_Build

thornado: $(thornado)

clean:
	rm -f *.o *.mod
