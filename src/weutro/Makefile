
#------------------------------------------------------------------------
#
#    Copyright (C) 2024-2026  Georg Umgiesser
#
#    This file is part of GESTAS.
#
#------------------------------------------------------------------------

DIR = weutro

LIBNAME = $(DIR)
LIBDIR  = .

LIB = $(LIBDIR)/lib$(LIBNAME).a

FGNU_BOUNDS = -fcheck=all
TRAP_LIST = zero,invalid,overflow
FGNU_NOOPT = -g -fbacktrace -ffpe-trap=$(TRAP_LIST) $(FGNU_BOUNDS)

FFLAGS = -g -fimplicit-none $(FGNU_NOOPT)
LFLAGS =

#FFLAGS = -cpp 

MOD  = \

OBJ  = \
        shympi_fake.o \
        writevars.o \
        weutro.o \
        weutro_param.o \

#------------------------------------------------------

default: lib

nothing:
	@true

lib:    $(LIB) nothing

$(LIB): $(MOD) $(OBJ)
	ar rvs $(LIB) $(MOD) $(OBJ)

clean:
	-rm -f *.o *.mod a.out

cleanall: clean
	-rm -f $(LIB)

#------------------------------------------------------

.SUFFIXES: .F .f90

.f.o:
	$(F77) -c $(FFLAGS) $<

.f90.o:
	$(F77) -c $(FFLAGS) $<

.F.o:
	$(F77) -c -cpp $(FFLAGS) $<

#------------------------------------------------------

