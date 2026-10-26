
DIR = aquabc

LIBNAME = $(DIR)
LIBDIR  = .

LIB = $(LIBDIR)/lib$(LIBNAME).a

FFLAGS = -cpp 

MOD  = \

OBJ  = \
        aquabc_fake.o \

#------------------------------------------------------

default: lib

nothing:
	@true

lib:    compile nothing

compile: $(MOD) $(OBJ)

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

