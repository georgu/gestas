
BFMDIR = /home/georg/work/shyfem_repo/bfm-for-shyfem-private
INCL_BFM     = -I$(BFMDIR)/include
INCL_BFM     = 

OBJS = bfm.f90 bfm_statevars.f90 bfm_fake.f90
OBJS = bfm_fake.o

F77 = gfortran
FFLAGS = -cpp $(INCL_BFM)

default: compile

nothing:
	@true

compile: $(OBJS) nothing

clean:
	-rm -f *.o *.mod a.lout

cleanall: clean

.SUFFIXES: .F .f90

.f.o:
	$(F77) -c $(FFLAGS) $<

.f90.o:
	$(F77) -c $(FFLAGS) $<

