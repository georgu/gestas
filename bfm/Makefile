
BFMDIR = /home/georg/work/shyfem_repo/bfm-for-shyfem-private
INCL_BFM     = -I$(BFMDIR)/include

FFLAGS = -cpp $(INCL_BFM)

default:

bfm:
	gfortran -c bfm.f90
	gfortran -c $(FFLAGS) bfm_statevars.f90

clean:
	-rm -f *.o *.mod a.lout

cleanall: clean

