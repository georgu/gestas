
#-----------------------------------------------------------

DIR = estas2

#-----------------------------------------------------------

FEMBIN = $(HOME)/shyfem/bin
MAKEDEPEND = $(FEMBIN)/mkdp.pl
MODDIR = .

F77 = gfortran
FFLAGS = -g -fimplicit-none
FFLAGS = -g -fimplicit-none -fcheck=all
LFLAGS =

AQUABCDIR = /home/georg/work/aquabc/ESTAS_AQUABC/modelCore/build
AQUABCLIB = -L$(AQUABCDIR) -laquabc

BFMDIR = /home/georg/work/bfm/bfm
BFMOBJSDIR = $(BFMDIR)/build/tmp/BFMLIB
BFMLIB = -L$(BFMDIR)/lib -lbfm
#MODULES= -I$(BFMDIR)/include
#include $(BFMOBJSDIR)/make_list.txt

NETCDFLIB = -L/usr/lib/x86_64-linux-gnu  -lnetcdf -lnetcdff

LIB = $(AQUABCLIB)

OBJS1 = subdate.o subscn.o sublintp.o
OBJS2 = estas2.o estas2_util.o \
		estas2_scalar.o estas2_wrt.o estas2_conz.o \
		estas2_bfm.o estas2_steady.o estas2_aquabc.o
OBJS = $(OBJS1) $(OBJS2)
EXES = boxinf fileinf boxindex estas2 boxscan simplediff intp
EXES_BFM = bfm_util test_bfm boxcompare

#-----------------------------------------------------------

default: intp

intp: $(OBJS1)
	gfortran -o $@ $(OBJS1)

boxinf: subclo.o boxutil.o boxinf.o
	gfortran -o $@ subclo.o boxutil.o $@.o

fileinf: subclo.o fileinf.o
	gfortran -o $@ subclo.o $@.o

boxindex: boxindex.f
	gfortran -o $@ $@.f

boxscan: boxscan.f
	gfortran -o $@ $@.f

estas2: $(OBJS)
	$(F77) $(LFLAGS) $(OBJS) $(LIB) -o $@

test_bfm: test_bfm.o
	gfortran -o $@ $(MODULES) $(BFMOBJS) $(NETCDFLIB) $@.f

bfm_util: bfm_util.o
	gfortran -o $@ $(MODULES) $(BFMOBJS) $(NETCDFLIB) $@.f

simplediff: simplediff.o
	gfortran -o $@ $@.f

#-----------------------------------------------------------

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.f

#-----------------------------------------------------------

clean:
	-rm -f a.out
	-rm -f $(EXES) $(EXES_BFM)
	-rm -f *.o
	-rm -f *.mod
	-rm -f ggg hhh

cleanall: clean
	-rm -f $(DIR).zip

zip: cleanall
	zip $(DIR).zip *

#-----------------------------------------------------------

.f.o:
	$(F77) -c $(FFLAGS) $(MODULES) $<

#---------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

./bfm_estas.mod: bfm_util.o
./boxes.mod: estas2.o
./boxes_aqua_bc.mod: estas2_aquabc.o
./boxes_bfm.mod: estas2_bfm.o
./boxes_conz.mod: estas2_conz.o
./boxes_scalar.mod: estas2_scalar.o
./boxes_steady.mod: estas2_steady.o
./boxes_wrt.mod: estas2_wrt.o
./clo.mod: subclo.o
./mod_boxinf.mod: boxinf.o
./mod_boxutil.mod: boxutil.o
./mod_date.mod: subdate.o
./mod_fileinf.mod: fileinf.o
./scan_string.mod: subscn.o
bfm_util.o: ./api_bfm.mod ./constants.mod ./init_var_bfm_local.mod \
		./mem.mod ./mem_par.mod ./time.mod 
boxinf.o: ./clo.mod ./mod_boxutil.mod
estas2.o: ./mod_date.mod
estas2_aquabc.o: ./boxes.mod
estas2_bfm.o: ./boxes.mod
estas2_conz.o: ./boxes.mod
estas2_scalar.o: ./boxes.mod
estas2_steady.o: ./boxes.mod
estas2_wrt.o: ./boxes.mod
fileinf.o: ./clo.mod
test_bfm.o: ./standalone.mod

