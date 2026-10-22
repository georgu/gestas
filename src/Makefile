
#-----------------------------------------------------------

DIR = gestas

#-----------------------------------------------------------

MAKEFLAGS += --no-print-directory

FEMBIN = $(HOME)/shyfem/bin
MAKEDEPEND = $(FEMBIN)/mkdp.pl
MODDIR = .

SUBDIRS = bfm test data weutro
SUBDIRS = bfm weutro

DATE=`date +%Y-%m-%d`

FGNU_BOUNDS = -fcheck=all
TRAP_LIST = zero,invalid,overflow
FGNU_NOOPT = -g -fbacktrace -ffpe-trap=$(TRAP_LIST) $(FGNU_BOUNDS)

F77 = gfortran
FFLAGS = -g -fimplicit-none
FFLAGS = -g -fimplicit-none $(FGNU_NOOPT)
LFLAGS =

WHAT = weutro

AQUABCDIR = /home/georg/work/aquabc/ESTAS_AQUABC/modelCore/build
AQUABCLIB = -L$(AQUABCDIR) -laquabc

WEUTRO = weutro
WEUTRODIR = $(WEUTRO)
WEUTRONAME = lib$(WEUTRO).a
WEUTROLIB = -L$(WEUTRODIR) -l$(WEUTRO)

BFMDIR = /home/georg/work/bfm/bfm
BFMOBJSDIR = $(BFMDIR)/build/tmp/BFMLIB
BFMLIB = -L$(BFMDIR)/lib -lbfm
#MODULES= -I$(BFMDIR)/include
#include $(BFMOBJSDIR)/make_list.txt

NETCDFLIB = -L/usr/lib/x86_64-linux-gnu  -lnetcdf -lnetcdff

LIB = $(AQUABCLIB)
LIB = $(WEUTROLIB)

OBJS1 = subdate.o subscn.o sublintp.o
OBJS2 = gestas.o gestas_util.o \
		gestas_scalar.o gestas_wrt.o gestas_conz.o \
		gestas_steady.o \
		gestas_weutro.o \
		gestas_perpetual.o \
		gestas_read.o \
		gestas_out.o \

		#gestas_aquabc.o \
		#gestas_bfm.o \

OBJS = $(OBJS1) $(OBJS2)
EXES = boxinf boxindex boxscan gestas fileinf simplediff
EXES_BFM = bfm_util test_bfm boxcompare

#-----------------------------------------------------------

default: weutro $(EXES)
#default: $(LIB) $(EXES)

all: default

$(AQUABCLIB):
weutro:
	cd $(WEUTRODIR); make

boxinf: subclo.o boxutil.o boxinf.o
	gfortran -o $@ subclo.o boxutil.o $@.o

fileinf: subclo.o fileinf.o
	gfortran -o $@ subclo.o $@.o

boxindex: boxindex.f
	gfortran -o $@ $@.f

boxscan: boxscan.f
	gfortran -o $@ $@.f

gestas: $(WHAT) $(OBJS)
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
	-rm -f *.tmp *.bak
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleanall: clean
	-rm -f $(DIR)*.zip
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

zip: cleanall
	zip $(DIR).zip *
	zip $(DIR).zip bfm/*
	zip $(DIR).zip data/*
	zip $(DIR).zip test/*
	zip $(DIR).zip weutro/*
	mv $(DIR).zip $(DIR)_$(DATE).zip

#-----------------------------------------------------------

.PHONY: clean cleanall weutro

.f.o:
	$(F77) -c $(FFLAGS) $(MODULES) $<

#---------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

./bfm_estas.mod: bfm_util.o
./boxes.mod: gestas.o
./boxes_aqua_bc.mod: gestas_aquabc.o
./boxes_bfm.mod: gestas_bfm.o
./boxes_conz.mod: gestas_conz.o
./boxes_scalar.mod: gestas_scalar.o
./boxes_steady.mod: gestas_steady.o
./boxes_wrt.mod: gestas_wrt.o
./clo.mod: subclo.o
./mod_boxinf.mod: boxinf.o
./mod_boxutil.mod: boxutil.o
./mod_date.mod: subdate.o
./mod_fileinf.mod: fileinf.o
./mod_lintp.mod: sublintp.o
./scan_string.mod: subscn.o
bfm_util.o: ./api_bfm.mod ./constants.mod ./init_var_bfm_local.mod \
		./mem.mod ./mem_par.mod ./time.mod 
boxinf.o: ./clo.mod ./mod_boxutil.mod
fileinf.o: ./clo.mod
gestas.o: ./mod_date.mod ./mod_lintp.mod
gestas_aquabc.o: ./boxes.mod ./mod_lintp.mod
gestas_bfm.o: ./boxes.mod
gestas_conz.o: ./boxes.mod
gestas_scalar.o: ./boxes.mod
gestas_steady.o: ./boxes.mod
gestas_wrt.o: ./boxes.mod
sublintp.o: ./mod_date.mod
test_bfm.o: ./standalone.mod

