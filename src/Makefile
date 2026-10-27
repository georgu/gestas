
#------------------------------------------------------------------------
#
#    Copyright (C) 2024-2026  Georg Umgiesser
#
#    This file is part of GESTAS.
#
#------------------------------------------------------------------------

#-----------------------------------------------------------

DIR = gestas

#-----------------------------------------------------------

MAKEFLAGS += --no-print-directory

FEMBIN = $(HOME)/shyfem/bin
MAKEDEPEND = $(FEMBIN)/mkdp.pl
MODDIR = .

SUBDIRS = bfm test data weutro
SUBDIRS = bfm weutro aquabc

DATE=`date +%Y-%m-%d`

FGNU_BOUNDS = -fcheck=all
TRAP_LIST = zero,invalid,overflow
FGNU_NOOPT = -g -fbacktrace -ffpe-trap=$(TRAP_LIST) $(FGNU_BOUNDS)

F77 = gfortran
FFLAGS = -g -fimplicit-none
FFLAGS = -g -fimplicit-none $(FGNU_NOOPT)
LFLAGS =

WHAT = weutro aquabc bfm

AQUABC = aquabc
AQUABCOBJ = aquabc/aquabc_fake.o
AQUABCDIR = /home/georg/work/aquabc/ESTAS_AQUABC/modelCore/build
AQUABCLIB = -L$(AQUABCDIR) -laquabc

WEUTRO = weutro
WEUTRODIR = $(WEUTRO)
WEUTRONAME = lib$(WEUTRO).a
WEUTROLIB = -L$(WEUTRODIR) -l$(WEUTRO)

BFM = bfm
BFMDIR = $(BFM)
BFMOBJ = bfm/bfm_fake.o

BFMLIBDIR = /home/georg/work/bfm/bfm
BFMOBJSDIR = $(BFMLIBDIR)/build/tmp/BFMLIB
BFMLIB = -L$(BFMLIBDIR)/lib -lbfm
#MODULES= -I$(BFMDIR)/include
#include $(BFMOBJSDIR)/make_list.txt

NETCDFLIB = -L/usr/lib/x86_64-linux-gnu  -lnetcdf -lnetcdff

LIB = $(AQUABCLIB)
LIB = $(WEUTROLIB)

OBJS1 = subdate.o subscn.o sublintp.o subclo.o
OBJS1 = date.o convert.o lintp.o clo.o
OBJS2 = gestas.o gestas_util.o \
		gestas_scalar.o gestas_wrt.o gestas_conz.o \
		gestas_steady.o \
		gestas_weutro.o \
		gestas_aquabc.o \
		gestas_perpetual.o \
		gestas_read.o \
		gestas_out.o \

		#gestas_aquabc.o \
		#gestas_bfm.o \

OBJS = $(OBJS1) $(OBJS2) $(AQUABCOBJ) $(BFMOBJ)
EXES = boxinf boxindex boxscan gestas fileinf simplediff
EXES_BFM = bfm_util test_bfm boxcompare

#-----------------------------------------------------------

default: weutro aquabc bfm $(EXES)
#default: $(LIB) $(EXES)

all: default

$(AQUABC):
	cd $(AQUABC); make

weutro:
	cd $(WEUTRODIR); make

bfm:
	cd $(BFMDIR); make

boxinf: clo.o boxutil.o boxinf.o
	gfortran -o $@ clo.o boxutil.o $@.o

fileinf: clo.o fileinf.o
	gfortran -o $@ clo.o $@.o

boxindex: boxindex.o
	gfortran -o $@ $@.f90

boxscan: boxscan.o
	gfortran -o $@ $@.f90

gestas: $(WHAT) $(OBJS)
	$(F77) $(LFLAGS) $(OBJS) $(LIB) -o $@

test_bfm: test_bfm.o
	gfortran -o $@ $(MODULES) $(BFMOBJS) $(NETCDFLIB) $@.f90

bfm_util: bfm_util.o
	gfortran -o $@ $(MODULES) $(BFMOBJS) $(NETCDFLIB) $@.f90

simplediff: simplediff.o
	gfortran -o $@ $@.f90

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

.PHONY: clean cleanall weutro aquabc bfm

.f.o:
	$(F77) -c $(FFLAGS) $(MODULES) $<

.f90.o:
	$(F77) -c $(FFLAGS) $(MODULES) $<

%.o: %.f90
	$(F77) $(FFLAGS) $(MODULES) -c $<

#---------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

./bfm_estas.mod: bfm_util.o
./boxes.mod: gestas.o
./boxes_aqua_bc.mod: gestas_aquabc.o
./boxes_bfm.mod: gestas_bfm.o
./boxes_conz.mod: gestas_conz.o
./boxes_scalar.mod: gestas_scalar.o
./boxes_steady.mod: gestas_steady.o
./boxes_weutro.mod: gestas_weutro.o
./boxes_wrt.mod: gestas_wrt.o
./mod_boxinf.mod: boxinf.o
./mod_boxutil.mod: boxutil.o
./mod_fileinf.mod: fileinf.o
bfm_util.o: ./api_bfm.mod ./constants.mod ./init_var_bfm_local.mod \
		./mem.mod ./mem_par.mod ./time.mod 
boxinf.o: ./clo.mod ./mod_boxutil.mod
fileinf.o: ./clo.mod
gestas.o: ./mod_date.mod ./mod_lintp.mod
gestas_aquabc.o: ./boxes.mod ./mod_lintp.mod
gestas_bfm.o: ./boxes.mod
gestas_conz.o: ./boxes.mod
gestas_out.o: ./boxes.mod
gestas_perpetual.o: ./boxes.mod ./mod_date.mod
gestas_read.o: ./boxes.mod ./mod_date.mod
gestas_scalar.o: ./boxes.mod ./mod_date.mod
gestas_steady.o: ./boxes.mod
gestas_weutro.o: ./boxes.mod ./mod_date.mod ./mod_lintp.mod
gestas_wrt.o: ./boxes.mod
test_bfm.o: ./standalone.mod

