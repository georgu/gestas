
#------------------------------------------------------------------------
#
#    Copyright (C) 2024-2026  Georg Umgiesser
#
#    This file is part of GESTAS.
#
#------------------------------------------------------------------------

MAKEFLAGS += --no-print-directory

default: compile

compile:
	cd src; make

clean:
	cd src; make clean

