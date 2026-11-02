
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

help:
	@echo "available targets"
	@echo "  help           this screen"
	@echo "  compile        compiles program"
	@echo "  clean          cleans directories"
	@echo "  cleanall       thoroughly cleans directories"

clean:
	cd src; make clean
	-rm -f *.tmp
	-rm -f *.bak

cleanall:
	cd src; make cleanall

