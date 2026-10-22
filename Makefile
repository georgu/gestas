
MAKEFLAGS += --no-print-directory

default: compile

compile:
	cd src; make

clean:
	cd src; make clean

