## Makefile for downloading the source code and building
## a binary package of R-INLA
##
## Run 'make' with no arguments for usage instructions
##
## This Makefile has been tested under Debian GNU/Linux, but should
## also work under Ubuntu.
## Requirements:
##   mercurial
##   A full gcc/g++/gfortran installation
##   lib*-dev packages, see target lib-install below
##
## If you already have a cloned inla code repository in a subdirectory 'inla',
## you can skip all 'download-' steps, and instead run
##   make -f inla/build-user/linux/Makefile
##
## Finn Lindgren, 2014-08-24, finn.lindgren@gmail.com
##
## Small changes for "Constructing Priors that Penalize the Complexity of
## Gaussian Random Fields
## Geir-Arne Fuglstad, 2017-03-03, geirarne.fuglstad@gmail.com

SHELL = /bin/bash

## This Makefile is only supported for platform 'linux':
PLATFORM = linux
## Must match the machine architecture (32 or 64):
BITS = 64
## Compiler flags
FLAGS = -m64 -std=gnu99 -Wall -O2 -mtune=generic -mfpmath=sse \
	-msse2 -funroll-loops -ftracer -fopenmp -pipe -DINLA_EXPERIMENTAL -I/usr/share/R/include

## You probably don't need to edit anything below this line

## For Debian or Ubuntu package installation in lib-install:
APT=apt-get

## Path configuration:
MAINDIR=$(PWD)
PREFIX = $(MAINDIR)/local
BINPREFIX = $(PREFIX)/bin/$(PLATFORM)
LIBPREFIX = $(PREFIX)/lib
INCPREFIX = $(PREFIX)/include

## Compiler options:
CC=gcc
CXX=g++
FC=gfortran
LD=$(CXX)
ARGS = CC="$(CC)" CXX="$(CXX)" FC="$(FC)" LD="$(LD)"

all :;
	@(echo "";\
	echo "Edit this Makefile, especially the 'FLAGS' variable,";\
	echo "check if you need to run 'sudo $(MAKE) lib-install', then do";\
	echo "";\
	echo "	$(MAKE) download-PC2017";\
	echo "	$(MAKE) bin";\
	echo "	$(MAKE) package";\
	echo "";\
	echo "If all goes well, you should have binaries for inla and fmesher at";\
	echo "";\
	echo "	$(BINPREFIX)";\
	echo "";\
	echo "Note: Not all lapack versions contain all the needed routines,";\
	echo "so make sure that /usr/lib/lapack/liblapack.so is selected in";\
	echo "";\
	echo "	sudo update-alternatives --config liblapack.so";\
	echo "";\
	echo "It's possible to link to both lapack and lapack_atlas.";\
	echo "";)

## Package list updated 2014-08-23
LIBPKG = liblapack-dev libgsl0-dev zlib1g-dev libsuitesparse-dev \
	libmetis-dev libxdmcp-dev libx11-dev libc6-dev libatlas-dev \
	libblas-dev libmuparser-dev r-base-core
lib-install :
	$(APT) install $(LIBPKG)

init :
	@mkdir -p $(LIBPREFIX)
	@mkdir -p $(BINPREFIX)
	@mkdir -p $(INCPREFIX)

download : init
	@( test -d inla && $(MAKE) download-update ) || \
	( rm -rf inla && hg clone https://bitbucket.org/hrue/r-inla inla )
	@ln -sTf $(PREFIX) inla/local

download-PC2017 : 
	@$(MAKE) download
	@cd inla ; hg update -r spatialPC2017
	@ln -sTf $(PREFIX) inla/local

bin :
	$(MAKE) taucs
	$(MAKE) fmesher
	$(MAKE) GMRFLib
	$(MAKE) inla

update :
	$(MAKE) download-update
	$(MAKE) bin
	$(MAKE) package


fmesher :
	$(MAKE) -C inla/fmesher PREFIX=$(PREFIX) $(ARGS) -k clean
	$(MAKE) -C inla/fmesher PREFIX=$(PREFIX) FLAGS="$(FLAGS)" $(ARGS) EXTLIBS="-lX11 -lgsl -lgslcblas -lxcb -lpthread -lXau -lXdmcp"
	$(MAKE) -C inla/fmesher PREFIX=$(PREFIX) $(ARGS) install
	@rsync -a $(PREFIX)/bin/fmesher $(BINPREFIX)/fmesher$(BITS)
	@rm $(PREFIX)/bin/fmesher

inla :
	ln -sTf /usr/include $(INCPREFIX)/muParser
	$(MAKE) -C inla/inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS) -k clean
	$(MAKE) -C inla/inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS)
	$(MAKE) -C inla/inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS) install
	@rsync -a $(PREFIX)/bin/inla $(BINPREFIX)/inla$(BITS)
	@rm $(PREFIX)/bin/inla $(PREFIX)/bin/inla-snapshot

GMRFLib :
	$(MAKE) -C inla/gmrflib PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS) -k clean
	$(MAKE) -C inla/gmrflib PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS)
	$(MAKE) -C inla/gmrflib PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS) install

sync-taucs :
	mkdir -p $(MAINDIR)/tmp
	@rsync -avz --delete inla/extlibs/taucs-2.2--my-fix/ tmp/taucs-2.2--my-fix/

taucs : sync-taucs
	$(MAKE) -C $(MAINDIR)/tmp/taucs-2.2--my-fix CFLAGS="$(FLAGS)" FFLAGS="$(FLAGS)" $(ARGS) -k clean
	$(MAKE) -C $(MAINDIR)/tmp/taucs-2.2--my-fix CFLAGS="$(FLAGS)" FFLAGS="$(FLAGS)" $(ARGS)
	@cp -v -f $(MAINDIR)/tmp/taucs-2.2--my-fix/lib/linux/libtaucs.a $(PREFIX)/lib

package :
	@rm -f inla/INLA_*.tgz
	cd inla ; utils/build-package-bin $(PREFIX)/bin $(MAINDIR)/tmp
#	@mv inla/INLA_*.tgz .

.PHONY: all inla GMRFLib taucs update fmesher bin init download download-PC2017 sync-taucs package
