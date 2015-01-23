# Compiler
CC=g++
CFLAGS=-O3

DATE		= $(shell date)
CURRENT		= $(shell git tag | tail -n 1)
PREVIOUS	= $(shell git tag | tail -n 2 | head -n 1)

TRM_SOFTWARE = /usr/local/trm_code

LDFLAGS = -L/usr/local/lib -L$(TRM_SOFTWARE)/lib -L/usr/X11/lib
CPPFLAGS = -I/usr/local/include -I$(TRM_SOFTWARE)/include -I/usr/X11/include -I./include

TRM_LIBS   =  -lcsla -lsubs -lroche
PGPLOT_LIBS =  -lpgplot -lcpgplot
EXTRA_LIBS = -lreadline -lgsl -lblas

SOURCES = src/BrightSpot.cc src/Disc.cc src/Donor.cc src/WhiteDwarf.cc \
src/lultracam.cc src/finddeg.cc src/lfit_io.cc src/lfit_mask.cc src/lfit_plot.cc src/lfit_rebin.cc \
src/lfit_scale.cc src/lfit_synthdata.cc src/lfit_tweak.cc src/lfit_bootstrap.cc \
src/lfit_errors.cc src/lfit_levmarq.cc src/cholesky.cc src/lfit_mcmc.cc  

INCLUDES = include/BrightSpot.h include/Donor.h include/cholesky.h include/lfit.h \
include/lultracam.h include/Disc.h include/WhiteDwarf.h include/geometry.h include/lfit_params.h

all: lfit lfit_fake

lfit: src/lfit.cc $(SOURCES) $(INCLUDES)
	$(CC) $(CFLAGS) src/lfit.cc $(SOURCES) $(LDFLAGS) $(CPPFLAGS) $(TRM_LIBS) $(PGPLOT_LIBS) $(EXTRA_LIBS) -o lfit
	
lfit_fake: src/lfit_fake.cc $(SOURCES) $(INCLUDES)
	$(CC) $(CFLAGS) src/lfit_fake.cc $(SOURCES) $(LDFLAGS) $(CPPFLAGS) $(TRM_LIBS) $(PGPLOT_LIBS) $(EXTRA_LIBS) -o lfit_fake