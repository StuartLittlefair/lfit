## Process this file with automake to generate Makefile.in

## This to avoid . being in the path of includes which causes problems for
## functions looking for the 'list' C++ header

AUTOMAKE_OPTIONS = nostdinc

## Programs

progdir = @bindir@/@PACKAGE@

prog_PROGRAMS  = lfit lfit_fake

lfit_SOURCES          = BrightSpot.cc Disc.cc Donor.cc WhiteDwarf.cc \
lultracam.cc finddeg.cc lfit_io.cc lfit_mask.cc lfit_plot.cc lfit_rebin.cc \
lfit_scale.cc lfit_synthdata.cc lfit_tweak.cc lfit_bootstrap.cc lfit_errors.cc \
lfit_levmarq.cc cholesky.cc lfit_mcmc.cc lfit.cc 


DATE = $(shell date)

install-data-hook: $(progdir)/lfit 
	echo "This is $(PACKAGE)-$(VERSION), built on $(DATE)" > $(progdir)/VERSION
	@echo "Installed link targets"



INCLUDES = -I../include -I../. 

# force re-build of reduce if configure.ac changes to allow for
# version number changes
lfit_DEPENDENCIES = ../configure.ac

#LDADD = libultracam.la @LIBCURL@

## Library

#lib_LTLIBRARIES = libultracam.la 

#libultracam_la_SOURCES = window.cc windata.cc ccd.cc frame.cc target.cc \
#mccd.cc skyline.cc defect.cc shift_and_add.cc WriteMemoryCallback.cc \
#parseXML.cc de_multiplex.cc read_header.cc aperture.cc get_server_frame.cc \
#loadXML.cc fdisk.cc lllccd.cc plot_images.cc findpos.cc plot_apers.cc \
#fitgaussian.cc ultracam.cc gauss_reject.cc profit_init.cc moffat_reject.cc \
#fitmoffat.cc pos_tweak.cc fit_plot_profile.cc covsrt.cc extract_flux.cc \
#sky_estimate.cc badInput.cc plot_defects.cc plot_setupwins.cc spectrum.cc \
#make_profile.cc specap.cc sky_move.cc sky_fit.cc ext_nor.cc plot_trail.cc \
#plot_spectrum.cc signal.cc
