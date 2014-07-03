#!/bin/csh 

g++ -arch x86_64 -O3  -fast BrightSpot.cc Disc.cc Donor.cc WhiteDwarf.cc lultracam.cc finddeg.cc lfit_io.cc lfit_mask.cc lfit_plot.cc lfit_rebin.cc lfit_scale.cc lfit_synthdata.cc lfit_tweak.cc lfit_bootstrap.cc lfit_errors.cc lfit_levmarq.cc cholesky.cc lfit_mcmc.cc lfit.cc `trm_link` -lpgplot -lcpgplot -lcsla -lroche -lreadline -lgsl -lblas -o lfit
