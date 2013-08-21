/*
 *  lfit_plot.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 19/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */
#include <cfloat>
#include "trm_subs.h"
#include "trm_array1d.h"
#include "trm_array2d.h"
#include "trm_buffer2d.h"
#include "trm_input.h"
#include "lfit.h"

using namespace std;

template <class X> X get_Y (const Subs::Array1D<X>& x, 
	     const Subs::Array1D<X>& y, X goal) { 

  X tmp;
  // check x array is ordered  
  if(! x.monotonic() ){
    Subs::Array1D<X> store(x);
    Subs::Buffer1D<unsigned long int> key = store.sort();
    int minpos = Subs::locate(store.ptr(),store.size(),goal);
    // linearly interpolate to find y at that x
    tmp = Subs::linterp(x[key[minpos-1]],y[key[minpos-1]],			    
			      x[key[minpos]],y[key[minpos]],goal);
  }else{
    int minpos = Subs::locate(x.ptr(),x.size(),goal);
    // linearly interpolate to find y at that x
    tmp = Subs::linterp(x[minpos-1],y[minpos-1],			    
			x[minpos],y[minpos],goal);
  }
  return tmp;
}

void plotLight(const LFIT::SynthData& model, 
			   const LFIT::Data& data, Subs::Plot& plot,
			   string autoscale, double& ymin, double& ymax){
	
    Subs::Array1D<float> time, actflux;
	int nregions=0, thisSize=1;
	Subs::Buffer1D<int> regionSizes;
	Subs::Buffer1D<bool> regionBad;

	// find how many separate good and bad regions exist,
	// and their size
	for(size_t i=0; i<data.size(); i++){
		// are we at end or region change?
		if( (i == data.size()-1) || (data[i].bad != data[i+1].bad) ){
			nregions++;
			regionSizes.push_back(thisSize);
			regionBad.push_back(data[i].bad);
			thisSize=1;
		}else{
			// no change
			thisSize++;
		}
	}
	
	// allocate buffers for good and bad regions
	Subs::Buffer2D<float> times(nregions,data.size()), fluxes(nregions,data.size());
	Subs::Array2D<float> residuals(nregions,data.size());
	int icount=0;
	for (int i=0; i<nregions; i++){
		for(int j=0; j<regionSizes[i]; j++){
			times[i][j]  = data[icount].time;
			fluxes[i][j] = data[icount].flux;
			if(model.flux.size() != 0){
				residuals[i][j] = data[icount].flux - 
					get_Y(model.phase,model.flux,float(data[icount].time));
			}
			icount++;
		}
	}
	
    for(size_t i=0; i<data.size(); i++){
		time.push_back(data[i].time);
		actflux.push_back(data[i].flux);
    }

    cpgslw(2);
	
    if(Subs::toupper(autoscale) == "Y"){
		Subs::Array1D<double> Reslim,Fluxlim;
		for(int i=0; i < nregions; i++){
			if(!regionBad[i]){
				for(int j=0; j < regionSizes[i]; j++){
					Reslim.push_back(residuals[i][j]);
					Fluxlim.push_back(fluxes[i][j]);
				}
			}
		}
		if(model.flux.size() != 0){
			ymin=min(float(0.0),float(Reslim.min()));
			ymax=max(model.flux.max(),float(Fluxlim.max()));
		}else{
			ymin=0.0;
			ymax=Fluxlim.max();
		}
    }	
	cpgenv(time.min(),time.max(),ymin,ymax,0,0);
    
    
    string xlab,ylab,plab;
    xlab = "Orbital Phase";
    ylab = "Flux (mJy)";
    plab = "";
    cpglab(xlab.c_str(),ylab.c_str(),plab.c_str());
    int i=1;
	if(model.flux.size() != 0){
		cpgslw(3);
		cpgsci(++i);
		cpgline(model.flux.size(),model.phase.ptr(),model.flux.ptr());
		cpgsci(++i);
		cpgline(model.dflux.size(),model.phase.ptr(),model.dflux.ptr());
		cpgsci(++i);
		cpgline(model.wflux.size(),model.phase.ptr(),model.wflux.ptr());
		cpgsci(++i);
		cpgline(model.bflux.size(),model.phase.ptr(),model.bflux.ptr());
		cpgsci(++i);
		cpgline(model.rflux.size(),model.phase.ptr(),model.rflux.ptr());
	}
	cpgsci(1);
	cpgslw(1);
	
	for(int i=0; i < nregions; i++){
		if(regionBad[i]){
			cpgsci(2);
		}else{
			cpgsci(1);
		}
		Subs::Array1D<float> x,y,z;
		for(int j=0; j < regionSizes[i]; j++){
			x.push_back(times[i][j]);
			y.push_back(fluxes[i][j]);
			if(model.flux.size() != 0){z.push_back(residuals[i][j]);}
		}
		cpgbin(x.size(),x.ptr(),y.ptr(),true);
		if(model.flux.size() != 0){cpgbin(x.size(),x.ptr(),z.ptr(),true);}
		x.clear(); y.clear(); z.clear();
	}
}

void light(const vector<string>& args, const LFIT::SynthData& model, const LFIT::Data& data){
	Subs::Input input(args,"LFIT_ENV",".lfit");
    //sign in input variables
    input.sign_in("autoscale",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("ymin",            Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("ymax",            Subs::Input::LOCAL, Subs::Input::PROMPT);
    
    //get input
    string autoscale;
    double ymin=0, ymax=0;
    input.get_value("autoscale", autoscale, "y", "autoscale (y/n)");
    if(Subs::toupper(autoscale) != "Y"){
        input.get_value("ymin",ymin,0.0,-DBL_MAX,DBL_MAX,"minimum y-value for plot");
        input.get_value("ymax",ymax,0.2,-DBL_MAX,DBL_MAX,"maximum y-value for plot");
    }
    Subs::Plot plot;
    plot.open("?");
	plotLight(model,data,plot,autoscale,ymin,ymax);
    plot.close();
}
