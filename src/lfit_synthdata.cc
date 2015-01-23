/*
 *  lfit_synthdata.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 19/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */
#include "trm/subs.h"
#include "lfit.h"

using namespace std;

void calc_synth_lightcurve(const LFIT::Params& params, const LFIT::Data& data,
                           LFIT::CV& cv, LFIT::SynthData& model){
    if(model.phase.size() != (int)data.size()){
        model.phase.resize(data.size());
        model.dflux.resize(data.size());
        model.bflux.resize(data.size());
        model.wflux.resize(data.size());
        model.rflux.resize(data.size());
        model.flux.resize(data.size());
    }
    
    for(size_t i=0; i<data.size(); ++i){
        double phi = data[i].time-params.phi0;
        double width = data[i].expose;
        model.phase[i]=(phi+params.phi0);
        model.dflux[i]=(params.dFlux*cv.disc.calcFlux(params.q,phi,width,params.incl));
        model.bflux[i]=(params.bFlux*cv.bspot.calcFlux(params.q,phi,width,params.incl));
        model.wflux[i]=(params.wFlux*cv.prim.calcFlux(params.q,phi,width,params.incl));
        model.rflux[i]=(params.rFlux*cv.donor.calcFlux(phi,width,params.incl));           
        model.flux[i]=(model.dflux[i]+model.bflux[i]+
                       model.wflux[i]+model.rflux[i]);
    }
}

int DOF(const LFIT::Data& data, const LFIT::Params& params){
	int nbad = 0;
	for(int i=0; i< data.size(); i++){
		if(data[i].bad)nbad++;
	}
	return (int)data.size() -(int)params.nvary() - nbad;
}

double chisq(const LFIT::Data& data, const LFIT::SynthData& model){
	int ndata = data.size();
	double sum=0.0;
	for(int i=0; i< ndata; i++){
		if(!data[i].bad){
			sum += Subs::sqr(data[i].flux - model.flux[i])/Subs::sqr(data[i].ferr);
		}
	}
	return sum;
}
