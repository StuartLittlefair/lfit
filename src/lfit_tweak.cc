/*
 *  lfit_tweak.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 19/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */

#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/input.h"
#include "lfit.h"
#include <cfloat>
#include <climits>


std::pair<int,double> run_amoeba(LFIT::Params& params, amFunc& func, double ftol, int nmax){
	
	
	// Each pair contains an array with the params which vary and
	// the value of chisq with those params.
	std::vector<std::pair<Subs::Array1D<double>, double> > paramsIn;
	int ndim = func.nvary();
	
	// get variable parameters and make one corner in param space
	Subs::Array1D<double> corner = func.get_variables();
	func.setParamRanges();
	Subs::Array1D<double> deltas = func.get_ranges();
	Subs::Array1D<double> ncorner=corner;
	
	// find chisq at this corner
	double chisq = func(corner);
	paramsIn.push_back(std::make_pair(corner,chisq));
	
	// now loop over all other corners in parameter space
	for(int i=0; i<ndim; i++){
		// set corner back to original parameters
		for(int j=0; j<ndim; j++) ncorner[j] = corner[j];
		// perturb one parameter
		ncorner[i] += deltas[i];
		chisq = func(ncorner);
		paramsIn.push_back(std::make_pair(ncorner, chisq));
	} 
	
	int nfunc;
	
	// Run amoeba
	Subs::amoeba(paramsIn,ftol,nmax,func,nfunc);
	
	// retrieve minimum chisq
	int nbest = 0;
	double cmin = paramsIn[0].second; // minimum chi**2
	for(size_t i=1; i<paramsIn.size(); i++){
		if(paramsIn[i].second < cmin){
			nbest = i;
			cmin = paramsIn[i].second;
		}
	}
	
	// update parameters
	params.updateWith(paramsIn[nbest].first);
		
	std::pair<int,double> results(nfunc,cmin);
	return results;
}

void tweak(const std::vector<std::string>& args, LFIT::Params& params, 
		   const LFIT::ThreeColour& three, const LFIT::CV& myCV, 
		   const LFIT::SynthData& myModel, 
		   int& colFit,const bool& pos){
	
	// which colour to fit?
	Subs::Input input(args,"LFIT_ENV",".lfit");
	input.sign_in("colour",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("ftol",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nmax",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("qvar",        Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("dpvar",       Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("rdvar",       Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("rwdvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("ulvar",       Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("scvar",       Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("azvar",       Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("fracvar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("exp1var",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("exp2var",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("tiltvar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("yawvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("expvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("phivar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	
	input.get_value("colour", colFit, 2, 1, 3, "lightcurve colour to fit");
    colFit -= 1;

    double ftol;
    input.get_value("ftol", ftol, 1.e-5, DBL_MIN, 0.1, "fractional tolerance for convergence");    
	
    int nmax;
    input.get_value("nmax", nmax, 1000, 1, INT_MAX, "maximum number of model evaluations"); 
	
	input.get_value("qvar", params.qvar, true,     "Mass ratio, q:      vary (y/n)?");
	input.get_value("dpvar", params.dpvar, true,   "Delta Phi:          vary (y/n)?");
	input.get_value("rdvar", params.rdvar, true,   "Disc Radius:        vary (y/n)?");
	input.get_value("rwdvar", params.rwdvar, true, "WD Radius:          vary (y/n)?");
	input.get_value("ulvar", params.ulvar, true,   "Limb darkening:     vary (y/n)?");
	input.get_value("scvar", params.scvar, true,   "BS scale:           vary (y/n)?");
	input.get_value("azvar", params.azvar, true,   "BS Azimuth:         vary (y/n)?");
	input.get_value("fracvar", params.fracvar, true, "Isotropic Fraction: vary (y/n)?");
	if(!params.complexSpot){
		params.exp1var=false; params.exp2var=false; params.tiltvar=false; params.yawvar=false;
	}else{
		input.get_value("exp1var", params.exp1var, true,   "BS Exponent1:           vary (y/n)?");
		input.get_value("exp2var", params.exp2var, true,   "BS Exponent1:           vary (y/n)?");
		input.get_value("tiltvar", params.tiltvar, true,   "BS Tilt:                vary (y/n)?");
		input.get_value("yawvar",  params.yawvar,  true,   "BS Yaw:                 vary (y/n)?");
	}
	input.get_value("expvar", params.expvar, true, "Disc Exponent:      vary (y/n)?");
	input.get_value("phivar", params.phivar, true, "Phase Offset:       vary (y/n)?");
	std::cout << "Number of variables: " << params.nvary() << std::endl;

	// this will render any previously calculated errors invalid
	params.errSet = false;


	// create function object
	// when called it scales model CV to the data (colour colFit)
	// and returns chisq
	amFunc func = amFunc(three,myCV,myModel,params,colFit,pos);
	
	std::pair<int,double> amoebaResult;
	amoebaResult = run_amoeba(params,func,ftol,nmax);
	
	Subs::Format form;
	form.precision(8); form.general();
	
	show(params);
	
	int dof = DOF(three[colFit],params);
	std::cout << "Minimum chi**2 = " << form(amoebaResult.second) 
	<< " DOF =  " << dof << std::endl
	<< " number of iterations: " << amoebaResult.first << std::endl;
	
	// force LFIT to scale lightcurve again before plotting
	colFit = -1;
}

