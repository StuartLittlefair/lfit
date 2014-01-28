/*
 *  lfit_errors.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 25/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */
#include<cfloat>
#include<climits>

#include "lfit.h"
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/input.h"

using namespace std;

void errors_sem(const std::vector<std::string>& args, LFIT::Params& params, LFIT::ThreeColour& three, LFIT::CV& myCV, LFIT::SynthData& myModel, 
int& colFit,const bool& pos){
	
	// which colour to fit?
	Subs::Input input(args,"LFIT_ENV",".lfit");
	input.sign_in("colour",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("ftol",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("nmax",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("qvar",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("dpvar",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("rdvar",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("rwdvar",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("ulvar",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("scvar",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("azvar",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("fracvar",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("expvar",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("phivar",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("opt",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("patient",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
	
	// are we ready to go?
	bool keepGoing=false;
	input.get_value("opt", keepGoing, true,     "Have you optimised parameters with TWEAK (y/n)?");
	if (!keepGoing) return;
	input.get_value("patient", keepGoing, true,     "This WILL take a while. Continue (y/n)?");
	if (!keepGoing) return;
	
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
	input.get_value("expvar", params.expvar, true, "Disc Exponent:      vary (y/n)?");
	input.get_value("phivar", params.phivar, true, "Phase Offset:       vary (y/n)?");
	std::cout << "Number of variables: " << params.nvary() << std::endl;
	
	Subs::Format form;
	form.precision(8); form.general();
	
	// run amoeba with nmax=0 to get chisq at optimal parameters
	amFunc func = amFunc(three,myCV,myModel,params,colFit,pos);
	
	pair<int,double> amoebaResult;
	amoebaResult = run_amoeba(params,func,ftol,0);
	
	double DOF = three[colFit].size()-(double)params.nvary();
	double minChi = amoebaResult.second;
	std::cout << "Best Chisq = " << form(minChi) << std::endl;
	std::cout << "Entering RTBIS loop to find Delta Chisq = 1.0" << std::endl;
	
	// RTBIS loop - solves for f = 0
	// where f = Chisq - OptimalChisq - 1.0;

	// save original params
	LFIT::Params masterParams = params;

	// start with Q
	params.qvar = false;
	
	
	double f; // difference in chisq between current position and best position
	double xacc = params.q/1000.0; // whats the minimum change we are prepared to accept?
	DOF = three[colFit].size()-1.0;

	
	
	// first find stepsize which brackets root
	// initial guess for this stepsize is 5% of parameter
	double stepsize = 0.05;
	int counter=0;
	bool stepTooSmall=true;
	while(stepTooSmall){
		// evaluate (Chisq-MinChi) - 1.0 at param + stepsize
		params.q = masterParams.q*(1.0+stepsize); 
		
		/* update params (this call syntax will force the function to recalculate all
			stored arrays the next time it is called. This ensures it deals correctly
			with the change in q)
		*/
		func.updateParams(params);

		// run amoeba
		amoebaResult = run_amoeba(params,func,ftol,nmax);
		double thisChi = amoebaResult.second;
		
		// is stepsize enough?
		f = thisChi-minChi-1.0;
		if(-1.0*f >= 0.0){
			// error, root not bracketed
			std::cout << "Chisq = " << form(thisChi) << std::endl;
			std::cout << "Delta Chisq = " << form(f) << std::endl;
			std::cout << "Root not bracketed, increasing step size" << std::endl;
			stepsize*=4.0;
			counter++;
		}else{
			stepTooSmall=false;
		}
		if(counter > 20) {
			std::string err = "Skip this Parameter";
			throw err;
		}
	}
	// new location in parameter space brackets root.
	// calculate RTBIS and dx
	double x2 = masterParams.q;
	double x1 = params.q;
	double dx;
	double rtbis = f<0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	
	int j;
	double oldErr=0.0;
	const int JMAX=40;
	for(j=0;j<JMAX;j++){
		
		// refine parameter closer to solution
		dx*=0.5;
		double xmid = rtbis + dx;
		params.q = xmid;
		
		// check change in error is still significant
		double error = fabs(masterParams.q - params.q);
		if( fabs(oldErr-error)/error < 0.01 ) break;
		oldErr = error;
		std::cout << "Current parameter error = " << form(error) << std::endl;
		
		// optimise again at new parameter value
		func.updateParams(params);
		amoebaResult = run_amoeba(params,func,ftol,nmax);
		double thisChi = amoebaResult.second;
		
		f= thisChi - minChi - 1.0;
		std::cout << "Delta Chisq = " << form(fabs(f)+1.0) << std::endl;
		
		if(f<0.0) rtbis=xmid;
		if(fabs(dx)<xacc || f == 0.0) break; // we are done
	
	}
	// did we exceed our limit?
	if(j == JMAX-1){
		string error = "Too many bisections in RTBIS";
		throw error;
	}
	
	// if not, accept RTBIS as parameter which gives delta_Chisq = 1.0
	double error = fabs(params.q - masterParams.q);
	std::cout << "Q = " << form(masterParams.q) << " +/- " << form(error); 
	

}


void errors_bootstrap(const std::vector<std::string>& args, LFIT::Params& params, 
					  LFIT::ThreeColour& three, LFIT::CV& myCV, LFIT::SynthData& myModel, 
					  int& colFit,const bool& pos){
		
	// which colour to fit?
	Subs::Input input(args,"LFIT_ENV",".lfit");
	input.sign_in("colour",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("ftol",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nmax",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("qvar",       Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("dpvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("rdvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("rwdvar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("ulvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("scvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("azvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("fracvar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("exp1var",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("exp2var",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("tiltvar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("yawvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("expvar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("phivar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
	input.sign_in("nboot",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("opt",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("patient",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("outfile",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

	// are we ready to go?
	bool keepGoing=false;
	input.get_value("opt", keepGoing, true,     "Have you optimised parameters with TWEAK (y/n)?");
	if (!keepGoing) return;
	input.get_value("patient", keepGoing, true,     "This WILL take a while. Continue (y/n)?");
	if (!keepGoing) return;
	
	input.get_value("colour", colFit, 2, 1, 3, "lightcurve colour to fit");
    colFit -= 1;

	int nboot;
	input.get_value("nboot", nboot, 200, 1, INT_MAX, "number of boostrap evaluations");
	
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
	
	
	/*
	 ** Enter bootstrapping loop
	 */
	LFIT::Params masterParams=params;
	vector<pair<double,LFIT::Params> > bootStrapResults;
	
	// create progress bar
	int nhashes=20;
	std::string message = "Bootstrapping to find errors";
	LFIT::ProgressBar bar = LFIT::ProgressBar(message,nhashes,nboot);
	bar.initialise();
	
	for(int iboot=0; iboot<nboot; iboot++){
		
		// print extra hash to status bar if necessary
		bar.update(iboot);
		
		// bootstrap resample lightcurve
		LFIT::ThreeColour booted = bootstrap(three, colFit);
		
		// reset parameters
		params = masterParams;
		
		// Make corners in parameter space for amoeba
		// paramsIn is a vector of N+1 pairs, where N is the number
		// of variable parameters
	
		// first, perturb params from optimal solution by 5%
		params.jitter(5);
		
		// create function object
		// when called it scales model CV to the data (colour colFit)
		// and returns chisq
		amFunc func = amFunc(booted,myCV,myModel,params,colFit,pos);
		
		pair<int,double> amoebaResult;
		amoebaResult = run_amoeba(params,func,ftol,nmax);

		// store params and chisq for this bootstrap
		pair<double,LFIT::Params> bestFit(amoebaResult.second,params);
		//std::cout << params.toString() << ": " << Subs::str(amoebaResult.second) << std::endl;
		bootStrapResults.push_back( bestFit );

	}
	// end progress bar
	bar.end();
	
	// Bootstrapping done.
	// now write out the results of the bootstrapping process...
	std::string ofname;
	input.get_value("outfile", ofname, "bootstrap.results", "Output file for bootstrap results");
	std::cout << "Bootstrapping done. Have a look in " << ofname << " for, well, you guessed it" << std::endl;

	std::ofstream fout;
	fout.open(ofname.c_str());
	fout << "#Chisq  Q  dPhi  Rd  Rwd  Ulimb  bsScale  bsAz  bsFrac  dExp  Incl  phi0\n";  
	for(size_t i=0; i<bootStrapResults.size(); i++){
		pair<double,LFIT::Params> thisResult=bootStrapResults[i];
		fout << thisResult.first << " " << thisResult.second.toString() << endl;
	}
	fout.close();
}
