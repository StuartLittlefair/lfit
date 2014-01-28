/*
 *  lfit.h
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 11/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <unistd.h>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/array2d.h"
#include "trm/plot.h"
#include "lultracam.h"
#include "trm/format.h"
#include "WhiteDwarf.h"
#include "BrightSpot.h"
#include "Disc.h"
#include "Donor.h"

namespace LFIT{
    
	
    // convenience structure that holds buffers containing total flux, wd flux, bs flux, disc flux and donor
    // flux of a synthetic lightcurve as well as phase
    struct SynthData{
        Subs::Array1D<float> phase, flux, wflux, bflux, dflux, rflux;
    };

  
    class ThreeColour : public std::vector<Data> {
    public:
      //!default constructor
    ThreeColour() : std::vector<Data>(3){}
      
      //! constructor from 3 Data objects
      ThreeColour(const Data& d1, const Data& d2, const Data& d3){
	this->push_back(d1);
	this->push_back(d2);
	this->push_back(d3);
      }
      
      ~ThreeColour(){
	//std::cout << "Three colours shouldn't get destroyed that often"<<std::endl;
      }
      
    };
    
    struct CV{
        LFIT::Donor donor;
        LFIT::WhiteDwarf prim;
        LFIT::BrightSpot bspot;
        LFIT::Disc disc;

		~CV(){
			//std::cout << "CV object going out of scope" << std::endl;
		}
    };
        
    int findcomm(std::string& trial);
	
	class ProgressBar{
		
	public:
		
		// Default constructor 
		ProgressBar() {
			message="";
			messageLength=0;
			barLength=0;
		}
		
		// constructor from string and length
		ProgressBar(std::string message_, int length_, int duration_){
			this->duration = duration_;
			this->message = message_ + ": completion 0% -";
			this->secondLineLength = this->message.size();
			this->barLength = length_;
			for(size_t i=0; i< this->barLength; i++) this->message += "#";
			this->message += "- 100%\n";
			this->messageLength = this->message.size();
			this->completion = 0;
		}
		
		// destructor
		~ProgressBar(){
			//std::cout <<"Destroying progress bar"<<std::endl;
		}
		
		// print out initial message Bar
		bool initialise(){
			if (this->barLength == 0){
				std::string err = "Error: message bar not initialised";
				throw err;
			}
			std::cout << this->message;
			std::string secondLine;
			for(size_t i=0; i<this->secondLineLength; i++) secondLine += " ";
			std::cout << secondLine;
			return true;
		}
		
		bool update(unsigned int progressCount){
			if (( (this->barLength+1)*progressCount/duration - this->completion) > 0 ){
				int toPrint = 
				(this->barLength+1)*progressCount/duration - this->completion;
				for(int i=0; i< toPrint; i++) {
					std::cout << "-";
					std::cout.flush();
				}
				this->completion++;
			}
			return true;
		}
		
		bool end(){
			std::cout << std::endl;
			return true;
		}
		
	private:
		std::string message;
		unsigned int messageLength;
		unsigned int secondLineLength;
		unsigned int barLength;
		unsigned int duration;
		unsigned int completion;
	};
	
}

// define function calls which will be used in helper routines such as tweak.cc
// display current parameters
void show(LFIT::Params& params);

void write(const std::vector<std::string>& args, const LFIT::Data& data, const LFIT::SynthData& model);

// calculate a synthetic lightcurve
void calc_synth_lightcurve(const LFIT::Params& params, const LFIT::Data& data,
                           LFIT::CV& cv, LFIT::SynthData& model);

// chisq of model w.r.t data
double chisq(const LFIT::Data& data, const LFIT::SynthData& model);

int DOF(const LFIT::Data& data, const LFIT::Params& params);

// scale components of model to data using linear least-squares
double scale(const std::vector<std::string>& args,LFIT::Params& params, 
			 const LFIT::ThreeColour& three, LFIT::CV& myCV, LFIT::SynthData& model, 
			 int& colFit, const bool& batch=false, const bool& pos=false);
			
double fastscale(const std::vector<std::string>& args,LFIT::Params& params, 
	const LFIT::ThreeColour& three, LFIT::CV& myCV, LFIT::SynthData& model, 
	int& colFit, const bool& batch=false, const bool& pos=false, const bool& wdCalc=true,
	const bool& donorCalc=true, const bool& discCalc=true, const bool& bsCalc=true);

// function object to pass to amoeba when tweaking model
class amFunc : public Subs::Afunc{
	
public: 
	
	//! Constructor
	amFunc(const LFIT::ThreeColour& three_, const LFIT::CV& myCV_, const LFIT::SynthData& myModel_, 
		const LFIT::Params& params_, int& colFit_, const bool& pos_) : 
	three(three_),CV(myCV_),model(myModel_),params(params_),
	old_params(params_),colFit(colFit_),pos(pos_),forceRecalc(false){
		//std::cout << "Making amFunc object with " << params.toString() << std::endl;
	}
	
	// destructor
	~amFunc(){
		//std::cout << "Destroying amFunc object" << std::endl;
	}
	
	void updateParams(const LFIT::Params& paramsIn){
		params = paramsIn;
		if(params != old_params){forceRecalc=true;}
	}
	
	void updateParams(const Subs::Array1D<double>& paramsIn, bool& wdCalc,
	bool& discCalc, bool& donorCalc, bool& bsCalc){
		
		LFIT::Params old_params = params;
		params.updateWith(paramsIn);
		
		if (params.q != old_params.q or params.dphi != old_params.dphi
		or params.phi0 != old_params.phi0){
			wdCalc=true;
			donorCalc=true;
			discCalc=true;
			bsCalc=true;
			return;
		}
		
		if(params.rwd != old_params.rwd
			or params.ulimb != old_params.ulimb){wdCalc=true;}

		if (params.rwd != old_params.rwd
			or params.rd != old_params.rd
			or params.dExp != old_params.dExp){discCalc=true;}
			
		if (params.bsAz != old_params.bsAz
			or params.bsFrac != old_params.bsFrac
			or params.bsScale != old_params.bsScale
			or params.bsExp1 != old_params.bsExp1
			or params.bsExp2 != old_params.bsExp2
			or params.bsTilt != old_params.bsTilt
			or params.bsYaw != old_params.bsYaw
			or params.rd != old_params.rd){bsCalc=true;}

	}
	
	std::string asString(){
		return params.toString();
	}
	
	Subs::Array1D<double> get_ranges(){
		return params.get_ranges();
	}
	
	Subs::Array1D<double> get_variables(){
		return params.get_variables();
	}
	
	std::vector< std::pair<double,double> > get_limits(){
	  return params.get_limits(CV.bspot.getTangent());
	}
	
	void setParamRanges(){
		params.setParamRanges();
	}
	
	int nvary(){
		return params.nvary();
	}
	
	double operator()(const Subs::Array1D<double>& paramsIn){
		
		/* this routine needs to be smarter. If colfit hasn't changed
		then the data being fitted (hence x points and x bin widths)
		also hasn't changed. So, we probably don't need to re-fill all the
		model arrays which scale uses
		*/
		bool wdCalc = false;
		bool discCalc = false;	
		bool donorCalc = false;
		bool bsCalc = false;
				
		// update params with new values, and set any boolean's we can
		try{
			updateParams(paramsIn,wdCalc,discCalc,donorCalc,bsCalc);
		}catch(const std::string& err){
			//std::cout << "Warning (non-critical): q, dphi pair was unphysical" << std::endl;
			// since we don't want amoeba or MCMC chain to go here, return very high chisq
			return 1.0e34;
		}
		
		// or have parameters been changed on top of changes implied by updateParams (as in error routine)?
		if (forceRecalc){
			wdCalc=true;
			discCalc=true;
			donorCalc=true;
			bsCalc=true;
			forceRecalc=false;
		}
		
		// use these params to scale model
		std::vector<std::string> dummy;
		dummy.push_back("");
		float chisq = fastscale(dummy,params,three,CV,model,colFit,true,pos,
			wdCalc,donorCalc,discCalc,bsCalc);

		Subs::Format form;
		form.precision(8); form.general();	
		//std::cout << params.toString() << " " << form(chisq) << " - " << wdCalc << discCalc << donorCalc << bsCalc << std::endl;	
		return chisq;	
	}
	
private:
	LFIT::ThreeColour three;
	LFIT::CV CV;
	LFIT::SynthData model;
	LFIT::Params params, old_params;
	int colFit;
	int lastFit;
	bool pos;
	bool forceRecalc;
};


class LFunc{
  
 public:
  
  int neval;
  double chisq_min;

  //! Constructor
 LFunc(const LFIT::ThreeColour& three_, const LFIT::CV& myCV_, const LFIT::SynthData& myModel_, 
       const LFIT::Params& params_, int& colFit_, const bool& pos_,Subs::Buffer1D<double> dSteps_) : 
  three(three_),CV(myCV_),model(myModel_),params(params_),
    old_params(params_),colFit(colFit_),pos(pos_),paramSteps(dSteps_){
    neval = 0;
    wdCalc=true;
    donorCalc=true;
    discCalc=true;
    bsCalc=true;
    newBorn = true;
  }
  void force_recalc();
  void check_recalc();
  void updateParams(const LFIT::Params& paramsIn);
  void updateParams(const Subs::Array1D<double>& paramsIn);
  std::string asString();
  Subs::Array1D<double> get_ranges();
  Subs::Array1D<double> get_variables();
  std::vector< std::pair<double,double> > get_limits();
  void setParamRanges();
  int nvary();
  double fitModel();
  void lmcomp(Subs::Buffer2D<double>& alpha, Subs::Buffer1D<double>& beta, double& chisq);
  
 private:
    
  Subs::Buffer1D<double> paramSteps;
  LFIT::ThreeColour three;
  LFIT::CV CV;
  LFIT::SynthData model;
  LFIT::Params params, old_params;
  int colFit;
  int lastFit;
  bool pos, newBorn;
  bool wdCalc, discCalc, donorCalc, bsCalc;
};

int lmfit(LFunc& func, double& chisq, double& lambda, Subs::Buffer2D<double>& covar);

// adjust parameters to fit model using downhill simplex (amoeba)
void tweak(const std::vector<std::string>& args, LFIT::Params& params, 
		   const LFIT::ThreeColour& three, const LFIT::CV& myCV, 
		   const LFIT::SynthData& myModel, 
		   int& colFit,const bool& pos);


std::pair<int,double> run_amoeba(LFIT::Params& params, amFunc& func, double ftol, int nmax);

// convert parameters to units of binary separation
void convert(LFIT::Params params);

// print commands to terminal
void printcomms(const bool& pos, const bool& steep);

// rebin ultracam three-colour lightcurve
LFIT::ThreeColour rebin(const std::vector<std::string>& args, const LFIT::ThreeColour& three);

// terminal entry of model parameters
void edit(const std::vector<std::string>& args, LFIT::Params& params);

// plot lightcurve with option to overplot model (and residuals)
void plotLight(const LFIT::SynthData& model, 
			   const LFIT::Data& data, Subs::Plot& plot,
			   std::string autoscale, double& ymin, double& ymax);

// wrapper function to call plotLight
void light(const std::vector<std::string>& args, const LFIT::SynthData& model, const LFIT::Data& data);

// load a model from file
void load(const std::vector<std::string>& args, LFIT::Params& params);

// dump a model to file
void dump(const std::vector<std::string>& args, const LFIT::Params& params);

// write a mask of bad data points to file
void wmask(const std::vector<std::string>& args, const LFIT::ThreeColour& three);

// load a mask from file
void lmask(const std::vector<std::string>& args, LFIT::ThreeColour& three);

// interactively set mask
void smask(const std::vector<std::string>& args, LFIT::ThreeColour& three);

// bootstrap resample lightcurve
LFIT::ThreeColour bootstrap(const LFIT::ThreeColour& inData, const int& colour);

void errors_bootstrap(const std::vector<std::string>& args, LFIT::Params& params, 
					  LFIT::ThreeColour& three, LFIT::CV& myCV, LFIT::SynthData& myModel, 
					  int& colFit,const bool& pos);

void errors_sem(const std::vector<std::string>& args, LFIT::Params& params, 
					  LFIT::ThreeColour& three, LFIT::CV& myCV, LFIT::SynthData& myModel, 
					  int& colFit,const bool& pos);

void mcmc(const std::vector<std::string>& args, LFIT::Params& params,
	const LFIT::ThreeColour& three, const LFIT::CV& myCV,
	const LFIT::SynthData& myModel,int& colFit, const bool& pos);
