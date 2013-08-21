/*
 *  lfit_scale.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 19/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */
#include<string>
#include<cfloat>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_buffer2d.h"
#include "lfit.h"

using namespace std;

enum components {
  WD,
  DISC,
  BSPOT,
  DONOR
};

string enumToString(const int& enumIn){
  switch(enumIn){
  case WD:
    return "White Dwarf";
    break;
  case DISC:
    return "Disc";
    break;
  case DONOR:
    return "Donor";
    break;
  case BSPOT:
    return "Bright Spot";
    break;
  default:
    string err = "Enum not valid";
    throw err;
  }
}

void scale_model(LFIT::SynthData& model, const LFIT::Params& params, 
		 const Subs::Buffer1D<double>& prim,
		 const Subs::Buffer1D<double>& disc, 
		 const Subs::Buffer1D<double>& bspot, 
		 const Subs::Buffer1D<double>& donor){
  for(int i=0; i<model.phase.size(); ++i){
    model.wflux[i] = prim[i] * params.wFlux;
    model.bflux[i] = bspot[i] * params.bFlux;
    model.dflux[i] = disc[i] * params.dFlux;
    model.rflux[i] = params.rFlux * donor[i];
    model.flux[i]  = model.wflux[i]+model.bflux[i]+model.dflux[i]+model.rflux[i];
  }
}

double fastscale(const vector<string>& args,LFIT::Params& params, 
		 const LFIT::ThreeColour& three, LFIT::CV& myCV, LFIT::SynthData& model, 
		 int& colFit, const bool& batch, const bool& pos,
		 const bool& wdCalc, const bool& donorCalc, const bool& discCalc, 
		 const bool& bsCalc){
    
  // argument pos dictates whether we try and hold it positive or not
  // batch dictates whether we need to interactively prompt for colour
  // both arguments are optional and default to to FALSE
	
  // the xxCalc arguments are also optional, and they tell fastscale whether to
  // recalculate the respective arrays or not, e.g if wdCalc=FALSE then the
  // white dwarf has not changed, nor has the times at which is should be calculated,
  // so the old white dwarf model can be used for speed
	
  // scales current model by least squares fit, returns chisq;
	
  if(!batch){
    Subs::Input input(args,"LFIT_ENV",".lfit");
    input.sign_in("colour",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.get_value("colour", colFit, 10, 1, 3, "lightcurve colour to fit");
    colFit -= 1;
  }
  LFIT::Data data = three[colFit];
    
  // calculate components (objects used will prevent unnecesary calculations)
  myCV.donor.tweak(params.q);
  myCV.prim.tweak(params.rwd,params.ulimb);
  myCV.bspot.tweak(params);
  myCV.disc.tweak(params.q,params.rwd,params.rd,params.dExp);    

  int ndata = data.size();
  int nfunc = 4; // one for each CV component we wish to scale
  const int maxfunc=4; 
  Subs::Buffer1D<double> y(ndata);	
  // buffers holding the model fluxes for each component remain from previous call,
  // in case they don't need recalculating
  bool static first=true;
  static Subs::Buffer1D<double> prim(ndata),disc(ndata),bspot(ndata),donor(ndata);
  Subs::Buffer1D<float> e(ndata);
  Subs::Buffer1D<double> coeff(nfunc);
  Subs::Buffer2D<double> tofit(ndata,nfunc), covar(nfunc,nfunc);
  double wFlux, wFluxErr, dFlux, dFluxErr, bFlux, bFluxErr, rFlux, rFluxErr;
	
  // resize model if necessary
  if(model.phase.size() != ndata){
    model.phase.resize(data.size());
    model.dflux.resize(data.size());
    model.bflux.resize(data.size());
    model.wflux.resize(data.size());
    model.rflux.resize(data.size());
    model.flux.resize(data.size());
  }
	
  // calculate model lightcurves for these parameters
  // this is done once, and the results stored in memory, for speed
  for(int i=0; i<ndata; i++){
    model.phase[i] = data[i].time; 
  }
	
  bool doWd, doDisc, doBS, doDonor;
  doWd = wdCalc; doDisc = discCalc; doBS = bsCalc; doDonor=donorCalc;
  if(first){
    doWd=true;
    doDisc=true;
    doBS=true;
    doDonor=true;
    first=false;
  }
	
  //#pragma omp parallel for shared(prim,disc,bspot,donor,data)
  for(int i=0; i< ndata; i++){
    // set obervation time
    double x = data[i].time - params.phi0;
    double width = data[i].expose;
		
    // calculate lightcurves of components once at start of routine
    if(doWd){prim[i]  = myCV.prim.calcFlux(params.q,x,width,params.incl);}
    if(doDisc){disc[i]  = myCV.disc.calcFlux(params.q,x,width,params.incl);}
    if(doBS){bspot[i] = myCV.bspot.calcFlux(params.q,x,width,params.incl);}
    if(doDonor){donor[i] = myCV.donor.calcFlux(x,width,params.incl);}
  }

  for(int i=0; i<ndata; i++){
    // put into matrix for LLSQR fitting
    tofit[i][0]=prim[i];
    tofit[i][1]=disc[i];
    tofit[i][2]=bspot[i];
    tofit[i][3]=donor[i];
		
    // set up data and error arrays for LLSQR fitting
    y[i]= data[i].flux;
    if (data[i].bad){
      e[i]=-fabs(data[i].ferr);
    }else{
      e[i]=data[i].ferr;
    }
  }
  //cout<< "components filled" << endl;		
  // establish map of CV components
  map<int,Subs::Buffer1D<double> > componentMap;
  componentMap[0]=prim;
  componentMap[1]=disc;
  componentMap[2]=bspot;
  componentMap[3]=donor;
  
  bool success;
  // attempt fit letting all components vary freely
  try{
    Subs::llsqr(ndata,y,e,nfunc,tofit,coeff,covar);
    success = true;
  }catch(const string err){
    //cout << "Error in LLSQR: " << err << endl;
    success = false;
  }
	
  // save results
  wFlux = coeff[0]; wFluxErr = sqrt(covar[0][0]);
  dFlux = coeff[1]; dFluxErr = sqrt(covar[1][1]);
  bFlux = coeff[2]; bFluxErr = sqrt(covar[2][2]);
  rFlux = coeff[3]; rFluxErr = sqrt(covar[3][3]);
	
  // set scaling of model
  params.wFlux = wFlux; 
  params.dFlux = dFlux; 
  params.bFlux = bFlux; 
  params.rFlux = rFlux; 
	
  // recalculate model with new scaling
  scale_model(model,params,prim,disc,bspot,donor);
	
  // are any of the parameters negative?
  bool negCoeffs=false;
  for(int i=0; i<4; i++){
    negCoeffs = negCoeffs || (coeff[i]<0.0);
  }

  // have we insisted on positivity, and do we need to enforce it?
  // alternatively, did we fail to fit with all four parameters?
  int barfed = 0;
  if((pos && negCoeffs) || !success){    

    if(!batch) { 
      std::cout << "Trying to enforce positivity..."<<std::endl;
      cout << "   Chisq from free fit: " << chisq(data,model) << endl;
    }
		
    // set each parameter to zero in turn and choose the combination which
    // minimises chi-squared
    nfunc=3;
    float chimin=FLT_MAX;
    tofit.resize(ndata,nfunc); covar.resize(nfunc,nfunc);
    coeff.resize(nfunc);
		
    for (int heldFixed=0; heldFixed<maxfunc; heldFixed++){

      // iterate over components, setting each one to zero in turn
      map<int,Subs::Buffer1D<double> >::iterator iter;
      int toFill=0;
      for(iter = componentMap.begin(); iter != componentMap.end(); iter++){
	if(iter->first != heldFixed){
	  for(int i=0; i<ndata; i++) tofit[i][toFill]=iter->second[i];
	  toFill++;
	}else{
	  if(!batch)
	    cout << "   " << enumToString(iter->first) << " set to zero" << endl;
	}
      }
			
      // fit with this component set to zero
      try{
	Subs::llsqr(ndata,y,e,nfunc,tofit,coeff,covar);
      }catch(const string err){
	//cout << "Error in LLSQR (force pos): " << err << endl;
	barfed += 1;
      }
			
      // set scaling of model from LLSQR results
      int fromFill=0;
      double wErr=0.0, dErr=0.0, bErr=0.0, rErr=0.0;
      for(iter = componentMap.begin(); iter != componentMap.end(); iter++){
	components comp = components(iter->first);
	switch(comp){
	case WD:
	  if(iter->first != heldFixed){
	    params.wFlux = coeff[fromFill];
	    wErr         = covar[fromFill][fromFill];
	    fromFill++;
	  }else{
	    wErr         = 0.0;
	    params.wFlux = 0.0;
	  }
	  break;
	case DISC:
	  if(iter->first != heldFixed){
	    params.dFlux = coeff[fromFill];
	    dErr         = covar[fromFill][fromFill];
	    fromFill++;
	  }else{
	    dErr         = 0.0;
	    params.dFlux = 0.0;
	  }
	  break;
	case DONOR:
	  if(iter->first != heldFixed){
	    params.rFlux = coeff[fromFill];
	    rErr         = covar[fromFill][fromFill];
	    fromFill++;
	  }else{
	    rErr         = 0.0;
	    params.rFlux = 0.0;
	  }
	  break;
	case BSPOT:
	  if(iter->first != heldFixed){
	    params.bFlux = coeff[fromFill];
	    bErr         = covar[fromFill][fromFill];
	    fromFill++;
	  }else{
	    bErr         = 0.0;
	    params.bFlux = 0.0;
	  }
	  break;
	default:
	  string err = "enum not a valid CV component";
	  throw err;
	}
      }
      
      scale_model(model,params,prim,disc,bspot,donor);	
      float thisChi = chisq(data,model);
      if(!batch) cout << "   Chisq: " << thisChi << endl << endl;
      if(thisChi<chimin){
	chimin=thisChi;
	// update fluxes and errors
	wFlux = params.wFlux; wFluxErr=wErr;
	rFlux = params.rFlux; rFluxErr=rErr;
	bFlux = params.bFlux; bFluxErr=bErr;
	dFlux = params.dFlux; dFluxErr=dErr;
      }
			 
    }
		
  }
  if(barfed != 4) success=true;

  if(!success){
    string err = "Scaling failed: LLSQR failed";
    throw err;
  }

  // write out results
  params.wFlux = wFlux; 
  params.dFlux = dFlux; 
  params.bFlux = bFlux; 
  params.rFlux = rFlux; 
  if(!batch){
    Subs::Format form;
    form.precision(8); form.general();
    cout << "White Dwarf flux = " << form(params.wFlux) << " +/- " << form(wFluxErr) << endl;
    cout << "Disc flux        = " << form(params.dFlux) << " +/- " << form(dFluxErr) << endl;
    cout << "Bright Spot flux = " << form(params.bFlux) << " +/- " << form(bFluxErr) << endl;
    cout << "Donor Flux       = " << form(params.rFlux) << " +/- " << form(rFluxErr) << endl << endl;
  }
	
  // scale best model and return
  scale_model(model,params,prim,disc,bspot,donor);
  return chisq(data,model);
}

double scale(const vector<string>& args,LFIT::Params& params, 
	     const LFIT::ThreeColour& three, LFIT::CV& myCV, LFIT::SynthData& model, 
	     int& colFit, const bool& batch, const bool& pos){
    
  // argument pos dictates whether we try and hold it positive or not
  // batch dictates whether we need to interactively prompt for colour
  // both arguments are optional and default to to FALSE
	
		
  // scales current model by least squares fit, returns chisq;
  if(!batch){
    Subs::Input input(args,"LFIT_ENV",".lfit");
    input.sign_in("colour",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.get_value("colour", colFit, 10, 1, 3, "lightcurve colour to fit");
    colFit -= 1;
  }
  LFIT::Data data = three[colFit];
    
  // calculate components (objects used should prevent unnecesary calculations)
  myCV.donor.tweak(params.q);
  myCV.prim.tweak(params.rwd,params.ulimb);
  myCV.bspot.tweak(params);
  myCV.disc.tweak(params.q,params.rwd,params.rd,params.dExp);    

  int ndata = data.size();
  int nfunc = 4; // one for each CV component we wish to scale
  const int maxfunc=4; 
  Subs::Buffer1D<double> y(ndata);	
  Subs::Buffer1D<double> prim(ndata),disc(ndata),bspot(ndata),donor(ndata);
  Subs::Buffer1D<float> e(ndata);
  Subs::Buffer1D<double> coeff(nfunc);
  Subs::Buffer2D<double> tofit(ndata,nfunc), covar(nfunc,nfunc);
  double wFlux, wFluxErr, dFlux, dFluxErr, bFlux, bFluxErr, rFlux, rFluxErr;
	
  // resize model if necessary
  if(model.phase.size() != ndata){
    model.phase.resize(data.size());
    model.dflux.resize(data.size());
    model.bflux.resize(data.size());
    model.wflux.resize(data.size());
    model.rflux.resize(data.size());
    model.flux.resize(data.size());
  }
	
  // calculate model lightcurves for these parameters
  // this is done once, and the results stored in memory, for speed
  for(int i=0; i<ndata; i++){
    model.phase[i] = data[i].time;
  }	
  for(int i=0; i< ndata; i++){
    // set obervation time
    double x = data[i].time - params.phi0;
    double width = data[i].expose;
		
    // calculate lightcurves of components once at start of routine
    prim[i]  = myCV.prim.calcFlux(params.q,x,width,params.incl);
    disc[i]  = myCV.disc.calcFlux(params.q,x,width,params.incl);
    bspot[i] = myCV.bspot.calcFlux(params.q,x,width,params.incl);
    donor[i] = myCV.donor.calcFlux(x,width,params.incl);
		
    // put into matrix for LLSQR fitting
    tofit[i][0]=prim[i];
    tofit[i][1]=disc[i];
    tofit[i][2]=bspot[i];
    tofit[i][3]=donor[i];
		
    // set up data and error arrays for LLSQR fitting
    y[i]= data[i].flux;
    if (data[i].bad){
      e[i]=-fabs(data[i].ferr);
    }else{
      e[i]=data[i].ferr;
    }
  }
  //cout<< "components filled" << endl;		
  // establish map of CV components
  map<int,Subs::Buffer1D<double> > componentMap;
  componentMap[0]=prim;
  componentMap[1]=disc;
  componentMap[2]=bspot;
  componentMap[3]=donor;
	
  // attempt fit letting all components vary freely
  try{
    Subs::llsqr(ndata,y,e,nfunc,tofit,coeff,covar);
  }catch(const string err){
    cout << "Error in LLSQR: " << err << endl;
  }
	
  // save results
  wFlux = coeff[0]; wFluxErr = sqrt(covar[0][0]);
  dFlux = coeff[1]; dFluxErr = sqrt(covar[1][1]);
  bFlux = coeff[2]; bFluxErr = sqrt(covar[2][2]);
  rFlux = coeff[3]; rFluxErr = sqrt(covar[3][3]);
	
  // set scaling of model
  params.wFlux = wFlux; 
  params.dFlux = dFlux; 
  params.bFlux = bFlux; 
  params.rFlux = rFlux; 
	
  // recalculate model with new scaling
  scale_model(model,params,prim,disc,bspot,donor);
	
  // are any of the parameters negative?
  bool negCoeffs=false;
  for(int i=0; i<4; i++){
    negCoeffs = negCoeffs || (coeff[i]<0.0);
  }

  // have we insisted on positivity, and do we need to enforce it?
  if(pos && negCoeffs){
		
    if(!batch) { 
      std::cout << "Trying to enforce positivity..."<<std::endl;
      cout << "   Chisq from free fit: " << chisq(data,model) << endl;
    }
		
    // set each parameter to zero in turn and choose the combination which
    // minimises chi-squared
    nfunc=3;
    float chimin=FLT_MAX;
    tofit.resize(ndata,nfunc); covar.resize(nfunc,nfunc);
    coeff.resize(nfunc);
		
    for (int heldFixed=0; heldFixed<maxfunc; heldFixed++){

      // iterate over components, setting each one to zero in turn
      map<int,Subs::Buffer1D<double> >::iterator iter;
      int toFill=0;
      for(iter = componentMap.begin(); iter != componentMap.end(); iter++){
	if(iter->first != heldFixed){
	  for(int i=0; i<ndata; i++) tofit[i][toFill]=iter->second[i];
	  toFill++;
	}else{
	  if(!batch)
	    cout << "   " << enumToString(iter->first) << " set to zero" << endl;
	}
      }
			
      // fit with this component set to zero
      try{
	Subs::llsqr(ndata,y,e,nfunc,tofit,coeff,covar);
      }catch(const string err){
	cout << "Error in LLSQR (force pos): " << err << endl;
      }
			
      // set scaling of model from LLSQR results
      int fromFill=0;
      double wErr=0.0, dErr=0.0, bErr=0.0, rErr=0.0;
      for(iter = componentMap.begin(); iter != componentMap.end(); iter++){
	components comp = components(iter->first);
	switch(comp){
	case WD:
	  if(iter->first != heldFixed){
	    params.wFlux = coeff[fromFill];
	    wErr         = covar[fromFill][fromFill];
	    fromFill++;
	  }else{
	    wErr         = 0.0;
	    params.wFlux = 0.0;
	  }
	  break;
	case DISC:
	  if(iter->first != heldFixed){
	    params.dFlux = coeff[fromFill];
	    dErr         = covar[fromFill][fromFill];
	    fromFill++;
	  }else{
	    dErr         = 0.0;
	    params.dFlux = 0.0;
	  }
	  break;
	case DONOR:
	  if(iter->first != heldFixed){
	    params.rFlux = coeff[fromFill];
	    rErr         = covar[fromFill][fromFill];
	    fromFill++;
	  }else{
	    rErr         = 0.0;
	    params.rFlux = 0.0;
	  }
	  break;
	case BSPOT:
	  if(iter->first != heldFixed){
	    params.bFlux = coeff[fromFill];
	    bErr         = covar[fromFill][fromFill];
	    fromFill++;
	  }else{
	    bErr         = 0.0;
	    params.bFlux = 0.0;
	  }
	  break;
	default:
	  string err = "enum not a valid CV component";
	  throw err;
	}
      }
			
      scale_model(model,params,prim,disc,bspot,donor);	
      float thisChi = chisq(data,model);
      if(!batch) cout << "   Chisq: " << thisChi << endl << endl;
      if(thisChi<chimin){
	chimin=thisChi;
	// update fluxes and errors
	wFlux = params.wFlux; wFluxErr=wErr;
	rFlux = params.rFlux; rFluxErr=rErr;
	bFlux = params.bFlux; bFluxErr=bErr;
	dFlux = params.dFlux; dFluxErr=dErr;
      }
			 
    }
		
  }
	
  // write out results
  params.wFlux = wFlux; 
  params.dFlux = dFlux; 
  params.bFlux = bFlux; 
  params.rFlux = rFlux; 
  if(!batch){
    Subs::Format form;
    form.precision(8); form.general();
    cout << "White Dwarf flux = " << form(params.wFlux) << " +/- " << form(wFluxErr) << endl;
    cout << "Disc flux        = " << form(params.dFlux) << " +/- " << form(dFluxErr) << endl;
    cout << "Bright Spot flux = " << form(params.bFlux) << " +/- " << form(bFluxErr) << endl;
    cout << "Donor Flux       = " << form(params.rFlux) << " +/- " << form(rFluxErr) << endl << endl;
  }
	
  // scale best model and return
  scale_model(model,params,prim,disc,bspot,donor);
  return chisq(data,model);
}
