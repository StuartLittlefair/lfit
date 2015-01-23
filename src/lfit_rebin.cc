/*
 *  lfit_rebin.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 19/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */
#include <cfloat>
#include<climits>
#include "trm/subs.h"
#include "trm/input.h"
#include "lfit.h"


LFIT::Data binUp(const LFIT::Data& dIn, const double& xstart, const double& xend, const int& nbin, const char& weight);

LFIT::ThreeColour rebin(const std::vector<std::string>& args, const LFIT::ThreeColour& three){

  Subs::Input input(args,"LFIT_ENV",".lfit");
  input.sign_in("nbin",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
  input.sign_in("xstart",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
  input.sign_in("xend",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
  input.sign_in("weight",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
  input.sign_in("yerrors",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
  int nbin;
  input.get_value("nbin", nbin, 400, 1, INT_MAX, "number of bins");
  double xstart, xend;
  input.get_value("xstart",xstart,0.9,0.5,0.999, "start phase");
  input.get_value("xend",xend    ,1.15,1.0001,1.5, "start phase");
  char weight;
  input.get_value("weight", weight, 'u', "uUvV", "U(niform) or inverse V(ariance) weighting?");
  weight = std::toupper(weight);
  if(weight != 'U' && weight != 'V'){
    string error = "Only valid replies are 'U' for uniform and 'V' for variance.";
    throw error;
  }
  
  LFIT::ThreeColour theMagicNumber;
  for(int i=0; i<3; i++){
    LFIT::Data data = three[i];
    if(xstart < data.get_xmin() || xend > data.get_xmax()){
      string error = "Requested start or end outside data limits of " + 
	Subs::str(data.get_xmin()) + " to " + Subs::str(data.get_xmax());
      throw error;
    }
    // OK, carry out work
    theMagicNumber[i] = binUp(data,xstart,xend,nbin,weight);
    data.resetLimits();
  }
  std::cout << "Red lcurve computed with " << theMagicNumber[0].size() << " data points\n";
  std::cout << "Green lcurve computed with " << theMagicNumber[1].size() << " data points\n";
  std::cout << "Blue lcurve computed with " << theMagicNumber[2].size() << " data points\n";
  return theMagicNumber;
}

LFIT::Data binUp(const LFIT::Data& data, const double& xstart, const double& xend, const int& nbin, const char& weight){
  
  
  LFIT::Datum datum;
  LFIT::Data newdata;
  newdata.clear();
  
  int start = 0;
  for(int j=0; j<nbin; j++){
    double x1 = xstart + (xend-xstart)*j/std::max(nbin-1,1);
    double x2 = xstart + (xend-xstart)*(j+1)/std::max(nbin-1,1);
    if(x1 > x2) std::swap(x1, x2);
    
    double sumx = 0., sumy = 0., sumex = 0., sumey = 0., sumw = 0., sumwsq = 0., w, sumysq = 0.;
    int k, np = 0;
    
    // assumes data is monotonic and increases in time
    if(j == 0) start = 0;
    for(k=start; k<data.size() && data[k].time<x2; k++){
      if(data[k].time >= x1){
		if(! data[k].bad){
	  
		  if(weight == 'V'){
			w = 1.0/Subs::sqr(data[k].ferr);
		  }else{
			w = 1.0;
		  }

		  sumey += Subs::sqr(w*data[k].ferr);
		  sumx  += w*data[k].time;
		  sumy  += w*data[k].flux;
		  sumysq += w*data[k].flux*data[k].flux;
		  sumw  += w;
		  sumwsq+= w*w;
		  np++;
		}
      }
    }
    start = k;
    
    if(np > 0){
		// store point
		datum.time = sumx/sumw;
		datum.flux = sumy/sumw;
		datum.expose = x2-x1;
		//set yerr from scatter in bin
		if(np>10){
			datum.ferr = sqrt(sumwsq*(sumysq-sumy*sumy/sumw)/sumw)/sumw;
		}else{
			datum.ferr   = sqrt(sumey)/sumw;
		}
		if (datum.ferr != datum.ferr){
			std::cout << datum.time << " " << datum.flux << " " << datum.ferr << " " << np << std::endl;
			std::cout << sumey << " " << sumw << std::endl;
		}
		newdata.push_back(datum);
    }    
  }
  return newdata;
}

