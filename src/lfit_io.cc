/*
 *  lfit_io.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 19/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */
#include <cfloat>
#include <iostream>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/input.h"
#include "lfit.h"

using namespace std;

void convert(LFIT::Params params){

  // convert from units of L1 to a
  Subs::Format form;
  form.precision(8); form.general();
	
  double xl1 = Roche::xl1(params.q);
  cout << "L1 distance/a = " << form(xl1) << endl;
  if(params.hasErrors()){
    cout << endl << "Current fit parameters are:" << endl <<
      "Mass ratio Q = " << form(params.q) << " +/- " << form(params.qErr) << endl <<
      "Delta Phi    = " << form(params.dphi) << " +/- " << form(params.dphiErr) << endl <<
      "Disc radius  = " << form(params.rd*xl1) << " +/- " << form(params.rdErr*xl1) << endl <<
      "Limb darkening = " << form(params.ulimb) << " +/- " << form(params.ulimbErr) << endl <<
      "White dwarf radius = " << form(params.rwd*xl1) << " +/- " << form(params.rwdErr*xl1) << endl <<
      "Bright spot scale = " << form(params.bsScale*xl1) << " +/- " << form(params.bsScaleErr*xl1) << endl;
    if (!params.complexSpot){
      cout << "Azimuth of spot = " << form(params.bsAz) << " +/- " << form(params.bsAzErr) << endl <<
	"Isotropic Fraction = " << form(params.bsFrac) << " +/- " << form(params.bsFracErr) << endl <<
	"Disc Exponent = " << form(params.dExp) << " +/- " << form(params.dExpErr) << endl <<
	"Phase offset = " << form(params.phi0) << " +/- " << form(params.phi0Err) << endl << endl;
    }else{
      cout << "Azimuth of spot = " << form(params.bsAz) << " +/- " << form(params.bsAzErr) << endl <<
	"Isotropic Fraction = " << form(params.bsFrac) << " +/- " << form(params.bsFracErr) << endl <<
	"BS Exponent1 = " << form(params.bsExp1) << " +/- " << form(params.bsExp1Err) << endl <<
	"BS Exponent2 = " << form(params.bsExp2) << " +/- " << form(params.bsExp2Err) << endl <<
	"BS Tilt = " << form(params.bsTilt) << " +/- " << form(params.bsTiltErr) << endl <<
	"BS Yaw = " << form(params.bsYaw) << " +/- " << form(params.bsYawErr) << endl <<
	"Disc Exponent = " << form(params.dExp) << " +/- " << form(params.dExpErr) << endl <<
	"Phase offset = " << form(params.phi0) << " +/- " << form(params.phi0Err) << endl << endl;				
    }		
  }else{
    cout << endl << "Current fit parameters are:" << endl <<
      "Mass ratio Q = " << form(params.q) << endl <<
      "Delta Phi    = " << form(params.dphi) << endl <<
      "Disc radius  = " << form(params.rd*xl1) << endl <<
      "Limb darkening = " << form(params.ulimb) << endl <<
      "White dwarf radius = " << form(params.rwd*xl1) << endl <<
      "Bright spot scale = " << form(params.bsScale*xl1) << endl;
    if (!params.complexSpot){
      cout << "Azimuth of spot = " << form(params.bsAz)  << endl <<
	"Isotropic Fraction = " << form(params.bsFrac)  << endl <<
	"Disc Exponent = " << form(params.dExp)  << endl <<
	"Phase offset = " << form(params.phi0)  << endl << endl;
    }else{
      cout << "Azimuth of spot = " << form(params.bsAz)  << endl <<
	"Isotropic Fraction = " << form(params.bsFrac)  << endl <<
	"BS Exponent1 = " << form(params.bsExp1)  << endl <<
	"BS Exponent2 = " << form(params.bsExp2)  << endl <<
	"BS Tilt = " << form(params.bsTilt)  << endl <<
	"BS Yaw = " << form(params.bsYaw)  << endl <<
	"Disc Exponent = " << form(params.dExp)  << endl <<
	"Phase offset = " << form(params.phi0)  << endl << endl;				
    }
  }
}

void printcomms(const bool& pos, const bool& complex){
  std::vector<std::string> commands;
  commands.push_back("REBIN   - rebin data");
	
  string onoff=" (Off)";
  if (pos) onoff=" (On)";
  commands.push_back("POS     - enforce positivity"+onoff);
  onoff = " (Off)";
  if(complex) onoff=" (On)";
  commands.push_back("COMPLEX   - use complex brightspot model"+onoff);
  commands.push_back("LMASK   - load mask of bad phases");
  commands.push_back("WMASK   - write mask of bad phases");
  commands.push_back("SMASK   - set mask interactively");
  commands.push_back("LOAD    - load file of fit parameters");
  commands.push_back("DUMP    - dump file of fit parameters");
  commands.push_back("WRITE   - write out current data and model");
  commands.push_back("EDIT    - edit current model parameters");
  commands.push_back("SHOW    - show current model parameters");
  commands.push_back("CONVERT - convert parameters to units of separation");
  commands.push_back("SCALE   - rescale model to fit data with current params");
  commands.push_back("TWEAK   - adjust model parameters by fitting to data");
  commands.push_back("MCMC    - find best parameters and errors using MCMC (slow)");
  commands.push_back("BOOT    - find errors by bootstrapping method (slow)");
  commands.push_back("PLOT    - plot current model & data");
  commands.push_back("QUIT    - give up");
  cout << endl;
  for(size_t i=0; i<commands.size(); i++){
    cout << commands[i] << endl;
  }
  cout << endl;
}

int LFIT::findcomm(std::string& trial){

  vector<string> commands;
  commands.push_back("quit");
  commands.push_back("lmask");
  commands.push_back("wmask");
  commands.push_back("smask");
  commands.push_back("plot");
  commands.push_back("pos");
  commands.push_back("write");
  commands.push_back("complex");
  commands.push_back("convert");
  commands.push_back("show");
  commands.push_back("tweak");
  commands.push_back("edit");
  commands.push_back("load");
  commands.push_back("rebin");
  commands.push_back("dump");
  commands.push_back("scale");
  commands.push_back("fscale");
  commands.push_back("mcmc");
  commands.push_back("boot");

  trial = Subs::tolower(trial);
  int nmatch=0;
  std::string save;
  for(size_t i=0; i<commands.size();i++){
    if(commands[i].find(trial) == 0){
      save = commands[i];
      nmatch++;
    }
  } 	
  if(nmatch==1) trial=save;
  return nmatch;
}

void edit(const vector<string>& args, LFIT::Params& params){
  Subs::Input input(args,"LFIT_ENV",".lfit");
  input.sign_in("q",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  input.sign_in("dphi",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  input.sign_in("rd",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  input.sign_in("rwd",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  input.sign_in("ulimb",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  input.sign_in("scale",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  input.sign_in("az",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  input.sign_in("frac",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  if(params.complexSpot){
    input.sign_in("exp1",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("exp2",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("tilt",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("yaw",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  }
  input.sign_in("dexp",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  input.sign_in("phi0",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    
  try{
    double qtmp, dphitmp;
    input.get_value("q",qtmp,0.1,0.0,DBL_MAX,"mass ratio q");
    input.get_value("dphi",dphitmp,0.04,0.0,DBL_MAX,"phase width delta phi");
    double incltmp = finddeg(qtmp,dphitmp);
    params.q = qtmp;
    params.dphi = dphitmp;
    params.incl = incltmp;
      
    input.get_value("rd",params.rd,0.3,0.0,DBL_MAX,"radius of disc (L1 units)");
    input.get_value("rwd",params.rwd,0.02,0.0,DBL_MAX,"radius of white dwarf (L1 units)");
    input.get_value("ulimb",params.ulimb,0.4,0.0,DBL_MAX,"limb darkening of white dwarf");
    input.get_value("scale",params.bsScale,0.02,0.0,DBL_MAX,"bright spot scale");
    input.get_value("az",params.bsAz,150.,0.0,180.0,"azimuth of spot");
    if(!params.complexSpot){
      input.get_value("frac",params.bsFrac,0.2,0.0,DBL_MAX,"isotropic fraction of spot");
    }else{
      input.get_value("frac",params.bsFrac,0.2,0.0,DBL_MAX,"fraction of light in planar bright spot strip");
      input.get_value("exp1",params.bsExp1,2.0,0.0,DBL_MAX,"Exponent1 of bright spot brightness profile");
      input.get_value("exp2",params.bsExp2,2.0,0.0,DBL_MAX,"Exponent2 of bright spot brightness profile");
      input.get_value("tilt",params.bsTilt,10.,0.0,180.0,"Tilt of secondary bright spot strip");
      input.get_value("yaw",params.bsYaw,10.,-90.0,90.0,"Yaw of secondary bright spot strip");
    }
    input.get_value("dexp",params.dExp,1.0,-DBL_MAX,DBL_MAX,"disc exponent");
    input.get_value("phi0",params.phi0,0.0,-DBL_MAX,DBL_MAX,"phase offset");
  }catch(const string& err){
    std::cout << err << std::endl;
  }
}


void show(LFIT::Params& params){
  
  Subs::Format form;
  
  if(params.hasErrors()){
    cout << endl << "Current fit parameters are:" << endl <<
      "Mass ratio Q = " << form(params.q) << " +/- " << form(params.qErr) << endl <<
      "Delta Phi    = " << form(params.dphi) << " +/- " << form(params.dphiErr) << endl <<
      "Disc radius  = " << form(params.rd) << " +/- " << form(params.rdErr) << endl <<
      "Limb darkening = " << form(params.ulimb) << " +/- " << form(params.ulimbErr) << endl <<
      "White dwarf radius = " << form(params.rwd) << " +/- " << form(params.rwdErr) << endl <<
      "Bright spot scale = " << form(params.bsScale) << " +/- " << form(params.bsScaleErr) << endl;
    if (!params.complexSpot){
      cout << "Azimuth of spot = " << form(params.bsAz) << " +/- " << form(params.bsAzErr) << endl <<
	"Isotropic Fraction = " << form(params.bsFrac) << " +/- " << form(params.bsFracErr) << endl <<
	"Disc Exponent = " << form(params.dExp) << " +/- " << form(params.dExpErr) << endl <<
	"Phase offset = " << form(params.phi0) << " +/- " << form(params.phi0Err) << endl << endl;
    }else{
      cout << "Azimuth of spot = " << form(params.bsAz) << " +/- " << form(params.bsAzErr) << endl <<
	"Isotropic Fraction = " << form(params.bsFrac) << " +/- " << form(params.bsFracErr) << endl <<
	"BS Exponent1 = " << form(params.bsExp1) << " +/- " << form(params.bsExp1Err) << endl <<
	"BS Exponent2 = " << form(params.bsExp2) << " +/- " << form(params.bsExp2Err) << endl <<
	"BS Tilt = " << form(params.bsTilt) << " +/- " << form(params.bsTiltErr) << endl <<
	"BS Yaw = " << form(params.bsYaw) << " +/- " << form(params.bsYawErr) << endl <<
	"Disc Exponent = " << form(params.dExp) << " +/- " << form(params.dExpErr) << endl <<
	"Phase offset = " << form(params.phi0) << " +/- " << form(params.phi0Err) << endl << endl;				
    }		
  }else{
    cout << " NO ERRORS: " << endl;
    cout << endl << "Current fit parameters are:" << endl <<
      "Mass ratio Q = " << form(params.q) << endl <<
      "Delta Phi    = " << form(params.dphi) << endl <<
      "Disc radius  = " << form(params.rd) << endl <<
      "Limb darkening = " << form(params.ulimb) << endl <<
      "White dwarf radius = " << form(params.rwd) << endl <<
      "Bright spot scale = " << form(params.bsScale) << endl;
    if (!params.complexSpot){
      cout << "Azimuth of spot = " << form(params.bsAz)  << endl <<
	"Isotropic Fraction = " << form(params.bsFrac)  << endl <<
	"Disc Exponent = " << form(params.dExp)  << endl <<
	"Phase offset = " << form(params.phi0)  << endl << endl;
    }else{
      cout << "Azimuth of spot = " << form(params.bsAz)  << endl <<
	"Isotropic Fraction = " << form(params.bsFrac)  << endl <<
	"BS Exponent1 = " << form(params.bsExp1)  << endl <<
	"BS Exponent2 = " << form(params.bsExp2)  << endl <<
	"BS Tilt = " << form(params.bsTilt)  << endl <<
	"BS Yaw = " << form(params.bsYaw)  << endl <<
	"Disc Exponent = " << form(params.dExp)  << endl <<
	"Phase offset = " << form(params.phi0)  << endl << endl;				
    }
    cout << "White Dwarf flux = " << form(params.wFlux)  << endl <<
      "Disc flux        = " << form(params.dFlux) << endl <<
      "Bright Spot flux = " << form(params.bFlux) << endl <<
      "Donor Flux       = " << form(params.rFlux) << endl << endl;
  }
  
}

void load(const vector<string>& args, LFIT::Params& params){
  
  Subs::Input input(args,"LFIT_ENV",".lfit");
  //sign in input variables
  input.sign_in("filename",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  
  // get input variables
  string filename;
  input.get_value("filename", filename, "params.dat", "file to read from");
  
  ifstream fin(filename.c_str());
  if(!fin){
    string error("Failed to open file for reading: ");
    error.append(filename);
    throw error;
  }
  string line;
  Subs::Buffer1D<float> pIn;
  while(fin){
    getline(fin, line);
    string::size_type loc = line.find("=",0);
    if(loc != string::npos ){
      string substr = line.substr(++loc);
      pIn.push_back( Subs::string_to_double(substr));
    }
  }
  if(pIn.size() > 10){
    params.complexSpot=true;
  }
  
  int pnum = 0;
  params.q       = pIn[pnum++];
  params.dphi    = pIn[pnum++];
  params.rd      = pIn[pnum++];
  params.ulimb   = pIn[pnum++];
  params.rwd     = pIn[pnum++];
  params.bsScale = pIn[pnum++];
  params.bsAz    = pIn[pnum++];
  params.bsFrac   = pIn[pnum++];
  if(params.complexSpot){
    params.bsExp1 = pIn[pnum++];
    params.bsExp2 = pIn[pnum++];
    params.bsTilt = pIn[pnum++];
    params.bsYaw = pIn[pnum++];
  }
  params.dExp    = pIn[pnum++];
  params.phi0    = pIn[pnum++];
  params.incl = finddeg(params.q,params.dphi);
  show(params);
}	

void dump(const vector<string>& args, const LFIT::Params& params){
  
  Subs::Input input(args,"LFIT_ENV",".lfit");
  //sign in input variables
  input.sign_in("filename",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  
  // get input variables
  string filename;
  input.get_value("filename", filename, "params.dat", "file to write to");
  
  ofstream fout(filename.c_str());
  if(!fout){
    string error("Failed to open file for writing: ");
    error.append(filename);
    throw error;
  }
  
  Subs::Format form;
  form.precision(8);
  form.general();
  fout << endl << "Current fit parameters are:" << endl <<
    "Mass ratio Q = " << form(params.q) << endl <<
    "Delta Phi    = " << form(params.dphi) << endl <<
    "Disc radius  = " << form(params.rd) << endl <<
    "Limb darkening = " << form(params.ulimb) << endl <<
    "White dwarf radius = " << form(params.rwd) << endl <<
    "Bright spot scale = " << form(params.bsScale) << endl <<
    "Azimuth of spot = " << form(params.bsAz) << endl;
  if(!params.complexSpot){
    fout << "Isotropic Fraction = " << form(params.bsFrac) << endl;
  }else{
    fout << "Isotropic Fraction = " << form(params.bsFrac)  << endl <<
      "BS Exponent1 = " << form(params.bsExp1)  << endl <<
      "BS Exponent2 = " << form(params.bsExp2)  << endl <<
      "BS Tilt = " << form(params.bsTilt)  << endl <<
      "BS Yaw = " << form(params.bsYaw)  << endl;
  }
  fout << "Disc Exponent = " << form(params.dExp) << endl <<
    "Phase offset = " << form(params.phi0) << endl << endl;
  fout.close();
}


void write(const vector<string>& args,
	   const LFIT::Data& data,
	   const LFIT::SynthData& model){
  
  Subs::Input input(args,"LFIT_ENV",".lfit");
  //sign in input variables
  input.sign_in("filename",       Subs::Input::LOCAL, Subs::Input::PROMPT);
  //get input
  string filename;
  input.get_value("filename", filename, "out.dat", "filename for output");
  
  if(model.flux.size() != data.size()){
    std::string err = "Misatched sizes in data and model: you need to scale\n";
    throw err;
  }
  
  std::ofstream fout(filename.c_str());
  for(int i=0; i<data.size(); i++){
    fout << data[i].time << " " << data[i].expose << " " << data[i].flux << " " << data[i].ferr <<  " " << model.flux[i] << std::endl;
  }
}
