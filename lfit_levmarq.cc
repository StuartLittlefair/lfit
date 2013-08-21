
#include "trm_subs.h"
#include "trm_array1d.h"
#include "trm_array2d.h"
#include "lfit.h"
#include <cfloat>
#include <climits>
#include <ctime>
#include <iostream>
#include <fstream>

void LFunc::force_recalc(){
  wdCalc=true;
  donorCalc=true;
  discCalc=true;
  bsCalc=true;
}

void LFunc::check_recalc(){

  wdCalc=false;
  donorCalc=false;
  discCalc=false;
  bsCalc=false;

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

void LFunc::updateParams(const LFIT::Params& paramsIn){
  old_params = params;
  params = paramsIn;
  check_recalc();
}

void LFunc::updateParams(const Subs::Array1D<double>& paramsIn){  
  LFIT::Params tmpParams = params;
  params.updateWith(paramsIn);
  old_params = tmpParams;
  check_recalc();
}

std::string LFunc::asString(){
  return params.toString();
}
    
Subs::Array1D<double> LFunc::get_ranges(){
  return params.get_ranges();
}
    
Subs::Array1D<double> LFunc::get_variables(){
  return params.get_variables();
}
    
std::vector< std::pair<double,double> > LFunc::get_limits(){
  return params.get_limits(CV.bspot.getTangent());
}
    
void LFunc::setParamRanges(){
  params.setParamRanges();
}
    
int LFunc::nvary(){
  return params.nvary();
}
    
double LFunc::fitModel(){
  
  /* evaluates function
   */
  // have we been called yet? if so, force full calculation
  if (newBorn){
    wdCalc=true;
    discCalc=true;
    donorCalc=true;
    bsCalc=true;
    newBorn=false;
  }
  
  // use these params to scale model
  std::vector<std::string> dummy;
  dummy.push_back("");
  float chisq = fastscale(dummy,params,three,CV,model,colFit,true,pos,
			  wdCalc,donorCalc,discCalc,bsCalc);
  Subs::Format form(8);
  //std::ofstream fout("levmarq.log",std::ios::ate);
  //std::cout << "Calling LFunc: " << wdCalc << ", " <<  donorCalc << ", " << discCalc << ", " << bsCalc << " - " << form(chisq) << std::endl;
  return chisq;	
}
    
    
void LFunc::lmcomp(Subs::Buffer2D<double>& alpha, Subs::Buffer1D<double>& beta, double& chisq){
      
  static bool first = true;
  static double cmin;
  bool start = (neval == 0);
  
  Subs::Array1D<double> centre = get_variables();
  chisq = fitModel();
  LFIT::SynthData savedFit = model;
  LFIT::CV savedCV = CV;
  neval++;
  
  if(first){
    cmin=chisq;
    first=false;
  }else if(chisq<cmin){
    cmin=chisq;
  }else{
    // return early, we are not progressing so derivatives will not be used
    return;
  }
  
  // Compute derivatives, using finite differences centred on current point
  LFIT::Data data = three[colFit];
  int ndata = data.size();
  Subs::Array1D<float> buff(ndata);
  Subs::Array1D<double> tparam;
  std::vector< Subs::Array1D<float> > deriv;
  double ignore;
  for(int i=0; i<centre.size(); i++){
    tparam = centre;
    tparam[i] -= paramSteps[i];
    
    // attempt parameter update. What to do if unphysical?
    try{
      updateParams(tparam);
    }catch(const std::string err){
      std::cout << "LMCOMP ERR: " << err << std::endl;
    }

    ignore = fitModel();
    neval++;
    deriv.push_back(model.flux);
    
    // attempt to reduce round off error
    volatile double temp = tparam[i];
    tparam[i] -= 2.0*paramSteps[i];
    double h = tparam[i]-temp;
    tparam[i] = temp+h;
    
    // attempt parameter update. What to do if unphysical?
    try{
      updateParams(tparam);
    }catch(const std::string err){
      std::cout << "LMCOMP ERR: " << err << std::endl;
    }

    ignore = fitModel();
    neval++;
    buff = model.flux;
    
    deriv[i] -= buff;
    deriv[i] /= h;
  }
  
  // restore original parameters
  // it shouldn't barf here, as we've already done this
  try{
    updateParams(centre);
  }catch(const std::string err){
    std::cout << "LMCOMP ERR: " << err << std::endl;
  }
  // restore original fit and CV model
  model = savedFit;
  CV = savedCV;
  
  // on with standard levmarq stuff
  Subs::Buffer1D<double> dyda(nvary());
  alpha = 0.0;
  beta  = 0.0;
  
  double wgt, dy, wt;
  for(int i=0; i<data.size(); i++){
    if(! data[i].bad){
      wgt = 1.0/Subs::sqr(data[i].ferr);
      dy  = data[i].flux - savedFit.flux[i];
      
      for(int l=0,j=0; l < nvary(); l++){
	wt = wgt*deriv[l][i];
	for(int m=0,k=0; m<=l; m++)
	  alpha[j][k++] += wt*deriv[m][i];
	beta[j++] += wt*dy;
      }
    }
  }
  
  for(int j=1; j<nvary(); j++)
    for(int k=0;k<j;k++) alpha[k][j] = alpha[j][k];
  
  if(start || chisq < chisq_min){
    chisq_min = chisq;
    Subs::Format form(8);
    std::cout << "Weighted chisq = " << form(chisq) << ", neval = " << neval << std::endl;
  }
}



int lmfit(LFunc& func, double& chisq, double& lambda, Subs::Buffer2D<double>& covar){

  static Subs::Buffer1D<double> beta;
  static Subs::Array1D<double> da;
  static Subs::Buffer2D<double> oneda, alpha;
  static double ochisq;
  
  // initialise
  int nvar = func.nvary();
  if(lambda < 0.){
    beta.resize(nvar);
    da.resize(nvar);
    oneda.resize(nvar,1);
    alpha.resize(nvar,nvar);
    covar.resize(nvar,nvar);
    lambda = 0.001;
    func.lmcomp(alpha, beta, chisq);
    ochisq = chisq;
  }
  
  //Alter linearised fitting matrix by augmenting diagonal elements
  for(int j=0; j<nvar; j++){
    for(int k=0; k<nvar; k++) covar[j][k]=alpha[j][k];
    covar[j][j] = alpha[j][j]*(1.0+lambda);
    oneda[j][0] = beta[j];
  }
  
  // Matrix solution to (hopefully) go to a better solution
  try{
    Subs::gaussj(covar,oneda);
  }catch(const std::string& err){
    std::cout << "Matrix solution failed: " << err << std::endl;
    Subs::Format form(7);
    std::cout << "Covariance Matrix: " << std::endl;
    for(int i=0; i<nvar; i++){
      for(int j=0;j<nvar;j++){
	std::cout << std::setw(13) <<  form(covar[j][i]) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "oneda: " << oneda << std::endl;
    throw err;
  }

  for(int j=0; j<nvar; j++) da[j] = oneda[j][0];
  
  // lambda set to 0 indicates convergence
  if (lambda == 0.) return 0;
  
  // Try new solution
  Subs::Array1D<double> asave = func.get_variables();
  //std::ofstream fout("levmarq.log",std::ios::ate);
  //std::cout << "Da = " << da << std::endl;
  try{
    func.updateParams(asave+da);
    func.lmcomp(covar,da,chisq);
  }catch(const std::string& err){
    // invalid q,dphi pair here, so we need smaller changes, set chisq high and lambda will increase
    std::cout << "invalid params - setting high chisq" << std::endl;
    func.updateParams(asave);
    func.force_recalc();
    func.lmcomp(covar,da,chisq);
    chisq = ochisq+10.0;
  }
    
  if(chisq<ochisq){
    //chisq reduce = joy. reduce lambda
    lambda *= 0.1;
    ochisq=chisq;
    alpha=covar;
    beta=da;
  }else{
    //chisq increase = poo. retrieve old parameters, increase lambda
    lambda *= 10.0;
    func.updateParams(asave);
    chisq=ochisq;
  }
  return 0;
}

  
