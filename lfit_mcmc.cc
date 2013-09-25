/*
 *  lfit_mcmc.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 06/01/2009.
 *  Copyright 2009 University of Sheffield. All rights reserved.
 *
 */
#include "trm_subs.h"
#include "trm_array1d.h"
#include "trm_array2d.h"
#include "trm_input.h"
#include "lfit.h"
#include <cfloat>
#include <climits>
#include <ctime>
#include <iostream>
#include <fstream>
#include "cholesky.h"

template <class X>
Subs::Array1D<X> matMul(const Subs::Array2D<X>& A, const Subs::Array1D<X>& b){

  int nx = A.get_nx();
  int ny = A.get_ny();
  if( b.size() != ny){
    throw Subs::Array2D_Error("Cannot multiply matrix of dimensions [" + Subs::str(nx) + "][" + Subs::str(ny) + "], with vector of size " + Subs::str(b.size()));
  }
  Subs::Array1D<X> prod(ny);
  prod = 0.0;
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++) {
      prod[i] += A[j][i]*b[j];
    }
  } 
  return prod;
}

bool stepAccept(const double& chiNew, const double& chiOld, int& seed);

double sdev(const Subs::Array1D<double>& arr){
	
  double mean = arr.mean();
  double sd = 0.0;
  for(int i=0; i<arr.size(); i++){
    sd += Subs::sqr(arr[i]-mean);
  }
  return sqrt(sd/double(arr.size()));
}

class mcmcParams{
	
public:	
  int seed;
  //! constructor from params object
  mcmcParams(const Subs::Array1D<double> paramsIn, 
	     std::vector< std::pair<double,double> > limitsIn,
	     Subs::Array2D<double> covarsIn, double factorIn,
	     int nstepsIn){
    seed = (int) time(NULL);
    gaussSeed = (int) time(NULL);
    startParams = paramsIn;
    limits      = limitsIn;
    covar       = covarsIn;
    lowTri      = getLowTri(covar);
    nsteps      = nstepsIn;
    stepFac     = factorIn; 
    currParams  = paramsIn;
    bestParams  = paramsIn;
    stepNum     = 0;
    successes   = 0;
    targetSuccessRate = 0.234; // 'ideal' success rate from Roberts & Rosenthal 1998)
    npars       = paramsIn.size();
    allSteps.resize(nsteps,paramsIn.size());
  }

  Subs::Array2D<double> getLowTri(const Subs::Array2D<double>& arrIn){
    // given a positive, definite matrix use cholesky decomposition to return the lower triangular matrix
    // this can be used to select variables from a normal multivariate distribution, given the method 
    // presented on the wikipedia article (how's that for a reference!)

    // we are supplied a covariance matrix. we need also to calculate a lower Triangular matrix using cholesky
    // decomposition, which we use to generate the next step in the chain by sampling from the multivariate normal
    // distribution associated with the covariance matrix
    Subs::Array2D<double> tmp;
    Subs::Array2D<double> sampleCovar = arrIn;

    bool notDone=true;
    int ntries=0;
    float ridge = 0.001;
    while(notDone){
      try{
	ntries++;
	tmp = choleskyDecomp(sampleCovar);
	notDone=false;
      }catch(const std::string& err){
	// the covariance matrix is not positive definite, so we add a small diagonal matrix until it is.
	// the size of the values in this diagonal matrix is multiplied by 10 each time we enter here.
	for(int j=0; j<nvary(); j++){sampleCovar[j][j] = arrIn[j][j] + arrIn[j][j]*ridge*10.0*(ntries-1);} 
      }
      if (ntries>50) {
	std::string err = "Can't update covariance matrix, even after changes. Call Stu!\n";
      }
    }
    return tmp;
  }
	
  int nvary(){return npars;}

  void tryStep(){	

    //int seed = (int) time(NULL);
    Subs::Array1D<double> gaussRans(nvary());
    int nAttempts=0;
    Subs::Array1D<int> badParams(0);

    //nAttempts++;
    for (int i=0; i<nvary(); i++){
      gaussRans[i] = Subs::gauss1(gaussSeed);
    }
    // sample new parameters from multivariate normal distribution with covariance given by covars
    // and apply current scaling
    Subs::Array1D<double> stepVec = matMul(lowTri,gaussRans);
    currParams = startParams + stepFac*stepVec;
  }

  bool outSideLimits(Subs::Array1D<int>& badParams){
    if(badParams.size() > 0) badParams.resize(0);
    bool outside = false;
    for(int i=0; i<nvary(); i++){
      if(currParams[i] > limits[i].second || currParams[i] < limits[i].first) {
	outside = true;
	badParams.push_back(i);
      }
    }
    return outside;
  }

  bool breaksPriors(){
    bool outside = false;
    for(int i=0; i<nvary(); i++){
      if(currParams[i] > limits[i].second || currParams[i] < limits[i].first) {
	outside = true;
      }
    }
    return outside;
  }
	
  void acceptStep(const double& chisqIn){		
    if(stepNum >= nsteps){
      std::string err = "Exceeded maximum number of steps for mcmc Object";
      throw err;
    }
    chisq.push_back(chisqIn);
    for (int i=0; i<nvary(); i++){
      allSteps[stepNum][i] = currParams[i];
    }
    stepNum++;
    successes++;
    // assign start parameters for next step
    startParams = currParams;
  }

  void rejectStep(){
    // WHAT DO WE DO HERE? PUT IT IN CHAIN OR NOT?
    if(stepNum >= nsteps){
      std::string err = "Exceeded maximum number of steps for mcmc Object";
      throw err;
    }
    chisq.push_back(chisq[chisq.size()-1]);
    for (int i=0; i<nvary(); i++){
      allSteps[stepNum][i] = startParams[i];
    }
    currParams = startParams;
    stepNum++;
  }
	
  int length(){
    return stepNum;
  }

  bool converged(){
    double thisChi = chisq[chisq.size()-1];
    Subs::Array1D<double> buff(chisq);
    double medChi = buff.median();
    if(thisChi > medChi){
      return true;
    }else{
      return false;
    }
  }

  Subs::Array1D<double> currPars(){
    Subs::Array1D<double> retArr(npars);
    for(int i=0; i<npars; i++){
      retArr[i] = currParams[i];
    }
    return retArr;
  }

  Subs::Array1D<double> bestPars(){
    Subs::Array1D<double> retArr(npars);
    for(int i=0; i<npars; i++){
      retArr[i] = bestParams[i];
    }
    return retArr;
  }
  
  double stepFactor(){
    return stepFac;
  }

  double covarEl(int i, int j, const int& nsteps){
    
    int ny = stepNum;
    Subs::Array1D<double> col1(nsteps), col2(nsteps);
    int l=0;
    for(int k=ny-nsteps; k<ny; k++){
      col1[l] = allSteps[k][i];
      col2[l] = allSteps[k][j];
      l++;
    }
    col1 -= col1.mean();
    col2 -= col2.mean();
    
    col1 *= col2;
    double sum=0.0;
    for(int k=0; k<col1.size(); k++) sum += col1[k];
    return sum/float(col1.size()-1);
  }

  void updateCovar(const int& nsteps){
    // calculate covariance from most recent part of chain
    // uses last nsteps steps
    int ny = nvary();
    for(int j=0; j<ny; j++){
      for(int i=0; i<ny; i++){
	covar[i][j] = covarEl(i,j,nsteps);
      }
    }
    lowTri = getLowTri(covar);
  }

  Subs::Array2D<double> currCovar(){
    return covar;
  }
  	
  double conflims(Subs::Array1D<double>& vals, Subs::Array1D<double>& chi){
		
    Subs::Array1D<double> chiSave = chi;
    Subs::Array1D<double> valSave = vals;

    // sort by chisq and data value
    // routine sorts array and returns key to ORIGINAL order
    Subs::Buffer1D<unsigned long int> chiKey = chi.sort();
    Subs::Buffer1D<unsigned long int> dataKey = vals.sort();
    int np = vals.size();

    // pick out maximum likelihood value (lowest chiSq)
    double xml = valSave[chiKey[0]];
		
    //Initialise the upper and lower confidence limits, setting them equal to xml
    double xmax = xml;
    double xmin = xml;
    double dxl=0.0;
    double dxh=0.0;
    int jlo=0; 
    int jhi=0;
    for(int i=1; i<np; i++){
      /*
	Compare the x value corresponding to next highest value of 
	chisq with the current upper and lower confidence limits.
      */
      xmax = std::max(xmax,valSave[chiKey[i]]);
      xmin = std::min(xmin,valSave[chiKey[i]]);
			
      // find index of xmax and xmin in the sorted x array
      vals.hunt(xmax,jhi);
      vals.hunt(xmin,jlo);
			
      // update confidence levels relative to maximum likelihood value
      dxl = xml-xmin;
      dxh = xmax-xml;
			
      /* As soon as xmax and xmin enclose 68.3% of the x values,
	 terminate the search */
      if(float(jhi-jlo) >= 0.683*float(np)){
	break;
      }
    }
		
    //std::cout << "Lower confidence limit is " << dxl << std::endl;
    //std::cout << "Upper confidence limit is " << dxh << std::endl;
    double avgErr = 0.5*(dxl+dxh);
    //std::cout << "Average confidence limit is " << avgErr << std::endl;
    return avgErr;
  }
	
  Subs::Array1D<double> conflims(){
    Subs::Array1D<double> retArr(npars);
    for(int i=0; i <npars; i++){
      Subs::Array1D<double> parVals(stepNum);
      for(int j=0; j<stepNum; j++){parVals[j]=allSteps[j][i];}
      retArr[i] = conflims(parVals,chisq);
    }
    return retArr;
  }

  void setCovar(const Subs::Array2D<double>& covarIn){
    if(covarIn.get_nx() != covar.get_nx() || covarIn.get_ny() != covar.get_ny()){
      std::string err = "mcmcParams: Conflicting sizes in covariance array";
      throw err;
    }
    covar = covarIn;
    Subs::Array2D<double> tmp;
    lowTri = getLowTri(covar);
  }
	
  void updateBestPars(const Subs::Array1D<double>& parsIn){
    if(parsIn.size() != npars){
      std::string err = "mcmcParams: Conflicting sizes in parameter array";
      throw err;
    }
    bestParams = parsIn;
  }
	
  void reportSuccessRate(std::ofstream& fout){
    double successRate = double(successes)/double(stepNum);      
    fout << "Step num: " << stepNum << ", successRate = " << successRate << ", step scaling factor = " << stepFac 
	 << std::endl; 
  }
    
  void setStepFacFromChain(std::ofstream& fout){
    double successRate = double(successes)/double(stepNum);      
    double delta = successRate / targetSuccessRate;
    double deltaStep = std::min(0.01,pow(float(stepNum),-0.5));
    if(delta > 1.0){
      // too successful, increase stepsize
      stepFac *= exp(deltaStep);
    }else{
      stepFac /= exp(deltaStep);
    }
    fout << "Step num: " << stepNum << ", successRate = " << successRate << ", step scaling factor = " << stepFac 
	 << std::endl; 
  }

  void setStepFacFromLocalChain(const int& nsteps, std::ofstream& fout){
    double successRate = double(successes)/double(nsteps);
    double delta = successRate / targetSuccessRate;
    stepFac *= delta;
    fout << "Step num: " << stepNum << ", successRate = " << successRate << ", step scaling factor = " << stepFac 
	 << std::endl; 
    successes = 0;
  }

  double correl(Subs::Array1D<double>& chain){
    double xbar = chain.mean();
		
    double aii=0.0;
    for(int i=0;i<chain.size();i++){
      aii += Subs::sqr(chain[i]-xbar);
    }
    double denom = aii/double(chain.size());
		
    // calculate auto-correlation function
    int acNp=100;
    Subs::Array1D<double> cfunc(acNp);
    if(denom > 1.0e-20){
      for(int j=1; j<=acNp; j++){
	double sij = 0.0;
	int ioff = j-1;
	for(int i=0; i<chain.size()-ioff; i++){
	  sij += (chain[i]-xbar)*(chain[i+ioff]-xbar);
	}
	double aij = sij/double(chain.size()-ioff);
	cfunc[ioff] = aij/denom;
      }
			
      int jlo=0;
      // find where correlation function drops below 0.5
      cfunc.hunt(0.5,jlo);

      // interpolate to get correlation length
      if(jlo >= 99){
	std::string err = "Correl: correlation length insufficient to determine correlation length";
	throw err;
      }else if(jlo == 0){
	return 0.0;
      }else{
	int jhi = jlo + 1;
	double clo = cfunc[jlo];
	double chi = cfunc[jhi];
	double f = (0.5-clo)/(chi-clo);
	return double(jlo) + f;
      }
    }
    // variance is zero, return stupid correlation length.
    return -99.0;
  }

  Subs::Array1D<double> correl(){
    // returns correlation length for MCMC chain parameters
    // recipe is Tegmark et al 2004.
    Subs::Array1D<double> retArr(npars);
    for(int i=0; i <npars; i++){
      Subs::Array1D<double> parVals(stepNum);
      for(int j=0; j<stepNum; j++){parVals[j]=allSteps[j][i];}
      try{
	retArr[i] = correl(parVals);
      }catch(std::string err){
	std::cout << err << std::endl;
      }
    }
    return retArr;
  }
		
private:
  Subs::Array1D<double> startParams;
  Subs::Array1D<double> currParams;
  Subs::Array1D<double> bestParams;
  Subs::Array2D<double> covar,lowTri;
  Subs::Array2D<double> allSteps;
  Subs::Array1D<double> chisq;
  std::vector< std::pair<double,double> > limits;
  double targetSuccessRate, stepFac;
  int nsteps, stepNum, npars, successes, gaussSeed;
};

void mcmc(const std::vector<std::string>& args, LFIT::Params& params,
	  const LFIT::ThreeColour& three, const LFIT::CV& myCV,
	  const LFIT::SynthData& myModel,int& colFit, const bool& pos){

  int nburn, njump;
  int ntweak=1500;
  double scaleFactor;
  bool reScale=true;
  { // get input values
    Subs::Input input(args,"LFIT_ENV",".lfit");
    input.sign_in("patient",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("colour",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nburn",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("njump",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("scaleFac",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("reScale",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("qvar",       Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("dpvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("rdvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("rwdvar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("ulvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("scvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("azvar",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("fracvar",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("exp1var",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("exp2var",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("tiltvar",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("yawvar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);		
    input.sign_in("expvar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("phivar",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
		
    // are we ready to go?
    bool keepGoing=false;
    input.get_value("patient", keepGoing, true,     "This WILL take a while. Continue (y/n)?");
    if (!keepGoing) return;
		
    input.get_value("colour", colFit, 10, 1, 3, "lightcurve colour to fit");
    input.get_value("nburn", nburn, 15000, 50, 10000000, "no of steps for burn-in phase (>=15,000 recommended)");
    input.get_value("njump", njump, 15000, 50, 10000000, "no. of steps for MCMC chain (>=15,000 recommended)");
    input.get_value("scaleFac", scaleFactor, 0.2, 0.01, 50.0, "scale factor for step length");
    input.get_value("reScale", reScale, true, "Re-scale step size as we proceed? (recommended)");
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
    colFit -= 1;
  }

  /* before we do anything, we must get a good estimate of what step sizes to use for
     each MCMC parameter. The best way of doing this is to run a levmarq optimisation, and
     use the covariance array for the parameters. Each step in the MCMC chain then draws a value from a multivariate
     normal distribution defined by the means and covariance array from the levmarq step.
     Optionally, we can add a scaling to the covariance array which is tuned to keep step
     acceptance rate near 0.25
     
     ...slight complication - have made q,dphi change in levarq much smaller to avoid unphysical changes in q,dphi 
     
     ... also, if we make dphi -ve then dphi always gets smaller, so a physical parameter set will always get through lmcomp
  */
  Subs::Array1D<double> origParams = params.get_variables(),dStep=params.get_variables();
  double fact = 0.01;
  dStep *= fact;
  if(params.qvar) dStep[0] = origParams[0]*0.001;
  if(params.qvar && params.dpvar) {
    dStep[1] = -origParams[1]*0.0005;
  }else if(params.dpvar && ! params.qvar){
    dStep[0] = -origParams[0]*0.0005;
  }
  LFunc lfunc = LFunc(three,myCV,myModel,params,colFit,pos,dStep);
  Subs::Array2D<double> covar;
  double lambda=-1.,lambda_old,thisChi=0.,lastChi;
  int ncount=0, nfail=0;
  const int NFMAX=3;
  const int NMAX=50;
  const double LAMBDAMAX=2.0e11;
  while(ncount<NMAX && nfail<NFMAX && lambda<LAMBDAMAX){
    lambda_old = lambda;
    lastChi = thisChi;
    int status = lmfit(lfunc,thisChi,lambda,covar);
    
    // check whether we are making enough progress
    if(ncount){
      if(thisChi<lastChi){
	if(lastChi-thisChi < 0.01) 
	  nfail++;
	else
	  nfail=0;
      }
    }
    
    if (status==0) ncount++;
  }
  
  std::string message;
  if(ncount == NMAX)  message = " on reaching maximum number = " + Subs::str(NMAX);
  if(nfail  == NFMAX) message = " after chi**2 decreased by less than 0.01 " + Subs::str(NFMAX) + " times";
  if(lambda >= LAMBDAMAX) message = " when lambda exceeded maximum = " + Subs::str(LAMBDAMAX);
  std::cout << "Initial levmarq iterations done: halted "  << message << std::endl;
  
  // run again to get solution and covariance array
  lambda = 0.0;
  lmfit(lfunc,thisChi,lambda,covar);
  
  // we can either start from originally supplied parameters, or from end point of levmarq...
  // on one hand, it doesn't make sense to start from a poorer fit, as the covariance array
  // we have calculated is for the region around best fit. On the other hand you want to 
  // be able to start the MCMC chain from anywhere to test for convergence. 
  // currently starting from where levmarq left off
  origParams = lfunc.get_variables();
  params.updateWith(origParams);
  
  Subs::Format form(8);
  std::cout << "Levmarq fit returned:" << std::endl;
  for(int i=0; i<params.nvary(); i++){
    std::cout << "Param #" << i << " = " << form(origParams[i]) << " +/- " << 
      form(sqrt(covar[i][i])) << std::endl; 
  }
  std::cout << "chi**2 = " << form(lfunc.chisq_min) << " and neval = " << lfunc.neval << std::endl;

  // create ameoba/mcmc function object
  // when called it scales model CV to the data (colour colFit)
  // and returns chisq
  amFunc func = amFunc(three,myCV,myModel,params,colFit,pos);

  // start at initial parameter values and % errors
  std::vector< std::pair<double,double> > limits = func.get_limits();
  //scaleFactor = 2.38*2.38/double(params.nvary()); // optimal variance from Roberts & Rosenthal 2001

  mcmcParams burnIn = mcmcParams(origParams,limits,covar,scaleFactor,nburn); 

  // take a random step from initial values, acknowledging the priors
  do{
    burnIn.tryStep();
  }while(burnIn.breaksPriors());

  
  // evaluate chisq
  Subs::Buffer1D<double> chisq(nburn);
  chisq[0] = func(burnIn.currPars());
  double chiMin = chisq[0];
    
  //accept step without test
  burnIn.acceptStep(chiMin);

  /*
    During burn-in, the algorithm should converge to the optimum solution.
    c	MCMC allows us to explore the covariance matrix around the optimum 
    c	solution. The random step sizes in each parameter are determined
    c	by user input. We can optionally evaluate the optimal step size
    c	empirically from variance in parameters after a few hundred successful jumps
  */	
  std::ofstream myfile;
  myfile.open("burnInLog.txt");
  std::ofstream scalefile;
  scalefile.open("burnInScale.txt");

  // create progress bar
  int nhashes=20;
  message = "Starting burn-in phase";
  LFIT::ProgressBar bar = LFIT::ProgressBar(message,nhashes,nburn);
  bar.initialise();

  // run burn-in
  int itry=0, isucceed=0, itotal=0;
  int iburn;
  bool converged=false;
  for (iburn=1; iburn < nburn; iburn++){
		
    // print extra hash to status bar if necessary
    bar.update(iburn);
		
    bool accept = false;
    double chisqNow;

    // get new parameter set
    // make a step, until that step is accepted
    // note: you will not end up with a formal MCMC chain if you rescale step size based on local acceptance rate

    itry++;
    itotal++;
    burnIn.tryStep();
    try{
      // this will return very high chisq if parameters are unphysical
      chisqNow = func(burnIn.currPars());
    }catch(const string& err){
      // if something has gone wrong, don't accept this step!!!
      chisqNow = 1.0e35;
    }			
    accept = stepAccept(chisqNow,chisq[iburn-1],burnIn.seed) && !burnIn.breaksPriors();
    
    if(accept){
      isucceed++;
      burnIn.acceptStep(chisqNow);
      chisq[iburn] = chisqNow;
      myfile << "BurnIn: Step " << iburn << " accepted after "<< itry << " attempts" << std::endl; 
      myfile << burnIn.currPars() << " " << chisq[iburn] << std::endl;
      itry=0;
    }else{
      burnIn.rejectStep();
      chisq[iburn] = chisq[iburn-1];
      myfile << "BurnIn: Step " << iburn << " rejected after "<< itry << " attempts" << std::endl; 
      myfile << burnIn.currPars() << " " << chisq[iburn] << std::endl;
    }
		
    // after 30 successful steps, rescale step sizes 
    if(reScale){
      if(isucceed == 30){
	//rescale step sizes
	burnIn.setStepFacFromLocalChain(itotal,scalefile);
	itotal=0;
	isucceed=0;
      }
    }

    // update covariance array from whole chain
    if(iburn>=5000 && iburn%5000 == 0){burnIn.updateCovar(iburn-1);}

    // update best parameters if this is the lowest chisq so far?
    if(chisq[iburn] < chiMin){
      burnIn.updateBestPars(burnIn.currPars());
      chiMin = chisq[iburn];
    }

    // have we converged? If so, stop the burn in
    if(accept && iburn > 1000){
      if (burnIn.converged()){
	converged=true;
      }
    }
  }
  // set final covariance from last 15000 steps of burnIn chains (or full length)
  int nsteps = std::max(14999,iburn-1);
  burnIn.updateCovar(nsteps);
  myfile << "Finished burn-in phase" << std::endl;
  // end progress bar
  bar.end();
	
  // close logfile
  myfile.close();
  scalefile.close();


  /* Now let's check to see if we converged */
  if(converged){
    // oh no, we didn't. throw a warning, and prompt user to retry from here!
    std::cout << "Warning! burn-in phase failed to meet convergence criterion!" << std::endl;
    std::cout << "Parameters have been updated to reflect best parameters found during burn-in phase" <<std::endl;
    std::cout << "I recommend starting another MCMC from here, perhaps with a longer burn-in phase..." << std::endl;
  }

  /* Production run 
     start from wherever the burn-in phase took us.
  */
  // new logfile
  myfile.open("mcmcLog.txt");
  scalefile.open("mcmcScale.txt");

  // new objects to perform monte carlo	
  chisq.resize(njump);
  mcmcParams mcmc = mcmcParams(burnIn.bestPars(),limits,burnIn.currCovar(),burnIn.stepFactor(),njump);
		
  // take a random step from initial values, without violating priors
  do{
    mcmc.tryStep();
  }while(mcmc.breaksPriors());

	
  // evaluate chisq
  chisq[0] = func(mcmc.currPars());
  chiMin = chisq[0];
  mcmc.acceptStep(chiMin);
	
  // re-create progress bar
  nhashes=40;
  message = "Starting production phase";
  bar = LFIT::ProgressBar(message,nhashes,njump);
  bar.initialise();
	
  // run production chain
  // add step to chain whether accepted or rejected
  itry=0; isucceed=0; itotal=0;
  for (int istep=1; istep < njump; istep++){
    bool accept = false;
    double chisqNow;
		
    // print extra hash to status bar if necessary
    bar.update(istep);
		
    // get new parameter set

    // attempt a step
    mcmc.tryStep();
    itry++;
    itotal++;
    try{
      // this will go to catch block if parameters are invalid (i.e no eclipse)
      chisqNow = func(mcmc.currPars());
    }catch(const string& err){
      chisqNow = 1.0e35;
    }
    accept = stepAccept(chisqNow,chisq[istep-1],mcmc.seed) && !mcmc.breaksPriors();

    if(accept){
      isucceed++;
      mcmc.acceptStep(chisqNow);
      chisq[istep] = chisqNow;
      myfile << "Production: Step " << istep << " accepted after "<< itry << " attempts" << std::endl; 
      myfile << mcmc.currPars() << " " << chisq[istep] << std::endl;
      itry=0;
    }else{
      mcmc.rejectStep();
      chisq[istep] = chisq[istep-1];
      myfile << "Production: Step " << istep << " rejected: currently "<< itry << " attempts" << std::endl; 
      myfile << mcmc.currPars() << " " << chisq[istep] << std::endl;
    }

    // update if this is the best so far
    if(chisq[istep] < chiMin){
      mcmc.updateBestPars(mcmc.currPars());
      chiMin = chisq[istep];
    }

    // after 30 successful steps, rescale step sizes 
    // scaling from local chain is right out. You can scale from whole chain, but I might remove this
    if(reScale){
      if(isucceed == 30){
	//rescale step sizes
	//mcmc.setStepFacFromLocalChain(itotal,scalefile);
	//mcmc.setStepFacFromChain(scalefile);
	mcmc.reportSuccessRate(scalefile);
	itotal=0;
	isucceed=0;
      }
    }
    // update covariance array from whole chain every 1000 steps
    //if(istep>=5000 && istep%1000 == 0){mcmc.updateCovar(istep-1);}
    
  }
	
  // done!
  // close progress bar
  bar.end();
	
  // close files
  myfile.close();
  // get two-tailed, one sigma confidence limits on parameters
  Subs::Array1D<double> vals = mcmc.bestPars();
  Subs::Array1D<double> errs = mcmc.conflims();
	
  // update parameters to reflect values and errors
  params.updateWith(vals);
  params.updateErrsWith(errs);	

  // output of values
  show(params);

  // sanity-check Markov-Chain correlation lengths
  try{
    Subs::Buffer1D<double> clengths = mcmc.correl();
    std::cout << "Sanity-check: correlation lengths of MCMC chains are: " << clengths << std::endl;
  }catch(std::string err){
    std::cout << "WARNING: Correlation length calculation failed " << std::endl;
    std::cout << err << std::endl;
  }
	
}

bool stepAccept(const double& chiNew, const double& chiOld, int& seed){
	
  // Metropolis-Hastings decision maker
  double dChi = chiNew - chiOld;
	
  // accept if chisq is lower
  if(dChi < 0.0){return true;}
	
  double pjump = exp(-dChi/2.0);
	
  // pick uniform random deviate in range 0 - 1
  double prob = Subs::ran1(seed);
  if(prob < pjump){return true;}
	
  return false;
}
