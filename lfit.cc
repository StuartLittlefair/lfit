#include <iostream>
#include <string>
#include <readline/readline.h>
#include <readline/history.h>
#include "trm_subs.h"
#include "lfit.h"

using namespace std;

// common variables
int nphi = 5; // number of bins to divide each point into for trapezoidal integration

void updatewithCheck(const LFIT::CV& myCV, LFIT::Params& params, const LFIT::Params& pin){

  Subs::Array1D<double> currParams=pin.get_variables();
  std::vector< std::pair<double,double> > limits = params.get_limits(myCV.bspot.getTangent());
  Subs::Format form;
  form.precision(3); form.general();
  bool outside = false;
  Subs::Array1D<int> badPars;
  for(int i=0; i<params.nvary(); i++){
    if(currParams[i] > limits[i].second || currParams[i] < limits[i].first) {
      badPars.push_back(i);
      outside = true;
    }
  }
  if(outside){
    std::cout << "Bright spot angle should be within 40 degrees of " << form(180.*myCV.bspot.getTangent()/Constants::PI) << std::endl;
    std::cout << "Problem values are: " << std::endl;
    for(int j=0; j<badPars.size(); j++){
      int i = badPars[j];
      cout << form(limits[i].first) << ", " << form(currParams[i]) << ", " << form(limits[i].second) << std::endl;
    }      
    string error = "Suggested parameters violate limits\n";
    throw error;
  }
  params=pin;
}

int main (int argc, char* argv[]) {
    
  try{
       
       
    int lastFit=-1; // which colour was last to be fit
    // set to -1 on entry. 0=red, 1=grn, 2=blu
    bool pos=true; // force fluxes of all components to be positive?
	   
    LFIT::Params params;
    LFIT::CV myCV;
    myCV.donor = LFIT::Donor(params.q,400);
    myCV.prim  = LFIT::WhiteDwarf(params.rwd,params.ulimb);
    myCV.bspot = LFIT::BrightSpot(params);
    myCV.disc  = LFIT::Disc(params.q,params.rwd,params.rd,params.dExp,1000);
       
    LFIT::Data red, grn, blu;
    LFIT::ucam_read(red,grn,blu);
    LFIT::ThreeColour tricolor     = LFIT::ThreeColour(red,grn,blu);
    LFIT::ThreeColour tricolorOrig = LFIT::ThreeColour(red,grn,blu);
    LFIT::SynthData model;
              
    // READY FOR USER INPUT!!!!
    printcomms(pos,params.complexSpot);

    // Batch mode control
    bool batch = (isatty(0) == 0 || isatty(1) == 0);
    
    // prompting structure is cribbed from Tom's PONTO program
    bool more = true;
    std::vector<std::vector<std::string> > args;
    while(more){

      if(batch){

	// batch mode input, avoid 'readline'
	while(1){
	  // Skip comment lines
	  char c = cin.peek();
	  while(std::cin && c == '#'){
	    cin.ignore(10000,'\n');
	    c = cin.peek();
	  }

	  args = Subs::read_multi_string_line(cin);
	  if(args.size() > 0) break;
	}

      }else{

	// terminal input
	char* entry = (char*)NULL;
	while(!entry){
	  entry = readline("lfit> ");
	  if(entry && *entry) add_history(entry);
	  if(entry && entry[0] == '!'){
	    system(entry + 1);
	    args.clear();

	  }else if(entry){
	    std::istringstream istr(entry);
	    args = Subs::read_multi_string_line(istr);
	  }
	}
	if(entry) free(entry);
      }
      // find a matching command and run it
      for(size_t ncom=0; ncom<args.size(); ncom++){
	if(args[ncom].size() > 0){
	  int nmatch = LFIT::findcomm(args[ncom][0]);
	  std::string clower = args[ncom][0];
	  if(nmatch ==1){
	    // Start the command block
	    if(clower  == "rebin"){
	      tricolor.clear();
	      Subs::Format form;
	      form.precision(3); form.general();
	      std::cout << "Red limits: " << form(tricolorOrig[0].get_xmin()) 
			<< " - " << form(tricolorOrig[0].get_xmax()) << std::endl;
	      std::cout << "Grn limits: " << form(tricolorOrig[1].get_xmin()) 
			<< " - " << form(tricolorOrig[1].get_xmax()) << std::endl;
	      std::cout << "Blu limits: " << form(tricolorOrig[2].get_xmin()) 
			<< " - " << form(tricolorOrig[2].get_xmax()) << std::endl;
	      try{
		tricolor = rebin(args[ncom],tricolorOrig);
	      }catch(std::string err){
		cout << err << endl;
	      }
	    }else if(clower == "plot"){
	      if (lastFit == -1){
		cout << "You must scale a model to the lightcurve before plotting\n";
	      }else{
		light(args[ncom],model,tricolor[lastFit]);
	      }
	    }else if(clower == "write"){
	      write(args[ncom],tricolor[lastFit],model);
	    }else if(clower == "pos"){
	      pos = !pos;
	      std::string onoff=" off";
	      if(pos)onoff=" on";
	      std::cout<<"Positivity is"<<onoff<<std::endl;
	    }else if(clower == "complex"){
	      params.complexSpot = !params.complexSpot;
	      std::string onoff=" off";
	      if(params.complexSpot)onoff=" on";
	      std::cout<<"Complex Spot model is"<<onoff<<std::endl;
	      myCV.bspot = LFIT::BrightSpot(params);
	    }else if(clower == "lmask"){
	      lmask(args[ncom],tricolor);
	    }else if(clower == "smask"){
	      smask(args[ncom],tricolor);
	    }else if(clower == "wmask"){
	      wmask(args[ncom],tricolor);
	    }else if(clower == "load"){
	      try{
		load(args[ncom],params);
		myCV.bspot = LFIT::BrightSpot(params);
	      }catch(string error){
		cout << error << endl;
	      }
	    }else if(clower == "dump"){
	      dump(args[ncom],params);
	    }else if(clower == "edit"){
	      LFIT::Params psave = params;
	      try{
		edit(args[ncom],psave);
		updatewithCheck(myCV,params,psave);
	      }catch(string error){
		cout << "Failed" << endl;
		cout << error << endl;
	      }
	    }else if(clower == "show"){
	      show(params);
	    }else if(clower == "convert"){
	      convert(params);
	    }else if(clower == "scale"){
	      Subs::Format form;
	      form.precision(8); form.general();
	      double chisq=fastscale(args[ncom],params,tricolor,myCV,model,
				     lastFit, false, pos, true, true, true, true);
	      int dof = DOF(tricolor[lastFit],params);
	      cout << "Chisq: " << form(chisq) << ", DOF: " << dof << endl;
	    }else if(clower == "fscale"){
	      Subs::Format form;
	      form.precision(8); form.general();
	      double chisq=fastscale(args[ncom],params,tricolor,myCV,model,
				     lastFit, false, pos,false,false,false,false);
	      int dof = DOF(tricolor[lastFit],params);
	      cout << "Chisq: " << form(chisq) << ", DOF: " << dof << endl;
	    }else if(clower == "tweak"){
	      try{
		tweak(args[ncom],params, tricolor, myCV, model, lastFit, pos);
	      }catch(std::string err){
		std::cout << "TWEAK failed: " << err << std::endl;
	      }
	    }else if(clower == "mcmc"){
	      mcmc(args[ncom],params, tricolor, myCV, model, lastFit, pos);
	    }else if(clower == "errors"){
	      errors_sem(args[ncom],params, tricolor, myCV, model, lastFit, pos);
	    }else if(clower == "boot"){
	      errors_bootstrap(args[ncom],params, tricolor, myCV, model, lastFit, pos);
	    }else if(clower  == "exit" || clower  == "quit"){
	      std::cout << "Are you sure you want to " << clower << " LFIT? ";
	      std::string reply;
	      getline(std::cin,reply);
	      if(reply == "y" || reply == "Y"){
		std::cout << "bye then." << std::endl;
		more = false;
		break;
	      }
	    }
	    // if we got here, the command was either ambiguous or not found
	  }else if(nmatch == 0){
	    cout << "Command: '" << args[ncom][0] << "' not recognised." << std::endl;
	    printcomms(pos,params.complexSpot);
	  }else{
	    cout << "Command: '" << args[ncom][0] << "' not unique." << std::endl;
	    printcomms(pos,params.complexSpot);
	  }
	}
      }
    }
  }catch(const string& err){
    cerr << "\nError occured inside LFIT:" << endl;
    cerr << err << endl;
  }

  return 0;   

}

