#include <string>
#include "trm/subs.h"
#include "trm/input.h"
#include "lfit.h"

using namespace std;

// common variables
int nphi = 5; // number of bins to divide each point into for trapezoidal integration

void load(const std::string filename, LFIT::Params& params){
  
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
}

void write(const std::string filename, const LFIT::SynthData& model){

    std::ofstream fout(filename.c_str());
    fout << "#Phase    Total (mJy)     WD(mJy)     disk(mjy)     HS(mJy)     donor(mJy)" << std::endl;
    for(int i=0; i<model.phase.size(); i++){
        fout << model.phase[i] << " " << model.flux[i] << " " << model.wflux[i] << " " << model.dflux[i] << " "<< model.bflux[i] << " "<< model.rflux[i] << std::endl;
    }
    fout.close();
}

int main (int argc, char* argv[]) {
    
  try{
  
    LFIT::Params params;
    LFIT::CV myCV;
    LFIT::SynthData model;
    
    Subs::Input input(argc,argv,"LFIT_ENV",".lfit");
    input.sign_in("infile",              Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("outfile",           Subs::Input::LOCAL, Subs::Input::NOPROMPT);

    input.sign_in("spotScale",        Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("discScale",        Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("wdScale",        Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("donorScale",        Subs::Input::LOCAL, Subs::Input::PROMPT);

    input.sign_in("startPhase",     Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("endPhase",       Subs::Input::LOCAL, Subs::Input::PROMPT);
    input.sign_in("numPhase",       Subs::Input::LOCAL, Subs::Input::PROMPT);

    std::string infile,outfile;
    double startPhase,endPhase;
    int numPhase;
    
    input.get_value("infile",infile,"gfit.in","file of parameters");
    load(infile, params);
    params.phi0 = 0.0;
    
    input.get_value("wdScale",params.wFlux,0.1,0.0,10.0,"white dwarf scaling (mJy)");
    input.get_value("discScale",params.dFlux,0.1,0.0,10.0,"disc scaling (mJy)");
    input.get_value("spotScale",params.bFlux,0.1,0.0,10.0,"bright spot scaling (mJy)");
    input.get_value("donorScale",params.rFlux,0.1,0.0,10.0,"donor scaling (mJy)");
    
    input.get_value("outfile",outfile,"model.txt","output file");
    input.get_value("startPhase",startPhase,0.0,0.0,1.0,"start phase");
    input.get_value("endPhase",endPhase,1.0,0.0,2.0,"end phase");
    input.get_value("numPhase",numPhase,200,1,10000000,"number of phases");

    myCV.donor = LFIT::Donor(params.q,400);
    myCV.prim  = LFIT::WhiteDwarf(params.rwd,params.ulimb);
    myCV.bspot = LFIT::BrightSpot(params);
    myCV.disc  = LFIT::Disc(params.q,params.rwd,params.rd,params.dExp,1000);
    
    model.phase.resize(numPhase);
    model.dflux.resize(numPhase);
    model.bflux.resize(numPhase);
    model.wflux.resize(numPhase);
    model.rflux.resize(numPhase);
    model.flux.resize(numPhase);
    
    
    for(size_t i=0; i<numPhase; ++i){
        double phi = startPhase + double(i)*(endPhase-startPhase)/double(numPhase);
        double width = (endPhase-startPhase)/double(numPhase);
        model.phase[i]=(phi+params.phi0);
        model.dflux[i]=(params.dFlux*myCV.disc.calcFlux(params.q,phi,width,params.incl));
        model.bflux[i]=(params.bFlux*myCV.bspot.calcFlux(params.q,phi,width,params.incl));
        model.wflux[i]=(params.wFlux*myCV.prim.calcFlux(params.q,phi,width,params.incl));
        model.rflux[i]=(params.rFlux*myCV.donor.calcFlux(phi,width,params.incl));           
        model.flux[i]=(model.dflux[i]+model.bflux[i]+
                       model.wflux[i]+model.rflux[i]);
    }
    
    write(outfile,model);
    
   }catch(const string& err){
    cerr << "\nError occured inside lfit_fake:" << endl;
    cerr << err << endl;
   }

  return 0;   

}