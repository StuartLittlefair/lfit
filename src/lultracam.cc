/*
 *  lultracam.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 10/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <cfloat>
#include <string>
#include "lultracam.h"
#include "trm/subs.h"
#include "trm/input.h"
#include "trm/array1d.h"
#include "trm/plot.h"
#include "trm/telescope.h"
#include "trm/position.h"
#include "trm/format.h"

void LFIT::split_on(const std::string& str, const std::string& delim,
              std::string& str1, std::string& str2){
    Subs::strip_trailing_whitespace(str);
    // split on equals
    std::string::size_type pos = str.find(delim);
    str1 = str.substr(0,pos-1);
    str2 = str.substr(pos+1);
    str2 = Subs::strip_leading_whitespace(str2);
    str2 = Subs::strip_trailing_whitespace(str2);
    str1 = Subs::strip_leading_whitespace(str1);
    str1 = Subs::strip_trailing_whitespace(str1);
}

bool LFIT::operator<(LFIT::Datum x, LFIT::Datum y){
    if(x.time < y.time) return true;
    return false;
}

void LFIT::ucam_read(LFIT::Data& red, LFIT::Data& grn, LFIT::Data& blu){
    
    LFIT::Input input;
    LFIT::read_infile(input);

    LFIT::Data targ1,targ2,targ3;
    LFIT::Data comp1,comp2,comp3;
    
    std::cout << "Reading data from " << input.files.size() << " files\n";

    /*
     now we have to read in data from each file,
     process each lightcurve to divide by standard,
     correct times to HJD, phasefold, sort in phase
     and cut down to requested phase range!
     */

    bool success=true;
    if(input.format == "ASCII") {
      if (input.files.size() != 1){
	std::string error = "Reading more than one ascii file is not supported!";
	throw error;
      }
      std::cout << "Reading " << input.files[0].fname << std::endl;
      success = LFIT::read_lcurve(targ1,input.files[0].fname) && success;
      std::cout << "...processing data" << std::endl;
      LFIT::process_ascii_data(input,targ1,red);
      LFIT::process_ascii_data(input,targ1,grn);
      LFIT::process_ascii_data(input,targ1,blu);
      targ1.clear();
    }else if(input.format == "LCURVE"){
      if (input.files.size() != 1){
	std::string error = "Reading more than one LCURVE format file is not supported!";
	throw error;
      }
      std::cout << "Reading " << input.files[0].fname << std::endl;
      success = LFIT::read_lcurve(targ1,input.files[0].fname) && success;
      red = targ1;
      grn = targ1;
      blu = targ1;
      targ1.clear();
    }else{
      for(int i=0; i<input.files.size();i++){
	std::cout << "Reading " << input.files[i].fname << std::endl;
	success = LFIT::loadultracam(targ1,1,
				     input.files[i].targap,
				     input.files[i].fname) && success;
	success = LFIT::loadultracam(targ2,2,
				     input.files[i].targap,
				     input.files[i].fname) && success;
	success = LFIT::loadultracam(targ3,3,
				     input.files[i].targap,
				     input.files[i].fname) && success;
	success = LFIT::loadultracam(comp1,1,
				     input.files[i].compap,
				     input.files[i].fname) && success;
	success = LFIT::loadultracam(comp2,2,
				     input.files[i].compap,
				     input.files[i].fname) && success;
	success = LFIT::loadultracam(comp3,3,
				     input.files[i].compap,
				     input.files[i].fname) && success;
	
	std::cout << "...processing data" << std::endl;
	// do all the data processing, and append to rgb lcurves
	LFIT::process_ucam_data(input,targ1,comp1,i,red);
	LFIT::process_ucam_data(input,targ2,comp2,i,grn);
	LFIT::process_ucam_data(input,targ3,comp3,i,blu);
	//clean up
	targ1.clear(); targ2.clear(); targ3.clear();
	comp1.clear(); comp2.clear(); comp3.clear();
      }
    }
    
    if(!success){
      std::string err = "Failed to read data from logfiles";
      throw err;
    }
    
    std::cout << "Red lcurve computed with " << red.size() << " data points" << std::endl;
    std::cout << "Green lcurve computed with " << grn.size() << " data points" << std::endl;
    std::cout << "Blue lcurve computed with " << blu.size() << " data points" << std::endl;
	
    std::cout << "...trimming data" << std::endl;
    //remove data which is outside the requested min, max phases.
    LFIT::trim(red,input);
    LFIT::trim(grn,input);
    LFIT::trim(blu,input);

    std::cout << "Red lcurve computed with " << red.size() << " data points" << std::endl;
    std::cout << "Green lcurve computed with " << grn.size() << " data points" << std::endl;
    std::cout << "Blue lcurve computed with " << blu.size() << " data points" << std::endl;

    std::cout << "...sorting data" << std::endl;
    // sort array in phase
    LFIT::psort(red);
    LFIT::psort(grn);
    LFIT::psort(blu);
    std::cout << "Red lcurve computed with " << red.size() << " data points" << std::endl;
    std::cout << "Green lcurve computed with " << grn.size() << " data points" << std::endl;
    std::cout << "Blue lcurve computed with " << blu.size() << " data points" << std::endl;

    std::cout << "...finding data limits" << std::endl;
    LFIT::setlims(red);
    LFIT::setlims(grn);
    LFIT::setlims(blu);
    
    // done!
    std::cout << "Red lcurve computed with " << red.size() << " data points" << std::endl;
    std::cout << "Green lcurve computed with " << grn.size() << " data points" << std::endl;
    std::cout << "Blue lcurve computed with " << blu.size() << " data points" << std::endl;
}

void LFIT::setlims(LFIT::Data& inData){
    double x1=1.0e30,x2=-1.0e30,y1=1.0e30,y2=-1.0e30;
    for (size_t i=0; i < inData.size(); i++ ){
        if(inData[i].time < x1) x1 = inData[i].time;
        if(inData[i].time > x2) x2 = inData[i].time;
        if(inData[i].flux < y1) y1 = inData[i].flux;
        if(inData[i].flux > y2) y2 = inData[i].flux;
    }
    inData.set_xmin(x1);
    inData.set_xmax(x2);
    inData.set_ymin(y1);
    inData.set_ymax(y2);
}

void LFIT::trim(LFIT::Data& inData, const LFIT::Input& input){
    
    LFIT::Data tmp;
    for (size_t i=0; i < inData.size(); i++ ){
    	//std::cout << inData[i].time << std::endl;
        if(inData[i].time < input.xmax && inData[i].time > input.xmin){
            tmp.push_back(inData[i]);
        }
    }
    inData = tmp;
}

void LFIT::psort(LFIT::Data& inData){
    // sorts a Data into increasing time (phase) order

    Subs::Buffer1D<double> tmp;
    Subs::Buffer1D<int> key;
    for(size_t i=0; i<inData.size(); i++){
        tmp.push_back(inData[i].time);
    }
    Subs::heaprank(tmp,key);
    LFIT::Data tmpdata;
    for(size_t i=0; i<inData.size(); i++){
        tmpdata.push_back(inData[key[i]]);
    }
    inData = tmpdata;
}

void LFIT::process_ucam_data(const LFIT::Input& input,const LFIT::Data& targ,
const LFIT::Data& comp,const int& fileNum,LFIT::Data& out){
    
    // takes the comparison star and target data objects and appends processed data
    // onto supplied LFIT::Data object out.
    if(targ.size() != comp.size()){
        std::string err = "differing lightcurve sizes for target and comparison\n";
        throw err;
    }
	
	// need to allocate size for array outside of parallel as dynamic filling breaks
	int sizeIn = out.size()-1;
	Subs::Buffer1D<int> deRef;
	for(size_t i=0; i<targ.size(); i++){
		if(! targ[i].bad && ! comp[i].bad){
			out.push_back(targ[i]);
			deRef.push_back(i);
		}
	}
	
	#pragma omp parallel for
    for(int i=0; i<deRef.size(); i++){

		int point = i+sizeIn+1;
		
		LFIT::Datum targtmp = targ[deRef[i]];
		LFIT::Datum comptmp = comp[deRef[i]];
		
		// divide by standard & convert to mags
		double magTarg = -2.5*log10(targtmp.flux/comptmp.flux);
		switch(targ.chip){
			case RED:
				magTarg += input.files[fileNum].comprMag;
				break;
			case GRN:
				magTarg += input.files[fileNum].compgMag;
				break;
			case BLU:
				magTarg += input.files[fileNum].compuMag;
				break;
		}
		double fracErr = sqrt(Subs::sqr(targtmp.ferr/targtmp.flux) + 
							  Subs::sqr(comptmp.ferr/comptmp.flux));
		
		// correct the times for light travel corrections etc
		Subs::Time time = Subs::Time(targtmp.time);
		Subs::Pinfo pinfo = input.objPos.pinfo(time,input.site);
		double tcorr = pinfo.tcor_hel/Constants::DAY;
		targtmp.time += tcorr;
		
		// no need to extinction correct as standard is on same frame
		// and is thus taken at same airmass...
		
		// convert back into fluxes!
		targtmp.flux = 3.631e6 * pow(10.0,-0.4*magTarg);
		
		// correct errors
		targtmp.ferr = targtmp.flux*fracErr;

		// now phase fold
		double MJD0 = input.hjd0 - 2400000.5;
		double phase = (targtmp.time - MJD0)/input.per;
		targtmp.time = phase - floor(phase);
		if(targtmp.time < 0.5) targtmp.time += 1.0;
		
		// correct exposure time to same units
		targtmp.expose /= input.per;
		
		out[point]=targtmp;
	}

}

void LFIT::process_ascii_data(const LFIT::Input& input, const LFIT::Data& targ,  LFIT::Data& out){
	
	for(int i=0; i<targ.size(); i++){

		LFIT::Datum targtmp = targ[i];
		// now phase fold
		double MJD0 = input.hjd0 - 2400000.5;
		double phase = (targtmp.time - MJD0)/input.per;
		targtmp.time = phase - floor(phase);
		if(targtmp.time < 0.5) targtmp.time += 1.0;
	
		// correct exposure time to same units
		targtmp.expose /= input.per;
		out.push_back(targtmp);
	}
}

void LFIT::read_infile(LFIT::Input& input){
    
    //open log file
    std::ifstream fin;

	// get log file name
	std::string file;
	std::vector<std::string> dummyArgs;
	dummyArgs.push_back("");
	Subs::Input fileIn(dummyArgs,"LFIT_ENV",".lfit");
	fileIn.sign_in("file", Subs::Input::LOCAL,  Subs::Input::PROMPT);
	fileIn.get_value("file",file,"input.dat","Input file name");
	

    fin.open(file.c_str());
    if(!fin) {
        string error("Failed to open file: ");
        error.append(file);
        throw error;
    }
    
    std::map <std::string, std::string> entries; 
    //read input file
    while(fin){
        char c = fin.peek();
        if(c == '#'){
            //skip comment and blank lines
            fin.ignore(100000,'\n');
        }else{
            // read in data and parse. 
            std::string str,str1,str2;
            getline(fin,str);
            if(str != ""){
                LFIT::split_on(str,"=",str1,str2);
                str1 = Subs::toupper(str1);
                entries[str1] = str2;
            }
        }
    }
    
    // now unpack entries into required information
    int nfiles=0;
    std::vector<std::string> params;
    params.push_back("TELESCOPE"); params.push_back("CO-ORDINATES");
    params.push_back("HJD0"); params.push_back("PERIOD");
    params.push_back("XMIN"); params.push_back("XMAX");
    params.push_back("NFILES");
	params.push_back("FORMAT");
    for(size_t i=0; i < params.size(); i++){
        if(entries.find(params[i]) == entries.end()){
            std::string err = params[i] + " entry not found in " + file + "\n";
            throw err;   
        } else {
            if(params[i] == "TELESCOPE"){
                input.site = Subs::Telescope("WHT");
            }
            if(params[i] == "CO-ORDINATES"){
                input.objPos = Subs::Position(entries[params[i]]);
            }
            if(params[i] == "HJD0")
                input.hjd0 = Subs::string_to_double(entries[params[i]]);
            if(params[i] == "FORMAT")
                input.format = entries[params[i]];
            if(params[i] == "PERIOD")
                input.per = Subs::string_to_double(entries[params[i]]);
            if(params[i] == "XMIN")
                input.xmin = Subs::string_to_double(entries[params[i]]);
            if(params[i] == "XMAX")
                input.xmax = Subs::string_to_double(entries[params[i]]);
            if(params[i] == "NFILES")
                nfiles = Subs::string_to_int(entries[params[i]]);
        }
    }
    
    // now have to loop through and extract information for each logfile
    params.clear();
    for(int i=1; i<=nfiles; i++){
        LFIT::FileInfo info;        
        // create a vector of appropriate keys to search
        params.push_back("FILE"+Subs::str(i)); 
        params.push_back("REDFILTER"+Subs::str(i));
        params.push_back("TARGAP"+Subs::str(i));
        params.push_back("COMPAP"+Subs::str(i));
        params.push_back("COMPIFLUX"+Subs::str(i));
        params.push_back("COMPGFLUX"+Subs::str(i));
        params.push_back("COMPUFLUX"+Subs::str(i));
        // loop over all entries and push into this file's info
        for(size_t j=0; j<params.size(); j++){
            if(entries.find(params[j]) == entries.end()){
                std::string err = params[j] + " entry not found in " + file + "\n";
                throw err;   
            } else {
                if(params[j]=="FILE"+Subs::str(i))
                    info.fname = entries[params[j]];
                if(params[j]=="REDFILTER"+Subs::str(i))
                    info.redfilt = entries[params[j]];
                if(params[j]=="TARGAP"+Subs::str(i))
                    info.targap = Subs::string_to_int(entries[params[j]]);
                if(params[j]=="COMPAP"+Subs::str(i))
                    info.compap = Subs::string_to_int(entries[params[j]]);
                if(params[j]=="COMPIFLUX"+Subs::str(i))
                    info.comprMag = Subs::string_to_double(entries[params[j]]);
                if(params[j]=="COMPGFLUX"+Subs::str(i))
                    info.compgMag = Subs::string_to_double(entries[params[j]]);
                if(params[j]=="COMPUFLUX"+Subs::str(i))
                    info.compuMag = Subs::string_to_double(entries[params[j]]);
            }
        }
        input.files.push_back(info);
        params.clear();
    }
}

bool LFIT::read_lcurve(LFIT::Data& lcurve, std::string file){
	
	//open log file
    std::ifstream fin;
    fin.open(file.c_str());
    if(!fin) {
        string error("Failed to open file: ");
        error.append(file);
        throw error;
    }
	
    lcurve.clear();
    while(fin){
      char c = fin.peek();
      if(c == '#' || c == ' '){
	//skip header lines
	fin.ignore(100000,'\n');
      }else{
	LFIT::Datum tmp;
	double weight;
	fin >> tmp.time >> tmp.expose >> tmp.flux >> tmp.ferr >> weight;
	tmp.bad = false;
	lcurve.push_back(tmp);
      }
    }
    return true;
}

bool LFIT::loadultracam(LFIT::Data& lcurve, int nccd, int naper, std::string file){
  
  
  /* this is actually the bit that reads in a log file 
   */
  
  int nbad=0, nfail=0, ngood=0;
  lcurve.clear();
  try{
    
    switch(nccd){
    case(1):
      lcurve.chip = RED;
      break;
    case(2):
      lcurve.chip = GRN;
      break;
    case(3):
      lcurve.chip = BLU;
      break;
    }
    
    //open log file
    std::ifstream fin;
    fin.open(file.c_str());
    if(!fin) {
      string error("Failed to open file: ");
      error.append(file);
      throw error;
    }

    string version;
    int intVer;
    //read log file
    while(fin){
      char c = fin.peek();
      if(c == '#'){
	//parse header lines
	string header;
	getline(fin,header);
	if(header.find("software version") != string::npos){
	  version = header.substr(37,header.length());
	  size_t j;
	  for ( ; (j = version.find(".")) != string::npos ; ) {
	    version.erase( j, 1);
	  }
	  intVer = Subs::string_to_int(version);
	}
      }else{
	// genuine data. read it in.
	std::string name;
	double time, expose;
	int flag, nsat, ccd;
	float fwhm, beta;
	if (intVer < 8136){
	  fin >> name >> time >> flag >> nsat >> expose >> ccd >> fwhm >> beta;
	}else{
	  fin >> name >> time >> flag >> expose >> ccd >> fwhm >> beta;
	}
	if(fin && ccd == nccd){
	  // nccd tells us which ccd to put in lcurve
	  int nape = 0, nrej, error_flag, nsky, worst = 0;
	  double counts, sigma, sky;
	  double xf, yf, xm, ym, ex, ey;
	  while(nape != naper && fin)
	    fin >> nape >> xf >> yf >> xm >> 
	      ym >> ex >> ey >> counts >> 
	      sigma >> sky >> nsky >> 
	      nrej >> worst >> error_flag;
	  if(!fin){
	    string err = 
	      "Error finding aperture number " + 
	      Subs::str(naper) + " in CCD " + Subs::str(ccd);
	    throw err;
	  }
          
	  fin.ignore(100000, '\n');
          
	  // store data, including duff data to keep it in synchronism with
	  // other reads of the same data file
	  LFIT::Datum tmp;
	  tmp.time = time;
	  tmp.expose = expose/Constants::DAY;
	  tmp.flux = counts/expose;
	  tmp.ferr = sigma/expose;
	  tmp.bad  = false;
	  if(worst > 0){
	    // this point ruined by bad pixels
	    tmp.bad=true;
	    nbad++;
	  }
	  if(tmp.flux <= 0.0){
	    tmp.bad = true;
	  }
	  if(sigma < 0){
	    // negative errors mean a fatal error code encountered
	    tmp.bad=true;
	    nfail++;
	  }else if(worst <= 0){
	    ngood++;
	  }
	  // push the datum onto appropriate Data object
	  lcurve.push_back(tmp);
	}else{
	  fin.ignore(100000, '\n');
	}
      }
    }
    return true;
  }catch (const std::string& err) {
    std::cerr << "\nError occurred inside lultracam. Detailed error message follows:" << std::endl;
    std::cerr << err << std::endl;
    return false;
  }
  return true;
}

std::istream& LFIT::operator>>(std::istream& s, const FileInfo& el){
    std::string err = "Error: FileInfo structures are not designed to be read in\n";
    throw err;
}
std::ostream& LFIT::operator<<(std::ostream& os, const FileInfo& el){
    std::string err = "Error: FileInfo structures are not designed to be written out\n";
    throw err;
}

std::istream& LFIT::operator>>(std::istream& s, const Datum& el){
    std::string err = "Error: Datum structures are not designed to be read in\n";
    throw err;
}
std::ostream& LFIT::operator<<(std::ostream& os, const Datum& el){
	std::cout << el.time << ", " << el.flux << ", " << el.ferr << std::endl;
}

std::istream& LFIT::operator>>(std::istream& s, const Data& el){
    std::string err = "Error: Data structures are not designed to be read in\n";
    throw err;
}
std::ostream& LFIT::operator<<(std::ostream& os, const Data& el){
    std::string err = "Error: Data structures are not designed to be written out\n";
    throw err;
}
    


