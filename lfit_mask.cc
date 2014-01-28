/*
 *  lfit_mask.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 19/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */
#include <iostream>
#include <cfloat>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/plot.h"
#include "trm/input.h"
#include "lfit.h"

using namespace std;

void wmask(const vector<string>& args, const LFIT::ThreeColour& three){
	
	try{
		
		
		Subs::Input input(args,"LFIT_ENV",".lfit");
		//sign in input variables
		input.sign_in("filename",       Subs::Input::LOCAL, Subs::Input::PROMPT);
		input.sign_in("colour",         Subs::Input::LOCAL,  Subs::Input::PROMPT);
		
		// get input variables
		int colour;
		input.get_value("colour", colour, 2, 1, 3, "lightcurve colour to use [1=r,2=g,3=b]");	
		// make zero offset
		colour--;
		string filename;
		input.get_value("filename", filename, "mask.msk", "mask file to write to");
		
		ofstream fout(filename.c_str());
		if(!fout){
			string error("Failed to open file for writing: ");
			error.append(filename);
			throw error;
		}
		
		LFIT::Data data = three[colour];
		for(size_t i=0; i<data.size(); i++){
			if(data[i].bad){
				fout << 1 << " ";
			}else{
				fout << 0 << " ";
			}
		}
		fout << endl;
		fout.close();
	}catch(string& err){
		cout << "Error occurred inside wmask: " << endl << err << endl;
	}
}

void lmask(const vector<string>& args, LFIT::ThreeColour& three){
	
	try{
		
		
		Subs::Input input(args,"LFIT_ENV",".lfit");
		//sign in input variables
		input.sign_in("filename",       Subs::Input::LOCAL, Subs::Input::PROMPT);
		input.sign_in("colour",         Subs::Input::LOCAL,  Subs::Input::PROMPT);
		
		// get input variables
		int colour;
		input.get_value("colour", colour, 2, 1, 3, "colour to mask [1=r,2=g,3=b]");	
		// make zero offset
		colour--;
		string filename;
		input.get_value("filename", filename, "mask.msk", "mask file to read from");
		
		ifstream fin(filename.c_str());
		if(!fin){
			string error("Failed to open mask file for reading: ");
			error.append(filename);
			throw error;
		}
		
		Subs::Array1D<int> mask;
		int tmp=0;
		while(fin>>tmp){
			if(tmp ==0 or tmp==1)
				mask.push_back(tmp);
		}
		
		if(mask.size() != (int)three[colour].size()) {
			string err = "Mask size of "+ Subs::str(mask.size()) + " is not equal to data size of "
			+ Subs::str(three[colour].size()) + "\n";
			throw err;
		}
		
		for(size_t i=0; i<three[colour].size(); i++){
			if(mask[i] == 1){
				three[colour][i].bad = true;
			}else{
				three[colour][i].bad = false;
			}
		}
		fin.close();
	}catch(string& err){
		cout << "Error occured in lmask: " << endl << err << endl;
	}
}

void smask(const vector<string>& args, LFIT::ThreeColour& three){
	
	try{
		
		Subs::Input input(args,"LFIT_ENV",".lfit");
		//sign in input variables
		input.sign_in("colour",          Subs::Input::LOCAL,  Subs::Input::PROMPT);
		input.sign_in("autoscale",       Subs::Input::LOCAL, Subs::Input::PROMPT);
		input.sign_in("ymin",            Subs::Input::LOCAL, Subs::Input::PROMPT);
		input.sign_in("ymax",            Subs::Input::LOCAL, Subs::Input::PROMPT);
		input.sign_in("clear",           Subs::Input::LOCAL, Subs::Input::PROMPT);
		
		//get input
		bool clear;
		input.get_value("clear", clear, "n", "clear existing mask? (y/n)");
		
		int colour;
		input.get_value("colour", colour, 2, 1, 3, "lightcurve colour to use [1=r,2=g,3=b]");	
		// make zero offset
		colour--;
		
		// clear existing mask if required
		if(clear){
			for(size_t i=0; i<three[colour].size();i++){
				three[colour][i].bad=false;
			}
		}
		
		// set new mask
		string autoscale;
		double ymin=0, ymax=0;
		input.get_value("autoscale", autoscale, "y", "autoscale (y/n)");
		if(Subs::toupper(autoscale) != "Y"){
			input.get_value("ymin",ymin,0.0,-DBL_MAX,DBL_MAX,"minimum y-value for plot");
			input.get_value("ymax",ymax,0.2,-DBL_MAX,DBL_MAX,"maximum y-value for plot");
		}
		
		Subs::Plot plot;
		plot.open("?");
		// pass in null model to plotLight
		LFIT::SynthData model;
		plotLight(model,three[colour],plot,autoscale,ymin,ymax);	

		// fill an array with times
		Subs::Array1D<float> times(three[colour].size());
		for(size_t i=0; i<three[colour].size(); i++){
			times[i]=three[colour][i].time;
		}
		
		vector<pair<int,int> >  masks;
		float xl,yl,xr,yr;
		char cret, reply;
		bool ok;
		cret='M';
		cpgsci(2);
		while(cret != 'Q'){
			ok = (cpgband(6, 0, 0.f, 0.f, &xl, &yl, &cret) == 1);
			cret = std::toupper(cret);
			if(!ok){
				string err = "Error calling cursor";
				throw err;
			}
			if(cret != 'Q'){
				cpgband(4, 0, xl, yl, &xr, &yr, &reply);
				reply = std::toupper(reply);
				if(reply == 'Q') cret=reply;
				cpgsls(2);
				cpgsci(2);
				cpgslw(3);
				cpgmove(xl, ymin);
				cpgdraw(xl, ymax);
				cpgmove(xr, ymin);
				cpgdraw(xr, ymax);
				cpgsls(1);
				cpgmove(xl, (ymin+ymax)/2.);
				cpgdraw(xr, (ymin+ymax)/2.);
				// find nearest index to cursor positions
				int ixl,ixr;
				ixl=Subs::locate(times.ptr(),times.size(),xl);
				ixr=Subs::locate(times.ptr(),times.size(),xr);
				pair<int,int> coord(ixl,ixr);
				masks.push_back(coord);
			}
		}
		
		// apply masks to data
		for(size_t i=0; i<masks.size(); i++){
			for(int j=masks[i].first; j<masks[i].second; j++){
				three[colour][j].bad=true;
			}
		}
		
		plot.close();
		
	}catch(string& err){
		cout << "Error occurred inside smask: " << endl << err << endl;
	}
	
}
