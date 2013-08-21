/*
 *  lfit_bootstrap.cc
 *  lfit_xcode
 *
 *  Bootstrap resamples a lightcurve by adjusting errors, sets errors on rejected points to 
 *  large value. Also, test setting to -1 or setting 'bad' flag, which should have same effect.
 *
 *  Created by Stuart Littlefair on 24/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */

#include<ctime>
#include "lfit.h"


using namespace std;

LFIT::ThreeColour bootstrap(const LFIT::ThreeColour& inThree, const int& colour){
	
	/* initialize random seed: */
	srand ( (unsigned)time(0) );
	
	LFIT::ThreeColour outThree=inThree;
	LFIT::Data outData = inThree[colour];
	int ndata = outData.size();
	
	
	Subs::Buffer1D<int> bootstrap_int(ndata);
	bootstrap_int = 0;
	
	for(int i=0; i<ndata; i++){
		
		// random number between 0 and NDATA
		int selected = rand()%(ndata+1);
		
		// check point is valid
		if(selected > ndata || selected < 0) {
			string err="Selected data point out of limits";
			throw err;
		}
		
		// increment select count for this point
		bootstrap_int[selected]++;
	}
	
	
	// now adjust error bars (or bad pixel mask) according to number of times a point is selected
	// if point not selected at all, then adjust error bar to large number, or set pixel as bad
	int nsel=0;
	for(int i=0; i<ndata; i++){
		if(bootstrap_int[i] == 0) {
			outData[i].bad = true;
		}else{
			nsel++;
			outData[i].ferr /= double(bootstrap_int[i]);
		}
	}
	//std::cout << "Selected " << nsel << " out of " << ndata << " points" << std::endl;
	outThree[colour]=outData;
	return outThree;
}
