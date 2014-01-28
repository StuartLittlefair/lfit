#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "trm/subs.h"
#include "trm/buffer2d.h"
#include "trm/array1d.h"
#include "trm/plot.h"
#include "trm/constants.h"
#include "trm/roche.h"

using namespace std;

// global variables
float q, eq,  inc, einc, dphi, edphi;

Subs::Vec3 getEarthVec(float phase, float angle){

  /* returns vector pointing to earth given binary phase
     and inclination angle */
  float rad_incl = angle*2.0*Constants::PI/360.0;
  float rad_phase = phase*2.0*Constants::PI;
  float x = cos(rad_phase)*sin(rad_incl);
  float y = -sin(rad_phase)*sin(rad_incl);
  float z = cos(rad_incl);
  Subs::Vec3 earth = Subs::Vec3(x,y,z);
  return earth;

}

int main(void) {

  try{

    cout << "Give mass ratio, q, and error: ";
    cin >> q >> eq;
    cout << "Give delta phi and error: ";
    cin >> dphi >> edphi;

    Subs::Array1D<float> istore;
    float qsave = q;
    float dphisave = dphi;
    for (int i=0; i<=1000; i++){
      
      // randomise params within errors
      q = qsave + Subs::gauss3(i)*eq;
      dphi = dphisave + Subs::gauss3(i)*edphi;
      
      // check that eclipse of desired width is possible
      // set phase to end of wd eclipse
      float phase = dphi/2.0;
      
      if(q < 0){
	string error("invalid q in incl: q<0\n");
	throw error;
      }
      
      // wd at 0,0,0
      Subs::Vec3 wd=Subs::Vec3();

      // check eclipsed at 90 deg
      float angle=90.0;
      Subs::Vec3 earth = getEarthVec(phase,angle);
      if (! Roche::blink(q,wd,earth,0.1)){
	istore.push_back(angle);
      }else{
	
	/* choose two angles (54 and 90 degrees) which
	   bracket the range. Then execute binary chop NLOOP times.
	   Final accuracy = 36./2**NLOOP
	   => 0.005 deg at worst. On average 1/2 of this.
	*/
	float angle1=54.0; float angle2=90.0;
	for (int j=0; j<=14; j++){
	  angle = (angle1+angle2)/2.0;
	  earth = getEarthVec(phase,angle);
	  if( Roche::blink(q,wd,earth,0.1) ) {
	    angle2=angle;
	  } else {
	    angle1=angle;
	  }
	}
	istore.push_back((angle1+angle2)/2.0);
      }
    }
    float mean = istore.mean();
    istore -= mean;
    float var = Subs::sqr(istore.length())/istore.size();
    cout << "Inclination: " << mean << " +/- " << sqrt(var) << endl;

  }
  catch(const string& e){
    
    cerr << "\nError occured inside incl:" << endl;
    cerr << e << endl;
    return 1;
  }
  
  
  return 0;
  
}
