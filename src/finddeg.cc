/*
 *  finddeg.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 19/06/2008.
 *  Copyright 2008 University of Sheffield. All rights reserved.
 *
 */
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

double finddeg(const double& q, const double& dphi){
	
    double phase = dphi/2.0;
    // white dwarf is a origin
    Subs::Vec3 p = Subs::Vec3(0.0,0.0,0.0);
    
    // if not eclipsed for 90degress, drop out
    double angle = 90.0;
    Subs::Vec3 earth;
    earth=Roche::set_earth(angle,phase);
    if(! Roche::blink(q,p,earth,0.01)){
        std::string err = "no valid solution in finddeg for q=";
		err += Subs::str(q) + " and dphi=" + Subs::str(dphi);
        throw err;
    }
    
    /*
     Choose two angles (54 and 90 degrees) which bracket the range
     Then execute binary chop NLOOP times.
     Final accuracy = 36./2**NLOOP
     = 0.005 deg at worst. On average 1/2 of this.
     */
    double angle1=54.0;
    double angle2=90.0;
    for(int i=0; i<=56; i++){
        angle = 0.5*(angle1+angle2);
        earth=Roche::set_earth(angle,phase);
        if(Roche::blink(q,p,earth,0.01)){
            angle2=angle;
        }else{
            angle1=angle;
        }
    }
    return (angle1+angle2)/2.0;
}
