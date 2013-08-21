/*
 *  Disc.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 09/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "Disc.h"
#include <cmath>
#include <string>
#include "trm_subs.h"
#include "trm_array1d.h"
#include "trm_vec3.h"
#include "trm_roche.h"
#include "trm_constants.h"

std::istream& LFIT::operator>>(std::istream& s, const DiscEl& el){
    std::string err = "Calling LFIT::operator >> with DonorEl is an error";
    throw err;
}
std::ostream& LFIT::operator<<(std::ostream& os, const DiscEl& el){
    std::string err = "Calling LFIT::operator << with DonorEl is an error";
    throw err;
}

bool LFIT::Disc::paramsChanged(const double& q, const double& rwd, const double& rdisc,
                   const double& exp_){

	double xl1 = Roche::xl1(q);
	double rinCheck = rwd*xl1; double routCheck = rdisc*xl1; 
	if(this->exp != exp_ or this->rin != rinCheck or this->rout != routCheck){
		this->exp = exp_;
		this->rin = rinCheck;
		this->rout = routCheck;
		return true;
	}	
	return false;
}

void LFIT::Disc::tweak(const double& q, const double& rwd, const double& rdisc,
                       const double& exp_){
    
	if (not paramsChanged(q,rwd,rdisc,exp_)){
		return;
	}
    
    int nrad, ntheta;
    // let's fill the grid
    nrad = int(sqrt(this->tiles.size()));
    ntheta = nrad;
    int icount=-1;
    double delta_r = (this->rout - this->rin)/double(nrad);
    for (int i=0; i<nrad; ++i){
        double r = this->rin + (this->rout - this->rin)*double(i+1)/double(nrad);
        for(int j=0; j<ntheta;++j){
            icount++;
            double theta = Constants::TWOPI*double(j)/double(ntheta);
            this->tiles[icount].pos.set(r*cos(theta),r*sin(theta),0.0);
            this->tiles[icount].norm.set(0.0,0.0,1.0);
            // area is r*dtheta*dr
            this->tiles[icount].area = r*delta_r*Constants::TWOPI/double(ntheta);
            
            // nominal reference flux is used here. The disc's light will later be normalised
            // to 1 out of eclipse
            this->tiles[icount].bright = pow(r/this->rout,-this->exp); 
            // visibility is set later
            this->tiles[icount].vis=0.0;
        }
    }
}

void LFIT::Disc::setVis(const double& q, const double& phi, const double& incl){
    
    Subs::Vec3 earth;
    earth=Roche::set_earth(incl,phi);
    //#pragma omp parallel for 
    for(int i=0; i<this->tiles.size(); ++i){
        double dotP = Subs::dot(earth,this->tiles[i].norm); // projected area
        if(dotP < 0.0){
            std::string err = "Shouldn't get invisible disc elements";
            throw err;
        }
        if(Roche::blink(q,this->tiles[i].pos,earth,0.05)){
            // eclipsed, so not visible
            this->tiles[i].vis = 0.0;
        }else{
            this->tiles[i].vis = fabs(dotP);
        }
    }
}

double LFIT::Disc::calcFlux(const double& q, const double& phi, 
                            const double& width, const double& incl){

    /* 
     computes flux of disc relative to flux outside
     eclipse.
     
     integrates over bin of finite phase width, width using trapezoidal int
     */
    
    double phi1 = phi - width/2.0;
    double dflux=0.0;
    int nphi = 5;
    for(int i=0; i<5; i++){
        double p = phi1 + width*double(i)/double(nphi-1);
        if(i == 0 || i == nphi-1){
            dflux += this->calcFlux(q,p,incl)/2.0;
        }else{
            dflux += this->calcFlux(q,p,incl);
        }
    }
    dflux /= double(nphi-1);
    
    return dflux;
    
    
}

double LFIT::Disc::calcFlux(const double& q, const double& phi, const double& incl){
    
    
    // maximum flux at phase 0.25 
    this->setVis(q,0.25,incl);
    double sum=0.0;
	//#pragma omp parallel for reduction(+:sum) 
    for(int i=0; i< this->tiles.size();i++){
        sum+=this->tiles[i].bright*this->tiles[i].area*this->tiles[i].vis;
    }
    double maxflux = sum;
    
    this->setVis(q,phi,incl);
    sum=0.0;
	//#pragma omp parallel for reduction(+:sum)
    for(int i=0; i< this->tiles.size();i++){
        sum+=this->tiles[i].bright*this->tiles[i].area*this->tiles[i].vis;
    }
    return sum/maxflux;
}

