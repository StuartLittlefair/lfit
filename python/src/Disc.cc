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
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "trm/constants.h"

bool LFIT::Disc::paramsChanged(const double& q_, const double& rwd, const double& rdisc,
                   const double& exp_){

	double xl1 = Roche::xl1(q_);
	double rinCheck = rwd*xl1; double routCheck = rdisc*xl1; 
	if(this->exp != exp_ or this->rin != rinCheck or this->rout != routCheck){
	    this->q  = q_;
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
    this->tiles.clear();
}

void LFIT::Disc::setup_grid(const double& incl){

    int nrad, ntheta;

    // let's create the grid
    // we need more elements at smaller radii, since that's where most flux comes
    // from, but total number of elements must equal ntiles.
    // if ntheta = 2*nrad - i
    // total number of tiles = 0.5*nrad*(3*nrad+1)
    // nrad = (-1 + sqrt(1 + 24*ntiles) ) / 6
    nrad = int((-1 + sqrt(1 +24*this->ntiles))/6);
    

    // now let's set it up
    int icount=-1;
    double delta_r = (this->rout - this->rin)/double(nrad);
    LFIT::Point::etype eclipses;
    for (int i=0; i<nrad; ++i){
        double r = this->rin + (this->rout - this->rin)*double(i+1)/double(nrad);
        
        ntheta = 2*nrad - i;
        for(int j=0; j<ntheta;++j){
            icount++;
            double theta = Constants::TWOPI*double(j)/double(ntheta);
                
            Subs::Vec3 posn = Subs::Vec3(r*cos(theta),r*sin(theta),0.0);
            Subs::Vec3 dirn = Subs::Vec3(0.0,0.0,1.0);

            // area is r*dtheta*dr
            double area = r*delta_r*Constants::TWOPI/double(ntheta);
            
            // nominal reference flux is used here. The disc's light will later be normalised
            // to 1 out of eclipse
            double flux = pow(r/this->rout,-this->exp); 

            // ingress, egress phases
            eclipses.clear();
            double ingress, egress;
            if (Roche::ingress_egress(this->q, Roche::SECONDARY, 1.0, 1.0, incl, 1.0e-5, posn, ingress, egress)){
                eclipses.push_back(std::make_pair(ingress,egress));
            }     
            
            this->tiles.push_back(LFIT::Point(posn,dirn,area,eclipses));
            this->tiles[icount].flux = flux;        
        }
    }
}

double LFIT::Disc::calcFlux(const double& q, const double& phi,const double& width, const double& incl){

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
    
    double static maxflux;
    bool static first=true;

    // have we been called without the tiles calculated
    if (this->tiles.size() == 0){
        std::cout << "LFIT::Disc::Calcflux shouldn't be called before calculating grid.\nThis is inefficient" << std::endl;
        this->setup_grid(incl);   
    }
        
    // maximum flux at phase 0.25
    if(first){
        double maxphi = 0.25;
        Subs::Vec3 earth = Roche::set_earth(incl,maxphi);
        double sum = 0.0; 
        for(int i=0; i< this->tiles.size();i++){
            double mu = Subs::dot(earth,this->tiles[i].dirn);
            if(mu > 0. && this->tiles[i].visible(maxphi)){
                sum+=this->tiles[i].flux*this->tiles[i].area*mu;
            }
        }
        first=false;
        maxflux = sum;
    }
    
    Subs::Vec3 earth = Roche::set_earth(incl,phi);
    double sum = 0.0; 
    for(int i=0; i< this->tiles.size();i++){
        double mu = Subs::dot(earth,this->tiles[i].dirn);
        if(mu > 0. && this->tiles[i].visible(phi)){
            sum+=this->tiles[i].flux*this->tiles[i].area*mu;
        }
    }
    return sum/maxflux;
}

