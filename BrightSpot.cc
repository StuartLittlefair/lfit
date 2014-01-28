/*
 *  brightspot.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 09/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "BrightSpot.h"
#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "trm/constants.h"

double LFIT::BrightSpot::calcFlux(const double& q, const double& phi,
const double& width, const double& incl){
    /* 
     computes flux of bright spot relative to flux outside
     eclipse.
     
     integrates over bin of finite phase width, width using trapezoidal int
     */
    
    double phi1 = phi - width/2.0;
    double bflux=0.0;
    int nphi = 5;
	//#pragma omp parallel for reduction(+:bflux)
    for(int i=0; i<5; i++){
        double p = phi1 + width*double(i)/double(nphi-1);
        if(i == 0 || i == nphi-1){
            bflux += LFIT::BrightSpot::calcFlux(q,p,incl)/2.0;
        }else{
            bflux += LFIT::BrightSpot::calcFlux(q,p,incl);
        }
    }
    bflux /= double(nphi-1);
    
    return bflux;
}

double LFIT::BrightSpot::getTangent() const{
  float alpha = atan2(this->y,this->x);
  if (alpha < 0){
    // alpha is returned between -pi and pi
    // if negative, spotis slightly behind disc (unlikely), i.e alpha > pi
    alpha = 90-alpha;
  }
  return alpha + Constants::TWOPI/4.;
}

double LFIT::BrightSpot::simpleFlux(const double& q,
const double& phi, const double& incl){
   /*
     * Computes light curve of spot assuming that it is a linear feature
     * with an intensity that varies as (X/LS)**N*EXP(-X/LS) where X
     * is measured from one end. The peak at X = N*LS is positioned on the
     * at XS, YS and the line is rotated relative to X,Y axes by AZ
     * Normalised to equal 1 at maximum if no eclipse.
     *
	 * updated 2009-02-18 to allow N to be 1, 2 or 3
     * PHI    = Orbital phase
     * Q      = Mass ratio = M2/M1
     * INCL   = Inclination angle (degrees)
     
     and from spot object itself we get
     * XS, YS, SCALE = position and scale of spot/RL1
     * frac    = Fraction of light which is isotropic.
     * AZ     = azimuth in degrees measured anti-clockwise from line of centres
     */

    const int ept=2; // ept cannot be changed at will but must be 1,2 or 3
    const double acc=1.0e-3;
    const double tot=8.0;
    double bflux;
    
    // compute some time saving variables
    Subs::Vec3 earth;
    earth=Roche::set_earth(incl,phi);
    double cphi = cos(Constants::TWOPI*phi);
    double sphi = sin(Constants::TWOPI*phi);
    double raz  = Constants::TWOPI*this->az/360.0;
    double caz  = cos(raz);
    double saz  = sin(raz);
    
    // ept is exponent in x**ept*exp(-x)
    // integration is carried out from x=0 to x=tot
	double out = 2.0-(tot*(tot+2.0)+2.0)*exp(-tot);
	if(ept == 1) out = 2.0*(1.0-exp(-tot));
	if(ept == 3) out = 3.0*out - pow(tot,3)*exp(-tot);
    double x1 = this->x - ept*this->scale*caz;
    double y1 = this->y - ept*this->scale*saz;
    double dx = tot*this->scale*caz;
    double dy = tot*this->scale*saz;
    double proj = std::max(0.0,saz*cphi+caz*sphi);
    
    // are ends of element occulted?
    Subs::Vec3 x = Subs::Vec3(x1,y1,0.0);
    bool p1,p2;
    if(Roche::blink(q,x,earth,0.05)){
        p1 = true;
    }else{
        p1 = false;
    }
    x.set(x1+dx,y1+dy,0.0);
    if(Roche::blink(q,x,earth,0.05)){
        p2 = true;
    }else{
        p2 = false;
    }
    if(p1 && p2){
        bflux = 0.0;
    }else if(!p1 && !p2){
        bflux = this->frac + (1.0-this->frac)*proj;
    }else{
        double a, a1,a2;
        // a1 is the occulted end
        if(p1){
            a1=0.0;
            a2=1.0;
        }else{
            a1=1.0;
            a2=0.0;
        }
        do {
            a=(a1+a2)/2.0;
            x.set(x1+a*dx,y1+a*dy,0.0);
            if(Roche::blink(q,x,earth,0.05)){
                a1=a;
            }else{
                a2=a;
            }
        }while(fabs(a2-a1) > acc);
        a = tot*(a1+a2)/2.0;
        double frac = 2.0-(a*(a+2.0)+2.0)*exp(-a);
		if(ept==1) frac = 2.0*(1.0-exp(-a));
		if(ept==3) frac = 3*frac - pow(a,3)*exp(-a);
        if(p1){
            bflux=(this->frac+(1.0-this->frac)*proj)*(1.0-frac/out);
        }else{
            bflux=(this->frac+(1.0-this->frac)*proj)*frac/out;
        }
    }
    return bflux;
}



double LFIT::BrightSpot::calcFlux(const double& q,
const double& phi, const double& incl){
	
	// SIMPLE: do old-style LFIT calculation with single bright spot strip hard-coded to exponent of 2
	if (!complex){return LFIT::BrightSpot::simpleFlux(q,phi,incl);}
	
	// NOT SIMPLE - shamelessly ripped off from Tom Marsh's LCURVE
	double theta = Constants::TWOPI*this->az/360.0;
	double alpha = Constants::TWOPI*this->yaw/360.0;
	double tilt_spot  = Constants::TWOPI*this->tilt/360.0;
	
	// the direction of the bright spot line is set by az, but the beaming direction adds in yaw as well
	Subs::Vec3 posn, dirn, bvec(cos(theta),sin(theta),0), earth;
	Subs::Vec3 pvec(0,0,1), tvec(sin(tilt_spot)*sin(theta+alpha), -sin(tilt_spot)*cos(theta+alpha), cos(tilt_spot));
	
	// position of bright spot
	Subs::Vec3 bspot(this->x,this->y,0);
	
	// length of bright spt in scale lengths
	const double BMAX = pow(this->exp1/this->exp2,1.0/this->exp2);
	const double SFAC = 20.0 + BMAX;
	
	// first find maximum bright spot flux phase - assumed to be out of eclipse, and when spot is facing towards us.
	double tanMaxPhi = 1.0/tan(theta+alpha);
	double MaxPhi = atan(tanMaxPhi)/Constants::PI/2.0; //lies between -1/2 and 1/2
	if(MaxPhi < 0.0){MaxPhi = 1.0+MaxPhi;}
	
	// calculate max flux
	earth=Roche::set_earth(incl,MaxPhi);
	double projTilted   = std::max(0.0,Subs::dot(tvec,earth));
	double maxProj = projTilted;
	int nspot=300;
	Subs::Vec3 xp;
	
	double maxflux=0.0;
	for(int i=0; i<nspot; i++){
		double dist = SFAC*i/(nspot-1);
		xp = bspot + this->scale*(dist-BMAX)*bvec;
		// assume Bspot out of eclipse, so no need to blink here
		double bright = pow(dist/BMAX,this->exp1)*exp(this->exp1/this->exp2 - pow(dist,this->exp2));
		maxflux += projTilted*bright*(1.0-this->frac);
		maxflux += maxProj*bright*this->frac;			
	}
	//std::cout << "Max flux is: " << maxflux << " , at phase " << MaxPhi << std::endl;
	
	// now find actual flux at this phase
	earth=Roche::set_earth(incl,phi);
	projTilted   = std::max(0.0,Subs::dot(tvec,earth));
	double projParallel = std::max(0.0,Subs::dot(pvec,earth));
	double bflux=0.0;

	for(int i=0; i<nspot; i++){
		double dist = SFAC*i/(nspot-1);
		xp = bspot + this->scale*(dist-BMAX)*bvec;
		
		if(!Roche::blink(q,xp,earth,0.05)){
			double bright = pow(dist/BMAX,this->exp1)*exp(this->exp1/this->exp2 - pow(dist,this->exp2));
			if(projTilted > 0.0) {bflux += projTilted*bright*(1.0-this->frac);}
			if(projParallel > 0.0) {bflux += maxProj*bright*this->frac;}		// this should be multiplied by something so frac is truly correct	
		}
	}
	return bflux/maxflux;
}

void LFIT::BrightSpot::spotPos(const double& q, const double& rd){
    /*     
     Determines location of bright-spot from mass ratio and radius
     of disc assuming that it is located on free-particle trajectory
     at edge of disc.
     
     a = Mass ratio = M2/M2
     RD = radius of disc/XL1
     Calculates X,Y position of spot in units of XL1
     */
    static double eps = 1.0e-8;
    static double acc = 1.0e-5;
    
    double smax,vel,ttry,tdid,time,tnext;
    
    double xl1 = Roche::xl1(q);
    double rtest = rd*xl1;
    
    Subs::Vec3 r,v;
    Roche::strinit(q,r,v);
    time = 0.0;
    
    Subs::Vec3 r0(r),v0(v); // store original pos and vel.
    
    //integrate stream until inside disc
    double rold = xl1;
    double rad;
    double step=1.0e-3;
    ttry = 1.0e-3;
    smax = std::min(1.0e-3,step/2.0);
    while(true){
        Roche::gsint(q,r,v,ttry,tdid,tnext,time,eps);
        vel = sqrt(Subs::sqr(v.x()) + Subs::sqr(v.y()));
        ttry = std::min(smax/vel,tnext);
        
        rad = sqrt(Subs::sqr(r.x()) + Subs::sqr(r.y()));
        if(rad > rold){
            std::string err="error in spotPos";
            throw err;
        }else if(rad>rtest){
            rold=rad;
            r0=r;
            v0=v;
        }else{
            break;
        }
    }
    
    // now we have one point with r > rd and one with r < rd
    // with a time step of tdid separating them
    double d1=0.0;
    double d2=tdid;
    double r1=rold;
    double r2=rad;
    ttry = tdid;
    while(true){
        r=r0;
        v=v0;
        double d=(d1+d2)/2.0;
        Roche::gsint(q,r,v,d,tdid,tnext,time,eps);
        vel = sqrt(Subs::sqr(v.x()) + Subs::sqr(v.y()));
        ttry = std::min(smax/vel,tnext);
    
        rad = sqrt(Subs::sqr(r.x()) + Subs::sqr(r.y()));
        if(rad > rtest){
            d1=tdid;
            r1=rad;
        }else{
            d2=tdid;
            r2=rad;
        }
        if(fabs(r1-r2) < acc){
            break;
        }
    }
    this->x = r.x();
    this->y = r.y();
}
