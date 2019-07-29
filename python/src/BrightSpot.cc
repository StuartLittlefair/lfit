/*
 *  brightspot.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 09/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "BrightSpot.h"
#include "Point.h"
#include <cmath>
#include <string>
#include <stdexcept>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "trm/constants.h"

void LFIT::BrightSpot::setNspot(int nelem){
    this->nspot = nelem;
}

int LFIT::BrightSpot::getNspot(){
    return this->nspot;
}

double LFIT::BrightSpot::calcFlux(const double& q, const double& phi, const double& width, const double& incl){
    /*
     computes flux of bright spot relative to flux outside
     eclipse.

     integrates over bin of finite phase width, width using trapezoidal int
     */

    double phi1 = phi - width/2.0;
    double bflux=0.0;
    int nphi = 5;
	//#pragma omp parallel for reduction(+:bflux)
    for(int i=0; i<nphi; i++){
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

double LFIT::BrightSpot::simpleFlux(const double& q,const double& phi, const double& incl){
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

void LFIT::BrightSpot::setup_grid( const double& incl ){
    if (!this->complex){
        return;
    }

	// NOT SIMPLE - shamelessly ripped off from Tom Marsh's LCURVE
	double theta = Constants::TWOPI*this->az/360.0;
	double alpha = Constants::TWOPI*this->yaw/360.0;
	double tilt_spot  = Constants::TWOPI*this->tilt/360.0;

	// the direction of the bright spot line is set by az, but the beaming direction adds in yaw as well
	Subs::Vec3 posn, bvec(cos(theta),sin(theta),0);
	Subs::Vec3 pvec(0,0,1), tvec(sin(tilt_spot)*sin(theta+alpha), -sin(tilt_spot)*cos(theta+alpha), cos(tilt_spot));

	// position of bright spot
	Subs::Vec3 bspot(this->x,this->y,0);

	// length of bright spt in scale lengths
	const double BMAX = pow(this->exp1/this->exp2,1.0/this->exp2);
	const double SFAC = 20.0 + BMAX;

    // the scaling of nspot below is tuned to give reasonable results
    this->nspot = std::max(200,int(200*this->scale/0.02));

	// create buffer of point objects
	int nspot = this->nspot;
	this->spot.resize(nspot*2);

	// setup point objects, including their ingress/egress phases
	const double AREA   = SFAC*this->scale*this->scale/(nspot-1);
    LFIT::Point::etype eclipses;
    std::cout << "setup_grid: initialisation done, looping over " << nspot << " elements" << std::endl;
	for(int i=0; i<nspot; i++){
        std::cout << "starting el " << i << "...";
        // spot position
	    double dist = SFAC*i/(nspot-1);
	    posn = bspot + this->scale*(dist-BMAX)*bvec;
        std::cout << "pos " << posn << "...";

	    // ingress, egress phases
	    eclipses.clear();
        double ingress, egress;
        std::cout << "eclipses cleared...";

        // sometimes dies with linmin error: allow this to bubble up
        if (Roche::ingress_egress(q, Roche::SECONDARY, 1.0, 1.0, incl, 1.0e-5, posn, ingress, egress)){
            eclipses.push_back(std::make_pair(ingress,egress));
            std::cout << "(" << ingress << ", " << egress << ")...";
        }else{
            std::cout << "no eclipses - " << "(" << ingress << ", " << egress << ")...";
        }

        // Factor here is adjusted to equal 1 at its peak
        double bright = pow(dist/BMAX,this->exp1)*exp(this->exp1/this->exp2 - pow(dist,this->exp2));

        // the tilted strip
        this->spot[i]      = LFIT::Point(posn,tvec,AREA,eclipses);
        this->spot[i].flux = bright*(1.0-this->frac)*this->spot[i].area;

        // the parallel strip
        this->spot[i+nspot]      = LFIT::Point(posn,pvec,AREA,eclipses);
        this->spot[i+nspot].flux = bright*this->frac*this->spot[i].area;
        std::cout << "done element " << i << " of " << nspot << std::endl;
    }
}

double LFIT::BrightSpot::calcFlux(const double& q,const double& phi, const double& incl){

	// SIMPLE: do old-style LFIT calculation with single bright spot strip hard-coded to exponent of 2
	if (!this->complex){return LFIT::BrightSpot::simpleFlux(q,phi,incl);}

    //  IS  CALLED WITH SPOT UNCALCULATED?
    if (this->spot.size() == 0){
        std::cout << "This shouldn't happen" <<std::endl;
        this->setup_grid(incl);
    }

    // first find maximum bright spot flux phase - assumed to be out of eclipse, and when spot is facing towards us.
    double theta = Constants::TWOPI*this->az/360.0;
	double alpha = Constants::TWOPI*this->yaw/360.0;
	double tanMaxPhi = 1.0/tan(theta+alpha);
	double MaxPhi = atan(tanMaxPhi)/Constants::PI/2.0; //lies between -1/2 and 1/2
	if(MaxPhi < 0.0){MaxPhi = 1.0+MaxPhi;}

	// setup variables
	Subs::Vec3 earthmax, earthact, tvec;
	double mumax, muact, maxProj;
	double bflux=0,maxflux=0;
	earthmax=Roche::set_earth(incl,MaxPhi);
	earthact=Roche::set_earth(incl,phi);
	// first element of spot is from tilted strip, so
	tvec = this->spot[0].dirn;
	maxProj = Subs::dot(tvec,earthmax);

    if(this->normalisation < 0.0){
        for(int i=0; i<this->spot.size(); i++){
            mumax = Subs::dot(this->spot[i].dirn,earthmax);
            if(i<nspot){
                // tilted strip
                if(mumax > 0.) maxflux += mumax*this->spot[i].flux;
            }else{
                // parallel strip
                if(mumax > 0.) maxflux += maxProj*this->spot[i].flux;
            }
        }
        this->normalisation = maxflux;
    }

    for(int i=0; i<this->spot.size(); i++){
        muact = Subs::dot(this->spot[i].dirn,earthact);

        if(i<nspot){
            // tilted strip
            if(muact > 0. && this->spot[i].visible(phi)) bflux += muact*this->spot[i].flux;
        }else{
            // parallel strip
            if(muact > 0. && this->spot[i].visible(phi)) bflux += maxProj*this->spot[i].flux;
        }
	}
	return bflux/this->normalisation;
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

    double xl1 = Roche::xl1(q);
    double rtest = rd*xl1;

    Subs::Vec3 r,v;
    Roche::strinit(q,r,v);
    Roche::stradv(q,r,v,rtest,1.0e-10,1.0e-3);
    this->x = r.x();
    this->y = r.y();
}

