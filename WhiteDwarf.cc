/*
 *  WhiteDwarf.cc
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 09/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#include "WhiteDwarf.h"
#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "trm/constants.h"
#include "trm/format.h"

/*
 fits circle through three critical points to approximate
 eclipse by red star
 returns 0 if OK
         1 if no eclipse of WD
         2 if total eclipse of white dwarf
*/

double findphi(const double& q, const double& incl){
    
    double angle = incl*Constants::PI/180.0;
    double sini = sin(angle);
    double cosi = cos(angle);
    // white dwarf is a origin
    Subs::Vec3 p = Subs::Vec3(0.0,0.0,0.0);

    Subs::Vec3 earth = Subs::Vec3(sini,0.0,cosi);
    
    if(! Roche::blink(q,p,earth,0.01)){
        return 0.0;
    }
    // next trial phase to make sure white dwarf out of eclipse
    double trial2 = 0.3*Constants::PI, trial1=0.0;
    // Binary chop. Worst case accuracy in value of DELTA = 0.3/2**NLOOP
    // = 0.00004 for NLOOP = 13. Average case is 1/2 of this.
    //
    for(int i=1; i<=56; i++){
        double trial = (trial1+trial2)/2.0;
        double sint = sin(trial);
        double cost = cos(trial);
        //  Unit vector toward earth in the rotating frame of the binary
        //  at egress of the disc centre
        earth = Subs::Vec3(sini*cost,-sini*sint,cosi);
        if(Roche::blink(q,p,earth,0.01)){
            //white dwarf is occulted at this phase, move lower limit up
            trial1 = trial;
        }else{
            trial2=trial;
        }
    }    
    return (trial1+trial2)/2.0/Constants::PI;
}

double LFIT::WhiteDwarf::calcFlux(const double& q,
                                  const double& phi,
                                  const double& width,
                                  const double& incl){
    /* 
     computes flux of white dwarf relative to flux outside
     eclipse.
     
     integrates over bin of finite phase width, width using trapezoidal int
     */

    Subs::Format form;	
    double wflux=0.0;
    static double err=1.0e-5,qold=0.0,iold=0.0;
    double plo=0.0, pup=0.0, delta;
    if(fabs(q-qold) > err || fabs(incl-iold) > err){
      delta = findphi(q,incl);
      //std::cout << "q = " << form(q) << std::endl <<
      //	" incl =  " << form(incl) << std::endl << " delta = " << form(delta) << std::endl;
      plo = delta/2.0 - this->radius;
      pup = delta/2.0 + this->radius;
    }
    
    double ptest = fabs(phi - 1.0);

    
    if(ptest-width/2.0 > pup) {
        wflux= 1.0;
    } else if(ptest+width/2.0 < plo){
        wflux= 0.0;
    } else {
      double phi1 = phi - width/2.0;
      int nphi = 5;
      //#pragma omp parallel for reduction(+:wflux)
      for(int i=0; i<nphi; i++){
	double p = phi1 + width*double(i)/double(nphi-1);
	if(i == 0 || i == nphi-1){
	  wflux += LFIT::WhiteDwarf::calcFlux(q,p,incl)/2.0;
	}else{
	  wflux += LFIT::WhiteDwarf::calcFlux(q,p,incl);
	}
      }
      wflux /= double(nphi-1);
    }
    return wflux;
}

double LFIT::WhiteDwarf::calcFlux(const double& q,
                                  const double& phi,
                                  const double& incl){
    /* 
     computes flux of white dwarf relative to flux outside
     eclipse.
     */
    static int maxrad=100;
    double rc,xc,yc;
    
    int n_ecl = LFIT::WhiteDwarf::circfit(q,phi,incl,rc,xc,yc);
    // convert from units of separation to xl1
    double xl1 = Roche::xl1(q);
    rc /= xl1; xc /= xl1; yc /= xl1;
    
    double wflux;
    
    if(n_ecl == 1){
        // no eclipse at this phase,incl
        wflux= 1.0;
    }else if(n_ecl ==2){
        // total eclipse at this phase,incl
        wflux= 0.0;
    }else{
        // integrate over annuli on white dwarf
        // first compute projected distance between centres of annuli
        double pd = sqrt(Subs::sqr(xc)+Subs::sqr(yc));
        if( rc >= pd+this->radius ){
            wflux=0.0;
        }else if( rc <= pd-this->radius ){
            wflux=1.0;
        }else{
            double sum=0.0;
            //sum over uneclipsed region
            if(rc > pd){
                double rlo=rc-pd;
                int nrad = int(double(maxrad)*(1.0-rlo/this->radius));
                nrad = std::max(10,nrad);
                double fac = rlo*(rc+pd);
                double xlo = Subs::sqr(rlo);
                double rw2 = Subs::sqr(this->radius);
                double range = rw2-xlo;
                double dx = range/double(nrad)/rw2;
                for(int i=1; i<= nrad; i++){
                    double x = xlo + range*(double(i)-0.5)/double(nrad);
                    double theta = acos((fac-x)/2.0/pd/sqrt(x));
                    sum += theta*(1.0 - this->ulimb*(1.0-sqrt(1.0-x/rw2)));
                }
                wflux = dx*sum/Constants::PI/(1.0-this->ulimb/3.0);
            }else{
                // sum over uneclipsed region
                double rlo = pd-rc;
                int nrad = int(double(maxrad)*(1.0-rlo/this->radius));
                nrad = std::max(10,nrad);
                double fac = -rlo*(rc+pd);
                double xlo = Subs::sqr(rlo);
                double rw2 = Subs::sqr(this->radius);
                double range = rw2-xlo;
                double dx = range/double(nrad)/rw2;
                for(int i=1; i<= nrad; i++){
                    double x = xlo + range*(double(i)-0.5)/double(nrad);
                    double theta = Constants::PI-acos((fac-x)/2.0/pd/sqrt(x));
                    sum += theta*(1.0 - this->ulimb*(1.0-sqrt(1.0-x/rw2)));
                }
                wflux = 1.0-dx*sum/Constants::PI/(1.0-this->ulimb/3.0);
            }
        }
    }
    return wflux;
}



int LFIT::WhiteDwarf::circfit(const double& q,
                              const double& phi, const double& incl,
                              double& rc, double& xc, double& yc){
    
   
    static const double acc = 5.0e-6; // default accuracy
    // create some speed saving variables....
    double phase = phi*Constants::TWOPI;
    double sphi  = sin(phase);
    double cphi  = cos(phase);
    double cosi  = cos(incl*Constants::PI/180.0);
    double sini  = sin(incl*Constants::PI/180.0);
    
    // first point is where donor cuts the line of centres. x0
    // is centre of wd. Let us see if there is any eclipse at all.
    // use something slightly larger than white dwarf as donor is not
    // tangential at this point.
    Subs::Vec3 earth;
    earth=Roche::set_earth(incl,phi);
    
    double alpha = 1.2*this->radius/sqrt(1.0-Subs::sqr(earth.x()));
    Subs::Vec3 x0 = Subs::Vec3(-alpha,0.0,0.0);
    
    // is x0 eclipsed?
    if(Roche::blink(q,x0,earth,0.05)){
        // total eclipse
        return 2;
    }
    x0.set(alpha,0.0,0.0);
    if(! Roche::blink(q,x0,earth,0.05)){
        // no eclipse at all
        return 1;
    }
    
    // ok. easy cases done. from here on in it's Tom Marsh magic.
    // don't ask me.....
    double a1 = -alpha;
    double a2 =  alpha;
    //std::cout << phi << " alpha in: " << alpha << " " << q << std::endl
    //<< earth.x() << " " << earth.y() << " " << earth.z() << std::endl;
    while(1){
        double x1 = (a1+a2)/2.0;
        x0.set(x1,0.0,0.0);
        if(Roche::blink(q,x0,earth,0.05)){
            a2 = x1;
        }else{
            a1 = x1;
        }
        if( fabs(a2-a1) <= acc) break;
    }
    alpha = (a1+a2)/2.0;
    //std::cout << alpha << " " << a1 << " " << a2 << std::endl;
    
    // we have now found position of 1st point
    double x1=alpha*sphi;
    double y1=-alpha*cosi*cphi;
    /*
     2nd and third points are on circle of twice radius of white dwarf
     existence of such points apparently guaranteed by existence of 1st point
     */
    double trw = 2.0*this->radius;
    double ax1 = trw*sphi;
    double ax2 = -trw*cphi*cosi;
    double ay1 = trw*cphi;
    double ay2 = trw*sphi*cosi;
    double az  = trw*sini;
    /*
     theta 1 is angle on circle where line of centres towards donor
     crosses and is thus in eclipse
     */
    double theta1 = atan2(-cosi*cphi,sphi);

    // binary chop to find next two points
    double x2,y2,x3,y3;
	x2=0.0; y2=0.0; x3=0.0; y3=0.0;
    for(int i=1; i<=2; i++){
        double t1,t2;
        if(i==1){
            t1 = theta1;
            t2 = theta1+Constants::PI;
        }else{
            t1 = theta1;
            t2 = theta1-Constants::PI;
        }
        
        double theta,cthet,sthet;
        do{
            theta=(t1+t2)/2.0;
            cthet = cos(theta);
            sthet = sin(theta);
            x0.set(ax1*cthet+ax2*sthet,
                   ay1*cthet+ay2*sthet,
                   az*sthet);
            if(Roche::blink(q,x0,earth,0.05)){
                t1 = theta;
            }else{
                t2 = theta;
            }
        }while(trw*fabs(t2-t1) > acc);
 
        theta=(t1+t2)/2.0;
        cthet=cos(theta);
        sthet=sin(theta);
        if(i==1){
            x2=trw*cthet;
            y2=trw*sthet;
        }else{
            x3=trw*cthet;
            y3=trw*sthet;
        }
    }
    
    // now we have our three points. Fit circle to them and return
    // parameters of circle (xc,yc,rc)
    double c11,c12,c13,c21,c22,c23,c31,c32,c33,b1,b2,b3,delta,d1,d2,d3;
    c11 = 2.*(y2-y3);
    c12 = 2.*(y3-y1);
    c13 = 2.*(y1-y2);
    c21 = 2.*(x3-x2);
    c22 = 2.*(x1-x3);
    c23 = 2.*(x2-x1);
    c31 = 4.*(x2*y3-x3*y2);
    c32 = 4.*(x3*y1-x1*y3);
    c33 = 4.*(x1*y2-x2*y1);
    b1  = Subs::sqr(x1)+Subs::sqr(y1);
    b2  = Subs::sqr(x2)+Subs::sqr(y2);
    b3  = Subs::sqr(x3)+Subs::sqr(y3);
    delta = c31+c32+c33;
    d1 = (c11*b1+c12*b2+c13*b3)/delta;
    d2 = (c21*b1+c22*b2+c23*b3)/delta;
    d3 = (c31*b1+c32*b2+c33*b3)/delta;
    xc = d1;
    yc = d2;
    rc = sqrt(d3+Subs::sqr(d1)+Subs::sqr(d2));
    
    return 0;
}

