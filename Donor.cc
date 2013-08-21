#include <cmath>
#include <string>
#include "trm_subs.h"
#include "trm_array1d.h"
#include "trm_vec3.h"
#include "trm_roche.h"
#include "Donor.h"
#include "trm_constants.h"
#include "geometry.h"

std::istream& LFIT::operator>>(std::istream& s, const DonorEl& el){
    std::string err = "Calling LFIT::operator >> with DonorEl is an error";
    throw err;
}
std::ostream& LFIT::operator<<(std::ostream& os, const DonorEl& el){
    std::string err = "Calling LFIT::operator << with DonorEl is an error";
    throw err;
}

//! Constructor
void LFIT::Donor::tweak(const double& q){
    
    if(q == this->q_old){
        //No need to recalculate donor star
        return;
    }
    this->q_old = q;
    // we've changed the roche lobe, so make sure next call to setVis does something
    this->phi_old = -1.0e30;
    
    // find inner lagrangian pt (distance from 1ry in units of separation)
    double xl1 = Roche::xl1(q);
    // set vector a to XL1 point (xl1,0,0)
    Subs::Vec3 a = Subs::Vec3(xl1,0.0,0.0);
    // calculate critical potential at xl1 point
    double cpot = Roche::rpot(q,a);
    
    // now find the back of 2ry
    a.set(1.0,0.0,0.0);
    Subs::Vec3 b(a);
    // set reasonable limits for length along vec b
    double t1 = 0.05*(1.0-xl1);
    double t2 = 1.10*(1.0-xl1);
    // solve to find critical potential along b....
    double t = Subs::rtsafe(Lfunc(a,b,q,cpot),t1,t2,1.0e-6);
    double xback = 1.0+t; // x-coord of back of donor
    
    //also find pole
    b.set(0.0,0.0,1.0);
    t2 = sqrt(Subs::sqr(1.01*(1.0-xl1))-Subs::sqr(a.x()-1.0));
    t1 = t2/1.0e2;
    // solve to find critical potential along b....
    t = Subs::rtsafe(Lfunc(a,b,q,cpot),t1,t2,1.0e-6);
    Subs::Vec3 pole(a+t*b);
    // calculate derivative of roche potential at pole
    Subs::Vec3 grad = Roche::drpot(q,pole);
    this->gmin = grad.length();
    
    // OK, let's make the roche surface
    int icount=-1;
    int np, nmer;
    Subs::Vec3 xHat = Subs::Vec3(1.0,0.0,0.0);
    np = int(sqrt(this->tiles.size()));
    nmer=np;
    double delta_x = (xback-xl1)*double(nmer-1)/(double(nmer))/double(nmer+1);
    for(int i=0; i< nmer; ++i){
        a.set(xl1 + (xback-xl1)*double(i+1)/double(nmer+1),0.0,0.0);
        for(int j=0; j<np; ++j){
            icount++;
            // b points to surface of star at various angles
            double theta = Constants::TWOPI*double(j)/double(np);
            b.set(0.0,cos(theta),sin(theta));
            // limit search along b to sensible values
            t2 = sqrt(Subs::sqr(1.01*(1.0-xl1))-Subs::sqr(a.x()-1.0));
            t1 = t2/1.0e2;
            // solve to find critical potential along b....
            t = Subs::rtsafe(Lfunc(a,b,q,cpot),t1,t2,1.0e-6);
            
            // so this is surface of star
            this->tiles[icount].pos = a + t*b;
            grad = Roche::drpot(q,this->tiles[icount].pos);
            
            // normal to tile is direction of derivative
            this->tiles[icount].norm = grad;
            this->tiles[icount].norm.unit();
            
            // we also need the element area, which is the cicumference
            // of the roche lobe at this point, divided by the number
            // of theta steps, and multiplied by delta_x/cos(alpha),
            // where alpha is angle between element and x-axis
            // cirumference is twopi*t
            // area needs to be multiplied by sep**2.0 to become physical
            double cos_alpha = cos( asin(Subs::dot(this->tiles[icount].norm,xHat) ) );
            this->tiles[icount].area = delta_x*Constants::TWOPI*t/double(np)/cos_alpha;

            // now set temperature of element, scaling for limb and gravity darkening
            this->tiles[icount].gravity = b.length();
            // temp should be multiplied by pow(grav/gmin,beta)
            
            // visibility set later
            this->tiles[icount].vis    =0.0;
            this->tiles[icount].bright =0.0;
        }
    }
}//matches constructor

void LFIT::Donor::setVis(const double& phi, const double& incl){
    
    // don't re-calculate unless absolutely necessary
    if(this->phi_old == phi && this->incl_old == incl){
        return;
    }
    this->phi_old = phi;
    this->incl_old = incl;
    
    Subs::Vec3 earth;
    earth = Roche::set_earth(incl,phi);
	//#pragma omp parallel for
    for(int i=0; i<this->tiles.size();i++){
        double dotP = Subs::dot(earth,this->tiles[i].norm);
        if(dotP >= -0.01){
            //element is visible. set vis to projected area
            this->tiles[i].vis=dotP;
            //set element intensity in units of W/m**2/Hz/sr
            //account for limb and gravity darkening
            // nominal values of temp and wavelenght are
            // used here. The donor's light will later be normalised
            // to 1.
            double temp = 3000.0*pow(this->tiles[i].gravity/this->gmin,this->beta);
            this->tiles[i].bright = (1.0-this->ulimb+this->ulimb*fabs(dotP)) *
            Subs::planck(6560.0,temp);
            
        }else{
            //element is invisible. 
            this->tiles[i].vis=0.0;
            //set element intensity in units of W/m**2/Hz/sr
            //account for limb and gravity darkening
            // nominal values of temp and wavelenght are
            // used here. The donor's light will later be normalised
            // to 1.
            double temp = 3000.0*pow(this->tiles[i].gravity/this->gmin,this->beta);
            tiles[i].bright = (1.0-this->ulimb+this->ulimb*fabs(dotP)) *
            Subs::planck(6560.0,temp);            
        }
    }
}

double LFIT::Donor::calcFlux(const double& phi, const double& width,
                             const double& incl){
    /*
    integrates over bin of finite phase width, width using trapezoidal int
    */
    
    double phi1 = phi - width/2.0;
    double rflux=0.0;
    int nphi = 5;
	
    for(int i=0; i<5; i++){
        double p = phi1 + width*double(i)/double(nphi-1);
        if(i == 0 || i == nphi-1){
            rflux += LFIT::Donor::calcFlux(p,incl)/2.0;
        }else{
            rflux += LFIT::Donor::calcFlux(p,incl);
        }
    }
    rflux /= double(nphi-1);
    
    return rflux;
}

double LFIT::Donor::calcFlux(const double& phi, const double& incl){
    
	double static maxflux;
	bool static first=true;

	if(first){
	    // maximum flux is at phi=0.75
		this->setVis(0.75,incl);
		double sum=0.0;
	    //#pragma omp parallel for reduction(+:sum)
		for(int i=0; i< this->tiles.size();i++){
			sum+=this->tiles[i].bright*this->tiles[i].area*this->tiles[i].vis;
		}
		first = false;
		maxflux = sum;
	}
	
    this->setVis(phi,incl);
    double sum=0.0;
	//#pragma omp parallel for reduction(+:sum)
    for(int i=0; i< this->tiles.size();i++){
        sum+=this->tiles[i].bright*this->tiles[i].area*this->tiles[i].vis;
    }

    return sum/maxflux;
}

// works in units of wd to XL1
bool LFIT::blink3(const double& q, const Subs::Vec3& p, const Subs::Vec3& phat){
    
    static double qq=-100.0;
    double xl1, crit, step, xcm, c1, c2, rsphere, pp;
    Subs::Vec3 a1;
    
    //if(qq != q){
    //first time in
    xcm = 1.0/(1.0+q);
    c1 = 2.0*xcm;
    xcm = q*xcm;
    c2 = 2.0*xcm;
    
    xl1 = Roche::xl1(q);
    a1.set(xl1,0.0,0.0);
    crit = Roche::rpot(q,a1);
    step = 0.05*(1.0-xl1);
    rsphere = 1.0-xl1;
    pp = Subs::sqr(rsphere);
    qq = q;
    //}
    double dx1 = phat.x();
    double dx2 = phat.y();
    double dx3 = phat.z();
    double z1  = p.x()*xl1;
    double z2  = p.y()*xl1;
    double z3  = p.z()*xl1;
    double b   = dx1*(z1-1.0) + dx2*z2 + dx3*z3;
    double c   = Subs::sqr(z1-1.0) + Subs::sqr(z2) + Subs::sqr(z3);
    
    double fac = b*b - c + pp;
    if(fac <= 0.0) return false;
    fac = sqrt(fac);
    double par1 = -b - fac;
    double par2 = -b + fac;
    
    if(par2 <= 0.0) return false;
    
    par1 = std::max(0.0,par1);
    double par = (par1+par2)/2.0;
    double x1  = z1 + par*dx1;
    double x2  = z2 + par*dx2;
    double x3  = z3 + par*dx3;
    
    double yy  = x2*x2;
    double rr  = yy + x3*x3;
    double rs1 = x1*x1 + rr;
    double r1  = sqrt(rs1);
    double xm  = x1-1.0;
    double rs2 = xm*xm + rr;
    
    if(rs2 <= 0.0) return true;
    double r2 = sqrt(rs2);
    double xc = x1 - xcm;
            c = c1/r1 + c2/r2 + xc*xc + yy;
    
    if (c > crit) return true;
    
    double a = x1*dx1 + x2*dx2;
    b = a + x3*dx3;
    double fderiv = -c1*b/rs1/r1-c2*(b-dx1)/rs2/r2+2.*(a-xcm*dx1);
    double p1,p2;
    if(fderiv > 0.0){
        p1 = par;
        p2 = par2;
    }else{
        p1 = par;
        p2 = par1;
    }
    int nstep = std::max(2,int(fabs(p2-p1)/step));
    double dp = (p2-p1)/double(nstep);
    for(int i=0; i< nstep; i++){
        double cmax=c;
        double pt = p1 + dp*double(i+1);
        x1 = z1 + pt*dx1;
        x2 = z2 + pt*dx2;
        x3 = z3 + pt*dx3;
        yy = x2*x2;
        rr = yy+x3*x3;
        r1 = sqrt(x1*x1 + rr);
        xm = x1-1.0;
        r2 = sqrt(xm*xm + rr);
        if(r2 <= 0.0) return true;
        xc = x1-xcm;
        c  = c1/r1 + c2/r2 + xc*xc + yy;
        if(c > crit) return true;
        if(c < cmax) return false;
    }
    return false;
}

// works in units of separation
bool LFIT::blink2(const double& q, const Subs::Vec3& p, const Subs::Vec3& phat){

    static double qq=-100.0, xl1, crit, step;
    Subs::Vec3 a,b;
    
    if(qq != q){
        //first time in
        xl1 = Roche::xl1(q);
        a.set(xl1,0.0,0.0);
        crit = Roche::rpot(q,a);
        step = 0.05*(1.0-xl1);
        qq = q;
    }
    
    // from here on in, everything is calculated every time
    
    // the red star lies entirely within a sphere centred on it's centre of
    // mass and with radius xl1
    a.set(1.0,0.0,0.0); // centre of sphere
    Subs::Buffer1D<double> tint = geometry::intRaySphere(p, phat, a, 1.0-xl1); 
    
    // ok, if there is no intersection with this sphere, why go further
    if(tint.size() != 2){
        return false;
    }
    
    // or perhaps there are intersections, but they are past sphere, or
    // inside sphere
    if(tint[0] < 0.0 || tint[1] < 0.0){
        return false;
    }
    
    // ok, now we need to test carefully. only check positive path...
    double t1=std::max(tint[0],0.0);
    double t2=std::max(tint[1],0.0);
    
    // start halfway between intersections
    double t = 0.5*(t1+t2);
    a = p + t*phat;
    
    // are we at centre of star?
    Subs::Vec3 xhat = Subs::Vec3(1.0,0.0,0.0);
    if( (a-xhat).length() == 0.0 ){
        return true;
    }
    
    // are we inside roche potential
    double pot = Roche::rpot(q,a);
    b = Roche::drpot(q,a);
    if(pot <= crit){
        return true;
    }
    
    // no? then we step in dir of decreasing roche potential
    // we need vector of derivative
    double dotP = Subs::dot(phat,b);
    
    // if dotP is +ve then >t implies >potential and vice versa
    double start,stop;
    if(dotP<-1.0e-30){
        start=std::min(t1,t2);
        stop =std::max(t1,t2);
    }else{
        start=std::max(t1,t2);
        stop =std::min(t1,t2);
    }
    
    // so we loop whilst roche potential is increasing
    int nsteps = std::max(2,int(fabs(start-stop)/step));
    double dt = (start-stop)/double(nsteps);

    for(int i=0; i< nsteps; ++i){
        double potmax = pot;
        t=start + dt*double(i);
        a = p + t*phat;
        
        //centre of donor?
        if( (a-xhat).length() == 0.0 ){
            return true;
        }    
        
        // test roche potential
        pot = Roche::rpot(q,a);
        if(pot > potmax){            
            // no longer increasing
            return false;
        }
        if(pot < crit){            
            return true;
        }
            
    }    
    return false;
}
   
   
    
