#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "Donor.h"
#include "trm/constants.h"
#include "geometry.h"

void LFIT::Donor::tweak(const double& q_){
    
    if(q_ == this->q){
        //No need to recalculate donor star
        return;
    }
    this->q = q_;
    // we've changed the roche lobe, so empty the tiles array
    this->tiles.clear();
    // also indicate we need to recalculate normalisation
    this->normalisation = -1.;
}

void LFIT::Donor::setup_grid(const double& incl){
    
    double q = this->q;
    
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
    np = int(sqrt(this->ntiles));
    this->tiles.resize(this->ntiles);

    // create tiles
    LFIT::Point::etype eclipses;
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
            Subs::Vec3 posn = a + t*b;
            // normal to tile is direction of derivative
            Subs::Vec3 dirn = Roche::drpot(q,posn);
            double gravity = dirn.length();
            dirn.unit();
            
            
            // ingress, egress phases
            //eclipses.clear();
            //double ingress, egress;
            //if (Roche::ingress_egress(q, Roche::SECONDARY, 1.0, 1.0, incl, 1.0e-5, posn, ingress, egress)){
            //    eclipses.push_back(std::make_pair(ingress,egress));
            //}            
                        
            // we also need the element area, which is the circumference
            // of the roche lobe at this point, divided by the number
            // of theta steps, and multiplied by delta_x/cos(alpha),
            // where alpha is angle between element and x-axis
            // cirumference is twopi*t
            // area needs to be multiplied by sep**2.0 to become physical
            double cos_alpha = cos( asin(Subs::dot(dirn,xHat) ) );
            double area = delta_x*Constants::TWOPI*t/double(np)/cos_alpha;

            // now set temperature of element, scaling for limb and gravity darkening
            // temp should be multiplied by pow(grav/gmin,beta)
            double temp = 3000.0*pow(gravity/this->gmin,this->beta);

            // flux, not accounting for limb darkening
            double flux = Subs::planck(6560.0,temp);

            // create point
            this->tiles[icount] = LFIT::Point(posn,dirn,area,gravity,eclipses);
            this->tiles[icount].flux = flux;
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
    
    // have we been called without the tiles calculated
    if (this->tiles.size() == 0){
        std::cout << "LFIT::Calcflux shouldn't be called before calculating grid.\nThis is inefficient" << std::endl;
        this->setup_grid(incl);   
    }

	if(this->normalisation < 0.0){
	    // maximum flux is at phi=0.75
		double maxphi = 0.75;
		Subs::Vec3 earth = Roche::set_earth(incl,maxphi);
		double sum=0.0;
	    //#pragma omp parallel for reduction(+:sum)
		for(int i=0; i< this->tiles.size();i++){
            double mu = Subs::dot(earth,this->tiles[i].dirn);
            if(mu > -0.01 && this->tiles[i].visible(maxphi)){
                double flux = this->tiles[i].flux * (1.-this->ulimb+fabs(mu)*this->ulimb);
                sum+=flux*this->tiles[i].area*mu;
            }
		}
		this->normalisation = sum;
	}
	
    double sum=0.0;
    Subs::Vec3 earth = Roche::set_earth(incl,phi);
    //#pragma omp parallel for reduction(+:sum)
    for(int i=0; i< this->tiles.size();i++){
        double mu = Subs::dot(earth,this->tiles[i].dirn);
        if(mu > -0.01 && this->tiles[i].visible(phi)){
            double flux = this->tiles[i].flux * (1.-this->ulimb+fabs(mu)*this->ulimb);
            sum+=flux*this->tiles[i].area*mu;
        }
    }

    return sum/this->normalisation;
}

   
   
    
