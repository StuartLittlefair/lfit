/*
 *  BrightSpot.h
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 09/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "lfit_params.h"
#include "Point.h"

namespace LFIT{
    
    class BrightSpot{
    public:
        //! Default Constructor
        BrightSpot() : complex(false),q(1.0e30),rd(1.0e30),az(),frac(),scale(),
			exp1(),exp2(),tilt(),yaw(),x(),y(),nspot(300),spot(),normalisation(-1.){
			}
        
        //! Constructor from all parameters needed, including
        // mass ratio, q and disc radius in units of XL1
        BrightSpot(const LFIT::Params& params) : 
		complex(params.complexSpot),q(params.q),rd(params.rd),az(params.bsAz),frac(params.bsFrac),
		scale(params.bsScale),nspot(300),spot(),normalisation(-1.){
			if (this->complex){
				this->exp1 = params.bsExp1;
				this->exp2 = params.bsExp2;
				this->tilt = params.bsTilt;
				this->yaw  = params.bsYaw;
			}
            this->spotPos(params.q,params.rd);
        }
        
        //! Constructor from all parameters needed, including
        // mass ratio, q and disc radius in units of XL1
        // allows custom sizing
        BrightSpot(const LFIT::Params& params, const size_t& nspot_) : 
		complex(params.complexSpot),q(params.q),rd(params.rd),az(params.bsAz),frac(params.bsFrac),
		scale(params.bsScale),nspot(nspot_),spot(),normalisation(-1.){
			if (this->complex){
				this->exp1 = params.bsExp1;
				this->exp2 = params.bsExp2;
				this->tilt = params.bsTilt;
				this->yaw  = params.bsYaw;
			}
            this->spotPos(params.q,params.rd);
        }        
        
        void tweak(const LFIT::Params& params){
            this->az = params.bsAz;
            this->frac = params.bsFrac;
            this->scale = params.bsScale;
			if(this->complex){
				this->exp1 = params.bsExp1;
				this->exp2 = params.bsExp2;
				this->tilt = params.bsTilt;
				this->yaw  = params.bsYaw;				
			}
			if (this->q != params.q or this->rd != params.rd){
				this->spotPos(params.q,params.rd);
				this->q = params.q; this->rd = params.rd;
			}
			// empty spot elements so a recompute is needed
			this->spot.clear();
			// indicate that normalisation needs recalculating
			this->normalisation = -1.;
        } 
        
        void spotPos(const double& q, const double& rd);
        
        double   calcFlux(const double& q, const double& phi, const double& incl);
		double simpleFlux(const double& q, const double& phi, const double& incl);
		void   setup_grid(const double& incl);
        double   calcFlux(const double& q, const double& phi, const double& width,
                        const double& incl);
	    double getTangent() const;
    private:
		bool  complex; //use the complex bright spot model or not
		double q; // last computed q 
		double rd; // last computed disc radius 
		double az; // azimuth of spot
		double frac; // fraction of light in strip which sits in plane of disc
		double scale; // scale of spot
		double exp1; // exponent governing brightness function of strip
		double exp2; // exponent governing brightness function of strip
		double tilt; // tilt of non-planar strip w.r.t. disc plane (90.0 to start)
		double yaw; // yaw of non-planar strip around planar strip
		double x; // x position of spot (z is zero cos spot is in plane)
		double y; // y position of spot
		int    nspot; // number of positions along spot
		Subs::Buffer1D<LFIT::Point> spot; // buffer of point objects
		double normalisation; //maximum flux for normalising curve
    };
}
