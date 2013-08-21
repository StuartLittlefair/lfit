/*
 *  WhiteDwarf.h
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 09/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <cmath>
#include <string>
#include "trm_subs.h"
#include "trm_array1d.h"
#include "trm_vec3.h"
#include "trm_roche.h"

namespace LFIT{

    class WhiteDwarf{
        
    public:
        //!Default constructor
        WhiteDwarf() : radius(),ulimb(){}

        //!Constructor
        WhiteDwarf(const double& radius_,
                   const double& ulimb_) : radius(radius_),ulimb(ulimb_){}

        /*
         fit circle through three critical points to approximate
         the eclipse by the donor star
         */
        
        void tweak(const double& radius_, const double& ulimb_){
            this->radius = radius_;
            this->ulimb = ulimb_;
        }
        
        int circfit(const double& q,
                    const double& phi, const double& incl,
                    double& rc,double& xc, double& yc);
        
        double calcFlux(const double& q,
                        const double& phi, 
                        const double& incl);
        
        double calcFlux(const double& q, const double& phi, 
                        const double& width,
                        const double& incl);
        
    private:
        double radius; // radius of wd in units of XL1
        double ulimb; // limb darkening coefficient
    };
        
}
