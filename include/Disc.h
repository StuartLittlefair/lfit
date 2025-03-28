/*
 *  Disc.h
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
#include "Point.h"

namespace LFIT{
    
    class Disc{
    public:
        //!Default Constructor
        Disc() : tiles(),ntiles(1000),q(),rin(),rout(),exp(),normalisation(-1.0){}
        //!Constructor
        // rdisc is in units of xl1
        Disc(const double& q_, const double& rwd_, 
             const double& rd_, const double& exp_,
             const size_t& size) : tiles(), ntiles(size),normalisation(-1.0){
            this->tweak(q_,rwd_,rd_,exp_);
        }
		
        bool paramsChanged(const double& q, const double& rwd, const double& rd, const double& exp_);
        void tweak(const double& q, const double& rwd, const double& rd, const double& exp_);
        void setup_grid(const double& incl);
        double calcFlux(const double& q, const double& phi, const double& incl);
        double calcFlux(const double& q, const double& phi, const double& width,
                        const double& incl);

    private:
        Subs::Buffer1D<LFIT::Point> tiles;
        int ntiles;
        double q;
        double rin;
        double rout;
        double exp; // flux scales as (r/rdisc)**-exp
        double normalisation;
    };
    
    
    
}
