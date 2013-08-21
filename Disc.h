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
#include "trm_subs.h"
#include "trm_array1d.h"
#include "trm_vec3.h"
#include "trm_roche.h"

namespace LFIT{
    
    struct DiscEl{
        //!Default Constructor
        DiscEl() : pos(),norm(),area(0.0),bright(0.0),vis(0.0){}
        
        //! Constructor
        DiscEl(const Subs::Vec3& pos_, const Subs::Vec3& norm_,
               const double& area_,const double& bright_,const double& vis_) : 
        pos(pos_),norm(norm_),area(area_),bright(bright_),vis(vis_){}
        
		// destructor
		~DiscEl(){}
		
        Subs::Vec3 pos, norm; // position vector of element (units of sep)
        double area,bright,vis;
    };
    // Dummy ASCII input operator to allow DonorEl's to be put in Buffer1D's
    std::istream& operator>>(std::istream& s, const DiscEl& p);
    // Dummy ASCII output operator to allow use of Buffer1D
    std::ostream& operator<<(std::ostream& s, const DiscEl& p);
    
    class Disc{
    public:
        //!Default Constructor
        Disc() : tiles(),rin(),rout(),exp(){}
        //!Constructor
        // rdisc is in units of xl1
        Disc(const double& q, const double& rwd_, 
             const double& rd_, const double& exp_,
             const size_t& size) : tiles(){
            tiles.resize(size);
            this->tweak(q,rwd_,rd_,exp_);
        }

		
        bool paramsChanged(const double& q, const double& rwd, const double& rd, const double& exp_);
        void tweak(const double& q, const double& rwd, const double& rd, const double& exp_);
        void setVis(const double& q,const double& phi, const double& incl);
        double calcFlux(const double& q, const double& phi, const double& incl);
        double calcFlux(const double& q, const double& phi, const double& width,
                        const double& incl);
    private:
        Subs::Buffer1D<DiscEl> tiles;
        double rin;
        double rout;
        double exp; // flux scales as (r/rdisc)**-exp
    };
    
    
    
}
