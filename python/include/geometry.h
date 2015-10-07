/*
 *  geometry.h
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 09/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <cmath>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/roche.h"

namespace geometry {
    
    void project(const Subs::Vec3& vec, 
                 const Subs::Vec3& earth,
                 double& x, double& y){
        // projects a vector onto plane of sky
        Subs::Vec3 yp = Subs::Vec3(-earth.z()*earth.x(),
                                   -earth.z()*earth.y(),
                                   1.0-earth.z()*earth.z());
        double um = yp.length();
        if(um < 0.0){
            um=1.0;
            yp.set(-earth.z()*earth.x(),
                   1.0,
                   1.0-earth.z()*earth.z());
        }
        
        Subs::Vec3 xp = Subs::cross(yp,earth);
        x = Subs::dot(xp,vec);
        y = Subs::dot(yp,vec);
    }
    
    Subs::Buffer1D<double> intRaySphere(const Subs::Vec3& p,
                                        const Subs::Vec3& pHat,
                                        const Subs::Vec3& cen,
                                        const double& rad){
        //returns vectors of intersection between ray and 
        // a sphere of radius rad, centred on cen

        Subs::Buffer1D<double> results;

        // line equation is supplied by user as P = p + t*phat
        // 1, 2, or no values of t are returned depending on the number
        // of intersections found
        Subs::Vec3 store = p - cen;
        double a=pHat.sqr();
        double b=2.0*Subs::dot(pHat,store);
        double c=store.sqr() - Subs::sqr(rad);
        double test = Subs::sqr(b) - 4.0*a*c;
        
        if(test > 0.0){
            // two roots
            double t1 = (-b + sqrt(test))/2.0/a;
            double t2 = (-b - sqrt(test))/2.0/a;
            results.push_back(t1);
            results.push_back(t2);
            return results;
        }else if(test<0.0){
            // no roots
            return NULL;
        }else{
            // one root
            results.push_back(-b/2.0/a);
            return results;
        }
    }
}
