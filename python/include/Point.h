#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"

#ifndef __POINT_H_INCLUDED__   // if x.h hasn't been included yet...
#define __POINT_H_INCLUDED__   //   #define this so the compiler knows it has been included

namespace LFIT {
    struct Point {

        typedef std::vector<std::pair<double,double> > etype;

        //! Default constructor
        Point() : posn(), dirn(), area(0.), gravity(1.), eclipse(), flux(0.) {}

        //! Constructor
        Point(const Subs::Vec3& posn_, const Subs::Vec3& dirn_, double area_, const etype& eclipse) :
            posn(posn_), dirn(dirn_), area(area_), gravity(1.0), eclipse(eclipse), flux(0.) {}

        //! Constructor with gravity
        Point(const Subs::Vec3& posn_, const Subs::Vec3& dirn_, double area_, double gravity_, 
            const etype& eclipse) :
            posn(posn_), dirn(dirn_), area(area_), gravity(gravity_), eclipse(eclipse), flux(0.) {}
            
        //! Position vector of element (units of binary separation)
        Subs::Vec3 posn;

        //! Outward facing direction of element (unit vector)
        Subs::Vec3 dirn;

        //! Area of element (units of binary separation**2)
        float area;

        //! Gravity of element
        float gravity;

        //! Ingress and egress phases of eclipses, if any
        etype eclipse;

        //! Brightness * area
        float flux;

        //! Computes whether a point is visible (else eclipsed)
        bool visible(double phase) const {
            double phi = phase - floor(phase);
            for(size_t i=0; i<eclipse.size(); i++){
                const std::pair<double,double> p = eclipse[i];
                if((phi >= p.first && phi <= p.second) || phi <= p.second-1.0)
                    return false;
            }
            return true;
        }
    };
    // Dummy ASCII input operator to allow Points to be put in Buffer1D's
    std::istream& operator>>(std::istream& s, const Point& p);
    
    
    // Dummy ASCII output operator to allow use of Buffer1D
    std::ostream& operator<<(std::ostream& s, const Point& p);
}

#endif