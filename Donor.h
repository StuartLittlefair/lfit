#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"

namespace LFIT{


    bool blink2(const double& q, const Subs::Vec3& p, const Subs::Vec3& phat);
    
    bool blink3(const double& q, const Subs::Vec3& p, const Subs::Vec3& phat);

    //! Structure defining a single element
    /** This defines the position, area, direction, gravity and brightness of an element
     * and also any phases during which it is eclipsed.
     */
    struct DonorEl {

        //! Default constructor
        DonorEl() : pos(), norm(), area(0.), gravity(1.), bright(0.), vis(0.) {}

        //! Constructor
        DonorEl(const Subs::Vec3& pos_, const Subs::Vec3& norm_, double area_, 
              double gravity_, double bright_, double vis_) : 
        pos(pos_), norm(norm_), area(area_), gravity(gravity_), 
        bright(bright_), vis(vis_) {}

        //! Position vector of element (units of binary separation)
        Subs::Vec3 pos;

        //! Outward facing direction of element (unit vector)
        Subs::Vec3 norm;

        //! Area of element (units of binary separation**2)
        double area;

        //! Gravity of element 
        double gravity;

        //! Brightness
        double bright;

        //! Visible? if so, projected area
        double vis;
    };
    
    // Dummy ASCII input operator to allow DonorEl's to be put in Buffer1D's
    std::istream& operator>>(std::istream& s, const DonorEl& p);
    
    
    // Dummy ASCII output operator to allow use of Buffer1D
    std::ostream& operator<<(std::ostream& s, const DonorEl& p);
    
    
    class Donor{
        
    public:
        //! Default constructor
        Donor() : tiles(), beta(0.08), gmin(), ulimb(0.9),
        q_old(1.0e30),phi_old(1.0e30),incl_old(1.0e30){}
        //! Constructor
        Donor(const double& q, const size_t& size) : tiles(),beta(0.08),
        gmin(),ulimb(0.9),q_old(1.0e30),phi_old(1.0e30),incl_old(1.0e30){
            tiles.resize(size);
            this->tweak(q);
        }
        void tweak(const double& q);
        void setVis(const double& phi, const double& incl);
        double calcFlux(const double& phi, const double& incl);
        double calcFlux(const double& phi, const double& width,
                        const double& incl);
        Subs::Buffer1D<double> getBright();
    private:
        Subs::Buffer1D<DonorEl> tiles;
        double beta; 
        double gmin; 
        double ulimb; 
        double q_old;
        double phi_old;
        double incl_old;
    };

    //! Function object to compute Roche potential and its derivative
    /**
     * This is needed for 'rtsafe' 
     */
    class Lfunc : public Subs::RTfunc {
        
        double qp, cpot;
        Subs::Vec3 x0, x1;
        
    public:
        
        //! Constructor storing fixed data
        Lfunc(Subs::Vec3 x0i, Subs::Vec3 x1i, double qpi, double cpoti) : 
        qp(qpi), cpot(cpoti), x0(x0i), x1(x1i) {}
        
        //! Function operator
        void operator()(double t, double& f, double& d) const { 
            Subs::Vec3 p, dp;
            p = x0 + t*x1;
            f   = Roche::rpot(qp, p)-cpot;
            dp  = Roche::drpot(qp, p);
            d   = Subs::dot(dp,x1);
        }
        
    };
    
}
