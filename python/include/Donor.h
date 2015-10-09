#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "Point.h"

namespace LFIT{

    class Donor{
        
    public:
        //! Default constructor
        Donor() : tiles(), beta(0.08), gmin(), ulimb(0.9),
        q(1.0e30), ntiles(400), normalisation(-1.){}
        //! Constructor
        Donor(const double& q_, const size_t& size) : tiles(),beta(0.08),
        gmin(),ulimb(0.9),q(1.0e30),ntiles(size),normalisation(-1.){
            this->tweak(q_);
        }
        void tweak(const double& q);
        double calcFlux(const double& phi, const double& incl);
        double calcFlux(const double& phi, const double& width,
                        const double& incl);
        void setup_grid(const double& incl);
        Subs::Buffer1D<double> getBright();
        double get_q(){
            return this->q;
        }
        void set_q(double q){
            this->tweak(q);
        }
        int get_size(){
            return this->tiles.size();
        }
        void set_size(int size){
            this->tiles.resize(size);
            this->ntiles = size;
            this->tweak(this->q);
        }
    private:
        Subs::Buffer1D<LFIT::Point> tiles;
        double beta; 
        double gmin; 
        double ulimb; 
        double q;
        int ntiles;
        double normalisation;
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
