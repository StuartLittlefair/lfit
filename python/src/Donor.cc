#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "Donor.h"
#include "trm/constants.h"
#include "geometry.h"

#include <iostream>
#include <fstream>

void LFIT::Donor::tweak(const double &q_)
{

    if (q_ == this->q)
    {
        // No need to recalculate donor star
        return;
    }
    this->q = q_;
    // we've changed the roche lobe, so empty the tiles array
    this->tiles.clear();
    // also indicate we need to recalculate normalisation
    this->normalisation = -1.;
}

void LFIT::Donor::setup_grid(const double &incl)
{

    double q = this->q;

    // Calculate a reference radius and potential
    double rref2, pref2;
    Roche::ref_sphere(q, Roche::SECONDARY, 1.0, 1.0, rref2, pref2);

    // variables for finding roche surface
    // position of surface, normal to surface and direction for search
    Subs::Vec3 posn, normvec, dirn;
    double rad, gravity;
    double ACC = 1.0e-6;

    // now find the back of 2ry
    dirn.set(1., 0., 0.);
    Roche::face(q, Roche::SECONDARY, 1.0, dirn, rref2, pref2, ACC, posn, normvec, rad, gravity);
    // Compute reference gravity value, from the side of the star opposite from the L1 point
    this->gmin = gravity;

    // OK, let's make the roche surface
    int icount = -1;
    int np, nmer;
    Subs::Vec3 xHat = Subs::Vec3(1.0, 0.0, 0.0);
    Subs::Vec3 yHat = Subs::Vec3(0.0, 1.0, 0.0);
    np = int(sqrt(this->ntiles));
    this->tiles.resize(this->ntiles);

    // create tiles
    LFIT::Point::etype eclipses;
    nmer = np;
    double dphi = Constants::TWOPI / double(np);
    double dtheta = Constants::PI / double(nmer);
    std::ofstream outfile;
    outfile.open("faces.txt");

    for (int i = 0; i < nmer; ++i)
    {
        // start at theta = 0, pointing to back of donor
        double theta = Constants::PI * (i + 0.5) / double(nmer);
        double sint = sin(theta);
        double cost = cos(theta);

        for (int j = 0; j < np; ++j)
        {
            icount++;
            // b points to surface of star at various angles
            double phi = Constants::TWOPI * double(j) / double(np);
            double sinp = sin(phi);
            double cosp = cos(phi);

            dirn.set(cost, sint * cosp, sint * sinp);
            Roche::face(q, Roche::SECONDARY, 1.0, dirn, rref2, pref2, ACC, posn, normvec, rad, gravity);

            // ingress, egress phases
            // eclipses.clear();
            // double ingress, egress;
            // if (Roche::ingress_egress(q, Roche::SECONDARY, 1.0, 1.0, incl, 1.0e-5, posn, ingress, egress)){
            //    eclipses.push_back(std::make_pair(ingress,egress));
            //}

            // we also need the element area, which is the circumference
            // of the roche lobe at this point, divided by the number
            // of theta steps, and multiplied by delta_x/cos(alpha),
            // where alpha is angle between element and x-axis
            // cirumference is twopi*t
            // area needs to be multiplied by sep**2.0 to become physical
            double area = rad * rad * sint * dphi * dtheta;
            area /= Subs::dot(dirn, normvec);

            // now set temperature of element, scaling for limb and gravity darkening
            // temp should be multiplied by pow(grav/gmin,beta)
            double temp = 3000.0 * pow(gravity / this->gmin, this->beta);

            // flux, not accounting for limb darkening
            double flux = Subs::planck(6560.0, temp);
            outfile << posn << " " << area << " " << gravity << " " << flux << " " << normvec << std::endl;
            // create point
            this->tiles[icount] = LFIT::Point(posn, normvec, area, gravity, eclipses);
            this->tiles[icount].flux = flux;
        }
    }
    outfile.close();
}

double LFIT::Donor::calcFlux(const double &phi, const double &width,
                             const double &incl)
{
    /*
    integrates over bin of finite phase width, width using trapezoidal int
    */

    double phi1 = phi - width / 2.0;
    double rflux = 0.0;
    int nphi = 5;

    for (int i = 0; i < 5; i++)
    {
        double p = phi1 + width * double(i) / double(nphi - 1);
        if (i == 0 || i == nphi - 1)
        {
            rflux += LFIT::Donor::calcFlux(p, incl) / 2.0;
        }
        else
        {
            rflux += LFIT::Donor::calcFlux(p, incl);
        }
    }
    rflux /= double(nphi - 1);

    return rflux;
}

double LFIT::Donor::calcFlux(const double &phi, const double &incl)
{
    // have we been called without the tiles calculated
    if (this->tiles.size() == 0)
    {
        std::cout << "LFIT::Calcflux shouldn't be called before calculating grid.\nThis is inefficient" << std::endl;
        this->setup_grid(incl);
    }

    if (this->normalisation < 0.0)
    {
        // maximum flux is at phi=0.75
        double maxphi = 0.75;
        Subs::Vec3 earth = Roche::set_earth(incl, maxphi);
        double sum = 0.0;
        // #pragma omp parallel for reduction(+:sum)
        for (int i = 0; i < this->tiles.size(); i++)
        {
            double mu = Subs::dot(earth, this->tiles[i].dirn);
            if (mu > 0.0 && this->tiles[i].visible(maxphi))
            {
                double flux = this->tiles[i].flux * (1. - this->ulimb + fabs(mu) * this->ulimb);
                sum += flux * this->tiles[i].area * mu;
            }
        }
        this->normalisation = sum;
    }

    double sum = 0.0;
    Subs::Vec3 earth = Roche::set_earth(incl, phi);
    // #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < this->tiles.size(); i++)
    {
        double mu = Subs::dot(earth, this->tiles[i].dirn);
        if (mu > 0.0 && this->tiles[i].visible(phi))
        {
            double flux = this->tiles[i].flux * (1. - this->ulimb + fabs(mu) * this->ulimb);
            sum += flux * this->tiles[i].area * mu;
        }
    }

    return sum / this->normalisation;
}
