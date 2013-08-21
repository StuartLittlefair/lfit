#include <iostream>
#include <string>
#include <readline/readline.h>
#include <readline/history.h>
#include "trm_subs.h"
#include "lfit.h"

double get_delta_phi(LFIT::Disc disc, double q, double incl){

	int np = 5000;
	double range=0.2;
	Subs::Array1D<double> phase(np);
	Subs::Array1D<double> flux(np);
	// first half
	#pragma omp parallel for
	for(int i=0; i<np; i++){
		phase[i] = -range + float(i)*range/float(np);
		flux[i]  = disc.calcFlux(q,phase[i],range/float(np),incl);
	}
	flux /= flux.max();
	int index = flux.locate(0.5);
	double step = (0.5-flux[index-1]) / (flux[index]-flux[index-1]);
	double dPhi = -2.0 * ( phase[index-1] + step*(phase[index]-phase[index-1]) );
	return dPhi;
}

int main( int argc, char* argv[]) {

	try{
		if(argc != 3){
			std::cout << "Need parameters" << std::endl;
			return 1;
		}
		
		double q = 2.35;
		double xl1 = Roche::xl1(q);
		double rDisc = Subs::string_to_double(argv[1])*xl1;
		double discExp = Subs::string_to_double(argv[2]);
		double incl = 69.2;
		
		LFIT::Disc disc = LFIT::Disc(q,0.01,rDisc,discExp,2500);

		double dPhi = get_delta_phi(disc,q,incl);
		
		Subs::Format form;
		form.precision(3); form.general();
		std::cout << dPhi << std::endl;
		return 0;
	}catch(const std::string& err){
		std::cerr << err << std::endl;
	}
}