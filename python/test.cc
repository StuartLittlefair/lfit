#include <string>
#include "trm/subs.h"
#include "trm/roche.h"
#include "WhiteDwarf.h"

int main(int argc, char* argv[]){
    try{
        double q = 0.1;
        double incl = 85.9;
        double r1_a = 0.01;
        double xl1 = Roche::xl1(q);
        double r1_xl1 = r1_a/xl1;
        LFIT::WhiteDwarf wd = LFIT::WhiteDwarf(r1_xl1,0.4);
        for(size_t i=0; i<1000; ++i){
            double phi = -0.05 + double(i)*0.1/1000.0;
            double width = 0.1/1000.0;
            std::cout << phi << " " << wd.calcFlux(q,phi,width,incl) << std::endl;
        }
    }catch(const std::string& err){
    std::cerr << "\nError occured inside lfit_fake:" << std::endl;
    std::cerr << err << std::endl;
   }

  return 0;   

}