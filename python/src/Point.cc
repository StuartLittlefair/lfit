#include "Point.h"
#include <string>

std::istream& LFIT::operator>>(std::istream& s, const Point& el){
    std::string err = "Calling LFIT::operator >> with Point is an error";
    throw err;
}
std::ostream& LFIT::operator<<(std::ostream& os, const Point& el){
    std::string err = "Calling LFIT::operator << with Point is an error";
    throw err;
}