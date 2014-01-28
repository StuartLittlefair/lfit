#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "trm/roche.h"

using namespace std;

// global variables
float rw, erw, q, p;

int main(void) {

  cout << "Enter Rw/XL1 and error: ";
  cin >> rw >> erw;
  cout << "Enter q and error: ";
  cin >> q;
  double xl1 = Roche::xl1(q);


  cout << "Rw/a = " << rw*xl1 << " +/- " << erw << endl; 
  return 0;

}

			     
