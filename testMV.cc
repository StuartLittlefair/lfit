#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/array2d.h"
#include "trm/input.h"
#include "trm/plot.h"
#include "trm/format.h"
#include "cholesky.h"
#include <cfloat>
#include <climits>
#include <ctime>
#include <iostream>
#include <fstream>
#include <stdlib.h> 

using namespace std;

template <class X>
Subs::Array1D<X> matMul(const Subs::Array2D<X>& A, const Subs::Array1D<X>& b){

  int nx = A.get_nx();
  int ny = A.get_ny();
  if( b.size() != ny){
    throw Subs::Array2D_Error("Cannot multiply matrix of dimensions [" + Subs::str(nx) + "][" + Subs::str(ny) + "], with vector of size " + Subs::str(b.size()));
  }
  Subs::Array1D<X> prod(ny);
  prod = 0.0;
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++) {
      prod[i] += A[j][i]*b[j];
    }
  } 
  return prod;
}

int main(int argc, char *argv[]){

  int npar = 13;
  Subs::Array2D<double> Sigma(npar,npar), A(npar,npar);
  Subs::Array1D<double> z(npar), samp(npar), mu(npar);
  Subs::Array1D<float> x,y;

  std::ifstream fin;
  fin.open("means");
  for(int i=0; i<npar; i++){
    fin >> mu[i];
  }
  fin.close();

  fin.open("covar");
  for(int j=0; j<npar; j++){
    for(int i=0; i<npar; i++) fin >> Sigma[i][j];
  }


  A = choleskyDecomp(Sigma);
  /*
  Subs::Format form(7);
  std::cout << "Covariance Matrix: " << std::endl;
  for(int i=0; i<npar; i++){
    for(int j=0;j<npar;j++){
      std::cout << std::setw(15) <<  form(Sigma[j][i]) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Decomp: " << std::endl;
  for(int i=0; i<npar; i++){
    for(int j=0;j<npar;j++){
      std::cout << std::setw(15) <<  form(A[j][i]) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  */
  if (argc != 3) {
    std::string err = "Command arguments needed: two parameters to plot against each other";
    throw err;
  } 
  int p1=atoi(argv[1])-1, p2=atoi(argv[2])-1;

  ofstream fout("xy.dat");
  int seed = (int) time(NULL);
  for(int count=0; count < 10000; count++){
    for(int i=0;i<npar;i++) z[i] = Subs::gauss1(seed);
    samp = mu+matMul(A,z);
    x.push_back(samp[p1]);
    y.push_back(samp[p2]);
    fout << Subs::str(samp[p1]) << " " << Subs::str(samp[p2]) << endl;
    samp.clear();
  }
  fout.close();
  
  Subs::Plot plot;
  plot.open("?");
  cpgenv(x.min(),x.max(),y.min(),y.max(),0,0);
  cpgpt(x.size(),x.ptr(),y.ptr(),1);
  plot.close();
    
}
