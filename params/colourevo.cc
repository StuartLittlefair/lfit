#include <iostream>
#include <iomanip>
#include <fstream>
#include "trm/subs.h"
#include "trm/buffer2d.h"
#include "trm/array1d.h"
#include "trm/plot.h"
#include "nr.h"

using namespace std;

int minloc(const Subs::Array1D<float>& a){
  int i;
  for (i=0; i<a.size(); i++){
    if (a[i]==a.min()) break;
  }
  return i;
}

int maxloc(const Subs::Array1D<float>& a){
  int i;
  for (i=0; i<a.size(); i++){
    if (a[i]==a.max()) break;
  }
  return i;
}

int nearest(const Subs::Buffer1D<float>& val, float aim){
  Subs::Array1D<float> diff(val.size());
  for(int j=0; j<val.size(); j++){
    diff[j] = abs(val[j]-aim);
  }
  return minloc(diff);
}

template <class X> X get_R (const Subs::Array1D<X>& x,
             const Subs::Array1D<X>& y, X goal) {

  X tmp;
  // check x array is ordered
  if(! x.monotonic() ){
    Subs::Array1D<X> store(x);
    Subs::Buffer1D<unsigned long int> key = store.sort();
    int minpos = Subs::locate(store.ptr(),store.size(),goal);
    // linearly interpolate to find y at that x
    tmp = Subs::linterp(x[key[minpos-1]],y[key[minpos-1]],
                              x[key[minpos]],y[key[minpos]],goal);
  }else{
    int minpos = Subs::locate(x.ptr(),x.size(),goal);
    // linearly interpolate to find y at that x
    tmp = Subs::linterp(x[minpos-1],y[minpos-1],
                        x[minpos],y[minpos],goal);
  }
  return tmp;
}


void spline(const Subs::Array1D<float>&  x, const Subs::Array1D<float>&  y,
	    unsigned long size, Subs::Array1D<float>& xint, 
	    Subs::Array1D<float>& yint){
  if(x.size() != y.size()) 
    throw "Error returned inside spline: x and y arrays differ in size\n";
  Vec_DP vecx(x.size()), vecy(x.size()), d2y(x.size());
  for(int i=0; i< x.size(); i++){vecx[i] = x[i]; vecy[i]=y[i];}
  DP yp1, ypn;
  yp1 = 1.0e30; ypn = 1.0e30;
  NR::spline(vecx,vecy,yp1,ypn,d2y);
  for(int i=0; i<size; i++){
    DP xtmp, ytmp;
    xtmp = x.min() + (x.max()-x.min())*float(i)/float(size);
    NR::splint(vecx,vecy,d2y,xtmp,ytmp);
    xint.push_back(xtmp);
    yint.push_back(ytmp);
  }      
}

float will_rad(const float& mass){
  
  float cons = 7.795e6;
  float rSol = 6.9550000e8;
  float mch = 1.44;
  float expo = 2.0/3.0;
  float tmp1 = pow(mch/mass,expo) - pow(mass/mch,expo);
  float val = cons*sqrt(tmp1)/rSol;
  return val;
}

float eggleton_rad(const float& mass){

  const float mp=0.00057;
  const float mch=1.44;
  const float a = -2.0/3.0;
  const float b = 2.0/3.0;
  const float c = 0.5;
  const float d = -1.0;
  float tmp1 = pow(mass/mch,a) -  pow(mass/mch,b);
  float tmp2 = 1.0 + 3.5*pow(mass/mp,a) + pow(mass/mp,d);
  float val = 0.0114*( pow(tmp1,c) * pow(tmp2,a) );
  return val;
}

float nau_rad(const float& mass){

  float m3=5.816;
  float mu=2.0;
  m3 = m3/Subs::sqr(mu);
  const float a = 4.0/3.0;
  const float b = 1.0/3.0;
  const float c = 0.5;
  float tmp1 = 1.0-pow(mass/m3,a);
  float tmp2 = pow(tmp1,c);
  float tmp3 = pow(mass/m3,b);
  
  float val = 0.0225*tmp2/tmp3/mu;
  return val;
}

void read_benvenuto(const string& file, Subs::Buffer1D<float>& m, Subs::Buffer1D<float>& r){
  
  ifstream fin;
  fin.open(file.c_str());
  if(!fin) {
    string error("Failed to open file: ");
    error.append(file);
    throw error;
  }
  while(fin){
    float skip,x, y;
    fin >> skip >> x >> y;
    m.push_back(x);
    r.push_back(y/100.0);
  }
}

void read_wood(const string& file, Subs::Buffer1D<float>& teff, 
	    Subs::Buffer1D<float>& radius){
  
  ifstream fin;
  fin.open(file.c_str());
  if(!fin) {
    string error("Failed to open file: ");
    error.append(file);
    throw error;
  }
  float skip;
  float rSol = 6.955e10;
  float ten = 10.0;
  while(fin){
    float x, y;
    fin >> skip >> skip >> skip >> skip >> skip >> x >> y >> skip
	>> skip >> skip;
    teff.push_back(pow(ten,y));
    radius.push_back(pow(ten,x)/rSol);
  }
}


void read_colours(const string& file, Subs::Buffer1D<float>& teff,
		  Subs::Buffer1D<float>& logg,
		  Subs::Buffer1D<float>& mass,
		  Subs::Buffer1D<float>& mbol,
		  Subs::Buffer1D<float>& u,
		  Subs::Buffer1D<float>& g,
		  Subs::Buffer1D<float>& r,
		  Subs::Buffer1D<float>& i,
		  Subs::Buffer1D<float>& z,
		  Subs::Buffer1D<float>& age){
  
  ifstream fin;
  fin.open(file.c_str());
  if(!fin) {
    string error("Failed to open file: ");
    error.append(file);
    throw error;
  }
  float skip;
  // skip header
  fin.ignore(100000, '\n');
  fin.ignore(100000, '\n');
  while(fin){
    float x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
    fin >> x1 >>  x2 >> x3 >> x4 >> x5 >> x6 >> x7 >>x8
	>> x9 >> x10;
    teff.push_back(x1); logg.push_back(x2); mass.push_back(x3);
    mbol.push_back(x4);u.push_back(x5);g.push_back(x6);
    r.push_back(x7);i.push_back(x8);z.push_back(x9);age.push_back(x10);
  }
}


void hamada(const string& file, Subs::Buffer1D<float>& mass, 
	    Subs::Buffer1D<float>& radius){
  
  ifstream fin;
  fin.open(file.c_str());
  if(!fin) {
    string error("Failed to open file: ");
    error.append(file);
    throw error;
  }
  float skip;
  // skip header
  fin.ignore(100000, '\n');
  fin.ignore(100000, '\n');
  while(fin){
    float x, y;
    fin >> x >> skip >> skip >> skip >> skip >> y >> skip >> skip
	>> skip;
    mass.push_back(x);
    radius.push_back(y);
  }
}
    
void get_data(const string& file, Subs::Buffer1D<float>& teff,
	      Subs::Buffer1D<float>& m70,
	      Subs::Buffer1D<float>& r70,
	      Subs::Buffer1D<float>& m75,
	      Subs::Buffer1D<float>& r75,
	      Subs::Buffer1D<float>& m80,
	      Subs::Buffer1D<float>& r80,
	      Subs::Buffer1D<float>& m85,
	      Subs::Buffer1D<float>& r85,
	      Subs::Buffer1D<float>& m90,
	      Subs::Buffer1D<float>& r90
	      ){

  ifstream fin;
  fin.open(file.c_str());
  if(!fin) {
    string error("Failed to open file: ");
    error.append(file);
    throw error;
  }
  float skip;
  // skip header
  fin.ignore(100000, '\n');
  while(fin){
    float x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11;
    fin >> x1 >>  x2 >> x3 >> x4 >> x5 >> x6 >> x7 >>x8
	>> x9 >> x10 >> x11;
    teff.push_back(x1);
    m70.push_back(x2);
    r70.push_back(x3);
    m75.push_back(x4);
    r75.push_back(x5);
    m80.push_back(x6);
    r80.push_back(x7);
    m85.push_back(x8);
    r85.push_back(x9);
    m90.push_back(x10);
    r90.push_back(x11);
  }
}

void read_panei(const string file, 
		Subs::Buffer1D<Subs::xyz<float,float,float> >& surface){
 ifstream fin;
  fin.open(file.c_str());
  if(!fin) {    string error("Failed to open file: ");
    error.append(file);
    throw error;
  }
  while(fin){
    char c = fin.peek();
    if (c == '#') {
      fin.ignore(100000, '\n'); // comment line
    } else {
      float x, y, z;
      fin >> x >> y >> z;
      Subs::xyz<float,float,float> tmp = Subs::xyz<float,float,float>(x,y,z);
      surface.push_back(tmp);
    }
  }

}

template <class X> X twodlinterp (const Subs::Array1D<X>& x,
				  const Subs::Array1D<X>& y,
				  const int& xpos,
				  const int& ypos,
				  const Subs::Buffer2D<X>& z,
				  const X& xgoal,
				  const X& ygoal){

  X val11, val12, val21, val22;
  val11 = z[xpos-1][ypos-1];
  val12 = z[xpos-1][ypos];
  val21 = z[xpos][ypos-1];
  val22 = z[xpos][ypos];
  cout << "quads: " << val11 << " " <<
    val12 << " " <<
    val21 << " " <<
    val22 << " " << endl;

  cout << "Goals: " << xgoal << " " << ygoal << endl;
 
  X term1 = val11*(x[xpos]-xgoal)*(y[ypos]-ygoal)/
    (x[xpos]-x[xpos-1])/(y[ypos]-y[ypos-1]);

  X term2 = val21*(xgoal-x[xpos-1])*(y[ypos]-ygoal)/
    (x[xpos]-x[xpos-1])/(y[ypos]-y[ypos-1]);

  X term3 = val12*(x[xpos]-xgoal)*(ygoal-y[ypos-1])/
    (x[xpos]-x[xpos-1])/(y[ypos]-y[ypos-1]);

  X term4 = val22*(xgoal-x[xpos-1])*(ygoal-y[ypos-1])/
    (x[xpos]-x[xpos-1])/(y[ypos]-y[ypos-1]);

 
  cout <<  term1 + term2 + term3 + term4 << endl << endl;
  return term1 + term2 + term3 + term4;
}

template <class X> X get_Y (const Subs::Array1D<X>& x, 
	     const Subs::Array1D<X>& y, X goal) { 

  X tmp;
  // check x array is ordered  
  if(! x.monotonic() ){
    Subs::Array1D<X> store(x);
    Subs::Buffer1D<unsigned long int> key = store.sort();
    int minpos = Subs::locate(store.ptr(),store.size(),goal);
    // linearly interpolate to find y at that x
    tmp = Subs::linterp(x[key[minpos-1]],y[key[minpos-1]],			    
			      x[key[minpos]],y[key[minpos]],goal);
  }else{
    int minpos = Subs::locate(x.ptr(),x.size(),goal);
    // linearly interpolate to find y at that x
    tmp = Subs::linterp(x[minpos-1],y[minpos-1],			    
			x[minpos],y[minpos],goal);
  }
  return tmp;
}

int main(void) {

  try{

    Subs::Array1D<float> teff1, m70, m80, m90, m75, m85,
      r70,r75,r80,r85,r90;    
    get_data("Bergeron/mr.dat",teff1,m70,r70,m75,r75,m80,r80,m85,r85,m90,r90);
    float rSol = 6.9550000e10;
    r70 /= rSol; r75 /= rSol; r80 /= rSol; r85 /= rSol; r90 /= rSol;
    
    Subs::Array1D<float> teff,mass,mbol,u,g,r,i,z,age,logg;
    read_colours("Bergeron/da2.ugriz",teff,logg,mass,mbol,u,g,r,i,z,age);
  
    // calc  colours
    Subs::Array1D<float> bergUminG(u), bergGminR(g), Mwd, 
      Twd, Rwd, Lgwd, UminGwd, sloanGwd;
    bergUminG -= g; bergGminR -= r;
    
     // set up to get radius from panei models
    Subs::Array1D<Subs::xyz<float,float,float> > panei;
    read_panei("Panei/12C-MR-H-He/12C-MR-H-He.dat",panei);  

    // calc Teff and Rwd for each mass, given Mdot and Townsley models
    float mdot = 3.0e-11;
    mdot /= 1.0e-10;
    int nip = 100;
    float pwr =0.25;
    for(int i=0; i < nip; i++){
      Mwd.push_back(0.405+0.75*float(i)/99.);
      Twd.push_back(17000*pow(mdot,pwr)*Mwd[i]/0.9);

      
      Subs::Array1D<float> tmpx, tmpy, splx, sply;
      
      for(int j=0; j<panei.size()-10; j+=5){
	float mp1 = panei[j].y;
	float mp2 = panei[j+5].y;
	if(Mwd[i] > mp1 && Mwd[i] < mp2){
	  for(int k=j; k<j+5; k++){
	    tmpx.push_back(panei[k].x);
	    tmpy.push_back(panei[k].z);
	  }
	  spline(tmpx,tmpy,500,splx,sply);
	  float div = 1000.0;
	  float rp1 = get_Y(splx,sply,Twd[i]/div)/100.0;
	  tmpx.clear(); tmpy.clear(); splx.clear(); sply.clear();

	  for(int k=j+5; k<j+10; k++){
	    tmpx.push_back(panei[k].x);
	    tmpy.push_back(panei[k].z);
	  }
	  spline(tmpx,tmpy,500,splx,sply);
	  float rp2 = get_Y(splx,sply,Twd[i]/div)/100.0;
	  tmpx.clear(); tmpy.clear(); splx.clear(); sply.clear();

	  // now interpolate to get radius at this mass and teff
	  Rwd.push_back( Subs::linterp(mp1,rp1,mp2,rp2,Mwd[i]) );
  
	  // calc log g
	  float Msun=1.989e30, G=6.67e-11, Rsun=6.9599e8;
	  Lgwd.push_back(log10(G*1.0e3*Mwd[i]*Msun*1000.0/
			       Rwd[i]/100.0/Rsun/Rwd[i]/100.0/Rsun));
	}
      }
    }
    
    // great. now we have a range of M, R, log g and T for our
    // accreting CVs. We need to get colours for each M, log g combination,
    // which is done by interpolating on a regular grid of Teff, log g from
    // the bergeron models.....
    int nlogg = 5; int nteff = ((teff.size()-1)/nlogg);
    Subs::Buffer2D<float> ugStore, sloangStore;
    Subs::Array1D<float> tstore,gstore;
    tstore.resize(nteff); gstore.resize(nlogg); 
    ugStore.resize(nteff,nlogg);
    sloangStore.resize(nteff,nlogg);

    for(int i=0; i < teff.size()-1; i++){
      int j = i-(nteff*(i/nteff));
      int k = (i/nteff);
      ugStore[j][k] = bergUminG[i];
      sloangStore[j][k] = g[i];
      tstore[j] = teff[i];
      gstore[k] = logg[i];
    }

    // loop over our teff, logg's and find UminG
    for(int i=0; i < nip; i++){
      int tpos, gpos;  
      if (!tstore.monotonic() || !gstore.monotonic()) {
       string err = "tstore/gstore not ordered";
       throw(err);
      }
      tpos = Subs::locate(tstore.ptr(),tstore.size(),Twd[i]);
      gpos = Subs::locate(gstore.ptr(),gstore.size(),Lgwd[i]);
      UminGwd.push_back( 
			twodlinterp(tstore,gstore,tpos,gpos,
				  ugStore,Twd[i],Lgwd[i])
			 );
      sloanGwd.push_back( 
			twodlinterp(tstore,gstore,tpos,gpos,
				  sloangStore,Twd[i],Lgwd[i])
			 );
    }

    Subs::Plot plot;
    plot.open("?");
    //cpgsubp(2,1);
    cpgsch(1.5);
    cpgscf(2);
    cpgenv(Mwd.min(), Mwd.max(), UminGwd.min(), UminGwd.max(), 0,0);
    cpgline(Mwd.size(), Mwd.ptr(), UminGwd.ptr());
    string xlab, ylab, plab;
    ylab = "u'-g'"; xlab = "M\\dwd\\u (M\\d\\(2281)\\u)"; plab="";
    cpglab(xlab.c_str(),ylab.c_str(),plab.c_str());
    cpgsci(2);
    cpgsfs(3);
    cpgrect(Mwd.min(),Mwd.max(),0.45,UminGwd.max());
    cpgsci(1);
    //cpgenv(Mwd.min(), Mwd.max(), sloanGwd.min(), sloanGwd.max(),0,0);
    //cpgline(Mwd.size(),Mwd.ptr(),sloanGwd.ptr());
    //ylab = "g'";
    //cpglab(xlab.c_str(),ylab.c_str(),plab.c_str());
    plot.close();

  }
  catch(const string& e){
    
    cerr << "\nError occured inside hist" << endl;
    cerr << e << endl;
    return 1;
  }


  return 0;

}
