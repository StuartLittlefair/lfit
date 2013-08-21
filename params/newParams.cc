#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include "trm_subs.h"
#include "trm_buffer2d.h"
#include "trm_array1d.h"
#include "trm_vec3.h"
#include "trm_roche.h"
#include "trm_plot.h"
#include "trm_constants.h"

using namespace std;

class ProgressBar{
		
public:
  
  // Default constructor 
  ProgressBar() {
    message="";
    messageLength=0;
    barLength=0;
  }
  
  // constructor from string and length
  ProgressBar(std::string message_, int length_, int duration_){
    this->duration = duration_;
    this->message = message_ + ": completion 0% -";
    this->secondLineLength = this->message.size();
    this->barLength = length_;
    for(size_t i=0; i< this->barLength; i++) this->message += "#";
    this->message += "- 100%\n";
    this->messageLength = this->message.size();
    this->completion = 0;
  }
  
  // destructor
  ~ProgressBar(){
    //std::cout <<"Destroying progress bar"<<std::endl;
  }
  
  // print out initial message Bar
  bool initialise(){
    if (this->barLength == 0){
      std::string err = "Error: message bar not initialised";
      throw err;
    }
    std::cout << this->message;
    std::string secondLine;
    for(size_t i=0; i<this->secondLineLength; i++) secondLine += " ";
    std::cout << secondLine;
    return true;
  }
  
  bool update(unsigned int progressCount){
    if (( (this->barLength+1)*progressCount/duration - this->completion) > 0 ){
      int toPrint = 
	(this->barLength+1)*progressCount/duration - this->completion;
      for(int i=0; i< toPrint; i++) {
	std::cout << "-";
	std::cout.flush();
      }
      this->completion++;
    }
    return true;
  }
  
  bool end(){
    std::cout << std::endl;
    return true;
  }
  
private:
  std::string message;
  unsigned int messageLength;
  unsigned int secondLineLength;
  unsigned int barLength;
  unsigned int duration;
  unsigned int completion;
};

Subs::Vec3 getEarthVec(float phase, float angle){

  /* returns vector pointing to earth given binary phase
     and inclination angle */
  float rad_incl = angle*2.0*Constants::PI/360.0;
  float rad_phase = phase*2.0*Constants::PI;
  float x = cos(rad_phase)*sin(rad_incl);
  float y = -sin(rad_phase)*sin(rad_incl);
  float z = cos(rad_incl);
  Subs::Vec3 earth = Subs::Vec3(x,y,z);
  return earth;

}

float calc_inc(const float q, const float dphi){

  Subs::Vec3 wd=Subs::Vec3();
  float phase = dphi/2.0;

  float angle1=54.0; float angle2=90.0;
  for (int j=0; j<=14; j++){
    float angle = (angle1+angle2)/2.0;
    Subs::Vec3 earth = getEarthVec(phase,angle);
    if( Roche::blink(q,wd,earth,0.001) ) {
      angle2=angle;
    } else {
      angle1=angle;
    }
  }
  return 0.5*(angle1+angle2);

}

void read_mcmc(const string& file, Subs::Buffer1D<float>& qBuf,
	       Subs::Buffer1D<float>& dphiBuf, Subs::Buffer1D<float>& rwBuf,
	       const int& qcol, const int& dphicol, const int& rwcol){

  ifstream fin;
  fin.open(file.c_str());
  if(!fin) {
    string error("Failed to open file: ");
    error.append(file);
    throw error;
  }
  qBuf.clear();
  dphiBuf.clear();
  rwBuf.clear();

  // skip header
  float skip;
  float qIn, dphiIn, rwIn;
  while(fin){

    char c = fin.peek();
    if(c == '#' || c == 'P' || c == 'B'){
      //skip comment and blank lines
      fin.ignore(100000,'\n');
    }else{
      float skip;
      int count=0;
      string str;
      getline(fin,str);
      istringstream iss(str);
      while(iss){
	iss >> skip;
	if(count == qcol)    qIn = skip;
	if(count == dphicol) dphiIn = skip;
	if(count == rwcol)   rwIn = skip;
	count++;
      }
      qBuf.push_back(qIn);
      dphiBuf.push_back(dphiIn);
      rwBuf.push_back(rwIn);		   
    }

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

void read_althaus(const string& file, Subs::Buffer1D<float>& teff, 
	    Subs::Buffer1D<float>& radius){
  
  ifstream fin;
  fin.open(file.c_str());
  if(!fin) {    string error("Failed to open file: ");
    error.append(file);
    throw error;
  }
  float skip;
  float rSol = 6.955e10;
  float ten = 10.0;
  while(fin){
    float x, y;
    fin >> skip >> x >> y >> skip >> skip;
    
    teff.push_back(pow(ten,x));
    radius.push_back(pow(ten,y)/rSol);
  }
}

void read_wood(const string& file, Subs::Buffer1D<float>& teff, 
	    Subs::Buffer1D<float>& radius){
  
  ifstream fin;
  fin.open(file.c_str());
  if(!fin) {    string error("Failed to open file: ");
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

void read_benv(const string& file, Subs::Buffer1D<float>& mass, 
	    Subs::Buffer1D<float>& radius){
  
  ifstream fin;
  fin.open(file.c_str());
  if(!fin) {    string error("Failed to open file: ");
    error.append(file);
    throw error;
  }
  float skip;
  float rSol = 6.955e10;
  float ten = 10.0;
  while(fin){
    float x, y;
    fin >> skip >> x >> y;
    mass.push_back(x);
    radius.push_back(y/100.0);
  }
}

template <class X> void splint(Subs::Buffer1D<X>& xa, Subs::Buffer1D<X>& ya, 
			       Subs::Buffer1D<X>& y2a, const X& x, X &y){
	int k;
        X h,b,a;

	int n=xa.size();
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) 
	  throw Subs::Subs_Error("splint(Subs::Buffer1D<X>& xa, Subs::Buffer1D<X>& ya, Subs::Buffer1D<X>& y2a, const X& x, X &y): Bad xa input");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
		+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

template <class X> void spline(Subs::Buffer1D<X>& x, Subs::Buffer1D<X>& y, 
			       const X& yp1, const X& ypn, Subs::Buffer1D<X>& y2){

	int i,k;
	X p,qn,sig,un;

	int n=y2.size();
	Subs::Buffer1D<X> u(n-1);
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}


template <class X> void spline(const Subs::Array1D<X>&  x, const Subs::Array1D<X>&  y,
	    unsigned long size, Subs::Array1D<X>& xint, 
	    Subs::Array1D<X>& yint){

  if(x.size() != y.size()) 
    throw Subs::Subs_Error("Error returned inside spline: x and y arrays differ in size\n");
  Subs::Buffer1D<X> vecx(x.size()), vecy(x.size()), d2y(x.size());
  for(int i=0; i< x.size(); i++){vecx[i] = x[i]; vecy[i]=y[i];}
  X yp1, ypn;
  yp1 = 1.0e30; ypn = 1.0e30;
  spline(vecx,vecy,yp1,ypn,d2y);
  for(int i=0; i<size; i++){
    X xtmp, ytmp;
    xtmp = x.min() + (x.max()-x.min())*i/size;
    splint(vecx,vecy,d2y,xtmp,ytmp);
    xint.push_back(xtmp);
    yint.push_back(ytmp);
  }      
}

template <class X> X get_Y (const Subs::Array1D<X> x, 
	     const Subs::Array1D<X> y, X goal) { 

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

template <class X> Subs::xy<X,X> meetPoint (const Subs::Array1D<X> x1, 
					    const Subs::Array1D<X> y1,
					    const Subs::Array1D<X> x2, 
					    const Subs::Array1D<X> y2){

  bool met=false;

  // only work on overlap
  int start, stop;
  start=0; stop=x1.size()-1;
  while(x1[stop] != Subs::clamp(x2.min(),x1[stop],x2.max())){
    stop--;
  }
  while(x1[start] != Subs::clamp(x2.min(),x1[start],x2.max())){
    start++;
  }
  
  // check they meet at all
  X sval = get_Y(x2,y2,x1[start]);
  X eval = get_Y(x2,y2,x1[stop]);
  X sdiff = y1[start]-sval; float ediff=y1[stop]-eval;
  if(ediff == Subs::sign(ediff,sdiff)) {
    return Subs::xy<X,X>(0.0,0.0);
    //throw Subs::Subs_Error("template <class X> Subs::xy<X,X> meetPoint (const Subs::Array1D<X> x1, const Subs::Array1D<X> y1,const Subs::Array1D<X> x2, const Subs::Array1D<X> y2): curves do not cross");
  }
  
  // binary chop search
  int klo=start; int khi=stop; int k;
  while(khi-klo>1){
    k = (khi+klo) >> 1; // split the differnce
    X tval=get_Y(x2,y2,x1[k]);
    X tdiff=y1[k]-tval;
    if(tdiff == Subs::sign(tdiff,sdiff)){
      // sign hasn't changed, we're still before the meet point
      klo=k;
    }else{
      khi=k;
    }
  }

  // loop over selected range in first array, and calculate difference
  // between curves
  Subs::Array1D<X> diff,xdiff,ydiff,diffsq;
  int begin, end;
#pragma omp critical
  {
    begin = max(klo-2,start);
    end = min(khi+2,stop);
  }
  for(int i=begin; i<=end; i++){
    // x lies within range of second array
    X val2 = get_Y(x2, y2, x1[i]);
    diff.push_back(y1[i]-val2);
    ydiff.push_back(y1[i]);
    xdiff.push_back(x1[i]);
    diffsq.push_back(abs(ydiff[i]));
  }
  X zero = 0.0;
  X xval = get_Y(diff,xdiff,zero);
  X yval = get_Y(diff,ydiff,zero);
  met = true;
  return Subs::xy<X,X>(xval,yval);
}

template <class X> Subs::Buffer1D<float> getParams (const Subs::Array1D<X> x1, 
						    const Subs::Array1D<X> y1,
						    const Subs::Array1D<X> x2, 
						    const Subs::Array1D<X> y2,
						    const float q, 
						    const float p,
						    const float rw,
						    const float inc) {
  
  /**
     
  This subroutine calculates the meeting point between two curves. Passed are
  a Subs:Array1D object representing the x and y values of each curve, and a 
  second set of arrays representing the "error curve" it then prints the binary
  parameters and their errors. Returns false and doesn't print anything if no meet point
 
   **/

  Subs::xy<float,float> meet = meetPoint(x1, y1, x2, y2);

  std::string err= "No Solution Found";
  if(meet.x == 0.0 && meet.y == 0.0){throw err;}

  float mw = meet.x; 

  // CHECK this is consistent
  float check = (Constants::G*meet.x*Constants::MSUN*(1+q)*
		 p*p/4.0/Constants::PI/Constants::PI);
  check = rw*pow(check,float(1.0/3.0));
  float m1 = meet.x;
  float r1 = meet.y;

  float m2 = q*meet.x;

  float a = (Constants::G*(1.0+q)*mw*Constants::MSUN*p*p
	     /4./Constants::PI/Constants::PI);
  a = pow(a,float(1.0/3.0))/Constants::RSUN;
  float sini = sin(inc*Constants::PI/180.0);
  float kw=sini*2.0*Constants::PI*q*a*
    Constants::RSUN/(1.0+q)/p/1.0e3;
  float kr=kw/q;
  // secondary radius
  // use from Warner 1995 (eggleton)
  // radius of sphere with same volume as Roche Lobe...
  float r2 = 0.49*a*pow(q,float(2.0/3.0));
  r2 /= 0.6 * pow(q,float(2.0/3.0)) + log(1.0+pow(q,float(1.0/3.0)));

  Subs::Buffer1D<float> results(7);
  results[0] = m1;
  results[1] = r1;
  results[2] = m2;
  results[3] = r2;
  results[4] = a;
  results[5] = kw;
  results[6] = kr;

  return results;  
}



int main(void) {

  try{

    Subs::Array1D<float> teff, r;
    float thi,tlo,tbest, terr;
    cout << "Enter Twd: ";
    cin >> tbest;
    Subs::Array1D<float> mass, radius;
    for(int i=40; i<=100; i+=10){
      mass.push_back(float(i)/100.0);
      // read in m-r relationship
      string file = "Wood95/co040.2.04dt1";
      stringstream append;
      append << setw(3) << setfill('0') << i;
      file.replace(9,3,append.str());
      read_wood(file,teff,r);

      // linearly interpolate to find radius at this teff.
      radius.push_back(get_Y(teff,r,tbest));
      teff.clear(); r.clear(); 
    }
    // spline interpolate function
    Subs::Array1D<float> wood_splineX, wood_splineY;
    spline(mass,radius,2000,wood_splineX,wood_splineY);

    // panei models
    Subs::Array1D<Subs::xyz<float,float,float> > panei;
    read_panei("Panei/12C-MR-H-He/12C-MR-H-He.dat",panei);
    Subs::Array1D<float> mp, rp;
    for(int i=0; i<panei.size()-5; i+=5){
      Subs::Array1D<float> tmpx, tmpy, splx, sply;
      mp.push_back(panei[i].y);
      for (int j=i; j<i+5; j++){
	tmpx.push_back(panei[j].x);
	tmpy.push_back(panei[j].z);
      }
      spline(tmpx,tmpy,500,splx,sply);
      float div = 1000.0;
      rp.push_back(get_Y(splx,sply,tbest/div)/100.0);
      tmpx.clear(); tmpy.clear();
      splx.clear(); sply.clear();
    }

    // hamada mass radius relationship
    Subs::Array1D<float> mh, rh;
    hamada("Hamada/mr_hamada.dat",mh,rh);
    rh = rh/100.0;

    Subs::Buffer1D<float> qBuf,dphiBuf,rwBuf;
    int qcol=1,dphicol=2,rwcol=4;
    string file;
    cout << "Enter MCMC file name: ";
    cin >> file;
    read_mcmc(file,qBuf,dphiBuf,rwBuf,qcol,dphicol,rwcol);
    int np = qBuf.size();
    cout << "Read in " << np << " lines" << endl;


    Subs::Array1D<float> m1Buf(np), r1Buf(np), m2Buf(np), r2Buf(np), aBuf(np), kwBuf(np), krBuf(np);

    cout << "Enter period and error (days): ";
    float p, ep;
    cin >> p >> ep;
    p*=Constants::DAY; ep*=Constants::DAY;
    cout << endl;

    ProgressBar bar = ProgressBar("Calculating masses from MCMC steps",50,np);
    bar.initialise();
    ofstream fout;
    fout.open("newparams.log");
    fout << "#q M1  R1  M2  R2  a  Kw  Kr incl\n";
    
    int icount=0;
#pragma omp parallel for shared(icount,bar)
    for(int k=0; k<np; k++){


      // trying pragma single instead of critical
#pragma omp critical
      {
	icount++;
	bar.update(icount);
      }

      
      float q          = qBuf[k];
      float dphi       = dphiBuf[k];
      float rw         = rwBuf[k]*Roche::xl1(q);
      float inc = calc_inc(q,dphi);
      float pThis = p + Subs::gauss3(k)*ep;      

      float t = Constants::G*(1.+q)*pThis*pThis/4.0/Constants::PI/Constants::PI;
      

      // for each value of a, we get a value of M1 from t and R1 from R1/a
      Subs::Array1D<float> mgeom(2000), rgeom(2000);
      for (int i=0; i<=2000; i++){
	float a = Constants::RSUN * (0.1 + 0.9*float(i)/1000.0);
	mgeom[i] = (pow(a,float(3.0))/t)/Constants::MSUN ;
	rgeom[i] = rw*a/Constants::RSUN;
      }

      Subs::Buffer1D<float> results(7);
      for(int bob=0; bob<results.size(); bob++){results[bob]=0.0;}
      bool fail = false;
      try{
	results = getParams(mgeom,rgeom,wood_splineX,wood_splineY,q,pThis,rw,inc);
      }catch(const std::string err){
	fail = true;
      }
      if(fail){
	fail=false;
	try{
	  // wood models failed. Try a panei model.
	  results = getParams(mgeom,rgeom,mp,rp,q,pThis,rw,inc);
	}catch(const std::string err){
	  fail = true;
	}
      }
      if(fail){
	fail = false;
	try{
	  //fout << "# " << q << " " << pThis << " " << rw << " " << inc << endl;
	  results = getParams(mgeom, rgeom, mh, rh, q, pThis,rw,inc); 
	  //fout << "# " << q << " ";
	  //for(int j=0; j<results.size(); j++) fout << results[j] << " ";
	  //fout << inc << endl;
	}catch(const std::string err){
	  fail = true;
	}
      }
      
      //if (fail) {cout << "No solution found" << endl;}
      
#pragma omp critical
      {
	if (!fail){
	  m1Buf[k] = results[0];
	  r1Buf[k] = results[1];
	  m2Buf[k] = results[2];
	  r2Buf[k] = results[3];
	  aBuf[k]  = results[4];
	  kwBuf[k] = results[5];
	  krBuf[k] = results[6];
	  fout << q << " ";
	  for(int j=0; j<results.size(); j++) fout << results[j] << " ";
	  fout << inc << endl;
	}
      }
      
      
    }
    bar.end();
    fout.close();
    Subs::Plot plot("/xs");
    cpgsch(1.5); cpgscf(2); 
    string xlab, ylab, plab;

    cpgenv(m2Buf.min(),m2Buf.max(),m1Buf.min(),m1Buf.max(),0,0);
    cpgpt(m2Buf.size(),m2Buf.ptr(),m1Buf.ptr(),1);
    xlab = "M\\dr\\u (M\\d\\(2281)\\u)"; ylab = "M\\dw\\u (M\\d\\(2281)\\u)"; plab="";
    cpglab(xlab.c_str(), ylab.c_str(), plab.c_str());
    plot.close();
  }
  catch(const string& e){
    
    cerr << "\nError occured inside parameters" << endl;
    cerr << e << endl;
    return 1;
  }


  return 0;

}

			     
