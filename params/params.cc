#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "trm_subs.h"
#include "trm_buffer2d.h"
#include "trm_array1d.h"
#include "trm_plot.h"
#include "trm_constants.h"
#include <algorithm>

using namespace std;

// global variables
float rw, erw, q, eq,  inc, einc, p, ep;

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

template <class X> Subs::xy<X,X> meetPoint (const Subs::Array1D<X> x1, 
					    const Subs::Array1D<X> y1,
					    const Subs::Array1D<X> x2, 
					    const Subs::Array1D<X> y2,
					    bool& met){

  met=false;

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
  // between curves. Be sure not to step outside array here.
  Subs::Array1D<X> diff,xdiff,ydiff,diffsq;
  int begin = max(klo-2,start);
  int end = min(khi+2,stop);
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

template <class X> bool getParams (const Subs::Array1D<X> x1, 
				   const Subs::Array1D<X> y1,
				   const Subs::Array1D<X> ex1, 
				   const Subs::Array1D<X> ey1,
				   const Subs::Array1D<X> x2, 
				   const Subs::Array1D<X> y2) {
  
  /**
     
  This subroutine calculates the meeting point between two curves. Passed are
  a Subs:Array1D object representing the x and y values of each curve, and a 
  second set of arrays representing the "error curve" it then prints the binary
  parameters and their errors. Returns false and doesn't print anything if no meet point
 
   **/

  bool status1=false,status2=false;
  Subs::xy<float,float> meet = meetPoint(x1, y1, x2, y2, status1);
  Subs::xy<float,float> meethi = meetPoint(ex1, ey1, x2, y2, status2);
  if (status1 == false || status2 == false) {
    cout << "No Solution Found" << endl;
    return false;
  }
  
  float mw = meet.x; 
  float mwerr = abs(meet.x-meethi.x);

  // CHECK this is consistent
  float check = (Constants::G*meet.x*Constants::MSUN*(1+q)*
		 Subs::sqr(p)/4.0/Constants::PI/Constants::PI);
  check = rw*pow(check,float(1.0/3.0));
  float m2 = q*meet.x;
  float em2 = m2*sqrt(Subs::sqr(eq/q)+ Subs::sqr(mwerr/mw));
  float a = (Constants::G*(1.0+q)*mw*Constants::MSUN*Subs::sqr(p)
	     /4/Constants::PI/Constants::PI);
  a = pow(a,float(1.0/3.0))/Constants::RSUN;
  float aerr = a*sqrt(Subs::sqr(eq/(1.0+q)) + Subs::sqr(2.0*ep/p) +
		      Subs::sqr(mwerr/mw))/3.0;
  float sini = sin(inc*Constants::PI/180.0);
  float errsini=0.5*abs( sin((inc+einc)*Constants::PI/180.0) 
			 - sin((inc-einc)*Constants::PI/180.0) );
  float kw=sini*2.0*Constants::PI*q*a*
    Constants::RSUN/(1.0+q)/p/1.0e3;
  float kr=kw/q;
  float kwerr = kw*sqrt(Subs::sqr(errsini/sini) + Subs::sqr(eq/q) + 
			Subs::sqr(aerr/a) + Subs::sqr(eq/(1.0+q)) +
			Subs::sqr(ep/p) );
  float krerr = kr*sqrt(Subs::sqr(errsini/sini) + 
			Subs::sqr(aerr/a) + Subs::sqr(eq/(1.0+q)) +
			Subs::sqr(ep/p) );
  // secondary radius
  // use from Warner 1995 (eggleton)
  // radius of sphere with same volume as Roche Lobe...
  float r2 = 0.49*a*pow(q,float(2.0/3.0));
  r2 /= 0.6 * pow(q,float(2.0/3.0)) + log(1.0+pow(q,float(1.0/3.0)));
  float r2err = r2*sqrt( Subs::sqr(aerr/a) + Subs::sqr(2.0*eq/q/3.0) +
			 Subs::sqr(2.0*eq/3.0/pow(q,float(1.0/3.0)))  + 
			 Subs::sqr(eq/(3.0*pow(q,float(2.0/3.0))*(1.0+q)))/
			 (0.6 * pow(q,float(2.0/3.0)) + log(1.0+pow(q,float(1.0/3.0)))));
  
  cout << endl << endl << "M1    =   " << meet.x << " +/-  "  << 
    mwerr << "  Msun" << endl;
  cout  << "R1    =   " << meet.y << " +/-  " 
	<<abs(meet.y-meethi.y) << "  Rsun" << endl;
  cout << "R1(2) =   " << check/Constants::RSUN << endl;;
  cout << "M2    =   " << m2 <<  " +/- " << em2 << " Msun" << endl;
  cout << "R2    =   " << r2 <<  " +/- " << r2err << " Rsun" << endl << endl;
  cout << "a     =   " << a <<  " +/- " << aerr << " Rsun" << endl;
  cout << "Kw    =   " << kw <<  " +/- " << kwerr << " km/s" << endl;
  cout << "Kr    =   " << kr <<  " +/- " << krerr << " km/s" << endl<< endl;
  return true;  
}





int main(void) {

  try{

    Subs::Array1D<float> teff, r;

    Subs::Plot plot("?");
    cpgsch(1.5); cpgscf(2); 
    string xlab, ylab, plab;

    // plot mass against radius for allowed teff
    /* for each mass track we need to
       a) find radius at selected teff.
       b) find radius at lower end
       c) find radius at upper end.

       we do this using the Wood95 CO models
    */
    Subs::Array1D<float> mass, radius, rhi, rlo, rhi2, rlo2;
    float thi,tlo,tbest, terr;
    cout << "Enter Twd and error: ";
    cin >> tbest >> terr;
    thi =tbest+terr; tlo=tbest-terr;
    float tlo2=tbest-5.0*terr,thi2=tbest+5.0*terr;
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
      rhi.push_back(get_Y(teff,r,thi));
      rlo.push_back(get_Y(teff,r,tlo));
      rhi2.push_back(get_Y(teff,r,thi2));
      rlo2.push_back(get_Y(teff,r,tlo2));
      teff.clear(); r.clear(); 
    }

    // spline interpolate function
    Subs::Array1D<float> spline_x, spline_y;
    spline(mass,radius,2000,spline_x,spline_y);
    // plot main relation
    //cpgenv(spline_x.min(),spline_x.max(),spline_y.min(),spline_y.max(),0,0);
    //cpgenv(0.6,1.0,0.007,0.013,0,0);
    cpgenv(0.1,1.44,0.0005,0.028,0,0);
    xlab="M\\dwd\\u/M\\d\\(2281)";
    ylab="R\\dwd\\u/R\\d\\(2281)";
    cpglab(xlab.c_str(),ylab.c_str(),plab.c_str());
    cpgsci(2);
    cpgline(spline_x.size(),spline_x.ptr(),spline_y.ptr());
    // plot errors
    cpgsci(2); cpgsls(2);
    spline_x.clear(); spline_y.clear();
    spline(mass,rhi,2000,spline_x,spline_y);
    cpgline(spline_x.size(),spline_x.ptr(),spline_y.ptr());
    spline_x.clear(); spline_y.clear();
    spline(mass,rlo,2000,spline_x,spline_y);
    cpgline(spline_x.size(),spline_x.ptr(),spline_y.ptr());

    /**
    // plot extreme errors
    cpgsci(4);
    spline_x.clear(); spline_y.clear();
    spline(mass,rhi2,2000,spline_x,spline_y);
    cpgline(spline_x.size(),spline_x.ptr(),spline_y.ptr());
    spline_x.clear(); spline_y.clear();
    spline(mass,rlo2,2000,spline_x,spline_y);
    cpgline(spline_x.size(),spline_x.ptr(),spline_y.ptr());
    **/

    // fix spline_x so it refers to best fit again...
    spline_x.clear(); spline_y.clear();
    spline(mass,radius,2000,spline_x,spline_y);    
    // plot non-interpolated points
    cpgsci(1); cpgsls(2);
    cpgpt(mass.size(),mass.ptr(),radius.ptr(),4);

    // hamada mass radius relationship
    Subs::Array1D<float> mh, rh;
    hamada("Hamada/mr_hamada.dat",mh,rh);
    rh = rh/100.0;
    cpgsci(3);
    cpgline(mh.size(),mh,rh);
    cpgsci(1);

    // benvenuto 15K mass radius relationship
    //Subs::Array1D<float> mb, rb;
    //read_benv("benvenuto/15000K.dat",mb,rb);
    //cpgsci(6);
    //cpgline(mb.size(),mb,rb);
    //cpgsci(1);

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
    // spline interpolate and plot
    Subs::Array1D<float> splMp, splRp;
    spline(mp,rp,2000,splMp,splRp);
    cpgsci(3); cpgsls(1);
    cpgline(splMp.size(),splMp.ptr(),splRp.ptr());
    cpgsci(1); cpgsls(2);


    //panei iron core
    panei.clear();
    read_panei("Panei/othercores/56Fe-MR-H-He.dat",panei);
    Subs::Array1D<float> mfe, rfe;
    for(int i=0; i<panei.size()-7; i+=7){
      Subs::Array1D<float> tmpx, tmpy, splx, sply;
      mfe.push_back(panei[i].y);
      for (int j=i; j<i+7; j++){
	tmpx.push_back(panei[j].x);
	tmpy.push_back(panei[j].z);
      }
      spline(tmpx,tmpy,500,splx,sply);
      float div = 1000.0;
      rfe.push_back(get_Y(splx,sply,tbest/div)/100.0);
      tmpx.clear(); tmpy.clear();
      splx.clear(); sply.clear();
    }
    // spline interpolate and plot
    Subs::Array1D<float> splMfe, splRfe;
    spline(mfe,rfe,2000,splMfe,splRfe);
    cpgsci(3); cpgsls(1);
    cpgline(splMfe.size(),splMfe.ptr(),splRfe.ptr());
    cpgsci(1); cpgsls(2);


    //panei He core
    panei.clear();
    read_panei("Panei/othercores/04He-MR-H.dat",panei);
    Subs::Array1D<float> mhe, rhe;
    for(int i=0; i<panei.size()-5; i+=5){
      Subs::Array1D<float> tmpx, tmpy, splx, sply;
      mhe.push_back(panei[i].y);
      for (int j=i; j<i+5; j++){
	tmpx.push_back(panei[j].x);
	tmpy.push_back(panei[j].z);
      }
      spline(tmpx,tmpy,500,splx,sply);
      float div = 1000.0;
      rhe.push_back(get_Y(splx,sply,tbest/div)/100.0);
      tmpx.clear(); tmpy.clear();
      splx.clear(); sply.clear();
    }
    // spline interpolate and plot
    Subs::Array1D<float> splMhe, splRhe;
    spline(mhe,rhe,2000,splMhe,splRhe);
    cpgsci(3); cpgsls(1);
    cpgline(splMhe.size(),splMhe.ptr(),splRhe.ptr());
    cpgsci(1); cpgsls(2);

    // althaus models for hi mass ONe WDs
    Subs::Array1D<float> ma, ra;
    teff.clear(); r.clear();
    for(int i=106; i<=130; i+=2){
      ma.push_back(float(i)/100.0);
      // read in m-r relationship
      string file = "althaus/t106.one";
      stringstream append;
      append << setw(3) << setfill('0') << i;
      file.replace(9,3,append.str());
      read_althaus(file,teff,r);
      // linearly interpolate to find radius at this teff.
      ra.push_back(get_Y(teff,r,tbest));
      teff.clear(); r.clear(); 
    }
    // spline interpolate and plot
    Subs::Array1D<float> spline_altx, spline_alty;
    spline(ma,ra,2000,spline_altx,spline_alty);
    cpgsci(5);
    cpgline(spline_altx.size(),spline_altx.ptr(),spline_alty.ptr());
    cpgsci(1);

    // now let's plot M-R from kepler's law and radius constraints.
    cout << "Enter Rw/a and error: ";
    cin >> rw >> erw;
    cout << "Enter q and error: ";
    cin >> q >> eq;
    cout << "Enter period (days) and error: ";
    cin >> p >> ep;
    p*=Constants::DAY; ep*=Constants::DAY;
    cout << "Enter inclination and error: ";
    cin >> inc >> einc;

    float t = Constants::G*(1+q)*Subs::sqr(p)/4.0/Subs::sqr(Constants::PI);
    float dt = sqrt(Subs::sqr(eq/(1.0+q)) + Subs::sqr(2.0*ep/p));
    // for each value of a, we get a value of M1 from t and R1 from R1/a
    Subs::Array1D<float> mgeom, emgeom, rgeom, ergeom;
    for (int i=0; i<=2000; i++){
      float a = Constants::RSUN * (0.1 + 0.9*float(i)/1000.0);
      mgeom.push_back( (pow(a,float(3.0))/t)/Constants::MSUN );
      emgeom.push_back(mgeom[i]*dt);
      rgeom.push_back(rw*a/Constants::RSUN);
      ergeom.push_back(erw*rgeom[i]/rw);
    }
    cpgsls(1);
    cpgline(mgeom.size(),mgeom.ptr(),rgeom.ptr());
    cpgsls(2);
    Subs::Array1D<float> xlo(mgeom.size()), ylo(mgeom.size());
    Subs::Array1D<float> xhi(mgeom.size()), yhi(mgeom.size());
    xlo = mgeom-emgeom;
    ylo = rgeom-ergeom;
    cpgline(mgeom.size(),xlo.ptr(),ylo.ptr());
    xhi = mgeom+emgeom;
    yhi = rgeom+ergeom;
    cpgline(mgeom.size(),xhi.ptr(),yhi.ptr());
    plot.close();
    
    // find solutions for various white dwarf models
    // Wood 95
    cout << endl << "===== Wood 95 (thick layer models) =====" << endl;
    bool status = getParams(mgeom,rgeom,xhi,yhi,spline_x,spline_y);
    cout << endl << "===== Panei (thin layer models) =====" << endl;
    status = getParams(mgeom,rgeom,xhi,yhi,mp,rp);
    cout << endl << "===== Panei Fe core models =====" << endl;
    status = getParams(mgeom,rgeom,xhi,yhi,mfe,rfe);
    cout << endl << "===== Althaus ONe models =====" << endl;
    status = getParams(mgeom,rgeom,xhi,yhi,spline_altx,spline_alty);
    cout << endl << "===== Hamada-Saltpeter (C) =====" << endl;
    status = getParams(mgeom,rgeom,xhi,yhi,mh,rh);
    
  }
  catch(const string& e){
    
    cerr << "\nError occured inside parameters" << endl;
    cerr << e << endl;
    return 1;
  }


  return 0;

}

			     
