#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
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


int main(void) {

  try{

    Subs::Array1D<float> teff1, m70, m80, m90, m75, m85,
      r70,r75,r80,r85,r90;    
    get_data("Bergeron/mr.dat",teff1,m70,r70,m75,r75,m80,r80,m85,r85,m90,r90);
    float rSol = 6.9550000e10;
    r70 /= rSol; r75 /= rSol; r80 /= rSol; r85 /= rSol; r90 /= rSol;
      
    Subs::Plot plot;
    
    Subs::Array1D<float> teff,mass,mbol,u,g,r,i,z,age,logg;
    read_colours("Bergeron/da2.ugriz",teff,logg,mass,mbol,u,g,r,i,z,age);
    map<string,Subs::Array1D<float> > red;
    red["r"] = r;
    red["i"] = i;
    red["z"] = z;

    float up, uperr, gp, gperr, rp, rperr, av;
    cout << "Give u' band flux and err: ";
    cin >> up >> uperr;
    cout << "Give g' band flux and err: ";
    cin >> gp >> gperr;
    string reply;
    bool needBand = true;
    string redKey;
    while(needBand){
      cout << "which band in red? ";
      cin >> reply;
      if(Subs::toupper(reply) == "R"){
        cout << "Give r' band flux and err: ";
        cin >> rp >> rperr;
        needBand = false;
        redKey = "r";
      }else if(Subs::toupper(reply) == "I"){
        cout << "Give i' band flux and err: ";
        cin >> rp >> rperr;
        needBand = false;
        redKey="i";
      }else if(Subs::toupper(reply) == "Z"){
        cout << "Give z' band flux and err: ";
        cin >> rp >> rperr;
        needBand = false;
        redKey="z";
      }
      if(needBand){cout << "please enter r, i or z... " << endl;}
    }
    cout << "Give extinction (Av): ";
    cin >> av;
    // convert to magnitudes
    uperr = abs(-log10((up+0.01*up+uperr)/3631.0/1000.0)/0.4 
		+ log10((up-0.01*up-uperr)/3631.0/1000.0)/0.4)/2.0; 
    up = -log10(up/3631.0/1000.0)/0.4;
    gperr = abs(-log10((gp+0.01*gp+gperr)/3631.0/1000.0)/0.4 
		+ log10((gp-0.01*gp-gperr)/3631.0/1000.0)/0.4)/2.0 ;
    gp = -log10(gp/3631.0/1000.0)/0.4;
    rperr = abs(-log10((rp+0.01*rp+rperr)/3631.0/1000.0)/0.4 
		+ log10((rp-0.01*rp-rperr)/3631.0/1000.0)/0.4)/2.0;
    rp = -log10(rp/3631.0/1000.0)/0.4;

    // calculate colours
    float uming, del_uming, gminr, del_gminr;
    uming = up-gp; 
    del_uming = sqrt( Subs::sqr(gperr) + Subs::sqr(uperr) );
    gminr = gp-rp;
    del_gminr = sqrt( Subs::sqr(rperr) + Subs::sqr(gperr) );
      

    // plot colours
    Subs::Array1D<float> bergUminG(u), bergGminR(g);
    bergUminG -= g; bergGminR -= red[redKey];
    plot.open("?");
    cpgsch(1.5);
    cpgslw(2);
    cpgscf(2);
    
    cpgenv(bergUminG.min(), bergUminG.max(), bergGminR.min(), bergGminR.max(),0,0);
    string axisOptions = "BCNST";
    cpgbox(axisOptions.c_str(),0.0,0,axisOptions.c_str(),0.0,0);

    cpgpt(bergUminG.size(),bergUminG.ptr(),bergGminR.ptr(),1);
    string xlab, ylab, plab;
    xlab = "u'-g'"; ylab = "g'-" + redKey + "'"; plab="";
    cpglab(xlab.c_str(),ylab.c_str(),plab.c_str());
    cpgsci(2);
    cpgpt(1,&uming,&gminr,4);
    cpgerr1(5,uming,gminr,del_uming,0.2);
    cpgerr1(6,uming,gminr,del_gminr,0.2);
    // plot de-reddening line
    float umingprime, gminrprime;
    float ebv = av/3.1;
    float au = ebv*5.155;
    float ag = ebv*3.793;
    float ar;
    if (redKey == "r"){ar = 2.751*ebv;}
    if (redKey == "i"){ar = 2.086*ebv;}
    if (redKey == "z"){ar = 1.479*ebv;}
    // de-reddened colours
    umingprime = uming + ag - au;
    gminrprime = gminr + ar - ag;
    cpgsch(0.5);
    cpgarro(uming,gminr,umingprime,gminrprime);
    plot.close();

    Subs::Array1D<float> chisq(teff.size()-1), tstore, gstore;
    Subs::Buffer2D<float> chisq2, uMagStore, rMagStore, gMagStore;
    int nlogg = 5;
    int nteff = ((teff.size()-1.0)/nlogg);
    chisq2.resize(nteff,nlogg);
    tstore.resize(nteff); gstore.resize(nlogg); 
    uMagStore.resize(nteff,nlogg); gMagStore.resize(nteff,nlogg); 
    rMagStore.resize(nteff,nlogg);
    for(int i=0; i < teff.size()-1; i++){
      // get colours for this teff and logg
      float col1 = u[i]-g[i]; float col2=g[i]-red[redKey][i];
      float tmp = Subs::sqr( (col1-uming)/del_uming ) +
            Subs::sqr( (col2-gminr)/del_gminr );
      int j = i-(nteff*(i/nteff));
      int k = (i/nteff);
      //cout << j << " " << k << endl;
      chisq[i] = log10(tmp);
      chisq2[j][k] = tmp;
      uMagStore[j][k] = u[i];
      gMagStore[j][k] = g[i];
      rMagStore[j][k] = red[redKey][i];
      tstore[j]=teff[i];
      gstore[k]=logg[i];
    }
    // Plot chisq surface
    plot.open("chisq.ps/cps");
    cpgsch(1.5);
    cpgscf(2);
    cpgsubp(2,2);
    Subs::Buffer1D<float> tr(6);
    tr[3]=6.5; tr[4]=0.0; tr[5]=0.5;
    tr[0]=5000.0;tr[2]=0.0;tr[1]=500.0;
    cpgenv(7000.0,15000.0,6.5,9.5,0,0);
    cpggray(chisq.ptr(),nteff,nlogg,1,nteff,1,nlogg,5.0,1.0,tr.ptr());
    string side = "T"; string label = "log\\d10\\u\\gx\\u2";
    cpgwedg(side.c_str(),0.0,3.0,5.0,1.0,label.c_str());
    xlab = "T\\deff";
    ylab = "log g";
    plab = "";
    cpglab(xlab.c_str(),ylab.c_str(),plab.c_str());

    // find error in teff
    //loop over each teff
    Subs::Array1D<float> chimin, umin, gmin, rmin;

    string fix="n";
    float logg_fix;
    cout << "Restrict log g (y/n)? ";
    cin >> fix;
    fix = Subs::tolower(fix);
    if(fix == "y"){
      cout << "\nEnter log g: ";
      cin >> logg_fix;
    }
    if(fix == "n"){
        // method one: no constraints on log g
        for(int i=0; i<nteff; i++){
            // find the logg which gives the minimum chisq at each teff
            // and record which chisq it gives for each teff
            Subs::Array1D<float> store(nlogg);
            for (int j=0; j<nlogg; j++){
              store[j] = chisq2[i][j];
            }
            chimin.push_back(store.min());
            int jpos = minloc(store);
            umin.push_back(uMagStore[i][jpos]);
            gmin.push_back(gMagStore[i][jpos]);
            rmin.push_back(rMagStore[i][jpos]);
        }
    } else {
        // method 2: constrain log g from mass
        for(int i=0; i<nteff; i++){
            for (int j=0; j<nlogg; j++){
              if(gstore[j]==logg_fix) {
                chimin.push_back(chisq2[i][j]);
                umin.push_back(uMagStore[i][j]);
                gmin.push_back(gMagStore[i][j]);
                rmin.push_back(rMagStore[i][j]);
              }
            }
        }
    }
    cpgenv(tstore.min(),tstore.max(),chimin.min(),chimin.max(),0,0);
    //cpgenv(8000.0,13000.0,0.0,200.0,0,0);
    xlab="T\\deff"; ylab="\\gx\\u2"; plab="";
    cpglab(xlab.c_str(),ylab.c_str(),plab.c_str());
    cpgline(tstore.size(),tstore.ptr(),chimin.ptr());
    // spline interpolate function
    {
        Vec_DP x(tstore.size()), y(tstore.size()), d2y(tstore.size());
        for(int i=0; i< tstore.size(); i++){
            x[i] = tstore[i]; 
            y[i] = chimin[i];
        }
        int minpos = minloc(chimin);
        float uAct = umin[minpos];
        float gAct = gmin[minpos];
        float rAct = rmin[minpos];
        DP yp1, ypn;
        yp1 = 1.0e30; ypn = 1.0e30;
        // spline interpolate Teff
        NR::spline(x,y,yp1,ypn,d2y);
        Subs::Array1D<float>  spline_x(5000),spline_y(5000); 
        for(int i=0; i<5000; i++){
            DP xtmp, ytmp;
            xtmp = tstore.min() + (tstore.max()-tstore.min())*float(i)/float(5000.0);
            NR::splint(x,y,d2y,xtmp,ytmp);
            spline_x[i] = xtmp;
            spline_y[i] = ytmp;
        }

        cpgsci(3);
        cpgline(spline_x.size(),spline_x.ptr(),spline_y.ptr());
        cpgsci(1);
        // minimum of chisq interpolation
        float min = spline_y.min();
        minpos = minloc(spline_y);
        // find errors in parameter by looking for delta chisq = 1 in increasing parameter direction
        float check=min; int loop=minpos;
        while(check <= min+1.0){
            loop++;
            check = spline_y[loop];
        }
        // linearly interpolate between last two values in order to find exact point
        float stretch = (min+1-spline_y[loop-1])/(spline_y[loop] - spline_y[loop-1]) ;
        float plusSigma = spline_x[loop-1] 
            + stretch*(spline_x[loop]-spline_x[loop-1]);
        plusSigma -= spline_x[minpos];
        // find lower error on parameter
        loop = minpos; check=min;
        while(check <= min+1.0){
            loop--;
            check = spline_y[loop];
        }
        // linearly interpolate between last two values in order to find exact point
        stretch = (min+1-spline_y[loop])/(spline_y[loop+1] - spline_y[loop]) ;
        float minusSigma = spline_x[loop] 
            + stretch*(spline_x[loop+1]-spline_x[loop]);
        minusSigma = spline_x[minpos]-minusSigma;
        cout << endl<< endl << "Teff  = " << spline_x[minpos] << 
            " + " << plusSigma << " - " << minusSigma << endl;

        float uDist = 10.0*pow(10.0,(up-uAct)/5.0);
        float gDist = 10.0*pow(10.0,(gp-gAct)/5.0);
        float rDist = 10.0*pow(10.0,(rp-rAct)/5.0);
        cout << "u,g,r distances are: " << uDist
            << " " << gDist << " " << rDist << endl;
    } 

       
    if(fix=="n"){  
      // find error in logg
      //loop over each logg
      chimin.clear();
      for(int i=0; i<nlogg; i++){
	// find the teff which gives the minimum chisq at each logg
	// and record which chisq it gives for each logg
	Subs::Array1D<float> store(nteff);
	for (int j=0; j<nteff; j++){
	  store[j] = chisq2[j][i];
	}
	chimin.push_back(store.min());
      }
      //cpgenv(tstore.min(),tstore.max(),chimin.min(),chimin.max(),0,0);
      cpgenv(7.0,9.0,0.0,100.0,0,0);
      xlab="log g";
      cpglab(xlab.c_str(),ylab.c_str(),plab.c_str());
      cpgline(gstore.size(),gstore.ptr(),chimin.ptr());
      // spline interpolate function
      {
	DP yp1, ypn;
	yp1 = 1.0e30; ypn = 1.0e30;
	Vec_DP x(gstore.size()), y(gstore.size()), d2y(gstore.size());
	for(int i=0; i< gstore.size(); i++){x[i] = gstore[i]; y[i]=chimin[i];}
	NR::spline(x,y,yp1,ypn,d2y);
	Subs::Array1D<float>  spline_x(1000),spline_y(1000); 
	for(int i=0; i<1000; i++){
	  DP xtmp, ytmp;
	  xtmp = gstore.min() + (gstore.max()-gstore.min())*float(i)/float(1000.0);
	  NR::splint(x,y,d2y,xtmp,ytmp);
	  spline_x[i] = xtmp;
	  spline_y[i] = ytmp;
	}
	cpgsci(3);
	cpgline(spline_x.size(),spline_x.ptr(),spline_y.ptr());
	cpgsci(1); 
	// minimum of chisq interpolation
	float min = spline_y.min();
	int minpos;
	minpos = minloc(spline_y);
	// find errors in parameter by looking for delta chisq = 1 in increasing parameter direction
	float check=min; int loop=minpos;
	while(check <= min+1.0){
	  loop++;
	  check = spline_y[loop];
	}
	// linearly interpolate between last two values in order to find exact point
	float stretch = (min+1-spline_y[loop-1])/(spline_y[loop] - spline_y[loop-1]) ;
	float plusSigma = spline_x[loop-1] 
	  + stretch*(spline_x[loop]-spline_x[loop-1]);
	plusSigma -= spline_x[minpos];
	// find lower error on parameter
	loop = minpos; check=min;
	while(check <= min+1.0){
	  loop--;
	  check = spline_y[loop];
	}
	// linearly interpolate between last two values in order to find exact point
	stretch = (min+1-spline_y[loop])/(spline_y[loop+1] - spline_y[loop]) ;
	float minusSigma = spline_x[loop] 
	  + stretch*(spline_x[loop+1]-spline_x[loop]);
	minusSigma = spline_x[minpos]-minusSigma;
	cout  << "log g  = " << spline_x[minpos] << 
	  " + " << plusSigma << " - " << minusSigma << endl;
      }
    }
    Subs::Plot plot2;
    plot.close();
    plot2.open("mr.ps/cps"); plot2.focus();
    cpgenv(0.1,1.2,0.006,0.04,0,0);
    //cpgenv(0.85,1.1,0.0065,0.013,0,0);
    // plot radius constraints
    cpgsfs(1); cpgsci(5); 
    cpgrect(0.0,1.4,0.00871,0.00819);
    cpgsfs(2); cpgsci(1);
    // plot log g and Teff constraints
    
    // draw mr plot at constant Log G
    cpgsls(2);
    cpgline(m70.size(),m70.ptr(),r70.ptr());
    cpgline(m75.size(),m75.ptr(),r75.ptr());
    cpgline(m80.size(),m80.ptr(),r80.ptr());
    cpgline(m85.size(),m85.ptr(),r85.ptr());
    cpgline(m90.size(),m90.ptr(),r90.ptr());
    // draw mr plot at constant Teff
    cpgsci(2);
    Subs::Array1D<float> xstore,ystore;
    for(int i=0;i<m70.size();i++){
      if(teff1[i] == 10000.0 ){
	xstore.push_back(m70[i]); xstore.push_back(m75[i]); 
	xstore.push_back(m80[i]); xstore.push_back(m85[i]);    
	xstore.push_back(m90[i]);
	ystore.push_back(r70[i]); ystore.push_back(r75[i]); 
	ystore.push_back(r80[i]); ystore.push_back(r85[i]);    
	ystore.push_back(r90[i]);
	cpgline(xstore.size(),xstore,ystore);
	xstore.clear(); ystore.clear();
      }
    }
    xstore.clear(); ystore.clear();
    for(int i=0; i<100; i++){
      xstore.push_back(0.1 + 1.3*i/99);
      ystore.push_back(eggleton_rad(xstore[i]));
    }
    cpgsci(1);
    cpgsls(1);
    cpgline(xstore.size(),xstore,ystore);

    xlab = "Mass/M\\d\\(2281)";
    ylab = "Radius/R\\d\\(2281)";
    plab = "";
    cpgsci(1);
    cpglab(xlab.c_str(),ylab.c_str(),plab.c_str());
    Subs::Array1D<float> mh, rh;
    hamada("Hamada/mr_hamada.dat",mh,rh);
    rh = rh/100.0;
    cpgsci(3);
    cpgline(mh.size(),mh,rh);
    cpgsci(1);

    // plot Matt Woods models
    cpgsci(6);
    Subs::Array1D<float> tmpT, tmpR, rwood, mwood;
    float tgoal = 10000.0;
    for(int i=40; i<=100; i+=10){
      string file = "Wood95/co040.2.04dt1";
      stringstream append;
      append << setw(3) << setfill('0') << i;
      file.replace(9,3,append.str());
      read_wood(file,tmpT,tmpR);
      rwood.push_back(get_R(tmpT,tmpR,tgoal));
      mwood.push_back(float(i)/100.0);
      tmpT.clear(); tmpR.clear();
    }
    // spline interpolate function
    Subs::Array1D<float> spline_x, spline_y;
    spline(mwood,rwood,2000,spline_x,spline_y);
    cpgline(spline_x.size(),spline_x.ptr(),spline_y.ptr());
    cpgsci(1);

    plot.close();
    
  }
  catch(const string& e){
    
    cerr << "\nError occured inside hist" << endl;
    cerr << e << endl;
    return 1;
  }


  return 0;

}
