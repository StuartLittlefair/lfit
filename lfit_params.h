#include "trm/subs.h"
#include "trm/array1d.h"
// find inclination from mass ratio and phase width of eclipse
double finddeg(const double& q, const double& dphi);

namespace LFIT{
	// another convenience structure to hold the parameters that make up a model
class Params{
	
public:
    double q,dphi,rd,rwd,ulimb,bsScale,bsAz,bsFrac,bsExp1,bsExp2,bsTilt,bsYaw,dExp,incl,
      dFlux,bFlux,wFlux,rFlux,phi0;
	double qErr,dphiErr,rdErr,rwdErr,ulimbErr,bsScaleErr,bsAzErr,
		bsFracErr,bsExp1Err,bsExp2Err,bsTiltErr,bsYawErr,dExpErr,phi0Err;
	bool qvar,dpvar,rdvar,rwdvar,ulvar,scvar,azvar,fracvar,
		exp1var,exp2var,tiltvar,yawvar,expvar,phivar,errSet,complexSpot;
	double qrange,dprange,rdrange,rwdrange,ulrange,scrange,azrange,
		fracrange,exp1range,exp2range,tiltrange,yawrange,exprange,phirange;

	//default constructor
	Params(){
		qErr = -1.0; dphiErr = -1.0; rdErr = -1.0; rwdErr = -1.0; ulimbErr = -1.0; bsScaleErr = -1.0;
		bsAzErr = -1.0; bsFracErr = -1.0; bsExp1Err = -1.0; bsExp2Err = -1.0;
		bsTiltErr = -1.0; bsYawErr = -1.0; dExpErr = -1.0; phi0Err = -1.0; errSet = false; complexSpot=false;
		q = 0.063;
		dphi=0.04054;
		rd= 0.433;
		rwd=0.02246;
		ulimb=0.346;
		bsScale = 0.039;
		bsAz = 157.0;
		bsFrac = 0.2;
		bsExp1 = 2.0;
		bsExp2 = 1.0;
		bsTilt = 90.0;
		bsYaw  = 1.0;
		dExp = 1.0;
		incl = finddeg(q,dphi);
		phi0 = -15.0e-5;

		qvar=true;dpvar=true;rdvar=true;rwdvar=true;ulvar=true;scvar=true;
		azvar=true;fracvar=true;exp1var=true;exp2var=true;tiltvar=true;
		yawvar=true;expvar=true;phivar=true;
		errSet=false;
		setParamRanges();

		dFlux=0.008;
		bFlux=0.07;
		wFlux=0.1;
		rFlux=0.005;
	}
	
	// destructor
	~Params(){
		//std::cout << "A params objects has gone out of scope" << std::endl;
	}
	
	void setParamRanges(){
		qrange=q/10.;
		dprange=dphi/10.;
		rdrange=rd/5.0;
		rwdrange=rwd/5.0;
		ulrange=0.1;
		scrange=bsScale/5.;
		azrange=10.0;
		fracrange=0.1;
		exp1range=0.05;
		exp2range=0.005;
		tiltrange=2.0;
		yawrange=2.0;
		exprange=0.2;
		phirange=5.0e-4;
	}
	
	Params& operator= (const Params& obj){
		if(this != &obj){
			bool newIncl = (q != obj.q || dphi != obj.dphi);
			if(newIncl) incl = finddeg(obj.q, obj.dphi);
			q = obj.q; qvar= obj.qvar;
			dphi = obj.dphi; dpvar= obj.dpvar;
			rd = obj.rd; rdvar= obj.rdvar;
			rwd = obj.rwd; rwdvar= obj.rwdvar;
			ulimb = obj.ulimb; ulvar= obj.ulvar;
			bsScale = obj.bsScale; scvar= obj.scvar;
			bsAz = obj.bsAz; azvar= obj.azvar;
			bsFrac = obj.bsFrac; fracvar= obj.fracvar;
			bsExp1 = obj.bsExp1; exp1var = obj.exp1var;
			bsExp2 = obj.bsExp2; exp2var = obj.exp2var;
			bsTilt = obj.bsTilt; tiltvar = obj.tiltvar;
			bsYaw = obj.bsYaw; yawvar = obj.yawvar;
			dExp = obj.dExp; expvar= obj.expvar;
			phi0 = obj.phi0; phivar= obj.phivar;
			errSet = obj.errSet;
			if(obj.errSet){
				qErr = obj.qErr; dphiErr = obj.dphiErr; rdErr = obj.rdErr;
				rwdErr = obj.rwdErr; ulimbErr = obj.ulimbErr; bsScaleErr = obj.bsScaleErr;
				bsAzErr = obj.bsAzErr; bsFracErr = obj.bsFracErr; 
				bsExp1Err = obj.bsExp1Err; bsExp2Err = obj.bsExp2Err;
				bsTiltErr = obj.bsTiltErr; bsYawErr = obj.bsYawErr;
				dExpErr = obj.dExpErr; phi0Err = obj.phi0Err; 
			}
			setParamRanges();
		}
		return *this;
	}
	
	bool operator== (const Params& obj) const{
		if(this == &obj){return true;}
		if(q != obj.q or dphi != obj.dphi or
			rd != obj.rd or rwd != obj.rwd or
			ulimb != obj.ulimb or bsScale != obj.bsScale or
			bsAz != obj.bsAz or bsFrac != obj.bsFrac or
			bsExp1 != obj.bsExp1 or bsExp2 != obj.bsExp2 or
			bsTilt != obj.bsTilt or bsYaw != obj.bsYaw or
		dExp != obj.dExp or phi0 != obj.phi0){
			return false;
		}else{
			return true;
		}
	}
	
	bool operator!= (const Params& obj) const{
		if(this == &obj){return false;}
		if(q != obj.q or dphi != obj.dphi or
			rd != obj.rd or rwd != obj.rwd or
			ulimb != obj.ulimb or bsScale != obj.bsScale or
			bsAz != obj.bsAz or bsFrac != obj.bsFrac or
			bsExp1 != obj.bsExp1 or bsExp2 != obj.bsExp2 or
			bsTilt != obj.bsTilt or bsYaw != obj.bsYaw or
		dExp != obj.dExp or phi0 != obj.phi0){
			return true;
		}else{
			return false;
		}
	}
	
	int nvary() const{
		int nout = 0;
		if(qvar) nout++;
		if(dpvar) nout++;
		if(rdvar) nout++;
		if(rwdvar) nout++;
		if(ulvar) nout++;
		if(scvar) nout++;
		if(azvar) nout++;
		if(fracvar) nout++;
		if(exp1var) nout++;
		if(exp2var) nout++;
		if(tiltvar) nout++;
		if(yawvar) nout++;
		if(expvar) nout++;
		if(phivar) nout++;
		return nout;
	}
	
	Subs::Array1D<double> get_variables() const{
		Subs::Array1D<double> outArr;
		if(qvar) outArr.push_back(q);
		if(dpvar) outArr.push_back(dphi);
		if(rdvar) outArr.push_back(rd);
		if(rwdvar) outArr.push_back(rwd);
		if(ulvar) outArr.push_back(ulimb);
		if(scvar) outArr.push_back(bsScale);
		if(azvar) outArr.push_back(bsAz);
		if(fracvar) outArr.push_back(bsFrac);
		if(exp1var) outArr.push_back(bsExp1);
		if(exp2var) outArr.push_back(bsExp2);
		if(tiltvar) outArr.push_back(bsTilt);
		if(yawvar) outArr.push_back(bsYaw);
		if(expvar) outArr.push_back(dExp);
		if(phivar) outArr.push_back(phi0);
		return outArr;
	}
    
	Subs::Array1D<double> get_errors() const{
		Subs::Array1D<double> outArr;
		if(qvar) outArr.push_back(qErr);
		if(dpvar) outArr.push_back(dphiErr);
		if(rdvar) outArr.push_back(rdErr);
		if(rwdvar) outArr.push_back(rwdErr);
		if(ulvar) outArr.push_back(ulimbErr);
		if(scvar) outArr.push_back(bsScaleErr);
		if(azvar) outArr.push_back(bsAzErr);
		if(fracvar) outArr.push_back(bsFracErr);
		if(exp1var) outArr.push_back(bsExp1Err);
		if(exp2var) outArr.push_back(bsExp2Err);
		if(tiltvar) outArr.push_back(bsTiltErr);
		if(yawvar) outArr.push_back(bsYawErr);
		if(expvar) outArr.push_back(dExpErr);
		if(phivar) outArr.push_back(phi0Err);
		return outArr;
	}
	
	Subs::Array1D<double> get_ranges() const{
		Subs::Array1D<double> outArr;
		if(qvar) outArr.push_back(qrange);
		if(dpvar) outArr.push_back(dprange);
		if(rdvar) outArr.push_back(rdrange);
		if(rwdvar) outArr.push_back(rwdrange);
		if(ulvar) outArr.push_back(ulrange);
		if(scvar) outArr.push_back(scrange);
		if(azvar) outArr.push_back(azrange);
		if(fracvar) outArr.push_back(fracrange);
		if(exp1var) outArr.push_back(exp1range);
		if(exp2var) outArr.push_back(exp2range);
		if(tiltvar) outArr.push_back(tiltrange);
		if(yawvar) outArr.push_back(yawrange);
		if(expvar) outArr.push_back(exprange);
		if(phivar) outArr.push_back(phirange);
		return outArr;
	}

	std::vector< std::pair<double,double> > get_limits(const float& discTangent) const{
		std::vector< std::pair<double,double> >outArr;
		if(qvar) outArr.push_back(std::pair<double,double>(0.0,3.5));
		if(dpvar) outArr.push_back(std::pair<double,double>(0.0,0.12));
		if(rdvar) outArr.push_back(std::pair<double,double>(0.0,1.0));
		if(rwdvar) outArr.push_back(std::pair<double,double>(0.0,0.1));
		if(ulvar) outArr.push_back(std::pair<double,double>(0.0,1.0));
		if(scvar) outArr.push_back(std::pair<double,double>(rwd/3.,0.5));
		float azSlop = 40.0; // allow a slack of 40 degrees around tangent to disc edge at impact
		float centre = discTangent*180/3.14159;
		if(azvar) outArr.push_back(std::pair<double,double>(centre-azSlop,centre+azSlop));
		if(fracvar) outArr.push_back(std::pair<double,double>(0.0,1.0));
		if(exp1var) outArr.push_back(std::pair<double,double>(0.01,4.0));
		if(exp2var) outArr.push_back(std::pair<double,double>(0.01,3.0));
		if(tiltvar) outArr.push_back(std::pair<double,double>(0.0,170.0));
		if(yawvar) outArr.push_back(std::pair<double,double>(-90.0,90.0));			
		if(expvar) outArr.push_back(std::pair<double,double>(0.0,1.5));
		if(phivar) outArr.push_back(std::pair<double,double>(-0.3,0.3));
		return outArr;
	}
	
	std::string toString() const{
		
		std::string line1 = Subs::str(q) + "  " + Subs::str(dphi) + "  " + Subs::str(rd) +
			"  " + Subs::str(rwd) + "  " + Subs::str(ulimb) + "  " + Subs::str(bsScale) +
			"  " + Subs::str(bsAz) + "  " + Subs::str(bsFrac) + "  " +
			Subs::str(bsExp1) + "  " + Subs::str(bsExp2) + "  " + Subs::str(bsTilt) +
			"  " + Subs::str(bsYaw) + "  " + Subs::str(dExp) +
			"  " + Subs::str(incl) + "  " + Subs::str(phi0);
		std::string line2 = Subs::str(qrange) + "  " + Subs::str(dprange) + "  " + Subs::str(rdrange) +
			"  " + Subs::str(rwdrange) + "  " + Subs::str(ulrange) + "  " + Subs::str(scrange) +
			"  " + Subs::str(azrange) + "  " + Subs::str(fracrange) + "  " + Subs::str(exprange) +
			"  " + "xxxxx" + "  " + Subs::str(phirange) + "\n";
		std::string line3 = Subs::str(qvar) + "  " + Subs::str(dpvar) + "  " + Subs::str(rdvar) +
			"  " + Subs::str(rwdvar) + "  " + Subs::str(ulvar) + "  " + Subs::str(scvar) +
			"  " + Subs::str(azvar) + "  " + Subs::str(fracvar) + "  " + Subs::str(expvar) +
			"  " + "xxxxx" + "  " + Subs::str(phivar);
		
		return line1;
	}
	
	void jitter(double percent){
		percent /= 100.0;
		// jitter parameters around current value by percentage
		int seed = (int)time(NULL); 
		
		double qtest, dphitest, incltest;
		if(qvar) {
			qtest = q + percent*q*Subs::gauss1(seed);
		} else {
			qtest = q;
		}
		if(dpvar) {
			dphitest = dphi + percent*dphi*Subs::gauss1(seed);
		} else{
			dphitest = dphi;
		}
		// recalc inclination
		if (qvar || dpvar) incltest = finddeg(qtest,dphitest);

		// no errors, let's update q, dphi and inclination
		if(qvar) q = qtest;
		if(dpvar) dphi = dphitest;
		if(qvar || dpvar) incl = incltest;
		
		if(rdvar) rd += percent*rd*Subs::gauss1(seed);
		if(rwdvar) rwd += percent*rwd*Subs::gauss1(seed);
		if(ulvar) ulimb += percent*ulimb*Subs::gauss1(seed);
		if(scvar) bsScale += percent*bsScale*Subs::gauss1(seed);
		if(azvar) bsAz += percent*bsAz*Subs::gauss1(seed);
		if(fracvar) bsFrac += percent*bsFrac*Subs::gauss1(seed);
		if(exp1var) bsFrac += percent*bsExp1*Subs::gauss1(seed);
		if(exp2var) bsFrac += percent*bsExp2*Subs::gauss1(seed);
		if(tiltvar) bsFrac += percent*bsTilt*Subs::gauss1(seed);
		if(yawvar) bsFrac += percent*bsYaw*Subs::gauss1(seed);
		if(expvar) dExp += percent*dExp*Subs::gauss1(seed);
		if(phivar) phi0 += percent*phi0*Subs::gauss1(seed);
	}
	
	void updateWith(Subs::Buffer1D<double> inArr){
		if (nvary() != inArr.size()) {
			std::string err = "params.updateWith: array size doesn't match num. variable params";
			err += "\n   " + Subs::str(inArr.size()) +  " v.s. " + Subs::str(nvary());
			throw err;
		}
		int np=0;
		
		double qTmp, dphiTmp, inclTmp;
		if(qvar) {
			qTmp = inArr[np++];
		}else{
			qTmp = q;
		}
		if(dpvar){
			dphiTmp = inArr[np++];
		}else{
			dphiTmp = dphi;
		}
		// recalc inclination
		if (qvar || dpvar) {
			inclTmp = finddeg(qTmp,dphiTmp);
		}else{
			inclTmp = incl;
		}
		
		// if we haven't thrown an error here, then accept q and/or dphi
		q = qTmp;
		dphi = dphiTmp;
		incl = inclTmp;
		
		if(rdvar) rd =inArr[np++];
		if(rwdvar) rwd = inArr[np++];
		if(ulvar) ulimb=inArr[np++];
		if(scvar) bsScale=inArr[np++];
		if(azvar) bsAz=inArr[np++];
		if(fracvar) bsFrac=inArr[np++];
		if(exp1var) bsExp1=inArr[np++];
		if(exp2var) bsExp2=inArr[np++];
		if(tiltvar) bsTilt=inArr[np++];
		if(yawvar) bsYaw=inArr[np++];
		if(expvar) dExp=inArr[np++];
		if(phivar) phi0=inArr[np++];
		
		setParamRanges();
	}

	void updateErrsWith(Subs::Buffer1D<double> inArr){
		if (nvary() != inArr.size()) {
			std::string err = "params.updateWith: array size doesn't match num. variable params";
			err += "\n   " + Subs::str(inArr.size()) +  " v.s. " + Subs::str(nvary());
			throw err;
		}
		errSet = true;
		int np=0;
		if(qvar) qErr = inArr[np++];
		if(dpvar) dphiErr = inArr[np++];
		if(rdvar) rdErr =inArr[np++];
		if(rwdvar) rwdErr = inArr[np++];
		if(ulvar) ulimbErr=inArr[np++];
		if(scvar) bsScaleErr=inArr[np++];
		if(azvar) bsAzErr=inArr[np++];
		if(fracvar) bsFracErr=inArr[np++];
		if(exp1var) bsExp1Err=inArr[np++];
		if(exp2var) bsExp2Err=inArr[np++];
		if(tiltvar) bsTiltErr=inArr[np++];
		if(yawvar) bsYawErr=inArr[np++];
		if(expvar) dExpErr=inArr[np++];
		if(phivar) phi0Err=inArr[np++];
	}
	
	bool hasErrors() const {return errSet;}
	
	/*
     mass ratio q = m1/m2
     dphi         = eclipse half width
     rwd          = radius of white dwarf in units of L1
     ulimb        = white dwarf limb darkening parameter
     rd           = disc radius in units of L1
     bsScale      = bright spot scale in units of L1
     bsAz         = azimuth of bright spot
     bsFrac        = iotropic fraction of bright spot
     dExp         = disc brightness exponent (flux propto r^-dExp)
     phi0         = phase offset
    */
};

}
