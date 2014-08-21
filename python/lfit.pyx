# distutils: language = c++
# distutils: sources = [../WhiteDwarf.cc, ../Disc.cc, ../BrightSpot.cc, ../Donor.cc, ../finddeg.cc]
cimport numpy as np
import numpy as np
from libcpp cimport bool
from cython.operator cimport dereference as deref
from trm import roche
cdef extern from "Donor.h" namespace "LFIT":
    cdef cppclass Donor:
        Donor(double, int)
        void tweak(double)
        double calcFlux(double,double)
        double calcFlux(double,double,double)
        
cdef class PyDonor:
    cdef Donor *thisptr
    def __cinit__(self,double q, int size):
        self.thisptr = new Donor(q,size)
    def __dealloc__(self):
        del self.thisptr
    def tweak(self,double q):
        self.thisptr.tweak(q)
    def calcFlux(self, double q, double incl, np.ndarray[np.double_t, ndim=1] phi, np.ndarray[np.double_t, ndim=1] width=None):
        cdef unsigned n = phi.shape[0]
        cdef unsigned int i
        phi = np.ascontiguousarray(phi)
        cdef np.ndarray[double, ndim=1] out = np.empty(n, dtype=np.double)
        if width is not None:
            width = np.ascontiguousarray(width)
            for i in range(n):
                out[i] = self.thisptr.calcFlux(phi[i],width[i],incl)
        else:
            for i in range(n):
                out[i] = self.thisptr.calcFlux(phi[i],incl)
        return out   
                    
cdef extern from "BrightSpot.h" namespace "LFIT":
    cdef cppclass Params:
        double q,dphi,rd,rwd,ulimb,bsScale,bsAz,bsFrac,bsExp1,bsExp2,bsTilt,bsYaw,dExp,incl,dFlux,bFlux,wFlux,rFlux,phi0,complexSpot

cdef extern from "BrightSpot.h" namespace "LFIT":
    cdef cppclass BrightSpot:
        BrightSpot(Params)
        void tweak(Params)
        void spotPos(double,double) # calculates x and y location of spot
        double calcFlux(double,double,double)
        double calcFlux(double,double,double,double)
        double getTangent()

cdef class PySpot:
    cdef BrightSpot *thisptr
    def __cinit__(self, double q, double rd, double az, double frac, double scale, bool complex=False, double exp1=2.0, double exp2=1.0, double tilt=90.0, double yaw=1.0):
        cdef Params *parptr = new Params()
        parptr.q = q
        parptr.rd = rd
        parptr.bsAz = az
        parptr.bsFrac = frac
        parptr.bsScale = scale
        parptr.bsExp1 = exp1
        parptr.bsExp2 = exp2
        parptr.complexSpot = complex
        parptr.bsTilt = tilt
        parptr.bsYaw = yaw
        self.thisptr = new BrightSpot(deref(parptr))
        del parptr  
            
    def tweak(self, double q, double rd, double az, double frac, double scale, bool complex=False,  double exp1=2.0, double exp2=1.0, double tilt=90.0, double yaw=1.0):
        cdef Params *parptr = new Params()
        parptr.q = q
        parptr.rd = rd
        parptr.bsAz = az
        parptr.bsFrac = frac
        parptr.bsScale = scale
        parptr.bsExp1 = exp1
        parptr.bsExp2 = exp2
        parptr.complexSpot = complex
        parptr.bsTilt = tilt
        parptr.bsYaw = yaw
        self.thisptr.tweak(deref(parptr))
        del parptr              
    def __dealloc__(self):
        del self.thisptr
    def getTangent(self):
        return self.thisptr.getTangent()
    def calcFlux(self, double q, double incl, np.ndarray[np.double_t, ndim=1] phi, np.ndarray[np.double_t, ndim=1] width=None):
        cdef unsigned n = phi.shape[0]
        cdef unsigned int i
        phi = np.ascontiguousarray(phi)
        cdef np.ndarray[double, ndim=1] out = np.empty(n, dtype=np.double)
        if width is not None:
            width = np.ascontiguousarray(width)
            for i in range(n):
                out[i] = self.thisptr.calcFlux(q,phi[i],width[i],incl)
        else:
            for i in range(n):
                out[i] = self.thisptr.calcFlux(q,phi[i],incl)
        return out
        
cdef extern from "Disc.h" namespace "LFIT":
    cdef cppclass Disc:
        Disc(double, double, double, double, int) except +
        void tweak(double, double, double, double)
        double calcFlux(double, double, double)
        double calcFlux(double, double, double, double)
               
cdef class PyDisc:
    cdef Disc *thisptr
    def __cinit__(self,double q, double rin, double rout, double exp, int nelem=1000):
        self.thisptr = new Disc(q,rin,rout,exp,nelem)
    def __dealloc__(self):
        del self.thisptr
    def calcFlux(self, double q, double incl, np.ndarray[np.double_t, ndim=1] phi, np.ndarray[np.double_t, ndim=1] width=None):
        cdef unsigned n = phi.shape[0]
        cdef unsigned int i
        phi = np.ascontiguousarray(phi)
        cdef np.ndarray[double, ndim=1] out = np.empty(n, dtype=np.double)
        if width is not None:
            width = np.ascontiguousarray(width)
            for i in range(n):
                out[i] = self.thisptr.calcFlux(q,phi[i],width[i],incl)
        else:
            for i in range(n):
                out[i] = self.thisptr.calcFlux(q,phi[i],incl)
        return out
    def tweak(self, double q, double rin, double rout, double exp):
        self.thisptr.tweak(q,rin,rout,exp)
              
cdef extern from "WhiteDwarf.h" namespace "LFIT":
    cdef cppclass WhiteDwarf:
        WhiteDwarf(double, double) except +
        double calcFlux(double, double, double)
        double calcFlux(double, double, double, double)
        void tweak(double, double)
        double get_radius()
        void   set_radius(double)
        double get_ulimb()
        void   set_ulimb(double)
                
cdef class PyWhiteDwarf:
    cdef WhiteDwarf *thisptr    # hold a C++ instance which we wrap
    def __cinit__(self, double radius, double ulimb):
        self.thisptr = new WhiteDwarf(radius,ulimb)
    def __dealloc__(self):
        del self.thisptr
    def tweak(self, double radius, double ulimb):
        self.thisptr.tweak(radius, ulimb)  
    def calcFlux(self, double q, double incl, np.ndarray[np.double_t, ndim=1] phi, np.ndarray[np.double_t, ndim=1] width=None):
        cdef unsigned n = phi.shape[0]
        cdef unsigned int i
        phi = np.ascontiguousarray(phi)
        cdef np.ndarray[double, ndim=1] out = np.empty(n, dtype=np.double)
        if width is not None:
            width = np.ascontiguousarray(width)
            for i in range(n):
                out[i] = self.thisptr.calcFlux(q,phi[i],width[i],incl)
        else:
            for i in range(n):
                out[i] = self.thisptr.calcFlux(q,phi[i],incl)
        return out
               
    property radius:
        def __get__(self): return self.thisptr.get_radius()
        def __set__(self,radius): self.thisptr.set_radius(radius)
    property ulimb:
        def __get__(self): return self.thisptr.get_ulimb()
        def __set__(self,ulimb): self.thisptr.set_ulimb(ulimb)
        
class CV(object):
    def __init__(self,pars,nel_disc=1000,nel_donor=400):
        print len(pars)
        assert (len(pars) == 18) or (len(pars) == 14)
        wdFlux,dFlux,sFlux,rsFlux,q,dphi,rdisc,ulimb,rwd,scale,az,fis,dexp,phi0 = pars[0:14]
        self.complex = False
        if len(pars) > 14:
            exp1, exp2, tilt, yaw = pars[14:]
            self.complex = True
        xl1 = roche.xl1(q)
        incl = roche.findi(q,dphi)
        
        self.wd = PyWhiteDwarf(rwd/xl1,ulimb)
        self.disc = PyDisc(q,rwd/xl1,rdisc,dexp,nel_disc)
        if complex:
            self.spot = PySpot(q,rdisc,az,fis,scale,exp1,exp2,tilt,yaw,self.complex)
        else:
            self.spot = PySpot(q,rdisc,az,fis,scale)
        self.donor = PyDonor(q,nel_donor)
        
    def calcFlux(self,pars,phi,width=None):
        assert (len(pars) == 18) or (len(pars) == 14)
        wdFlux,dFlux,sFlux,rsFlux,q,dphi,rdisc,ulimb,rwd,scale,az,fis,dexp,phi0 = pars[0:14]
        self.complex = False
        if len(pars) > 14:
            exp1, exp2, tilt, yaw = pars[14:]
            self.complex = True
        xl1 = roche.xl1(q)
        inc = roche.findi(q,dphi)
        self.wd.tweak(rwd/xl1,ulimb)
        self.disc.tweak(q,rwd/xl1,rdisc,dexp)
        self.donor.tweak(q)
        if self.complex:
            self.spot.tweak(q,rdisc,az,fis,scale,self.complex,exp1,exp2,tilt,yaw)
        else:
            self.spot.tweak(q,rdisc,az,fis,scale)
        ywd = self.wd.calcFlux(q,inc,phi,width)
        yd  = self.disc.calcFlux(q,inc,phi,width)
        ys  = self.spot.calcFlux(q,inc,phi,width)
        yrs = self.donor.calcFlux(q,inc,phi,width) 
        
        return wdFlux*ywd + dFlux*yd + sFlux*ys + rsFlux*yrs                       