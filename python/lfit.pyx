# distutils: language = c++
# distutils: sources = [../src/WhiteDwarf.cc, ../src/Disc.cc, ../src/BrightSpot.cc, ../src/Donor.cc, ../src/finddeg.cc]
cimport numpy as np
import numpy as np
from libcpp cimport bool
from cython.operator cimport dereference as deref
from trm import roche

############## Donor #######################   
cdef extern from "Donor.h" namespace "LFIT":
    cdef cppclass Donor:
        Donor(double, int)
        void tweak(double)
        double calcFlux(double,double)
        double calcFlux(double,double,double)
        double get_q()
        void set_q(double)
        int get_size()
        void set_size(int)
        
def rebuild_PyDonor(q,size):
    return PyDonor(q,size)
cdef class PyDonor:
    cdef Donor *thisptr
    def __cinit__(self,double q, int size):
        self.thisptr = new Donor(q,size)
    property q:
        def __get__(self): return self.thisptr.get_q()
        def __set__(self,q): self.thisptr.set_q(q)
    property nel:
        def __get__(self): return self.thisptr.get_size()
        def __set__(self,nel): self.thisptr.set_size(nel)      
    def __dealloc__(self):
        del self.thisptr
    def __reduce__(self):
        return (rebuild_PyDonor, (self.q, self.nel))        
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

############## BS #######################                       
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

def rebuild_PySpot(q,rd,az,frac,scale,complex=False,exp1=2.0,exp2=1.0,tilt=90.0,yaw=1.0):
    return PySpot(q,rd,az,frac,scale,complex,exp1,exp2,tilt,yaw)
cdef class PySpot:
    cdef BrightSpot *thisptr
    cdef Params *parptr
    def __cinit__(self, double q, double rd, double az, double frac, double scale, bool complex=False, double exp1=2.0, double exp2=1.0, double tilt=90.0, double yaw=1.0):
        self.parptr = new Params()
        self.parptr.q = q
        self.parptr.rd = rd
        self.parptr.bsAz = az
        self.parptr.bsFrac = frac
        self.parptr.bsScale = scale
        self.parptr.bsExp1 = exp1
        self.parptr.bsExp2 = exp2
        self.parptr.complexSpot = complex
        self.parptr.bsTilt = tilt
        self.parptr.bsYaw = yaw
        self.thisptr = new BrightSpot(deref(self.parptr))  
            
    def tweak(self, double q, double rd, double az, double frac, double scale, bool complex=False,  double exp1=2.0, double exp2=1.0, double tilt=90.0, double yaw=1.0):
        self.parptr.q = q
        self.parptr.rd = rd
        self.parptr.bsAz = az
        self.parptr.bsFrac = frac
        self.parptr.bsScale = scale
        self.parptr.bsExp1 = exp1
        self.parptr.bsExp2 = exp2
        self.parptr.complexSpot = complex
        self.parptr.bsTilt = tilt
        self.parptr.bsYaw = yaw
        self.thisptr.tweak(deref(self.parptr))   
        
    def __reduce__(self):
        q = self.parptr.q
        rd = self.parptr.rd
        bsAz = self.parptr.bsAz
        bsFrac = self.parptr.bsFrac
        bsScale = self.parptr.bsScale
        bsExp1 = self.parptr.bsExp1
        bsExp2 = self.parptr.bsExp2
        complexSpot = self.parptr.complexSpot
        bsTilt = self.parptr.bsTilt
        bsYaw = self.parptr.bsYaw
        return (rebuild_PySpot, (q,rd,bsAz,bsFrac,bsScale,bsExp1,bsExp2,complexSpot,bsTilt,bsYaw))
                   
    def __dealloc__(self):
        del self.thisptr
        del self.parptr
        
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
        
############## DISC #######################   
cdef extern from "Disc.h" namespace "LFIT":
    cdef cppclass Disc:
        Disc(double, double, double, double, int) except +
        void tweak(double, double, double, double)
        double calcFlux(double, double, double)
        double calcFlux(double, double, double, double)
        
def rebuild_PyDisc(q,rin,rout,exp,nelem=1000):
    return PyDisc(q,rin,rout,exp,nelem)    
cdef class PyDisc:
    cdef Disc *thisptr
    cdef double q, rin, rout, exp
    cdef int nelem
    def __cinit__(self,double q, double rin, double rout, double exp, int nelem=1000):
        self.thisptr = new Disc(q,rin,rout,exp,nelem)
        self.q = q
        self.rin = rin
        self.rout = rout
        self.exp = exp
        self.nelem = nelem
    def __reduce__(self):
        return (rebuild_PyDisc, (self.q,self.rin,self.rout,self.exp,self.nelem))
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
        self.q = q
        self.rin = rin
        self.exp = exp
                 
############## WD #######################   
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

def rebuild_PyWhiteDwarf(radius,ulimb):
    return PyWhiteDwarf(radius,ulimb)                
cdef class PyWhiteDwarf:
    cdef WhiteDwarf *thisptr    # hold a C++ instance which we wrap
    def __cinit__(self, double radius, double ulimb):
        self.thisptr = new WhiteDwarf(radius,ulimb)
        self.radius = radius
        self.ulimb = ulimb
    def __reduce__(self):
        return (rebuild_PyWhiteDwarf, (self.radius, self.ulimb))
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

############## CV #######################   
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
        if self.complex:
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