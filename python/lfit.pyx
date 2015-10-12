# cython: embedsignature=True
# cython: boundscheck=False
# distutils: language = c++
# distutils: sources = [./src/WhiteDwarf.cc, ./src/Disc.cc, ./src/BrightSpot.cc, ./src/Donor.cc, ./src/finddeg.cc, ./src/Point.cc]
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
        double calcFlux(double,double) except+
        double calcFlux(double,double,double) except+
        void setup_grid(double) except+
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
        cdef unsigned int n = phi.shape[0]
        cdef unsigned int i
        # setup grid
        self.thisptr.setup_grid(incl)
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
        void spotPos(double,double) except+# calculates x and y location of spot
        double calcFlux(double,double,double) except+
        double calcFlux(double,double,double,double) except+
        double getTangent() except+
        void setup_grid(double) except+

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
        self.thisptr.setup_grid(incl)
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
        void tweak(double, double, double, double) except+
        void setup_grid(double) except+
        double calcFlux(double, double, double) except+
        double calcFlux(double, double, double, double) except+
        
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
        self.thisptr.setup_grid(incl)
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
        double calcFlux(double, double, double) except+
        double calcFlux(double, double, double, double) except+
        void tweak(double, double) except+
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
    '''A Wrapper object that holds a disc, donor, bright spot and white dwarf
       and provides convenient routines for calculating the total flux.
       
       Access is provided to the underlying components for advanced use, but
       most users will only ever need to use the calcFlux method, and access
       the ywd, yd, ys and yrs properties which, when calculated provide arrays
       of the white dwarf, disc, bright spot, and donor star fluxes respectively''' 
    def __init__(self,pars,nel_disc=1000,nel_donor=400):
        '''initialiser for CV object. The parameters argument is a tuple, array or
        list which contains either 14 parameters, or 18 parameters for more complicated
        bright spot models. 
        
        The bright spot is modelled as a linear strip at an angle to the line of centres.
        A fraction of the bright spot strip radiates isotropically, whilst the remainder
        is beamed normal to the surface (in the simple model). In the more complex model
        the bright spot can be made to decay in brightness along it's length at different
        rates (by the two exponent parameters), and beam in a direction other than the normal
        to the direction (using the tilt and yaw parameters
        
        The CV parameters are (in order):
        wdFlux -  white dwarf flux at maximum light
        dFlux  -  disc flux at maximum light
        sFlux  -  bright spot flux at maximum light
        rsFlux -  donor flux at maximum light
        q      -  mass ratio
        dphi   -  full width of white dwarf at mid ingress/egress
        rdisc  -  radius of accretion disc (scaled by distance to inner lagrangian point XL1)
        ulimb  -  linear limb darkening parameter for white dwarf
        rwd    -  white dwarf radius (scaled to XL1)
        scale  -  bright spot scale (scaled to XL1)
        az     -  the azimuth of the bright spot strip (w.r.t to line of centres between stars)
        fis    -  the fraction of the bright spot's flux which radiates isotropically
        dexp   -  the exponent which governs how the brightness of the disc falls off with radius
        phi0   -  a phase offset
        
        the next four parameters are only used for complex bright spot models
        exp1, exp2, tilt, yaw. Their use is described above.
        
        The accretion disc and donor are broken into tiles covering their surface. You can
        override the defaults for these tiles by setting the nel_disc or nel_donor arguments.
        This can increase numerical accuracy at the expense of computing time'''
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
        self.ywd = None
        self.yd = None
        self.ys = None
        self.yrs = None
        self.computed = False
        
    def calcFlux(self,pars,phi,width=None):
        '''tweaks the parameters and calculates the flux from the CV as a whole, and from the
           components of the CV.
           
           the pars list is as described for creation of a CV, and you can switch between
           simple and complex bright spots on the fly just by providing different numbers
           of parameters.
           
           the flux at the phases given in phi is calculated and returned. If the optional
           width argument is provided, the flux is calculated in a bin of this width and the 
           average flux in this bin is returned'''
        assert (len(pars) == 18) or (len(pars) == 14)
        wdFlux,dFlux,sFlux,rsFlux,q,dphi,rdisc,ulimb,rwd,scale,az,fis,dexp,phi0 = pars[0:14]
        self.complex = False
        if len(pars) > 14:
            exp1, exp2, tilt, yaw = pars[14:]
            self.complex = True
            
        xl1 = roche.xl1(q)
        inc = roche.findi(q,dphi)
        if inc < 0:
            raise Exception('invalid combination of q and dphi: %f %f' % (q, dphi))   
            
        bsValid = True
        # check for valid BS parameters
        slop = 80.0
        rd_a = rdisc*xl1
        try:
            x,y,vx,vy = roche.bspot(q,rd_a)
            alpha = np.degrees(np.arctan2(y,x))
            # alpha is between -90 and 90. if negative spot lags disc ie alpha > 90
            if alpha < 0: alpha = 90-alpha
            tangent = alpha + 90 # disc tangent    
            
            # azimuth must be between 0 and 178, and less than 80 degrees from 
            # disc tangent
            if (az < 0) or (az > 178) or (np.fabs(tangent-az) > slop):
                raise Exception('invalid BS azimuth: %f' % az)
        except ValueError as e:
            # if roche.bspot raises error, we didn't hit disc
            raise ValueError('Gas stream does not hit disc for q, rw = %f, %f' % (q,rdisc))

        self.wd.tweak(rwd/xl1,ulimb)
        self.disc.tweak(q,rwd/xl1,rdisc,dexp)
        self.donor.tweak(q)
        if self.complex:
            self.spot.tweak(q,rdisc,az,fis,scale,self.complex,exp1,exp2,tilt,yaw)
        else:
            self.spot.tweak(q,rdisc,az,fis,scale,self.complex)
        self.ywd = wdFlux*self.wd.calcFlux(q,inc,phi-phi0,width)
        self.yd  = dFlux*self.disc.calcFlux(q,inc,phi-phi0,width)
        self.ys  = sFlux*self.spot.calcFlux(q,inc,phi-phi0,width)
        self.yrs = rsFlux*self.donor.calcFlux(q,inc,phi-phi0,width) 
        self.computed = True
        
        return self.ywd + self.yd + self.ys + self.yrs                       