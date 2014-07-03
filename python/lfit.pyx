# distutils: language = c++
# distutils: sources = ../WhiteDwarf.cc

cimport numpy as np
import numpy as np

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
                
cdef class Whitedwarf:
    cdef WhiteDwarf *thisptr    # hold a C++ instance which we wrap
    def __cinit__(self, double radius, double ulimb):
        self.thisptr = new WhiteDwarf(radius,ulimb)
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
                   
    property radius:
        def __get__(self): return self.thisptr.get_radius()
        def __set__(self,radius): self.thisptr.set_radius(radius)
    property ulimb:
        def __get__(self): return self.thisptr.get_ulimb()
        def __set__(self,ulimb): self.thisptr.set_ulimb(ulimb)
        