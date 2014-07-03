import matplotlib.pyplot as plt
import numpy as np
from trm import roche
import sys
sys.path.append('/Users/sl/data/SDSS1411/code/')
import whitedwarf as wd
import lfit
import time
import commands
import os
from trm import roche
q = 0.1
inc = 85.9


phi = np.linspace(-0.05,0.05,1000)
width = np.mean(np.diff(phi))*np.ones_like(phi)/2.

xl1 = roche.xl1(q) 
w = lfit.Whitedwarf(0.01/xl1,0.4)
print w.radius, w.ulimb
start = time.clock()
y = w.calcFlux(q,inc,phi,width)
stop = time.clock()
print 'LFIT version took %f' % (stop-start)

w2 = wd.WhiteDwarf(0.01,0.4)
r2 = 1.0-roche.xl1(q)
start = time.clock()
y2 = np.array([w2.instantFlux(q,x,inc,r2) for x in phi])
stop = time.clock()
print 'python version took %f' % (stop-start)

start = time.clock()
commands.getoutput('./test > test.out')
stop = time.clock()
print 'C++ version took %f' % (stop-start)
phi3,y3 = np.loadtxt('test.out').T
os.unlink('test.out')

phi4,_,y4,_,_,_ = np.loadtxt('lroche_data.txt').T
y4 /= y4.max()

plt.plot(phi,y2,'--k')
plt.plot(phi3,y3,'--r')
plt.plot(phi4,y4,'--b',lw=2)
plt.plot(phi,y,'-k')
plt.ylim((0.0,1.1))
plt.show()