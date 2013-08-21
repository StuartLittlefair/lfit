import pylab as p
import numpy as n
import sys
import tempfile
import os
from math import sqrt

f = open('burnInLog.txt')
lines = n.array(f.readlines())

olines = lines[n.arange(1,len(lines),2)]
ofName = tempfile.mktemp()
of = open(ofName,'w')
for line in olines:
	of.write(line)
	
f.close()
of.close()

data = n.loadtxt(ofName)
os.unlink(ofName)

nparams = int(data[0,0])
for i in range(nparams):
	p.subplot(nparams,1,i+1)
	x = data[:,i+1]
	label = 'parameter ' + str(i+1)
	p.plot(x,'r-',label=label)

p.show()
