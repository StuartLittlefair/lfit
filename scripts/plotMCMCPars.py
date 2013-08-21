import tempfile
import pylab as p
import numpy as n
import sys
import os

xcol = int(sys.argv[1]) or int(raw_input('Give column for x'))
ycol = int(sys.argv[2]) or int(raw_input('Give column for y'))
fileName = sys.argv[3] or raw_input('Give mcmc Log filename')

f = open(fileName)

lines = n.array(f.readlines())

olines = lines[n.arange(1,len(lines),2)]
ofName = tempfile.mktemp()
of = open(ofName,'w')
for line in olines:
	of.write(line)
	
of.close()
f.close()

data = n.loadtxt(ofName)
os.unlink(ofName)	


x = data[:,xcol-1]
y = data[:,ycol-1]
p.plot(x,y,'r.')
p.show()
