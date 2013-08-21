import tempfile
import pylab as p
import numpy as n
import sys
import os

xcol = int(sys.argv[1]) or int(raw_input('Give column for x'))
ycol = int(sys.argv[2]) or int(raw_input('Give column for y'))


def plotPars(col):
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


	x = data[:,xcol]
	y = data[:,ycol]
	#p.plot(x,y,col)
	H, xedges, yedges = n.histogram2d(x, y, bins=(20, 20))
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	levels = [0.33*H.max(),0.05*H.max()]
	p.contour(H,levels,extent=extent,origin='lower')

fileName = sys.argv[3] or raw_input('Give mcmc Log filename for red colour')
plotPars('r.')
fileName = sys.argv[4] or raw_input('Give mcmc Log filename for green colour')
plotPars('g.')
p.show()
