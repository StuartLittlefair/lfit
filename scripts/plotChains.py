import pylab as p
import numpy as n
import sys
import tempfile
import os
from math import sqrt

def conflims(vals,chi):
	chiKey  = chi.argsort()
	datakey = vals.argsort()
	chiSort = chi.copy()
	valsSort = vals.copy()
	valsSort.sort()
	chiSort.sort()
	xml     = vals[chiKey[0]]

	xmax    = xml
	xmin    = xml
	dxl     = 0.0
	dxh     = 0.0
	jlo     = 0
	jhi     = 0
	
	for i in range(chi.size):
		xmax = n.max([xmax,vals[chiKey[i]]])
		xmin = n.min([xmin,vals[chiKey[i]]])
		jhi = valsSort.searchsorted(xmax)
		jlo = valsSort.searchsorted(xmin)
		dxl = xml-xmin
		dxh = xmax-xml
		if jhi-jlo >= 0.683*chi.size: break
	return (xml,dxl,dxh)
		
def readLogFile(fname):
	f = open(fname)
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
	return data
	
if __name__ == '__main__':
	fname = raw_input('Give name of MCMC log file: ')
	data = readLogFile(fname)

	if len(sys.argv) > 1:
		paramNumList = []
		for arg in sys.argv[1:]:
			paramNumList.append(int(arg))
	else:
		paramNumList = None


	nparams = int(data[0,0])
	if paramNumList:
		subplot = 1
		for paramNum in paramNumList:
			p.subplot(len(paramNumList),1,subplot)
			subplot += 1
			x = data[:,paramNum]
			entryLength = len(data[0])
			if entryLength > nparams+1:
				chisq = data[:,nparams+1]
				mean,lowlim,uplim = conflims(x,chisq)
				print "Param #",paramNum, ": ", mean , " + ", uplim, " - ", lowlim
			label = 'parameter ' + str(paramNum)
			p.plot(x,'r-',label=label)
	else:
		entryLength = len(data[0])
		if entryLength> nparams+1:
			nplots=nparams+1
		else:	
			nplots=nparams
		for i in range(nparams):
			p.subplot(nplots,1,i+1)
			x = data[:,i+1]
			if entryLength > nparams+1:
				chisq = data[:,nparams+1]
				bestFit,lowlim,uplim = conflims(x,chisq)
				print "Param #",i+1, ": ", x.mean() , " (", bestFit, ")", " + ", uplim, " - ", lowlim
			label = 'parameter ' + str(i+1)
			p.plot(x,'r-',label=label)
		if entryLength> nparams+1:
			p.subplot(nplots,1,nparams+1)
			x = data[:,nparams+1]
			p.plot(x,'r-',label=label)
	p.show()
