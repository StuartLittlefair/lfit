import tempfile
#import ppgplot as p
import pylab as p
import numpy as n
import sys
import os

fileName = sys.argv[1] or raw_input('Give mcmc Log filename')
#dev = sys.argv[2] or raw_input('Give plot device')

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

npars = data.shape[1]-2
nsteps = data.shape[0]

iplot = range(1,npars+1,1)
if (npars == 13):
	labels = [r'$q$',r'$\Delta\phi$',r'$R_d$',r'$R_{wd}$',r'$BS_{Scale}$',r'$BS_{Az}$',r'$F_{is}$',r'$BS_{E1}$',r'$BS_{E2}$',r'$BS_{Tilt}$',r'$BS_{Yaw}$',r'$D_{exp}$',r'$\phi_0$']
if (npars == 11):
	labels = [r'$q$',r'$\Delta\phi$',r'$R_d$',r'$R_{wd}$',r'$BS_{Scale}$',r'$BS_{Az}$',r'$F_{is}$',r'$BS_{E1}$',r'$BS_{Tilt}$',r'$D_{exp}$',r'$\phi_0$']
if (npars == 9):
	labels = [r'$q$',r'$\Delta\phi$',r'$R_d$',r'$R_{wd}$',r'$BS_{Scale}$',r'$BS_{Az}$',r'$F_{is}$',r'$D_{exp}$',r'$\phi_0$']
nCols = 1
plotNum = npars*npars
for i in iplot:
	
	# fancy plot numbering to order subplots correctly
	plotNum = npars*npars - (i-1)*(npars)
	for j in range(npars,i,-1):
		x = data[:,i]
		y = data[:,j]

		#
		p.subplot(npars,npars,plotNum)
		# calculate next subplot number
		plotNum = plotNum-1

		# bin points into 2D histogram (density of points as a function of x and y)
		H, xedges, yedges = n.histogram2d(x, y, bins=(20, 20))
		# find x and y limits of edge of histogram
		extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
		# define 1sigma, 2sigma and 3sigma contours
		levels = [0.33*H.max(),0.05*H.max(),0.003*H.max()]
		# plot greyscale of 2d histogram
		p.imshow(H,extent=extent,aspect='auto',interpolation='nearest',cmap=p.cm.gray_r,origin='lower')
		# plot contours
		p.contour(H,levels,extent=extent,origin='lower')

		# disable numerical labelling
		a = p.gca()
		for tick in a.yaxis.get_major_ticks():
			tick.label1On = False
			tick.label2On = False
		for tick in a.xaxis.get_major_ticks():
			tick.label1On = False
			tick.label2On = False
			
		# label plots if along x or y axis
		#print i, j,labels[i-1], labels[j-1]
		if i == 1:
			p.xlabel(labels[j-1])
		if j == npars:
			rhAxis=p.twinx()
			for tick in a.yaxis.get_major_ticks():
				tick.label1On = False
				tick.label2On = False
			for tick in rhAxis.yaxis.get_major_ticks():
				tick.label1On = False
				tick.label2On = False
			p.ylabel(labels[i-1],rotation='horizontal')
		#	if i==1:
		#		rhAxis.set_xlabel(labels[j-1])



p.show()
