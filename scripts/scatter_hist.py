import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import load
import tempfile
import sys
import os
from scipy.optimize import leastsq as lsq

#skewed gaussian fitting function (Rusch & Lelieur, Anal Chem, 1973, 45, 1541)
fitfunc = lambda p, x: p[0]*np.exp( -np.log(2.0) * (np.log(1 + (2.0*p[1]*(x-p[2])/p[3]))/p[1])**2 )
errfunc = lambda p, x, y: y-fitfunc(p,x)

xcol = int(sys.argv[1]) or int(raw_input('Give column for x'))
ycol = int(sys.argv[2]) or int(raw_input('Give column for y'))
fileName = sys.argv[3] or raw_input('Give mcmc Log filename')

f = open(fileName)

lines = np.array(f.readlines())

olines = lines[np.arange(1,len(lines),2)]
ofName = tempfile.mktemp()
of = open(ofName,'w')
for line in olines:
	of.write(line)
	
of.close()
f.close()

data = np.loadtxt(ofName)
os.unlink(ofName)	


x = data[:,xcol-1]
y = data[:,ycol-1]


fig = plt.figure(1, figsize=(8.5,8.5))

from mpl_toolkits.axes_grid import make_axes_locatable

axScatter = plt.subplot(111)
divider = make_axes_locatable(axScatter)

# create a new axes with a height of 1.2 inch above the axScatter
axHistx = divider.new_vertical(1.2, pad=0.1, sharex=axScatter)

# create a new axes with a width of 1.2 inch on the right side of the
# axScatter
axHisty = divider.new_horizontal(1.2, pad=0.1, sharey=axScatter)

fig.add_axes(axHistx)
fig.add_axes(axHisty)


# make some labels invisible
plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
         visible=False)

# the scatter plot:
#axScatter.scatter(x, y)
H, xedges, yedges = np.histogram2d(x, y, bins=(20, 20))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
levels = [0.33*H.max(),0.05*H.max(),0.003*H.max()]
axScatter.imshow(H,extent=extent,aspect='auto',cmap=plt.cm.gray_r,origin='lower')
CS=axScatter.contour(H,levels,extent=extent,origin='lower')
axScatter.clabel(CS, inline=1, fontsize=9)
#axScatter.set_aspect(1.)

# now determine nice limits by hand:

xmax = np.max(x)
xmin = np.min(x)
ymin = np.min(y)
ymax = np.max(y) 


binwidth = (xmax-xmin)/100.
bins = np.arange(xmin, xmax + binwidth, binwidth)
nx, binx, patchx = axHistx.hist(x, bins=bins)
xp = (binx[:-1] + binx[1:])/2
binwidth = (ymax-ymin)/100.
bins = np.arange(ymin, ymax + binwidth, binwidth)
ny, biny, patchy = axHisty.hist(y, bins=bins, orientation='horizontal')
yp = (biny[:-1] + biny[1:])/2

# fit with skewed gaussian and plot
params = np.array([nx.max(),0.01,xp.mean(),np.sqrt(xp.var())])
out = lsq(errfunc,params,args=(xp,nx),full_output=1)
pfinal = out[0]
#print pfinal
fitX = np.linspace(xmin,xmax,200)
fitY = fitfunc(pfinal,fitX)
axHistx.plot(fitX,fitY)

params = np.array([ny.max(),0.01,yp.mean(),np.sqrt(yp.var())])
out = lsq(errfunc,params,args=(yp,ny),full_output=1)
pfinal = out[0]
#print pfinal
fitX = np.linspace(ymin,ymax,200)
fitY = fitfunc(pfinal,fitX)
axHisty.plot(fitX,fitY)

# the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
# thus there is no need to manually adjust the xlim and ylim of these
# axis.

#axHistx.axis["bottom"].major_ticklabels.set_visible(False)
for tl in axHistx.get_xticklabels():
    tl.set_visible(False)
#axHistx.set_yticks([0, 50, 100])

#axHisty.axis["left"].major_ticklabels.set_visible(False)
for tl in axHisty.get_yticklabels():
    tl.set_visible(False)
#axHisty.set_xticks([0, 50, 100])

plt.draw()
plt.show()
#plt.savefig("a.pdf")
