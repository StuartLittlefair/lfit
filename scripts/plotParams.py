from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import tempfile
import sys
import os

def readPars(fileName,xcol,ycol,zcol):

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
    z = data[:,zcol-1]
    return (x,y,z)

def plot(x,y,z):
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x,y,z)

xcol = int(sys.argv[1]) or int(raw_input('Give column for x'))
ycol = int(sys.argv[2]) or int(raw_input('Give column for y'))
zcol = int(sys.argv[3]) or int(raw_input('Give column for z'))
fileName = sys.argv[4] or raw_input('Give mcmc Log filename')
print fileName
x,y,z = readPars(fileName,xcol,ycol,zcol)

plot(x,y,z)
plt.show()
