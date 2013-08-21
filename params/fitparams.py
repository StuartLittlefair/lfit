from scipy.optimize import leastsq as lsq
from scipy.special import erf
from scipy.stats import skew
import scipy.stats
import numpy
import ppgplot as p
import warnings
warnings.filterwarnings("ignore")

class Param:
    def __init__(self,shortString,longString,index):
        self.shortString = shortString
        self.longString = longString
        self.index = index

def plotMult(x,parsList,total,label):

    fitList = []
    ymin = 1.0e32
    ymax = -1.0e32
    for par in parsList:
        fitList.append(fitfunc(par,x))
        fitList[-1] /= fitList[-1].sum()
        ymin = min(ymin,fitList[-1].min())
        ymax = max(ymax,fitList[-1].max())
    total   /= total.sum()
    ymin = min(ymin,total.min())
    ymax = max(ymax,total.max())

    p.pgenv(x.min(),x.max(),ymin,ymax)
    colInd = 4
    for fit in fitList:
        p.pgsci(colInd)
        p.pgline(x,fit)
        colInd -= 1
    p.pgsci(1)
    p.pgline(x,total)
    p.pglab(label,'PDF','')

def plot(array,label,params):
    (y,bins) = numpy.histogram(array,bins=50,normed=True)
    x = 0.5*(bins[:-1] + bins[1:])
    y /= float(len(array))
    maxloc = y.argmax()
    yFit = fitfunc(params,x)
    p.pgenv(x.min(),x.max(),y.min(),y.max())
    p.pgline(x,yFit)
    p.pgbin(x,y,True)
    p.pglab(label,'PDF','')

def fitSkewedGaussian(array):
    (y,bins) = numpy.histogram(array,bins=50,normed=True)
    x = 0.5*(bins[:-1] + bins[1:])
    y /= float(len(array))
    maxloc = y.argmax()
    mode = x[maxloc]
    # fit skewed Gaussian
    gamma = skew(array)
    delta = numpy.sqrt(numpy.pi*numpy.abs(gamma)**(2./3.)/2./(numpy.abs(gamma)**(2./3.) + ((4-numpy.pi)/2)**(2./3.)))
    if delta < 1:
        alpha = delta/numpy.sqrt(1-delta**2.0)
    else:
        alpha = 0.99 
    if gamma < 0:
        alpha *= -1
    params = numpy.array([mode,array.var(),alpha,y[maxloc]])
    out = lsq(errfunc,params,args=(x,y),full_output=1)
    pfinal = out[0]
    return pfinal


def percentile(x,y,perc):
    cdf = numpy.cumsum(y)
    cdf /= cdf.max()
    loc = numpy.abs(cdf-perc).argmin()
    x1 = x[loc-1]
    x2 = x[loc+1]
    y1 = cdf[loc-1]
    y2 = cdf[loc+1]
    return x2 - ( (y2-perc)*(x2-x1)/(y2-y1) )

def getStatsPDF(x,y,label):
    maxloc = y.argmax()
    mode = x[maxloc]
    # get 16th and 84th percentile (defines 1 sigma confidence range) 
    conflim = [percentile(x,y,0.16),percentile(x,y,0.84)]
    print "%s = %.8f + %.8f - %.8f" % (label, mode, conflim[1]-mode, mode-conflim[0])
    
def getStats(array,shortLabel):
    (y,bins) = numpy.histogram(array,bins=50,normed=True)
    x = 0.5*(bins[:-1] + bins[1:])
    y /= float(len(array))
    maxloc = y.argmax()
    mode = x[maxloc]
    # get 16th and 84th percentiles, which represent the upper and lower limits of the 68% confidence interval (1-sigma)
    conflim = [scipy.stats.scoreatpercentile(array,16),scipy.stats.scoreatpercentile(array,84)]
    print "%s = %.8f + %.8f - %.8f" % (shortLabel, mode, conflim[1]-mode, mode-conflim[0])

if __name__ == "__main__":
    
    fitfunc = lambda p, x: p[3]*numpy.exp( -(x-p[0])**2/2.0/p[1] ) * (1+ erf(p[2]*(x-p[0])/numpy.sqrt(p[1]*2)) )
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    p.pgopen('?')
    p.pgsubp(2,5)
    p.pgslw(2)
    p.pgsch(3)

    paramList = [Param('q','Mass Ratio (q)',0),
                 Param('m1','M\\d1\\u (M\\d\\(2281)\\u)',1),
                 Param('r1','R\\d1\\u (R\\d\\(2281)\\u)',2),
                 Param('m2','M\\d2\\u (M\\d\\(2281)\\u)',3),
                 Param('r2','R\\d2\\u (R\\d\\(2281)\\u)',4),
                 Param('i','Inclination',8),
                 Param('a','Binary Separation (R\\d\\(2281)\\u)',5),
                 Param('kw','K\\d1\\u (km s\\u-1\\d)',6),
                 Param('kr','K\\d2\\u (km s\\u-1\\d)',7)]

    while True:
        mode = raw_input('(S)ingle dataset or (M)ultiple datasets? ')
        if mode.upper() == 'M' or mode.upper() == 'S':
            break
        else:
            print "Please answer S or M "
            
            
    if mode.upper() == "S":
        asciiFile = raw_input('Give data file containing parameter samples: ')
        dataIn = numpy.loadtxt(asciiFile)
        for param in paramList:
            array=dataIn[:,param.index]
            plot(array,param.longString,fitSkewedGaussian(array))
            getStats(array,param.shortString)
    else:
        dataList = []
        colours = ['blu','grn','red']
        numSets = 0
        numSets = int(raw_input('How many datasets to combine? '))
        files = []
        for i in range(numSets):
            files.append( raw_input('Give data file containing parameter samples for ' + colours[i] + ' data: ') )

        for i in range(numSets):
            dataList.append(asciidata.open(files[i]))

        for param in paramList:
            parsList = []
            fitsList = []
            minX = 1.0e32
            maxX = -1.0e32
            for i in range(numSets):
                array = numpy.array(dataList[i][param.index].tonumarray(),dtype='float64')
                minX = min(minX,array.min())
                maxX = max(maxX,array.max())
                parsList.append(fitSkewedGaussian(array))
            x = numpy.linspace(minX,maxX,1000)
            if numSets == 2:
                result = fitfunc(parsList[0],x)*fitfunc(parsList[1],x)
            else:
                result = fitfunc(parsList[0],x)*fitfunc(parsList[1],x)*fitfunc(parsList[2],x)
            if numSets == 2:
                plotMult(x,[parsList[0],parsList[1]],result,param.longString)
            else:
                plotMult(x,[parsList[0],parsList[1],parsList[2]],result,param.longString)
            getStatsPDF(x,result,param.shortString)
    p.pgclos()



