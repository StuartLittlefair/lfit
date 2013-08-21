import numpy
import pylab
import sys

fileName = sys.argv[1] or raw_input('Give MCMC scale logfile name')

f = open(fileName)
lines = f.readlines()
x = []
succRate = []
scaleFac = []
for line in lines:
    arr = line.split()
    x.append(arr[2].replace(',',''))
    succRate.append(arr[5].replace(',',''))
    scaleFac.append(arr[10].replace(',',''))
x = numpy.array(x).astype('int')
succRate = numpy.array(succRate).astype('float')
scaleFac = numpy.array(scaleFac).astype('float')

pylab.subplot(2,1,1)
pylab.plot(x,succRate)
pylab.ylabel('Success Rate')
pylab.subplot(2,1,2)
pylab.plot(x,scaleFac)
pylab.ylabel('Scale Factor')
pylab.xlabel('Step Number')
pylab.show()
