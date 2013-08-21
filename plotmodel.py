import pylab
import glob

files = glob.glob('model_*.txt')
for file in files:
	data = pylab.loadtxt(file)
	pylab.plot(data[:,0],data[:,1],label=file)
pylab.legend()
pylab.xlabel('Phase')
pylab.ylabel('Flux')
pylab.show()
