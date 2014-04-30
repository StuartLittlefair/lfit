import numpy as np
from scipy.interpolate import interp2d, SmoothBivariateSpline

def ld (band,logg,teff):
    assert band in ['u','g','r','i','z']
    filename= 'Gianninas13/ld_coeffs_%s.txt' % band
    data=np.loadtxt(filename)
    x=data[:,0] #logg
    y=data[:,1] #teff
    z=data[:,2] #linear ld coefficient
    #func=interp2d(x,y,z,kind='linear', bounds_error=True)
    func = SmoothBivariateSpline(x,y,z)
    return func(logg,teff)[0]

def main():
    logg, gerr = raw_input('> Give log g and error: ').split()
    teff, terr = raw_input('> Give eff. temp. and error: ').split()
    logg = float(logg); gerr = float(gerr)
    teff = float(teff); terr = float(terr)

    gvals=np.random.normal(loc=logg,scale=gerr,size=100)
    tvals=np.random.normal(loc=teff,scale=terr,size=100)

    #ldvals = []
    #for g,t in zip(gvals,tvals):
    #    ldvals.extend( ld('i',g,t) )
    for band in ['u','g','r','i','z']:
        ldvals = [ld(band,g,t) for g,t in zip(gvals,tvals)]
        print '%s band LD coeff = %f +/- %f' % (band, np.median(ldvals),np.std(ldvals))


if __name__ == "__main__":
    main()
