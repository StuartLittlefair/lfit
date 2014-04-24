import numpy
import emcee
from mcmc_utils import *
import scipy.interpolate as interp
import matplotlib.pyplot as plt

def model(pars,rind=6):
    teff, logg, dist, ebv = pars
    
    loggs = np.array([7.0,7.5,8.0,8.5,9.0])
    teffs = np.array([   5500.,    6000.,    6500.,    7000.,    7500.,    8000., \
        8500.,    9000.,    9500.,   10000.,   10500.,   11000., \
        11500.,   12000.,   12500.,   13000.,   13500.,   14000., \
        14500.,   15000.,   15500.,   16000.,   16500.,   17000., \
        20000.,   25000.,   30000.,   35000.,   40000.,   45000., \
        50000.,   55000.,   60000.,   65000.,   70000.,   75000., \
        80000.,   85000.,   90000.,   95000.,  100000.])        
    nteff = len(teffs)
    nlogg = len(loggs)
    assert teff <= teffs.max()
    assert teff >= teffs.min()
    assert logg >= loggs.min()
    assert logg <= loggs.max()
        
    # find index of position just below
    try:
        teff_indx = np.argmax(np.nonzero(teffs-teff < 0)[0])
    except:
        #only fails when we want first element (teff = teffs[0])
        teff_indx=0
    try:
        logg_indx = np.argmax(np.nonzero(loggs-logg < 0)[0])    
    except:
        logg_indx=0
    
    # load bergeron models
    data = np.loadtxt('Bergeron/da2.ugriz')

    indx0 = logg_indx*nteff + teff_indx    
    indx1 = (1+logg_indx)*nteff + teff_indx    
    indx2 = logg_indx*nteff + teff_indx + 1     
    indx3 = (1+logg_indx)*nteff + teff_indx + 1


    x =  np.array([teffs[teff_indx],teffs[teff_indx+1]])
    y =  np.array([loggs[logg_indx],loggs[logg_indx+1]])

    abs_mags = []
    # u data in col 4, g in col 5, red in rind (r=6, i=7, z=8)
    for col_indx in [4,5,rind]:
        z = np.zeros((2,2))
        z[0,0] = data[indx0,col_indx:col_indx+1]
        z[0,1] = data[indx1,col_indx:col_indx+1]
        z[1,0] = data[indx2,col_indx:col_indx+1]
        z[1,1] = data[indx3,col_indx:col_indx+1]
        
        abs_mags.append(interp.RectBivariateSpline(x,y,z,kx=1,ky=1)(teff,logg)[0,0])
    abs_mags = np.array(abs_mags)

    # A_x/E(B-V) extinction from Cardelli (1998)
    r_ext_arr = [2.751, 2.086, 1.479]
    r_ext     = r_ext_arr[rind-6]
    ext       = ebv*np.array([5.155,3.793,r_ext])
    dmod      = 5.0*np.log10(dist/10.0)
    app_red_mags = abs_mags + ext + dmod
    #return app_red_mags
    return 3631e3*10**(-0.4*app_red_mags)

def ln_prior(pars):
    lnp = 0.0
    #teff, logg, dist, reddening

    #teff, uniform between allowed range (6 to 90,000)
    prior = Prior('log_uniform',6000.,60000.)
    lnp += prior.ln_prob(pars[0])

    #logg, uniform between allowed range (7.01 to 8.99)
    # or, in allowed range
    #prior = Prior('uniform',7.01,8.99)
    prior = Prior('gauss',8.0,0.2)
    lnp += prior.ln_prob(pars[1])
    
    # distance, uniform between 50 and 10,000 pc
    # (this is biassed against real distances vs actual prior)
    # so we scale by volume of thin radius step dr (prop. to r**2/50**2)
    prior = Prior('uniform',50,10000)
    #lnp += (pars[2]/50)**2 * prior.ln_prob(pars[2])
    lnp +=  prior.ln_prob(pars[2])

    # reddening, cannot exceed galactic value of 0.121
    prior = Prior('uniform',0.0,0.02)
    lnp += prior.ln_prob(pars[3])
    return lnp
    
def chisq(pars,y,yerr,rind):
    try:
        resids = (y - model(pars,rind))/ yerr
        return np.sum(resids*resids)
    except:
        return np.inf
        
def ln_likelihood(pars,y,yerr,rind):
    errs = yerr
    return -0.5*(np.sum( np.log( 2.0*np.pi*errs**2 ) ) + chisq(pars,y,errs,rind))
    
def ln_prob(pars,y,yerr,rind):
    lnp = ln_prior(pars)
    if np.isfinite(lnp):
        return lnp + ln_likelihood(pars,y,yerr,rind)
    else:
        return lnp

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Fit WD Fluxes')
    parser.add_argument('--fit','-f',action='store_true',dest='fit')
    parser.add_argument('--nwalkers',action='store',type=int,default=50)
    parser.add_argument('--nburn',action='store',type=int,default=100)
    parser.add_argument('--nprod',action='store',type=int,default=500)
    parser.add_argument('--nthread',action='store',type=int,default=4)
    args = parser.parse_args()
        
    u,ue = raw_input("Give u' band flux and err: ").split()
    g,ge = raw_input("Give g' band flux and err: ").split()
    needBand = True
    redband = 'r'
    redindex = 6
    while needBand:
        redband = raw_input("which band in red? ").lower()
        if redband in ['r','i','z']:
            needBand = False
            redindex = 6 + ['r','i','z'].index(redband)
    r,re = raw_input("Give %s' band flux and err: " % redband).split()

    u = float(u); ue = float(ue)
    g = float(g); ge = float(ge)
    r = float(r); re = float(re)

    # calculate mags and errors by brute force!
    #nump = 1e5
    #sys_err = 0.00
    #um = -2.5*np.log10( np.random.normal(loc=u,scale=ue,size=nump) / 3631.0 / 1000.0 ) 
    #um, ume = um.mean(), um.std() + sys_err

    #gm = -2.5*np.log10( np.random.normal(loc=g,scale=ge,size=nump) / 3631.0 / 1000.0 ) 
    #gm, gme = gm.mean(), gm.std() + sys_err

    #rm = -2.5*np.log10( np.random.normal(loc=r,scale=ue,size=nump) / 3631.0 / 1000.0 ) 
    #rm, rme = rm.mean(), rm.std() + sys_err

    um = u; ume = ue
    gm = g; gme = ge
    rm = r; rme = re
    print um, ume
    print gm, gme
    print rm, rme
    y = np.array([um,gm,rm])
    e = np.array([ume,gme,rme])
    #print model([23500,8.0,1650.0,0.01],redindex), chisq([23500,8.0,1650.0,0.01],y,e,redindex)
    #for teff in np.linspace(10000,30000,20):
    #    print teff, chisq([teff,8.0,1650.0,0.01],y,e,redindex)

    guessP = np.array([23500,8.0,1650.0,0.01])
    nameList = ['Teff','log g','d','egi']
    npars = len(guessP)
    nwalkers = args.nwalkers
    p0 = emcee.utils.sample_ball(guessP,0.01*guessP,size=nwalkers)
    sampler = emcee.EnsembleSampler(nwalkers,npars,ln_prob,args=[y,e,redindex],threads=args.nthread)
    
    #burnIn
    nburn = args.nburn
    pos, prob, state = sampler.run_mcmc(p0,nburn)
    
    #production
    sampler.reset()
    nprod = args.nprod
    sampler = run_mcmc_save(sampler,pos,nprod,state,"chain.txt")  
    chain = flatchain(sampler.chain,npars,thin=4)
    
    bestPars = []
    for i in range(npars):
        par = chain[:,i]
        lolim,best,uplim = np.percentile(par,[16,50,84])
        print "%s = %f +%f -%f" % (nameList[i],best,uplim-best,best-lolim)
        bestPars.append(best)
    fig = thumbPlot(chain,nameList)
    fig.savefig('cornerPlot.pdf')
    plt.close()
    
    
