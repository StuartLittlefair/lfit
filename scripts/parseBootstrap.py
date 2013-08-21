#!/usr/bin/env python

from scipy import io
from scipy import stats
from math import sqrt
from decimal import Decimal

file = raw_input('>Give name of bootstrap results file:  ')

data = io.read_array(file)
chisq = data[:,0]
q = data[:,1]
dphi = data[:,2]
rd = data[:,3]
rwd = data[:,4]
ulimb = data[:,5]
bsScale = data[:,6]
bsAz = data[:,7]
bsFis = data[:,8]
dExp = data[:,9]
incl = data[:,10]
phi0 = data[:,11]

def printParam(par,parString):
    mean = stats.mean(par)
    sd   = sqrt(stats.var(par))
    print parString,' = ', mean, ' +/- ', sd

printParam(q,'Q')
printParam(dphi,'DPHI')
printParam(rd,'RD/L1')
printParam(rwd,'RWD/L1')
printParam(ulimb,'ULIMB')
printParam(bsScale,'SCALE')
printParam(bsAz,'AZ')
printParam(bsFis,'FIS')
printParam(dExp,'EXP')
printParam(incl,'INCL')
printParam(phi0,'PHI0')
